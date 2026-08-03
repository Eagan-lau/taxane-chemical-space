from __future__ import annotations

import json
import pandas as pd

from .schema import SMARTS_LIBRARY_COLUMNS
from .utils import clean_text, sha256_text, truthy, write_json


NETWORK_READY_QC = {"ok"}


def _bool_mask(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([False] * len(df), index=df.index)
    return df[col].astype(str).str.lower().isin(["true", "1", "yes"])


def _base_smarts_mask(rules: pd.DataFrame, *, require_qc_ok: bool = True) -> pd.Series:
    if rules.empty:
        return pd.Series([], dtype=bool)
    mask = (
        (rules.get("template_scope", "").astype(str) == "generalized_template")
        & _bool_mask(rules, "predictive_rule_use")
        & rules.get("reaction_smarts", "").astype(str).str.strip().ne("")
    )
    if require_qc_ok:
        mask = mask & rules.get("template_qc_status", "").astype(str).isin(NETWORK_READY_QC)
    return mask


def _exclusive_tier_series(rules: pd.DataFrame) -> pd.Series:
    """Assign every predictive SMARTS to one contribution tier.

    The assignment is hierarchical and mutually exclusive: a rule that qualifies
    for T1 is counted only as T1; otherwise T2 is considered; otherwise T3.
    This keeps contribution analysis distinct from cumulative sensitivity
    releases.
    """
    if rules.empty:
        return pd.Series([], dtype=str)
    strict = _bool_mask(rules, "strict_core_use")
    expanded = _bool_mask(rules, "expanded_use")
    exploratory = _bool_mask(rules, "exploratory_use")
    tier = pd.Series([""] * len(rules), index=rules.index, dtype=object)
    tier.loc[strict] = "T1_only"
    tier.loc[~strict & expanded] = "T2_only"
    tier.loc[~strict & ~expanded & exploratory] = "T3_only"
    return tier


def export_reaction_smarts_library(
    rules: pd.DataFrame,
    *,
    tier: str = "all",
    require_qc_ok: bool = True,
) -> pd.DataFrame:
    """Return a network-ready reaction SMARTS library.

    This exporter intentionally excludes exact substrate-product anchors. It emits only
    generalized predictive reaction SMARTS records. The resulting table is the file
    intended for downstream network construction.
    """
    if rules.empty:
        return pd.DataFrame(columns=SMARTS_LIBRARY_COLUMNS)

    mask = _base_smarts_mask(rules, require_qc_ok=require_qc_ok)
    exclusive_tier = _exclusive_tier_series(rules)
    tier_low = clean_text(tier).lower()
    if tier_low in {"core", "t1", "t1_core", "t1_only"}:
        mask = mask & (exclusive_tier == "T1_only")
    elif tier_low in {"expanded", "t2", "t2_extended", "t2_only"}:
        mask = mask & (exclusive_tier == "T2_only")
    elif tier_low in {"exploratory", "t3", "t3_exploratory", "t3_only"}:
        mask = mask & (exclusive_tier == "T3_only")
    elif tier_low in {"t2_cumulative", "expanded_cumulative"}:
        mask = mask & _bool_mask(rules, "expanded_use")
    elif tier_low in {"t3_cumulative", "exploratory_cumulative"}:
        mask = mask & _bool_mask(rules, "exploratory_use")
    elif tier_low in {"all", "all_smarts", "network_ready_all"}:
        pass
    else:
        raise ValueError(f"Unknown SMARTS library tier: {tier!r}")

    out = rules.loc[mask].copy().reset_index(drop=True)
    if out.empty:
        return pd.DataFrame(columns=SMARTS_LIBRARY_COLUMNS)

    out["reaction_smarts_hash"] = out["reaction_smarts"].map(lambda x: sha256_text(str(x))[:20])
    out["smarts_library_tier"] = tier_low
    out["exclusive_release_tier"] = _exclusive_tier_series(out)
    out["smarts_rule_id"] = [f"SMRT_{tier_low.upper()}_{i+1:09d}" for i in range(len(out))]
    # A final safety invariant: anchors must never appear in SMARTS release files.
    out = out[(out["template_scope"] == "generalized_template") & out["reaction_smarts"].astype(str).str.strip().ne("")]
    missing_cols = [col for col in SMARTS_LIBRARY_COLUMNS if col not in out.columns]
    if missing_cols:
        out = pd.concat([out, pd.DataFrame("", index=out.index, columns=missing_cols)], axis=1)
    return out.loc[:, SMARTS_LIBRARY_COLUMNS].copy()



def validate_smarts_release(all_smarts: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    issues: list[dict] = []
    for _, row in all_smarts.iterrows():
        rid = str(row.get("rule_id", ""))
        if str(row.get("template_scope", "")) != "generalized_template":
            issues.append({"rule_id": rid, "validation_issue": "non_generalized_template_scope"})
        if not truthy(row.get("predictive_rule_use", "")):
            issues.append({"rule_id": rid, "validation_issue": "predictive_rule_use_false"})
        if not str(row.get("reaction_smarts", "")).strip():
            issues.append({"rule_id": rid, "validation_issue": "empty_reaction_smarts"})
        if truthy(row.get("abstracted_from_exact_reaction", "")):
            if not truthy(row.get("abstracted_smarts_applies_to_original_pair", "")):
                issues.append({"rule_id": rid, "validation_issue": "exact_derived_smarts_missing_replay_validation"})
            if str(row.get("exact_abstraction_qc_status", "")).lower() not in {"pass", "ok"}:
                issues.append({"rule_id": rid, "validation_issue": "exact_abstraction_qc_not_pass"})
    issues_df = pd.DataFrame(issues, columns=["rule_id", "validation_issue"])
    summary = {
        "smarts_rule_records": int(len(all_smarts)),
        "validation_issue_records": int(len(issues_df)),
        "is_valid_smarts_release": bool(len(issues_df) == 0),
        "contains_exact_anchors": bool(any(all_smarts.get("template_scope", pd.Series(dtype=str)).astype(str) != "generalized_template")) if not all_smarts.empty else False,
        "contains_empty_reaction_smarts": bool(any(all_smarts.get("reaction_smarts", pd.Series(dtype=str)).astype(str).str.strip() == "")) if not all_smarts.empty else False,
        "note": "Empty SMARTS releases are valid for anchor-only builds, but insufficient for downstream network construction." if len(all_smarts) == 0 else "SMARTS release contains only generalized predictive reaction SMARTS records; exact-derived SMARTS must pass replay validation."
    }
    return issues_df, summary

def reaction_smarts_library_summary(smarts_df: pd.DataFrame) -> dict:
    if smarts_df.empty:
        return {
            "reaction_smarts_rules": 0,
            "unique_reaction_smarts": 0,
            "unique_template_sources": 0,
        }
    return {
        "reaction_smarts_rules": int(len(smarts_df)),
        "unique_reaction_smarts": int(smarts_df["reaction_smarts"].nunique()),
        "unique_template_sources": int(smarts_df.get("template_sources", pd.Series(dtype=str)).astype(str).nunique()),
    }


def write_smarts_releases(rules: pd.DataFrame, release_dir: str | Path) -> dict:
    """Write final SMARTS-only libraries and return count summary.

    The canonical downstream files are `reaction_smarts_library.*.tsv`. Aliases
    named `reaction_smarts_rules.*.tsv` are also written for readability.
    """
    from .utils import ensure_dir
    release = ensure_dir(release_dir)
    all_smarts = export_reaction_smarts_library(rules, tier="all", require_qc_ok=True)
    core = export_reaction_smarts_library(rules, tier="T1_core", require_qc_ok=True)
    expanded = export_reaction_smarts_library(rules, tier="T2_extended", require_qc_ok=True)
    exploratory = export_reaction_smarts_library(rules, tier="T3_exploratory", require_qc_ok=True)
    t2_cumulative = export_reaction_smarts_library(rules, tier="T2_cumulative", require_qc_ok=True)
    t3_cumulative = export_reaction_smarts_library(rules, tier="T3_cumulative", require_qc_ok=True)

    # Canonical output names for the next graph-construction step.
    all_smarts.to_csv(release / "reaction_smarts_library.all.tsv", sep="\t", index=False)
    core.to_csv(release / "reaction_smarts_library.T1_core.tsv", sep="\t", index=False)
    expanded.to_csv(release / "reaction_smarts_library.T2_extended.tsv", sep="\t", index=False)
    exploratory.to_csv(release / "reaction_smarts_library.T3_exploratory.tsv", sep="\t", index=False)
    core.to_csv(release / "reaction_smarts_library.T1_only.tsv", sep="\t", index=False)
    expanded.to_csv(release / "reaction_smarts_library.T2_only.tsv", sep="\t", index=False)
    exploratory.to_csv(release / "reaction_smarts_library.T3_only.tsv", sep="\t", index=False)
    t2_cumulative.to_csv(release / "reaction_smarts_library.T2_cumulative.tsv", sep="\t", index=False)
    t3_cumulative.to_csv(release / "reaction_smarts_library.T3_cumulative.tsv", sep="\t", index=False)
    core.to_csv(release / "reaction_smarts_library.core.tsv", sep="\t", index=False)
    expanded.to_csv(release / "reaction_smarts_library.expanded.tsv", sep="\t", index=False)
    exploratory.to_csv(release / "reaction_smarts_library.exploratory.tsv", sep="\t", index=False)
    core.to_csv(release / "network_ready_reaction_smarts_library.tsv", sep="\t", index=False)

    # Human-readable aliases; same content as above.
    all_smarts.to_csv(release / "reaction_smarts_rules.all.tsv", sep="\t", index=False)
    core.to_csv(release / "reaction_smarts_rules.T1_core.tsv", sep="\t", index=False)
    expanded.to_csv(release / "reaction_smarts_rules.T2_extended.tsv", sep="\t", index=False)
    exploratory.to_csv(release / "reaction_smarts_rules.T3_exploratory.tsv", sep="\t", index=False)
    core.to_csv(release / "reaction_smarts_rules.T1_only.tsv", sep="\t", index=False)
    expanded.to_csv(release / "reaction_smarts_rules.T2_only.tsv", sep="\t", index=False)
    exploratory.to_csv(release / "reaction_smarts_rules.T3_only.tsv", sep="\t", index=False)
    t2_cumulative.to_csv(release / "reaction_smarts_rules.T2_cumulative.tsv", sep="\t", index=False)
    t3_cumulative.to_csv(release / "reaction_smarts_rules.T3_cumulative.tsv", sep="\t", index=False)
    core.to_csv(release / "reaction_smarts_rules.core.tsv", sep="\t", index=False)
    expanded.to_csv(release / "reaction_smarts_rules.expanded.tsv", sep="\t", index=False)
    exploratory.to_csv(release / "reaction_smarts_rules.exploratory.tsv", sep="\t", index=False)
    core.to_csv(release / "network_ready_reaction_smarts_rules.tsv", sep="\t", index=False)

    issues_df, validation_summary = validate_smarts_release(all_smarts)
    issues_df.to_csv(release / "reaction_smarts_rules.validation_issues.tsv", sep="\t", index=False)
    write_json(release / "reaction_smarts_rules.validation.json", validation_summary)

    # Backward-compatible names, but deliberately SMARTS-only. The mixed audit
    # table is only under 03_rules/.
    core.to_csv(release / "general_transformation_rules.core.tsv", sep="\t", index=False)
    expanded.to_csv(release / "general_transformation_rules.expanded.tsv", sep="\t", index=False)
    exploratory.to_csv(release / "general_transformation_rules.exploratory.tsv", sep="\t", index=False)

    summary = {
        "reaction_smarts_library_all": int(len(all_smarts)),
        "reaction_smarts_library_core": int(len(core)),
        "reaction_smarts_library_expanded": int(len(expanded)),
        "reaction_smarts_library_exploratory": int(len(exploratory)),
        "reaction_smarts_library_T1_core": int(len(core)),
        "reaction_smarts_library_T2_extended": int(len(expanded)),
        "reaction_smarts_library_T3_exploratory": int(len(exploratory)),
        "reaction_smarts_library_T1_only": int(len(core)),
        "reaction_smarts_library_T2_only": int(len(expanded)),
        "reaction_smarts_library_T3_only": int(len(exploratory)),
        "reaction_smarts_library_T2_cumulative": int(len(t2_cumulative)),
        "reaction_smarts_library_T3_cumulative": int(len(t3_cumulative)),
        "network_ready_reaction_smarts_library": int(len(core)),
        # aliases used by earlier drafts
        "reaction_smarts_rules_all": int(len(all_smarts)),
        "reaction_smarts_rules_core": int(len(core)),
        "reaction_smarts_rules_expanded": int(len(expanded)),
        "reaction_smarts_rules_exploratory": int(len(exploratory)),
        "network_ready_reaction_smarts_rules": int(len(core)),
        "smarts_release_contains_exact_anchors": False,
        "smarts_release_contains_empty_smarts": False,
        "reaction_smarts_release_valid": bool(validation_summary.get("is_valid_smarts_release", False)),
    }
    issues = []
    for name, df in {"all": all_smarts, "core": core, "expanded": expanded, "exploratory": exploratory}.items():
        if not df.empty:
            bad_scope = df[df["template_scope"].astype(str) != "generalized_template"]
            empty_smarts = df[df["reaction_smarts"].astype(str).str.strip() == ""]
            if len(bad_scope):
                summary["smarts_release_contains_exact_anchors"] = True
                for _, row in bad_scope.iterrows():
                    issues.append({"tier": name, "rule_id": row.get("rule_id", ""), "issue": "non_generalized_template_in_smarts_release"})
            if len(empty_smarts):
                summary["smarts_release_contains_empty_smarts"] = True
                for _, row in empty_smarts.iterrows():
                    issues.append({"tier": name, "rule_id": row.get("rule_id", ""), "issue": "empty_reaction_smarts"})
            exact_derived = df[df.get("abstracted_from_exact_reaction", pd.Series([""] * len(df), index=df.index)).astype(str).str.lower().isin(["true", "1", "yes"])]
            for _, row in exact_derived.iterrows():
                if not truthy(row.get("abstracted_smarts_applies_to_original_pair", "")):
                    issues.append({"tier": name, "rule_id": row.get("rule_id", ""), "issue": "exact_derived_smarts_missing_replay_validation"})
                if str(row.get("exact_abstraction_qc_status", "")).lower() not in {"pass", "ok"}:
                    issues.append({"tier": name, "rule_id": row.get("rule_id", ""), "issue": "exact_abstraction_qc_not_pass"})
    issues_df = pd.DataFrame(issues, columns=["tier", "rule_id", "issue"])
    issues_df.to_csv(release / "reaction_smarts_rules.validation_issues.tsv", sep="\t", index=False)
    validation = {**summary, "validation_issue_records": int(len(issues_df)), "validation_passed": len(issues_df) == 0}
    with open(release / "reaction_smarts_rules.validation.json", "w", encoding="utf-8") as fh:
        json.dump(validation, fh, indent=2, sort_keys=True)
    summary["reaction_smarts_release_valid"] = bool(validation["validation_passed"])
    return summary


# Compatibility with earlier pipeline drafts.
def release_smarts_libraries(rules: pd.DataFrame, release_dir: str | Path) -> dict:
    return write_smarts_releases(rules, release_dir)


def smarts_only_rules(rules: pd.DataFrame, tier: str = "all") -> pd.DataFrame:
    return export_reaction_smarts_library(rules, tier=tier, require_qc_ok=True)
