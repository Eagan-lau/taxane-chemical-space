from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .source_layers import best_layer, normalize_evidence_layers
from .utils import clean_text, ensure_dir, join_values, sha256_text, split_multi_value, truthy, write_json


CONSENSUS_GROUP_COLUMNS = [
    "consensus_group_id",
    "reaction_type",
    "main_functional_group_changes",
    "transferred_group_class",
    "acceptor_atom_class",
    "donor_class",
    "reaction_delta_fingerprint",
    "direction_qc_status",
    "evidence_rows",
    "source_database_count",
    "source_databases",
    "evidence_layer_support",
    "representative_rule_id",
    "representative_reaction_smarts",
    "representative_confidence",
    "promoted_to_rule",
    "skip_reason",
]


def _truthy_mask(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([False] * len(df), index=df.index)
    return df[col].astype(str).str.lower().isin(["true", "1", "yes"])


def _text(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([""] * len(df), index=df.index, dtype=str)
    return df[col].fillna("").astype(str)


def _support_tokens(series: pd.Series) -> list[str]:
    vals: list[str] = []
    for value in series.fillna("").astype(str).tolist():
        for token in split_multi_value(value):
            if token and token not in vals:
                vals.append(token)
    return vals


def _group_key(row: pd.Series) -> tuple[str, ...]:
    """Consensus grouping key.

    The key is intentionally reaction-data derived.  It avoids enzyme-name or
    Taxol-specific labels, and it keeps direction QC separate so unknown-direction
    transforms cannot silently promote a strict consensus.
    """
    reaction_type = clean_text(row.get("reaction_type", "")) or "unassigned"
    functional = clean_text(row.get("main_functional_group_changes", "")) or clean_text(row.get("reaction_subtype", ""))
    transferred = clean_text(row.get("transferred_group_class", "")) or clean_text(row.get("transferred_group", ""))
    acceptor = clean_text(row.get("acceptor_atom_class", ""))
    donor = clean_text(row.get("donor_class", "")) or clean_text(row.get("cofactor_or_donor_class", ""))
    delta = clean_text(row.get("reaction_delta_fingerprint", ""))
    direction = clean_text(row.get("direction_qc_status", "")) or clean_text(row.get("direction_evidence_type", ""))
    return (
        reaction_type,
        functional,
        transferred,
        acceptor,
        donor,
        delta,
        direction,
    )


def _representative(group: pd.DataFrame) -> pd.Series:
    g = group.copy()
    g["_rank_conf"] = pd.to_numeric(g.get("final_rule_confidence", 0), errors="coerce").fillna(0.0)
    g["_rank_templates"] = pd.to_numeric(g.get("template_count", 0), errors="coerce").fillna(0)
    g["_rank_sources"] = pd.to_numeric(g.get("source_record_count", 0), errors="coerce").fillna(0)
    return g.sort_values(["_rank_conf", "_rank_templates", "_rank_sources"], ascending=[False, False, False]).iloc[0]


def build_data_driven_consensus_rules(
    rules: pd.DataFrame,
    *,
    output_dir: str | Path | None = None,
    min_evidence_rows: int = 3,
    min_source_database_count: int = 1,
    min_template_count: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Promote high-support SMARTS clusters into auditable consensus rules.

    This flow is fully data-driven: every candidate starts from an existing,
    QC-passing generalized SMARTS row, then clusters across all databases by
    reaction metadata and structural delta.  It does not hard-code enzyme-family
    mappings or Taxol-specific chemistry.
    """
    if rules.empty:
        empty_report = pd.DataFrame(columns=CONSENSUS_GROUP_COLUMNS)
        summary = {"enabled": True, "input_rules": 0, "candidate_groups": 0, "promoted_consensus_rules": 0}
        if output_dir:
            out = ensure_dir(output_dir)
            empty_report.to_csv(out / "data_driven_consensus_candidate_groups.tsv", sep="\t", index=False)
            write_json(out / "data_driven_consensus_summary.json", summary)
        return pd.DataFrame(columns=rules.columns), empty_report, summary

    eligible = rules[
        (_text(rules, "template_scope") == "generalized_template")
        & _truthy_mask(rules, "predictive_rule_use")
        & (_text(rules, "template_qc_status") == "ok")
        & _text(rules, "reaction_smarts").str.strip().ne("")
    ].copy()
    if eligible.empty:
        report = pd.DataFrame(columns=CONSENSUS_GROUP_COLUMNS)
        summary = {
            "enabled": True,
            "input_rules": int(len(rules)),
            "eligible_rules": 0,
            "candidate_groups": 0,
            "promoted_consensus_rules": 0,
        }
        if output_dir:
            out = ensure_dir(output_dir)
            report.to_csv(out / "data_driven_consensus_candidate_groups.tsv", sep="\t", index=False)
            write_json(out / "data_driven_consensus_summary.json", summary)
        return pd.DataFrame(columns=rules.columns), report, summary

    eligible["_consensus_key"] = eligible.apply(_group_key, axis=1)
    rows: list[dict[str, Any]] = []
    consensus_rows: list[dict[str, Any]] = []
    for i, (key, group) in enumerate(eligible.groupby("_consensus_key", dropna=False), start=1):
        rep = _representative(group)
        source_tokens = _support_tokens(group.get("template_sources", group.get("source_database", pd.Series(dtype=str))))
        if not source_tokens:
            source_tokens = _support_tokens(group.get("source_database", pd.Series(dtype=str)))
        layer_support = normalize_evidence_layers(join_values(group.get("evidence_layers_all", group.get("evidence_layer_best", pd.Series(dtype=str))).astype(str).tolist()))
        group_id = "CONSENSUS_CLUSTER_" + sha256_text("||".join(key))[:20]
        evidence_rows = int(len(group))
        template_support = int(pd.to_numeric(group.get("template_count", 1), errors="coerce").fillna(1).sum())
        source_db_count = int(len(source_tokens))
        direction_status = clean_text(rep.get("direction_qc_status", ""))
        skip_reason = ""
        promoted = True
        if evidence_rows < int(min_evidence_rows):
            promoted = False
            skip_reason = "insufficient_cluster_rows"
        elif template_support < int(min_template_count):
            promoted = False
            skip_reason = "insufficient_template_support"
        elif source_db_count < int(min_source_database_count):
            promoted = False
            skip_reason = "insufficient_source_database_support"
        elif direction_status == "direction_qc_unknown_exploratory_only":
            promoted = False
            skip_reason = "unknown_direction_consensus_not_promoted"

        rows.append({
            "consensus_group_id": group_id,
            "reaction_type": key[0],
            "main_functional_group_changes": key[1],
            "transferred_group_class": key[2],
            "acceptor_atom_class": key[3],
            "donor_class": key[4],
            "reaction_delta_fingerprint": key[5],
            "direction_qc_status": key[6],
            "evidence_rows": evidence_rows,
            "source_database_count": source_db_count,
            "source_databases": join_values(source_tokens),
            "evidence_layer_support": layer_support,
            "representative_rule_id": clean_text(rep.get("rule_id", "")),
            "representative_reaction_smarts": clean_text(rep.get("reaction_smarts", "")),
            "representative_confidence": clean_text(rep.get("final_rule_confidence", "")),
            "promoted_to_rule": promoted,
            "skip_reason": skip_reason,
        })
        if not promoted:
            continue

        d = rep.to_dict()
        d["rule_id"] = f"CONSENSUS_DATA_DRIVEN_{len(consensus_rows)+1:09d}"
        d["template_hash"] = sha256_text(f"{group_id}::{d.get('reaction_smarts','')}")[:20]
        d["template_scope"] = "generalized_template"
        d["predictive_rule_use"] = "true"
        d["anchor_edge_use"] = "false"
        d["template_qc_status"] = "ok"
        d["template_qc_note"] = join_values([d.get("template_qc_note", ""), "data_driven_consensus_representative_qc_ok"])
        d["evidence_layer_best"] = best_layer(layer_support)
        d["evidence_layers_all"] = layer_support
        d["template_sources"] = join_values(source_tokens)
        d["template_count"] = template_support
        d["source_record_count"] = int(pd.to_numeric(group.get("source_record_count", 1), errors="coerce").fillna(1).sum())
        d["consensus_group_id"] = group_id
        d["consensus_generation_mode"] = "data_driven_cluster_representative_smarts"
        d["consensus_evidence_rows"] = evidence_rows
        d["consensus_source_database_count"] = source_db_count
        d["consensus_evidence_layer_support"] = layer_support
        d["consensus_qc_status"] = "passed_data_driven_thresholds"
        d["consensus_representative_rule_ids"] = clean_text(rep.get("rule_id", ""))
        d["consensus_supporting_reaction_types"] = join_values(group.get("reaction_type", pd.Series(dtype=str)).astype(str).tolist())
        d["curated_taxol_anchor"] = "false"
        d["curated_pathway_name"] = ""
        d["curated_pathway_step_ids"] = ""
        d["strict_core_use"] = bool(d["evidence_layer_best"] == "T1_Bio_Core" and truthy(d.get("strict_core_use", "")))
        d["expanded_use"] = bool(d["evidence_layer_best"] in {"T1_Bio_Core", "T2_Bio_Extended"} and truthy(d.get("expanded_use", "")))
        d["exploratory_use"] = True
        d["exclusive_release_tier"] = "T1_only" if truthy(d.get("strict_core_use", "")) else ("T2_only" if truthy(d.get("expanded_use", "")) else "T3_only")
        d["notes"] = join_values([d.get("notes", ""), f"data_driven_consensus_group={group_id}", f"consensus_evidence_rows={evidence_rows}"])
        consensus_rows.append(d)

    report = pd.DataFrame(rows, columns=CONSENSUS_GROUP_COLUMNS).fillna("")
    consensus = pd.DataFrame(consensus_rows).fillna("")
    if not consensus.empty:
        missing_cols = [col for col in rules.columns if col not in consensus.columns]
        if missing_cols:
            consensus = pd.concat(
                [consensus, pd.DataFrame("", index=consensus.index, columns=missing_cols)],
                axis=1,
            )
        consensus = consensus.loc[:, list(rules.columns)].copy()
    summary = {
        "enabled": True,
        "input_rules": int(len(rules)),
        "eligible_rules": int(len(eligible)),
        "candidate_groups": int(len(report)),
        "promoted_consensus_rules": int(len(consensus)),
        "min_evidence_rows": int(min_evidence_rows),
        "min_source_database_count": int(min_source_database_count),
        "min_template_count": int(min_template_count),
        "promotion_policy": "QC-passing generalized SMARTS clustered by reaction type, functional-group delta, transfer role, acceptor/donor class and direction QC.",
    }
    if output_dir:
        out = ensure_dir(output_dir)
        report.to_csv(out / "data_driven_consensus_candidate_groups.tsv", sep="\t", index=False)
        consensus.to_csv(out / "data_driven_consensus_rules.tsv", sep="\t", index=False)
        write_json(out / "data_driven_consensus_summary.json", summary)
    return consensus, report, summary
