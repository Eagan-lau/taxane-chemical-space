from __future__ import annotations

import math
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from .chemistry import as_float, as_int, parse_reaction_delta_fingerprint, require_rdkit
from .io import ensure_dir, read_table, write_json, write_table


def _text(row: dict[str, Any], key: str) -> str:
    value = row.get(key, "")
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() == "nan" else value


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "t"}


@lru_cache(maxsize=200_000)
def _product_pattern_metrics(smarts: str) -> dict[str, int]:
    _Chem, AllChem, *_rest = require_rdkit()
    reaction = AllChem.ReactionFromSmarts(smarts)
    if reaction is None or reaction.GetNumProductTemplates() != 1:
        return {
            "generic_product_atoms": -1,
            "product_pattern_heavy_atoms": -1,
        }
    product = reaction.GetProductTemplate(0)
    return {
        "generic_product_atoms": sum(
            1 for atom in product.GetAtoms() if atom.GetAtomicNum() == 0
        ),
        "product_pattern_heavy_atoms": sum(
            1 for atom in product.GetAtoms() if atom.GetAtomicNum() > 1
        ),
    }


def _structural_gate(
    row: dict[str, Any],
    *,
    max_heavy_atom_gain: int,
    max_product_pattern_growth: int,
) -> list[str]:
    reasons: list[str] = []
    metrics = _product_pattern_metrics(_text(row, "reaction_smarts"))
    row.update(metrics)
    if metrics["generic_product_atoms"] < 0:
        reasons.append("product_pattern_unreadable")
    if metrics["generic_product_atoms"] > 0:
        reasons.append("generic_product_atom")
    reactant_atoms = as_int(row.get("reactant_atoms"), 0) or 0
    product_atoms = metrics["product_pattern_heavy_atoms"]
    if product_atoms - reactant_atoms > max_product_pattern_growth:
        reasons.append("excessive_product_pattern_growth")

    delta = parse_reaction_delta_fingerprint(
        _text(row, "structural_element_delta")
        or _text(row, "effective_reaction_delta")
        or _text(row, "reaction_delta_fingerprint")
    )
    heavy_gain = sum(
        max(0, count) for element, count in delta.items() if element != "H"
    )
    row["expected_heavy_atom_gain"] = heavy_gain
    if heavy_gain > max_heavy_atom_gain:
        reasons.append("excessive_expected_heavy_atom_gain")
    if as_int(row.get("inferred_independent_reaction_centers"), 0) != 1:
        reasons.append("not_single_reaction_center")
    if as_float(row.get("mapped_atom_retention"), 0.0) < 0.70:
        reasons.append("template_mapping_retention_below_0.70")
    return reasons


def _evidence_gate(row: dict[str, Any], tier: str) -> list[str]:
    reasons: list[str] = []
    reaction_type = _text(row, "reaction_type").lower()
    if not reaction_type or reaction_type.startswith("unassigned"):
        reasons.append("reaction_semantics_unassigned")

    granularity = _text(row, "biochemical_step_granularity").lower()
    source_annotated = _text(row, "single_center_evidence_mode") == "source_annotation"
    consensus_mode = bool(_text(row, "consensus_generation_mode"))
    consensus_pass = _text(row, "consensus_qc_status").lower().startswith("pass")
    validated_consensus = consensus_mode and consensus_pass
    likely_single = (
        "single" in granularity
        and "multi" not in granularity
        and "composite" not in granularity
    )
    database_count = as_int(row.get("consensus_source_database_count"), 0) or 0
    sources = {
        token.strip()
        for token in _text(row, "template_sources").split(";")
        if token.strip()
    }
    non_retrorules_source = bool(sources - {"RetroRules"})
    biological_sources = {
        "Rhea",
        "BioNaviNP_BioChem",
        "MetaNetX",
        "KEGG",
        "TaxolKnownPathway_Curated",
    }
    biological_source_support = bool(sources & biological_sources)
    smarts_recomputed = (
        _text(row, "single_center_evidence_mode") == "smarts_recomputed"
    )
    explicitly_multistep = "multi" in granularity or "composite" in granularity
    semantics_assigned = bool(reaction_type) and not reaction_type.startswith(
        "unassigned"
    )
    primary_supported = (
        (source_annotated and likely_single)
        or (validated_consensus and likely_single)
        or (likely_single and database_count >= 2)
        or (
            smarts_recomputed
            and biological_source_support
            and not explicitly_multistep
            and semantics_assigned
        )
    )

    if tier == "primary":
        if not primary_supported:
            reasons.append("insufficient_single_step_evidence_for_primary")
    elif tier == "extended":
        if not (
            primary_supported
            or source_annotated
            or validated_consensus
            or likely_single
            or (non_retrorules_source and database_count >= 1)
        ):
            reasons.append("insufficient_biochemical_evidence_for_extended")
    else:
        raise ValueError(f"Unknown selection tier: {tier}")
    return reasons


def _taxane_selection_score(row: dict[str, Any]) -> float:
    base = as_float(row.get("grammar_selection_score"), 0.0)
    reactant_atoms = as_int(row.get("reactant_atoms"), 0) or 0
    product_atoms = as_int(row.get("product_pattern_heavy_atoms"), 0) or 0
    match_count = as_int(row.get("g0_match_count"), 0) or 0
    context_bonus = 0.015 * min(reactant_atoms, 20)
    breadth_penalty = 0.04 * math.log1p(match_count)
    growth_penalty = 0.01 * max(0, product_atoms - reactant_atoms - 16)
    return round(base + context_bonus - breadth_penalty - growth_penalty, 6)


def _select_tier(
    activated: pd.DataFrame,
    *,
    tier: str,
    representatives_per_group: int,
    max_heavy_atom_gain: int,
    max_product_pattern_growth: int,
) -> tuple[pd.DataFrame, pd.DataFrame, Counter[str]]:
    accepted: list[dict[str, Any]] = []
    audited: list[dict[str, Any]] = []
    reason_counts: Counter[str] = Counter()
    for row in activated.fillna("").to_dict("records"):
        reasons = _structural_gate(
            row,
            max_heavy_atom_gain=max_heavy_atom_gain,
            max_product_pattern_growth=max_product_pattern_growth,
        )
        reasons.extend(_evidence_gate(row, tier))
        reasons = sorted(set(reasons))
        row[f"{tier}_selection_status"] = "eligible" if not reasons else "excluded"
        row[f"{tier}_selection_reasons"] = ";".join(reasons)
        row["taxane_selection_score"] = _taxane_selection_score(row)
        audited.append(row)
        if reasons:
            reason_counts.update(reasons)
        else:
            accepted.append(row)

    candidates = pd.DataFrame(accepted)
    if candidates.empty:
        return candidates, pd.DataFrame(audited), reason_counts
    candidates = candidates.sort_values(
        [
            "semantic_group_id",
            "taxane_selection_score",
            "reactant_atoms",
            "reaction_smarts_hash",
        ],
        ascending=[True, False, False, True],
        kind="stable",
    )
    candidates[f"{tier}_semantic_group_rank"] = (
        candidates.groupby("semantic_group_id").cumcount() + 1
    )
    selected = candidates[
        candidates[f"{tier}_semantic_group_rank"] <= representatives_per_group
    ].copy()
    selected["taxane_grammar_tier"] = tier
    selected["taxane_grammar_rule_id"] = [
        f"TAXANE_{tier.upper()}_{index:06d}"
        for index in range(1, len(selected) + 1)
    ]
    return selected, pd.DataFrame(audited), reason_counts


def select_taxane_activated_grammar(
    activated_path: Path,
    output_dir: Path,
    *,
    representatives_per_group: int = 3,
    max_heavy_atom_gain: int = 24,
    max_product_pattern_growth: int = 32,
) -> dict[str, Path]:
    output_dir = ensure_dir(output_dir)
    activated = read_table(activated_path).fillna("")
    primary, primary_audit, primary_reasons = _select_tier(
        activated,
        tier="primary",
        representatives_per_group=representatives_per_group,
        max_heavy_atom_gain=max_heavy_atom_gain,
        max_product_pattern_growth=max_product_pattern_growth,
    )
    extended, extended_audit, extended_reasons = _select_tier(
        activated,
        tier="extended",
        representatives_per_group=representatives_per_group,
        max_heavy_atom_gain=max_heavy_atom_gain,
        max_product_pattern_growth=max_product_pattern_growth,
    )
    paths = {
        "primary": output_dir / "taxane_activated_grammar.primary.tsv",
        "extended": output_dir / "taxane_activated_grammar.extended.tsv",
        "primary_audit": output_dir / "taxane_activated_grammar.primary_audit.tsv",
        "extended_audit": output_dir / "taxane_activated_grammar.extended_audit.tsv",
        "summary": output_dir / "taxane_activated_grammar.summary.json",
    }
    write_table(primary, paths["primary"])
    write_table(extended, paths["extended"])
    write_table(primary_audit, paths["primary_audit"])
    write_table(extended_audit, paths["extended_audit"])
    summary = {
        "mode": "G0_activated_evidence_and_structure_gated_grammar",
        "activated_input": str(activated_path),
        "activated_rules": int(len(activated)),
        "activated_semantic_groups": int(activated["semantic_group_id"].nunique()),
        "primary_rules": int(len(primary)),
        "primary_semantic_groups": int(
            primary.get("semantic_group_id", pd.Series(dtype=str)).nunique()
        ),
        "extended_rules": int(len(extended)),
        "extended_semantic_groups": int(
            extended.get("semantic_group_id", pd.Series(dtype=str)).nunique()
        ),
        "representatives_per_group": representatives_per_group,
        "max_heavy_atom_gain": max_heavy_atom_gain,
        "max_product_pattern_growth": max_product_pattern_growth,
        "primary_exclusion_reason_counts": dict(sorted(primary_reasons.items())),
        "extended_exclusion_reason_counts": dict(sorted(extended_reasons.items())),
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths


def assemble_open_grammar(
    global_selected_path: Path,
    g0_selected_path: Path,
    output_dir: Path,
    *,
    representatives_per_group: int = 4,
) -> dict[str, Path]:
    """Combine global and G0-domain representatives for dynamic activation."""

    output_dir = ensure_dir(output_dir)
    global_selected = read_table(global_selected_path).fillna("")
    g0_selected = read_table(g0_selected_path).fillna("")
    required = {"reaction_smarts_hash", "semantic_group_id"}
    for label, frame in (
        ("global_selected", global_selected),
        ("g0_selected", g0_selected),
    ):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{label} is missing required columns: {sorted(missing)}")

    global_hashes = set(global_selected["reaction_smarts_hash"].astype(str))
    g0_hashes = set(g0_selected["reaction_smarts_hash"].astype(str))
    combined = pd.concat(
        [
            g0_selected.assign(_selection_source="G0_domain"),
            global_selected.assign(_selection_source="global"),
        ],
        ignore_index=True,
        sort=False,
    )
    combined["_source_priority"] = combined["_selection_source"].map(
        {"G0_domain": 0, "global": 1}
    )
    score_column = (
        "taxane_selection_score"
        if "taxane_selection_score" in combined.columns
        else "grammar_selection_score"
    )
    combined["_score"] = pd.to_numeric(
        combined.get(score_column, 0.0), errors="coerce"
    ).fillna(0.0)
    combined["_reactant_atoms"] = pd.to_numeric(
        combined.get("reactant_atoms", 0), errors="coerce"
    ).fillna(0)
    combined = combined.sort_values(
        [
            "semantic_group_id",
            "_source_priority",
            "_score",
            "_reactant_atoms",
            "reaction_smarts_hash",
        ],
        ascending=[True, True, False, False, True],
        kind="stable",
    )
    combined = combined.drop_duplicates("reaction_smarts_hash", keep="first")
    combined["open_grammar_group_rank"] = (
        combined.groupby("semantic_group_id").cumcount() + 1
    )
    combined = combined[
        combined["open_grammar_group_rank"] <= representatives_per_group
    ].copy()
    combined["open_grammar_scope"] = combined["reaction_smarts_hash"].map(
        lambda value: (
            "both"
            if value in global_hashes and value in g0_hashes
            else ("G0_domain" if value in g0_hashes else "global")
        )
    )
    combined["open_grammar_rule_id"] = [
        f"OPEN_{index:06d}" for index in range(1, len(combined) + 1)
    ]
    combined = combined.drop(
        columns=[
            "_selection_source",
            "_source_priority",
            "_score",
            "_reactant_atoms",
        ],
        errors="ignore",
    )

    group_summary = (
        combined.groupby(["semantic_group_id", "open_grammar_scope"], dropna=False)
        .size()
        .rename("rule_count")
        .reset_index()
    )
    paths = {
        "open_grammar": output_dir / "taxane_open_grammar.tsv",
        "group_summary": output_dir / "taxane_open_grammar.group_summary.tsv",
        "summary": output_dir / "taxane_open_grammar.summary.json",
    }
    write_table(combined, paths["open_grammar"])
    write_table(group_summary, paths["group_summary"])
    scope_counts = combined["open_grammar_scope"].value_counts().to_dict()
    summary = {
        "mode": "domain_seeded_open_reaction_grammar",
        "global_selected_input": str(global_selected_path),
        "g0_selected_input": str(g0_selected_path),
        "global_selected_rules": int(len(global_selected)),
        "g0_selected_rules": int(len(g0_selected)),
        "open_grammar_rules": int(len(combined)),
        "open_grammar_semantic_groups": int(
            combined["semantic_group_id"].nunique()
        ),
        "representatives_per_group": representatives_per_group,
        "scope_counts": {
            str(key): int(value) for key, value in sorted(scope_counts.items())
        },
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths


def augment_open_grammar_with_domain(
    open_grammar_path: Path,
    domain_grammar_path: Path,
    output_dir: Path,
) -> dict[str, Path]:
    """Add reviewed family-domain rules while preserving source separation."""

    output_dir = ensure_dir(output_dir)
    external = read_table(open_grammar_path).fillna("")
    domain = read_table(domain_grammar_path).fillna("")
    required = {"reaction_smarts_hash", "reaction_smarts"}
    for label, frame in (("open_grammar", external), ("domain_grammar", domain)):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{label} is missing required columns: {sorted(missing)}")

    external_hashes = set(external["reaction_smarts_hash"].astype(str))
    domain_hashes = set(domain["reaction_smarts_hash"].astype(str))
    combined = pd.concat(
        [
            domain.assign(_domain_priority=0),
            external.assign(_domain_priority=1),
        ],
        ignore_index=True,
        sort=False,
    )
    combined = combined.sort_values(
        ["_domain_priority", "reaction_smarts_hash"],
        kind="stable",
    ).drop_duplicates("reaction_smarts_hash", keep="first")
    combined["grammar_provenance_scope"] = combined["reaction_smarts_hash"].map(
        lambda value: (
            "external_and_taxane_domain"
            if value in external_hashes and value in domain_hashes
            else ("taxane_domain" if value in domain_hashes else "external")
        )
    )
    combined["final_grammar_rule_id"] = [
        f"FINAL_{index:06d}" for index in range(1, len(combined) + 1)
    ]
    combined = combined.drop(columns=["_domain_priority"], errors="ignore")
    paths = {
        "grammar": output_dir / "taxane_reaction_grammar.primary.tsv",
        "summary": output_dir / "taxane_reaction_grammar.primary.summary.json",
    }
    write_table(combined, paths["grammar"])
    scope_counts = combined["grammar_provenance_scope"].value_counts().to_dict()
    summary = {
        "mode": "external_open_plus_reviewed_taxane_domain_grammar",
        "external_open_rules": int(len(external)),
        "reviewed_domain_rules": int(len(domain)),
        "final_rules": int(len(combined)),
        "final_semantic_groups": int(
            combined.get("semantic_group_id", pd.Series(dtype=str)).nunique()
        ),
        "scope_counts": {
            str(key): int(value) for key, value in sorted(scope_counts.items())
        },
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths
