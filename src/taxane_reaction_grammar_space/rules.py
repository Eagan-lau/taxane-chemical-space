from __future__ import annotations

import json
import math
import time
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from .chemistry import (
    as_float,
    as_int,
    falsey,
    reaction_template_metrics,
    stable_hash,
    truthy,
)
from .io import ensure_dir, write_json, write_table


REQUIRED_INPUT_COLUMNS = {
    "reaction_smarts",
    "smarts_rule_id",
    "reaction_smarts_hash",
    "reaction_delta_fingerprint",
    "reaction_type",
    "exclusive_release_tier",
}

PASSTHROUGH_COLUMNS = [
    "smarts_rule_id",
    "rule_id",
    "reaction_smarts",
    "reaction_smarts_hash",
    "smarts_library_tier",
    "exclusive_release_tier",
    "template_hash",
    "template_scope",
    "predictive_rule_use",
    "template_qc_status",
    "abstracted_from_exact_reaction",
    "derived_from_exact_anchor",
    "abstracted_smarts_applies_to_original_pair",
    "exact_abstraction_qc_status",
    "benchmark_exclusion_flag",
    "reaction_type",
    "reaction_subtype",
    "rule_application_unit",
    "biochemical_step_granularity",
    "biochemical_step_granularity_confidence",
    "composite_rule_flag",
    "reaction_center_count",
    "independent_reaction_center_count",
    "functional_group_change_count",
    "main_functional_group_changes",
    "normalized_direction",
    "smarts_direction",
    "molecular_direction",
    "direction_evidence_type",
    "source_direction",
    "direction_handling",
    "direction_qc_status",
    "direction_qc_note",
    "cofactor_or_donor_class",
    "reaction_representation_scope",
    "transferred_group",
    "donor_class",
    "acceptor_atom_class",
    "transferred_group_class",
    "main_pair_projection_method",
    "main_pair_projection_note",
    "evidence_layer_best",
    "evidence_layers_all",
    "template_sources",
    "curated_taxol_anchor",
    "curated_pathway_name",
    "curated_pathway_step_ids",
    "source_reaction_ids",
    "rhea_ids",
    "kegg_ids",
    "metanetx_ids",
    "consensus_group_id",
    "consensus_generation_mode",
    "consensus_evidence_rows",
    "consensus_source_database_count",
    "consensus_evidence_layer_support",
    "consensus_qc_status",
    "final_rule_confidence",
    "strict_core_use",
    "expanded_use",
    "exploratory_use",
    "template_count",
    "source_record_count",
    "example_reaction_smiles",
    "example_substrate_smiles",
    "example_product_smiles",
    "reaction_delta_fingerprint",
    "notes",
]


def _text(row: dict[str, Any], column: str) -> str:
    value = row.get(column, "")
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() == "nan" else value


def _primary_exclusion_reasons(
    row: dict[str, Any],
    *,
    max_reactant_atoms: int,
    require_single_center: bool,
    release_tier: str,
) -> list[str]:
    reasons: list[str] = []
    if not truthy(row.get("predictive_rule_use", "")):
        reasons.append("not_predictive_rule")
    if _text(row, "template_qc_status").lower() != "ok":
        reasons.append("template_qc_not_ok")
    if release_tier == "T1":
        if not truthy(row.get("strict_core_use", "")):
            reasons.append("not_strict_core")
    elif release_tier == "T2":
        if not truthy(row.get("expanded_use", "")):
            reasons.append("not_expanded_use")
    elif release_tier == "T3":
        if not truthy(row.get("exploratory_use", "")):
            reasons.append("not_exploratory_use")
    else:
        raise ValueError(f"Unsupported release tier: {release_tier}")
    if truthy(row.get("composite_rule_flag", "")):
        reasons.append("composite_rule")

    granularity = _text(row, "biochemical_step_granularity").lower()
    if require_single_center:
        reaction_centers = as_int(row.get("reaction_center_count"), None)
        independent_centers = as_int(
            row.get("independent_reaction_center_count"), None
        )
        functional_changes = as_int(
            row.get("functional_group_change_count"), None
        )
        annotation_present = any(
            value is not None
            for value in (reaction_centers, independent_centers, functional_changes)
        )
        annotation_single = (
            reaction_centers == 1
            and independent_centers == 1
            and functional_changes == 1
            and "single" in granularity
            and not any(token in granularity for token in ("multi", "composite"))
        )
        annotation_explicit_multicenter = any(
            value is not None and value > 1
            for value in (reaction_centers, independent_centers, functional_changes)
        ) or any(token in granularity for token in ("multi", "composite"))
        derived_single = (
            as_int(row.get("inferred_independent_reaction_centers"), 0) == 1
            and as_int(row.get("changed_mapped_atom_count"), 0) > 0
            and as_float(row.get("mapped_atom_retention"), 0.0) >= 0.70
            and as_float(row.get("mapping_coverage"), 0.0) >= 0.30
        )
        if annotation_single:
            row["single_center_evidence_mode"] = "source_annotation"
        elif not annotation_explicit_multicenter and derived_single:
            row["single_center_evidence_mode"] = "smarts_recomputed"
        else:
            row["single_center_evidence_mode"] = (
                "annotation_conflict"
                if annotation_present or annotation_explicit_multicenter
                else "insufficient_evidence"
            )
            reasons.append("single_center_not_verified")

    direction_qc = _text(row, "direction_qc_status").lower()
    if not direction_qc.startswith("direction_qc_"):
        reasons.append("direction_qc_missing")
    elif any(token in direction_qc for token in ("fail", "conflict", "ambiguous")):
        reasons.append("direction_qc_failed")

    exact_status = _text(row, "exact_abstraction_qc_status").lower()
    exact_applies = _text(row, "abstracted_smarts_applies_to_original_pair")
    if exact_status and exact_status != "pass":
        reasons.append("exact_abstraction_not_replay_validated")
    if exact_applies and not truthy(exact_applies):
        reasons.append("abstracted_smarts_does_not_replay")

    if _text(row, "rule_application_unit") != "single_smarts_application":
        reasons.append("unsupported_rule_application_unit")

    if _text(row, "compile_status") != "ok":
        reasons.append("reaction_smarts_compile_failed")
    if as_int(row.get("n_reactants"), 0) != 1:
        reasons.append("not_single_reactant_smarts")
    if as_int(row.get("n_products"), 0) != 1:
        reasons.append("not_single_product_smarts")
    if as_int(row.get("reactant_atoms"), 0) <= 0:
        reasons.append("empty_reactant_template")
    if as_int(row.get("reactant_atoms"), 0) > max_reactant_atoms:
        reasons.append("reactant_context_too_large")
    if as_int(row.get("mapped_reactant_atoms"), 0) <= 0:
        reasons.append("no_mapped_reactant_atoms")
    if as_int(row.get("mapped_product_atoms"), 0) <= 0:
        reasons.append("no_mapped_product_atoms")
    if not _text(row, "reaction_edit_signature"):
        reasons.append("missing_reaction_edit_signature")
    return sorted(set(reasons))


def _semantic_group_key(row: dict[str, Any]) -> str:
    fields = [
        _text(row, "effective_reaction_delta"),
        _text(row, "reaction_edit_signature"),
        _text(row, "normalized_direction"),
        _text(row, "transferred_group_class"),
        _text(row, "donor_class"),
        _text(row, "acceptor_atom_class"),
    ]
    return stable_hash("\x1f".join(fields), length=24)


def _context_class(reactant_atoms: int) -> str:
    if reactant_atoms <= 8:
        return "local"
    if reactant_atoms <= 20:
        return "mesoscopic"
    return "contextual"


def _selection_score(row: dict[str, Any]) -> float:
    confidence = as_float(row.get("final_rule_confidence"), 0.0)
    source_count = as_int(row.get("consensus_source_database_count"), 0) or 0
    evidence_rows = as_int(row.get("consensus_evidence_rows"), 0) or 0
    template_count = as_int(row.get("template_count"), 0) or 0
    source_records = as_int(row.get("source_record_count"), 0) or 0
    reactant_atoms = as_int(row.get("reactant_atoms"), 0) or 0
    generic_atoms = as_int(row.get("generic_reactant_atoms"), 0) or 0
    consensus_bonus = (
        0.30 if _text(row, "consensus_qc_status").lower().startswith("pass") else 0.0
    )
    annotation_bonus = (
        0.08
        if _text(row, "single_center_evidence_mode") == "source_annotation"
        else 0.0
    )
    evidence_bonus = (
        0.08 * math.log1p(source_count)
        + 0.04 * math.log1p(evidence_rows)
        + 0.03 * math.log1p(template_count)
        + 0.02 * math.log1p(source_records)
    )
    excessive_context_penalty = max(0, reactant_atoms - 20) * 0.005
    generic_penalty = max(0, generic_atoms - 2) * 0.03
    return round(
        confidence
        + consensus_bonus
        + annotation_bonus
        + evidence_bonus
        - excessive_context_penalty
        - generic_penalty,
        6,
    )


def prepare_generative_grammar(
    rules_path: Path,
    output_dir: Path,
    *,
    chunk_size: int = 2000,
    representatives_per_group: int = 3,
    max_reactant_atoms: int = 48,
    require_single_center: bool = True,
    max_rules: int | None = None,
    release_tier: str = "T1",
) -> dict[str, Path]:
    start = time.time()
    output_dir = ensure_dir(output_dir)
    release_tier = str(release_tier).upper()
    if release_tier not in {"T1", "T2", "T3"}:
        raise ValueError("release_tier must be one of T1, T2, or T3")
    selected_candidates: list[dict[str, Any]] = []
    audit_records: list[dict[str, Any]] = []
    exclusion_counts: Counter[str] = Counter()
    rows_seen = 0

    header = pd.read_csv(rules_path, sep="\t", dtype=str, nrows=0)
    missing = sorted(REQUIRED_INPUT_COLUMNS - set(header.columns))
    if missing:
        raise ValueError(f"Rule library is missing required columns: {missing}")
    usecols = [column for column in PASSTHROUGH_COLUMNS if column in header.columns]

    iterator = pd.read_csv(
        rules_path,
        sep="\t",
        dtype=str,
        usecols=usecols,
        chunksize=chunk_size,
        keep_default_na=False,
    )
    for chunk_no, chunk in enumerate(iterator, start=1):
        if max_rules is not None and rows_seen >= max_rules:
            break
        if max_rules is not None and rows_seen + len(chunk) > max_rules:
            chunk = chunk.head(max_rules - rows_seen).copy()
        rows_seen += len(chunk)
        for row in chunk.to_dict("records"):
            metrics = reaction_template_metrics(_text(row, "reaction_smarts"))
            row.update(metrics.__dict__)
            row["effective_reaction_delta"] = (
                _text(row, "reaction_delta_fingerprint")
                or _text(row, "structural_element_delta")
                or f"EDIT:{_text(row, 'reaction_edit_signature')}"
            )
            reasons = _primary_exclusion_reasons(
                row,
                max_reactant_atoms=max_reactant_atoms,
                require_single_center=require_single_center,
                release_tier=release_tier,
            )
            semantic_group_id = _semantic_group_key(row)
            reactant_atoms = as_int(row.get("reactant_atoms"), 0) or 0
            row["semantic_group_id"] = semantic_group_id
            row["reactant_context_class"] = _context_class(reactant_atoms)
            row["grammar_selection_score"] = _selection_score(row)
            row["grammar_qc_status"] = "eligible" if not reasons else "excluded"
            row["grammar_qc_reasons"] = ";".join(reasons)
            audit_records.append(
                {
                    "smarts_rule_id": _text(row, "smarts_rule_id"),
                    "reaction_smarts_hash": _text(row, "reaction_smarts_hash"),
                    "semantic_group_id": semantic_group_id,
                    "reaction_type": _text(row, "reaction_type"),
                    "reaction_delta_fingerprint": _text(
                        row, "reaction_delta_fingerprint"
                    ),
                    "effective_reaction_delta": _text(
                        row, "effective_reaction_delta"
                    ),
                    "reaction_edit_signature": metrics.reaction_edit_signature,
                    "single_center_evidence_mode": _text(
                        row, "single_center_evidence_mode"
                    ),
                    "grammar_qc_status": row["grammar_qc_status"],
                    "grammar_qc_reasons": row["grammar_qc_reasons"],
                    "compile_status": metrics.compile_status,
                    "compile_error": metrics.compile_error,
                    "n_reactants": metrics.n_reactants,
                    "n_products": metrics.n_products,
                    "reactant_atoms": metrics.reactant_atoms,
                    "product_atoms": metrics.product_atoms,
                    "mapped_reactant_atoms": metrics.mapped_reactant_atoms,
                    "mapped_product_atoms": metrics.mapped_product_atoms,
                    "generic_reactant_atoms": metrics.generic_reactant_atoms,
                    "changed_mapped_atom_count": metrics.changed_mapped_atom_count,
                    "inferred_independent_reaction_centers": (
                        metrics.inferred_independent_reaction_centers
                    ),
                    "structural_element_delta": metrics.structural_element_delta,
                    "mapped_atom_retention": metrics.mapped_atom_retention,
                    "mapping_coverage": metrics.mapping_coverage,
                    "reactant_context_class": row["reactant_context_class"],
                    "grammar_selection_score": row["grammar_selection_score"],
                }
            )
            if reasons:
                exclusion_counts.update(reasons)
            else:
                selected_candidates.append(row)
        print(
            f"[prepare-rules] chunk={chunk_no} rows_seen={rows_seen} "
            f"eligible={len(selected_candidates)}",
            flush=True,
        )

    candidates = pd.DataFrame(selected_candidates)
    if candidates.empty:
        raise RuntimeError("No rules passed the primary generative grammar filters")
    candidates = candidates.sort_values(
        [
            "semantic_group_id",
            "grammar_selection_score",
            "reactant_atoms",
            "reaction_smarts_hash",
        ],
        ascending=[True, False, True, True],
        kind="stable",
    )
    candidates["semantic_group_rank"] = (
        candidates.groupby("semantic_group_id").cumcount() + 1
    )
    grammar = candidates[
        candidates["semantic_group_rank"] <= representatives_per_group
    ].copy()
    grammar["grammar_rule_id"] = [
        f"GRAMMAR_{release_tier}_{index:07d}"
        for index in range(1, len(grammar) + 1)
    ]

    group_summary = (
        candidates.groupby(
            [
                "semantic_group_id",
                "effective_reaction_delta",
                "reaction_edit_signature",
            ],
            dropna=False,
        )
        .agg(
            eligible_rule_count=("smarts_rule_id", "size"),
            reaction_types=(
                "reaction_type",
                lambda values: ";".join(
                    sorted(
                        {
                            str(value)
                            for value in values
                            if str(value).strip()
                            and str(value).strip().lower() != "nan"
                        }
                    )
                ),
            ),
            selected_rule_count=(
                "semantic_group_rank",
                lambda values: int((values <= representatives_per_group).sum()),
            ),
            max_selection_score=("grammar_selection_score", "max"),
            min_reactant_atoms=("reactant_atoms", "min"),
            max_reactant_atoms=("reactant_atoms", "max"),
        )
        .reset_index()
    )

    paths = {
        "grammar": output_dir / f"generative_grammar.{release_tier}_primary.tsv",
        "eligible": output_dir / f"generative_grammar.{release_tier}_eligible.tsv",
        "audit": output_dir / f"generative_grammar.{release_tier}_rule_audit.tsv",
        "groups": output_dir / f"generative_grammar.{release_tier}_semantic_groups.tsv",
        "summary": output_dir / f"generative_grammar.{release_tier}_summary.json",
    }
    write_table(grammar, paths["grammar"])
    write_table(candidates, paths["eligible"])
    write_table(pd.DataFrame(audit_records), paths["audit"])
    write_table(group_summary, paths["groups"])

    summary = {
        "mode": f"{release_tier}_single_center_semantic_representative_grammar",
        "release_tier": release_tier,
        "rules_input": str(rules_path),
        "rows_seen": rows_seen,
        "eligible_rules": int(len(candidates)),
        "semantic_groups": int(candidates["semantic_group_id"].nunique()),
        "selected_grammar_rules": int(len(grammar)),
        "representatives_per_group": representatives_per_group,
        "max_reactant_atoms": max_reactant_atoms,
        "require_single_center": require_single_center,
        "context_class_counts": {
            str(key): int(value)
            for key, value in grammar["reactant_context_class"]
            .value_counts()
            .to_dict()
            .items()
        },
        "single_center_evidence_mode_counts": {
            str(key): int(value)
            for key, value in grammar["single_center_evidence_mode"]
            .value_counts()
            .to_dict()
            .items()
        },
        "reaction_type_counts": {
            str(key): int(value)
            for key, value in grammar["reaction_type"].value_counts().to_dict().items()
        },
        "exclusion_reason_counts": dict(sorted(exclusion_counts.items())),
        "elapsed_seconds": round(time.time() - start, 3),
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths


def prepare_taxane_domain_grammar(
    rules_path: Path,
    output_dir: Path,
) -> dict[str, Path]:
    """Harmonize reviewed taxane-domain rules without carrying scaffold fields."""

    start = time.time()
    output_dir = ensure_dir(output_dir)
    raw = pd.read_csv(rules_path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "smarts_rule_id",
        "reaction_smarts",
        "reaction_smarts_hash",
        "reaction_type",
        "taxane_domain_core_use",
    }
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"Domain grammar is missing required columns: {missing}")

    accepted: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    exclusion_counts: Counter[str] = Counter()
    for source_row in raw.to_dict("records"):
        row = {
            column: source_row.get(column, "")
            for column in PASSTHROUGH_COLUMNS
            if column in source_row
        }
        row["exclusive_release_tier"] = "T1_taxane_domain"
        row["domain_rule_id"] = _text(source_row, "taxane_domain_rule_id")
        row["domain_consensus_type"] = _text(
            source_row, "taxane_consensus_type"
        )
        row["domain_release_layer"] = _text(source_row, "domain_release_layer")
        metrics = reaction_template_metrics(_text(row, "reaction_smarts"))
        row.update(metrics.__dict__)
        row["effective_reaction_delta"] = (
            _text(row, "reaction_delta_fingerprint")
            or metrics.structural_element_delta
            or f"EDIT:{metrics.reaction_edit_signature}"
        )
        reasons: list[str] = []
        if not truthy(source_row.get("taxane_domain_core_use", "")):
            reasons.append("not_domain_core")
        if not truthy(source_row.get("predictive_rule_use", "")):
            reasons.append("not_predictive_rule")
        if _text(source_row, "template_qc_status").lower() != "ok":
            reasons.append("template_qc_not_ok")
        granularity = _text(
            source_row, "biochemical_step_granularity"
        ).lower()
        if (
            "single" not in granularity
            or "multi" in granularity
            or "composite" in granularity
        ):
            reasons.append("not_likely_single_step")
        if truthy(source_row.get("composite_rule_flag", "")):
            reasons.append("composite_rule")
        if metrics.compile_status != "ok":
            reasons.append("reaction_smarts_compile_failed")
        if metrics.n_reactants != 1:
            reasons.append("not_single_reactant_smarts")
        if metrics.n_products != 1:
            reasons.append("not_single_product_smarts")
        if metrics.inferred_independent_reaction_centers != 1:
            reasons.append("single_center_not_verified")
        if metrics.mapped_atom_retention < 0.70:
            reasons.append("mapping_retention_below_0.70")
        if not metrics.reaction_edit_signature:
            reasons.append("missing_reaction_edit_signature")

        row["single_center_evidence_mode"] = "domain_source_plus_smarts_recomputed"
        row["semantic_group_id"] = _semantic_group_key(row)
        row["reactant_context_class"] = _context_class(metrics.reactant_atoms)
        row["grammar_selection_score"] = _selection_score(row)
        row["grammar_qc_status"] = "eligible" if not reasons else "excluded"
        row["grammar_qc_reasons"] = ";".join(sorted(set(reasons)))
        audit.append(
            {
                "smarts_rule_id": _text(row, "smarts_rule_id"),
                "domain_rule_id": row["domain_rule_id"],
                "domain_consensus_type": row["domain_consensus_type"],
                "reaction_type": _text(row, "reaction_type"),
                "reaction_smarts_hash": _text(row, "reaction_smarts_hash"),
                "semantic_group_id": row["semantic_group_id"],
                "grammar_qc_status": row["grammar_qc_status"],
                "grammar_qc_reasons": row["grammar_qc_reasons"],
                "compile_status": metrics.compile_status,
                "compile_error": metrics.compile_error,
                "n_reactants": metrics.n_reactants,
                "n_products": metrics.n_products,
                "reactant_atoms": metrics.reactant_atoms,
                "product_atoms": metrics.product_atoms,
                "mapped_atom_retention": metrics.mapped_atom_retention,
                "mapping_coverage": metrics.mapping_coverage,
                "changed_mapped_atom_count": metrics.changed_mapped_atom_count,
                "inferred_independent_reaction_centers": (
                    metrics.inferred_independent_reaction_centers
                ),
                "structural_element_delta": metrics.structural_element_delta,
                "reaction_edit_signature": metrics.reaction_edit_signature,
            }
        )
        if reasons:
            exclusion_counts.update(reasons)
        else:
            accepted.append(row)

    grammar = pd.DataFrame(accepted)
    if not grammar.empty:
        grammar = grammar.sort_values(
            ["semantic_group_id", "grammar_selection_score", "reaction_smarts_hash"],
            ascending=[True, False, True],
            kind="stable",
        ).drop_duplicates("reaction_smarts_hash", keep="first")
        grammar["semantic_group_rank"] = (
            grammar.groupby("semantic_group_id").cumcount() + 1
        )
    audit_frame = pd.DataFrame(audit)
    paths = {
        "grammar": output_dir / "taxane_domain_grammar.primary.tsv",
        "audit": output_dir / "taxane_domain_grammar.audit.tsv",
        "summary": output_dir / "taxane_domain_grammar.summary.json",
    }
    write_table(grammar, paths["grammar"])
    write_table(audit_frame, paths["audit"])
    summary = {
        "mode": "reviewed_taxane_domain_grammar_harmonization",
        "input_rules": int(len(raw)),
        "eligible_rules": int(len(grammar)),
        "eligible_semantic_groups": int(
            grammar.get("semantic_group_id", pd.Series(dtype=str)).nunique()
        ),
        "exclusion_reason_counts": dict(sorted(exclusion_counts.items())),
        "elapsed_seconds": round(time.time() - start, 3),
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    write_json(summary, paths["summary"])
    return paths
