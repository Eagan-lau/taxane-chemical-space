from __future__ import annotations

from collections import defaultdict

import pandas as pd

from .chem import reverse_reaction_smarts
from .ec import summarize_ec_evidence
from .family import annotate_families
from .granularity import annotate_biochemical_step_granularity
from .reaction_type import classify_reaction_type, direction_mode
from .schema import RULE_COLUMNS
from .scoring import DEFAULT_WEIGHTS, combine_scores, score_lookup
from .source_layers import best_layer, is_t3_only, normalize_evidence_layers
from .utils import clean_text, join_values, split_multi_value, truthy


def _aggregate_templates(templates: pd.DataFrame) -> pd.DataFrame:
    if templates.empty:
        return pd.DataFrame()
    df = templates.fillna("").copy()
    if "template_hash" not in df.columns:
        return pd.DataFrame()
    df = df[df["template_hash"].astype(str).str.strip().ne("")]
    if df.empty:
        return pd.DataFrame()
    if df["template_hash"].is_unique:
        if "template_count" not in df.columns:
            df["template_count"] = 1
        if "source_record_count" not in df.columns:
            df["source_record_count"] = 1
        return df

    rows = []
    merge_cols = [
        "source_database", "evidence_layer", "source_reaction_id", "record_id", "ec_numbers", "template_ec_candidates",
        "database_ec_candidates", "ec_prior_candidates", "rhea_ids", "kegg_ids", "metanetx_ids", "reaction_type_source",
        "reaction_subtype_source", "cofactor_or_donor_class", "external_participant_roles", "participant_role_confidence", "reaction_representation_scope", "transferred_group", "leaving_group_class", "enzyme_name", "protein_ids", "source_evidence_text", "source_file",
        "donor_class", "acceptor_atom_class", "transferred_group_class", "main_pair_projection_method", "main_pair_projection_note",
        "curated_taxol_anchor", "curated_pathway_name", "curated_pathway_step_id",
        "abstracted_from_exact_reaction", "derived_from_exact_anchor", "rxnmapper_confidence",
        "rdchiral_extraction_status", "abstracted_smarts_applies_to_original_pair",
        "exact_abstraction_qc_status", "benchmark_exclusion_flag", "direction_handling", "direction_variant",
        "normalized_direction", "direction_qc_status", "direction_qc_note", "reversible_group_id",
        "source_reaction_smarts", "reverse_template_hash",
    ]
    for h, group in df.groupby("template_hash", dropna=False):
        if not h:
            continue
        g2 = group.copy()
        g2["_rep_rank"] = g2["template_scope"].astype(str).map(lambda x: 0 if x == "generalized_template" else 1)
        rep = g2.sort_values("_rep_rank").iloc[0].to_dict()
        for col in merge_cols:
            if col in group.columns:
                rep[col] = join_values(group[col].astype(str).tolist())
        rep["template_count"] = int(len(group))
        rep["source_record_count"] = int(len(set(group.get("record_id", pd.Series(dtype=str)).astype(str))))
        rep["predictive_rule_use"] = "true" if any(group["predictive_rule_use"].astype(str).str.lower().isin(["true", "1", "yes"])) else "false"
        rep["anchor_edge_use"] = "true" if any(group["anchor_edge_use"].astype(str).str.lower().isin(["true", "1", "yes"])) else "false"
        rows.append(rep)
    return pd.DataFrame(rows).fillna("")


def _direction_evidence_type(layer_best: str, dir_mode: str, handling: str) -> str:
    h = clean_text(handling).lower()
    if h.startswith("reversed_from_source"):
        return "source_reverse_corrected"
    if h.startswith("split_reversible"):
        return "source_reversible_split"
    if dir_mode != "unknown":
        return "source_defined"
    return "template_left_to_right"


def _direction_qc_status_for_rule(d: dict, dir_mode: str) -> str:
    status = clean_text(d.get("direction_qc_status", ""))
    if status:
        return status
    handling = clean_text(d.get("direction_handling", "")).lower()
    if "unknown_direction" in handling:
        return "direction_qc_unknown_exploratory_only"
    if handling.startswith("reversed_from_source") or "source_reverse_corrected" in handling:
        return "direction_qc_reversed_to_substrate_product"
    if handling.startswith("split_reversible") or "source_reversible" in handling:
        return "direction_qc_reversible_split"
    if (
        handling.startswith("kept_forward")
        or "source_columns_as_curated_main_pair" in handling
        or "exact_pair_left_to_right" in handling
    ):
        return "direction_qc_ok"
    if dir_mode == "source_forward":
        return "direction_qc_ok"
    if dir_mode == "source_reverse_corrected":
        return "direction_qc_reversed_to_substrate_product"
    if dir_mode == "source_reversible" and handling.startswith("split_reversible"):
        return "direction_qc_reversible_split"
    return "direction_qc_unknown_exploratory_only"


def _direction_strict_ready(d: dict, dir_mode: str | None = None) -> bool:
    mode = dir_mode or direction_mode(d)
    status = _direction_qc_status_for_rule(d, mode)
    if status in {"direction_qc_ok", "direction_qc_reversed_to_substrate_product", "direction_qc_reversible_split"}:
        return True
    handling = clean_text(d.get("direction_handling", "")).lower()
    return handling.startswith("kept_forward") or handling.startswith("reversed_from_source") or handling.startswith("split_reversible")


def _append_direction_note(existing: str, handling: str) -> str:
    note = ""
    h = clean_text(handling).lower()
    if h.startswith("reversed_from_source"):
        note = "source_reverse_SMARTS_corrected_to_substrate_to_product"
    elif h.startswith("split_reversible"):
        note = "source_reversible_SMARTS_split_into_directional_rules"
    elif h.startswith("unknown_direction"):
        note = "unknown_source_direction_kept_left_to_right"
    return join_values([existing, note])


def _exact_abstraction_strict_ready(d: dict) -> bool:
    if not truthy(d.get("abstracted_from_exact_reaction", "")):
        return True
    qc_values = [x.lower() for x in split_multi_value(d.get("exact_abstraction_qc_status", ""))]
    replay_values = split_multi_value(d.get("abstracted_smarts_applies_to_original_pair", ""))
    if not qc_values or not replay_values:
        return False
    return all(x in {"pass", "ok"} for x in qc_values) and all(truthy(x) for x in replay_values)


def _exact_abstraction_notes(d: dict) -> list[str]:
    if not truthy(d.get("abstracted_from_exact_reaction", "")):
        return []
    notes: list[str] = []
    qc_status = clean_text(d.get("exact_abstraction_qc_status", ""))
    replay_status = clean_text(d.get("abstracted_smarts_applies_to_original_pair", ""))
    benchmark_exclusion = clean_text(d.get("benchmark_exclusion_flag", ""))
    if qc_status and not _exact_abstraction_strict_ready(d):
        notes.append(f"exact_anchor_generalization_qc={qc_status}")
    if replay_status and not all(truthy(x) for x in split_multi_value(replay_status)):
        notes.append(f"abstracted_smarts_applies_to_original_pair={replay_status}")
    if benchmark_exclusion and benchmark_exclusion != "none":
        notes.append(f"benchmark_exclusion_flag={benchmark_exclusion}")
    return notes


def build_rule_library(
    templates: pd.DataFrame,
    family_evidence: pd.DataFrame | None = None,
    *,
    evidence_weights: dict | None = None,
    strict_core_min_confidence: float = 0.80,
    expanded_min_confidence: float = 0.60,
    exploratory_min_confidence: float = 0.20,
) -> pd.DataFrame:
    weights = evidence_weights or DEFAULT_WEIGHTS
    agg = _aggregate_templates(templates)
    rows = []
    for _, row in agg.iterrows():
        d = row.to_dict()
        layers_all = normalize_evidence_layers(d.get("evidence_layer", ""))
        layer_best = best_layer(layers_all)
        reaction_type, reaction_subtype, rt_mode = classify_reaction_type(d)
        granularity = annotate_biochemical_step_granularity({**d, "reaction_type": reaction_type, "reaction_subtype": reaction_subtype})
        ec_summary = summarize_ec_evidence({**d, "reaction_type": reaction_type, "reaction_subtype": reaction_subtype}, layer_best)
        candidate_ecs = ec_summary.get("candidate_ec_numbers", "")
        broad = ec_summary.get("broad_ec_classes", "")
        primary_families, secondary_families, fam_mode, fam_evidence = annotate_families({**d, **ec_summary}, family_evidence)
        strict_ec_supported = truthy(ec_summary.get("strict_ec_annotation_use", ""))
        if (
            fam_mode == "none"
            and strict_ec_supported
            and layer_best == "T1_Bio_Core"
            and str(ec_summary.get("ec_conflict_level", "none")) != "high"
        ):
            fam_mode = "ec_supported_family_unassigned"
            fam_evidence = join_values([
                fam_evidence,
                "No specific enzyme-family evidence supplied; retained because high-confidence source/template-supported EC evidence supports enzymatic assignment. EC may be partial or full; granularity is reported separately.",
            ])
        dir_mode = direction_mode(d)
        direction_handling = clean_text(d.get("direction_handling", ""))
        scope = clean_text(d.get("template_scope", ""))
        qc = clean_text(d.get("template_qc_status", "")) or "invalid"
        exact_abstraction_ready = _exact_abstraction_strict_ready(d)
        core_excluded = "exclude_from_core" in clean_text(d.get("benchmark_exclusion_flag", "")).lower()
        predictive = truthy(d.get("predictive_rule_use", ""))
        anchor = truthy(d.get("anchor_edge_use", ""))
        source_s = score_lookup("source_score", layer_best, weights, 0.0)
        qc_s = score_lookup("qc_score", qc, weights, 0.0)
        fam_s = score_lookup("family_score", fam_mode, weights, 0.0)
        rt_s = score_lookup("reaction_type_score", rt_mode, weights, 0.0)
        direction_qc_status = _direction_qc_status_for_rule(d, dir_mode)
        dir_s = score_lookup("direction_score", direction_qc_status, weights, score_lookup("direction_score", dir_mode, weights, 0.10))
        scope_s = score_lookup("scope_score", scope, weights, 0.0)
        final = combine_scores({"source": source_s, "qc": qc_s, "family": 0.0, "reaction_type": rt_s, "direction": dir_s, "scope": scope_s}, weights)
        t3_only = is_t3_only(layers_all)
        # Family evidence is auxiliary and must not decide whether a SMARTS rule is valid.
        # Core inclusion is based on high-quality source/QC/reaction annotation evidence.
        direction_strict_ready = _direction_strict_ready(d, dir_mode)
        strict = bool(
            predictive and scope == "generalized_template" and layer_best == "T1_Bio_Core"
            and qc == "ok" and rt_mode in {"source", "structural_delta_specific", "structural_delta", "ec_broad"}
            and exact_abstraction_ready and not core_excluded and direction_strict_ready
            and str(ec_summary.get("ec_conflict_level", "none")) != "high"
            and final >= strict_core_min_confidence and not t3_only
        )
        expanded = bool(
            predictive and scope == "generalized_template" and layer_best in {"T1_Bio_Core", "T2_Bio_Extended"}
            and qc in {"ok", "rdkit_unavailable"} and exact_abstraction_ready and not core_excluded
            and direction_strict_ready and final >= expanded_min_confidence and not t3_only
        )
        exploratory = bool(predictive and final >= exploratory_min_confidence)
        exclusive_release_tier = "T1_only" if strict else ("T2_only" if expanded else ("T3_only" if exploratory else ""))
        curated_taxol_anchor = any(truthy(x) for x in split_multi_value(d.get("curated_taxol_anchor", "")))
        leakage = "curated_taxol_exact_anchor" if anchor and curated_taxol_anchor and scope == "exact_anchor" else "none"
        benchmark_exclusion = clean_text(d.get("benchmark_exclusion_flag", ""))
        if not benchmark_exclusion and curated_taxol_anchor and scope == "generalized_template":
            benchmark_exclusion = "curated_taxol_anchor_derived_exclude_from_external_recall"
        if not benchmark_exclusion:
            benchmark_exclusion = "none"
        direction_evidence_type = _direction_evidence_type(layer_best, dir_mode, direction_handling)
        base_note = "exact_anchor_not_generalized" if scope == "exact_anchor" else ("curated_taxol_anchor_derived_SMARTS_excluded_from_external_recall_benchmark" if benchmark_exclusion == "curated_taxol_anchor_derived_exclude_from_external_recall" else "")
        notes = join_values([_append_direction_note(base_note, direction_handling), *_exact_abstraction_notes(d)])
        rows.append({
            "rule_id": f"RULE_{len(rows)+1:09d}",
            "template_hash": d.get("template_hash", ""),
            "template_scope": scope,
            "predictive_rule_use": predictive,
            "anchor_edge_use": anchor,
            "template_qc_status": qc,
            "template_qc_note": d.get("template_qc_note", ""),
            "abstracted_from_exact_reaction": d.get("abstracted_from_exact_reaction", "false"),
            "derived_from_exact_anchor": d.get("derived_from_exact_anchor", "false"),
            "rxnmapper_confidence": d.get("rxnmapper_confidence", ""),
            "rdchiral_extraction_status": d.get("rdchiral_extraction_status", ""),
            "abstracted_smarts_applies_to_original_pair": d.get("abstracted_smarts_applies_to_original_pair", ""),
            "exact_abstraction_qc_status": d.get("exact_abstraction_qc_status", ""),
            "benchmark_exclusion_flag": benchmark_exclusion,
            "reaction_smarts": d.get("reaction_smarts", ""),
            "example_reaction_smiles": d.get("canonical_reaction_smiles", "") or d.get("reaction_smiles", ""),
            "example_substrate_smiles": d.get("main_substrate_smiles", ""),
            "example_product_smiles": d.get("main_product_smiles", ""),
            "reaction_delta_fingerprint": d.get("reaction_delta_fingerprint", ""),
            "reaction_type": reaction_type,
            "reaction_subtype": reaction_subtype,
            "reaction_type_assignment_mode": rt_mode,
            "rule_application_unit": granularity.get("rule_application_unit", "single_smarts_application"),
            "biochemical_step_granularity": granularity.get("biochemical_step_granularity", "uncertain"),
            "biochemical_step_granularity_confidence": granularity.get("biochemical_step_granularity_confidence", ""),
            "granularity_assignment_mode": granularity.get("granularity_assignment_mode", ""),
            "granularity_evidence_summary": granularity.get("granularity_evidence_summary", ""),
            "composite_rule_flag": granularity.get("composite_rule_flag", ""),
            "reaction_center_count": granularity.get("reaction_center_count", ""),
            "independent_reaction_center_count": granularity.get("independent_reaction_center_count", ""),
            "functional_group_change_count": granularity.get("functional_group_change_count", ""),
            "main_functional_group_changes": granularity.get("main_functional_group_changes", ""),
            "smarts_direction": "left_to_right",
            "molecular_direction": "substrate_to_product",
            "direction_evidence_type": direction_evidence_type,
            "source_direction": d.get("direction", ""),
            "direction_handling": direction_handling,
            "reversible_group_id": d.get("reversible_group_id", ""),
            "direction_variant": d.get("direction_variant", ""),
            "normalized_direction": d.get("normalized_direction", "substrate_to_product"),
            "direction_qc_status": direction_qc_status,
            "direction_qc_note": d.get("direction_qc_note", ""),
            "reverse_template_hash": d.get("reverse_template_hash", ""),
            "reverse_transform_available": "false",
            "reverse_rule_ids": "",
            "reverse_rule_relation": "not_evaluated_in_rule_builder",
            "cofactor_or_donor_class": d.get("cofactor_or_donor_class", ""),
            "external_participant_roles": d.get("external_participant_roles", ""),
            "participant_role_confidence": d.get("participant_role_confidence", ""),
            "reaction_representation_scope": d.get("reaction_representation_scope", ""),
            "transferred_group": d.get("transferred_group", ""),
            "leaving_group_class": d.get("leaving_group_class", ""),
            "donor_class": d.get("donor_class", ""),
            "acceptor_atom_class": d.get("acceptor_atom_class", ""),
            "transferred_group_class": d.get("transferred_group_class", ""),
            "main_pair_projection_method": d.get("main_pair_projection_method", ""),
            "main_pair_projection_note": d.get("main_pair_projection_note", ""),
            "source_ec_numbers": ec_summary.get("source_ec_numbers", ""),
            "candidate_ec_numbers": candidate_ecs,
            "template_ec_candidates": d.get("template_ec_candidates", ""),
            "database_ec_candidates": d.get("database_ec_candidates", ""),
            "ec_prior_candidates": ec_summary.get("ec_prior_candidates", d.get("ec_prior_candidates", "")),
            "full_ec_numbers": ec_summary.get("full_ec_numbers", ""),
            "partial_ec_numbers": ec_summary.get("partial_ec_numbers", ""),
            "supported_partial_ec_numbers": ec_summary.get("supported_partial_ec_numbers", ""),
            "prior_ec_numbers": ec_summary.get("prior_ec_numbers", ""),
            "ec_status": ec_summary.get("ec_status", "missing"),
            "ranked_ec_numbers": ec_summary.get("ranked_ec_numbers", ""),
            "top_ec_number": ec_summary.get("top_ec_number", ""),
            "top_ec_confidence": ec_summary.get("top_ec_confidence", ""),
            "top_ec_assignment_mode": ec_summary.get("top_ec_assignment_mode", ""),
            "top_ec_evidence_types": ec_summary.get("top_ec_evidence_types", ""),
            "top_ec_sources": ec_summary.get("top_ec_sources", ""),
            "top_ec_specificity": ec_summary.get("top_ec_specificity", ""),
            "top_ec_granularity": ec_summary.get("top_ec_granularity", ""),
            "top3_ec_numbers": ec_summary.get("top3_ec_numbers", ""),
            "top3_ec_confidences": ec_summary.get("top3_ec_confidences", ""),
            "top3_ec_assignment_modes": ec_summary.get("top3_ec_assignment_modes", ""),
            "top3_ec_sources": ec_summary.get("top3_ec_sources", ""),
            "top3_ec_specificities": ec_summary.get("top3_ec_specificities", ""),
            "top3_ec_granularities": ec_summary.get("top3_ec_granularities", ""),
            "ec_annotation_scope": ec_summary.get("ec_annotation_scope", ""),
            "broad_ec_classes": broad,
            "broad_ec_class_count": ec_summary.get("broad_ec_class_count", ""),
            "ec_candidate_count": ec_summary.get("ec_candidate_count", ""),
            "ec_conflict_flag": ec_summary.get("ec_conflict_flag", ""),
            "ec_conflict_level": ec_summary.get("ec_conflict_level", ""),
            "ec_reaction_type_consistency": ec_summary.get("ec_reaction_type_consistency", ""),
            "ec_reaction_type_expected_classes": ec_summary.get("ec_reaction_type_expected_classes", ""),
            "ec_reaction_type_observed_classes": ec_summary.get("ec_reaction_type_observed_classes", ""),
            "ec_reaction_type_top_class": ec_summary.get("ec_reaction_type_top_class", ""),
            "ec_reaction_type_consistency_note": ec_summary.get("ec_reaction_type_consistency_note", ""),
            "ec_reaction_type_consistency_mode": ec_summary.get("ec_reaction_type_consistency_mode", ""),
            "ec_directionality_scope": ec_summary.get("ec_directionality_scope", ""),
            "ec_directionality_warning": ec_summary.get("ec_directionality_warning", ""),
            "reverse_ec_inheritance_policy": ec_summary.get("reverse_ec_inheritance_policy", ""),
            "strict_ec_annotation_use": ec_summary.get("strict_ec_annotation_use", ""),
            "ec_evidence_summary_json": ec_summary.get("ec_evidence_summary_json", ""),
            "primary_candidate_families": primary_families,
            "secondary_candidate_families": secondary_families,
            "family_assignment_mode": fam_mode,
            "family_evidence": fam_evidence,
            "family_annotation_available": "true" if fam_mode not in {"none", "ec_supported_family_unassigned"} else "false",
            "family_annotation_confidence": "1.0" if fam_mode == "external_evidence" else ("0.0" if fam_mode == "none" else ""),
            "family_annotation_scope": "external_auxiliary_annotation" if fam_mode == "external_evidence" else ("not_assigned_at_rule_stage" if fam_mode in {"none", "ec_supported_family_unassigned"} else fam_mode),
            "family_evidence_sources": fam_evidence,
            "evidence_layer_best": layer_best,
            "evidence_layers_all": layers_all,
            "template_sources": d.get("source_database", ""),
            "source_reaction_ids": d.get("source_reaction_id", ""),
            "source_record_ids": d.get("record_id", ""),
            "curated_taxol_anchor": "true" if curated_taxol_anchor else "false",
            "curated_pathway_name": d.get("curated_pathway_name", ""),
            "curated_pathway_step_ids": d.get("curated_pathway_step_id", ""),
            "rhea_ids": d.get("rhea_ids", ""),
            "kegg_ids": d.get("kegg_ids", ""),
            "metanetx_ids": d.get("metanetx_ids", ""),
            "template_count": int(d.get("template_count", 1) or 1),
            "source_record_count": int(d.get("source_record_count", 1) or 1),
            "consensus_group_id": d.get("consensus_group_id", ""),
            "consensus_generation_mode": d.get("consensus_generation_mode", ""),
            "consensus_evidence_rows": d.get("consensus_evidence_rows", ""),
            "consensus_source_database_count": d.get("consensus_source_database_count", ""),
            "consensus_evidence_layer_support": d.get("consensus_evidence_layer_support", ""),
            "consensus_qc_status": d.get("consensus_qc_status", ""),
            "consensus_representative_rule_ids": d.get("consensus_representative_rule_ids", ""),
            "consensus_supporting_reaction_types": d.get("consensus_supporting_reaction_types", ""),
            "source_score": source_s,
            "qc_score": qc_s,
            "family_score": fam_s,
            "reaction_type_score": rt_s,
            "direction_score": dir_s,
            "scope_score": scope_s,
            "final_rule_confidence": final,
            "strict_core_use": strict,
            "expanded_use": expanded,
            "exploratory_use": exploratory,
            "exclusive_release_tier": exclusive_release_tier,
            "leakage_risk": leakage,
            "notes": notes,
        })
    df = pd.DataFrame(rows).fillna("")
    for col in RULE_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    df = df[RULE_COLUMNS]

    # Link paired directional rules after rule IDs are assigned. This does not
    # generate reverse chemistry; it only records that the opposite SMARTS is also
    # present as an explicit rule in the library.
    if not df.empty and "reaction_smarts" in df.columns:
        smarts_to_ids: dict[str, list[str]] = defaultdict(list)
        for _, row in df.iterrows():
            smarts = clean_text(row.get("reaction_smarts", ""))
            if smarts and row.get("template_scope", "") == "generalized_template":
                smarts_to_ids[smarts].append(str(row.get("rule_id", "")))
        for idx, row in df.iterrows():
            smarts = clean_text(row.get("reaction_smarts", ""))
            if not smarts or row.get("template_scope", "") != "generalized_template":
                continue
            rev = reverse_reaction_smarts(smarts)
            ids = [rid for rid in smarts_to_ids.get(rev, []) if rid and rid != row.get("rule_id", "")]
            if ids:
                df.at[idx, "reverse_transform_available"] = "true"
                df.at[idx, "reverse_rule_ids"] = join_values(ids)
                df.at[idx, "reverse_rule_relation"] = "explicit_directional_rule_in_library"
    return df[RULE_COLUMNS]
