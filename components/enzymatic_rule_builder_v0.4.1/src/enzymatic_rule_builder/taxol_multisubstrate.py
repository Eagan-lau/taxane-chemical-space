from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from .chem import canonical_reaction_smiles, canonical_smiles
from .utils import clean_text, ensure_dir, join_values, sha256_text, write_json


AUGMENT_COLUMNS = [
    "rule_form",
    "network_application_mode",
    "multi_reactant_included",
    "can_apply_single_substrate_network",
    "can_apply_projected_main_pair_network",
    "rdkit_compile_status",
    "n_reactants_rdkit",
    "n_products_rdkit",
    "taxol_rule_variant",
    "taxol_external_participant_inference",
    "required_external_participants_left",
    "required_external_participants_right",
    "external_participant_basis",
    "external_participant_confidence",
    "full_reaction_semantic",
    "main_pair_reaction_smiles",
    "projected_main_pair_reaction_smarts",
    "source_full_multisubstrate_reaction_smarts",
    "projected_from_multisubstrate_rule",
    "projection_qc_status",
    "projection_note",
    "source_taxol_pathway_file",
    "augmentation_source",
]


def _norm_bool(value: Any) -> str:
    return "true" if str(value).strip().lower() in {"1", "true", "yes", "y", "t"} else "false"


def _read_qc_index(rule_qc: str | Path) -> pd.DataFrame:
    cols = ["smarts_rule_id", "status", "n_reactants", "n_products"]
    df = pd.read_csv(rule_qc, sep="\t", dtype=str, usecols=lambda c: c in cols).fillna("")
    if "smarts_rule_id" not in df.columns:
        raise ValueError(f"rule QC file lacks smarts_rule_id: {rule_qc}")
    return df.drop_duplicates("smarts_rule_id", keep="first")


def _reaction_type_for_taxol(enzyme: str, ec: str) -> str:
    e = clean_text(enzyme).lower()
    ec_text = clean_text(ec)
    if ec_text.startswith(("1.14.13.", "1.14.14.", "1.14.20.")) or "h" in e or "ogd" in e or "ox" in e:
        return "hydroxylation_or_oxygenation"
    if ec_text.startswith("2.3.1."):
        if "bt" in e or "nbt" in e or "bapt" in e:
            return "benzoylation_or_aromatic_acyl_transfer"
        return "acetylation_or_deacetylation_like_acyl_transfer"
    if ec_text.startswith("3.1.1."):
        return "deacetylation_or_acetyl_ester_hydrolysis"
    if ec_text.startswith("3.3.2."):
        return "epoxide_hydrolase_or_ring_opening"
    if ec_text.startswith("4.2.3."):
        return "cyclization_or_ring_rearrangement"
    return "unassigned_structural_delta"


def _participant_profile(enzyme: str, ec: str) -> dict[str, Any]:
    """Infer missing Taxol pathway external participants from enzyme/EC annotations.

    The input taxol_pathway.csv contains curated main substrate/product pairs but not
    complete stoichiometric co-substrates. These profiles are intentionally semantic
    annotations; they are not claimed to be atom-balanced source records.
    """
    e = clean_text(enzyme)
    e_low = e.lower()
    ec_text = clean_text(ec)

    if ec_text.startswith("4.2.3."):
        return {
            "left": [],
            "right": ["diphosphate"],
            "basis": "EC 4.2.3 terpene cyclase: prenyl diphosphate elimination coproduct inferred",
            "confidence": "high",
        }

    if ec_text.startswith(("1.14.14.", "1.14.13.")):
        return {
            "left": ["O2", "NADPH", "H+"],
            "right": ["H2O", "NADP+"],
            "basis": f"EC {ec_text} monooxygenase/P450-like hydroxylation cofactor set inferred",
            "confidence": "high",
        }

    if ec_text.startswith("1.14.20.") or "ogd" in e_low:
        return {
            "left": ["O2", "2-oxoglutarate"],
            "right": ["succinate", "CO2"],
            "basis": "2-oxoglutarate-dependent oxygenase co-substrates inferred from EC/enzyme name",
            "confidence": "medium_high" if ec_text.startswith("1.14.20.") else "medium",
        }

    if ec_text.startswith("1.14.") or e_low.endswith("h") or "h-" in e_low or "ox" in e_low:
        return {
            "left": ["O2", "reduced electron donor"],
            "right": ["H2O", "oxidized electron donor"],
            "basis": "oxygenase-like external participants inferred from partial EC/enzyme name",
            "confidence": "medium",
        }

    if ec_text.startswith("2.3.1.") or e_low in {"tat", "dbat", "t7at", "tbt", "bapt", "t3'alpha nbt", "t3'\u03b1nbt"}:
        if "bapt" in e_low:
            donor = "beta-phenylalanoyl-CoA"
            basis = "BAPT acyltransferase side-chain donor inferred from Taxol pathway enzyme identity"
        elif "bt" in e_low or "nbt" in e_low:
            donor = "benzoyl-CoA"
            basis = "benzoyltransferase donor inferred from Taxol pathway enzyme identity"
        else:
            donor = "acetyl-CoA"
            basis = "acetyltransferase donor inferred from Taxol pathway enzyme identity"
        return {
            "left": [donor],
            "right": ["CoA"],
            "basis": basis,
            "confidence": "high",
        }

    if ec_text.startswith("3.1.1.") or "da" in e_low:
        return {
            "left": ["H2O"],
            "right": ["acetate"],
            "basis": "ester hydrolase/deacetylase water and acetate products inferred from EC/enzyme name",
            "confidence": "medium_high",
        }

    if ec_text.startswith("3.3.2."):
        return {
            "left": ["H2O"],
            "right": [],
            "basis": "epoxide hydrolase water participant inferred from EC",
            "confidence": "medium_high",
        }

    return {
        "left": [],
        "right": [],
        "basis": "no external participants inferred from available Taxol pathway columns",
        "confidence": "unknown",
    }


def _semantic_equation(substrate: str, product: str, left_parts: list[str], right_parts: list[str]) -> str:
    left = [substrate] + list(left_parts)
    right = [product] + list(right_parts)
    return " + ".join(left) + " >> " + " + ".join(right)


def _taxol_rows(taxol_pathway: str | Path, base_columns: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = Path(taxol_pathway)
    df = pd.read_csv(path, dtype=str).fillna("")
    required = {"Enzyme", "Substrate", "Product", "EC"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"taxol pathway file lacks required columns {sorted(missing)}: {path}")

    main_rows: list[dict[str, Any]] = []
    ext_rows: list[dict[str, Any]] = []
    mapping_rows: list[dict[str, Any]] = []
    all_cols = list(dict.fromkeys(list(base_columns) + AUGMENT_COLUMNS))

    for idx, row in df.iterrows():
        enzyme = clean_text(row.get("Enzyme", "")) or f"TAXOL_STEP_{idx+1:03d}"
        ec = clean_text(row.get("EC", ""))
        substrate = canonical_smiles(row.get("Substrate", "")) or clean_text(row.get("Substrate", ""))
        product = canonical_smiles(row.get("Product", "")) or clean_text(row.get("Product", ""))
        main_pair = canonical_reaction_smiles(f"{substrate}>>{product}") or f"{substrate}>>{product}"
        profile = _participant_profile(enzyme, ec)
        semantic = _semantic_equation(substrate, product, profile["left"], profile["right"])
        reaction_type = _reaction_type_for_taxol(enzyme, ec)

        common = {
            "reaction_smarts": main_pair,
            "reaction_smarts_hash": sha256_text(main_pair)[:24],
            "smarts_library_tier": "taxol_curated_augmented",
            "template_hash": sha256_text(f"taxol|{enzyme}|{main_pair}")[:24],
            "template_scope": "taxol_known_pathway",
            "predictive_rule_use": "false",
            "template_qc_status": "ok",
            "abstracted_from_exact_reaction": "false",
            "derived_from_exact_anchor": "true",
            "rxnmapper_confidence": "",
            "rdchiral_extraction_status": "not_applicable_exact_curated_pair",
            "abstracted_smarts_applies_to_original_pair": "true",
            "exact_abstraction_qc_status": "taxol_curated_main_pair",
            "benchmark_exclusion_flag": "none",
            "reaction_type": reaction_type,
            "reaction_subtype": enzyme,
            "reaction_type_assignment_mode": "taxol_ec_enzyme_curated",
            "rule_application_unit": "taxol_main_pair_exact_anchor",
            "biochemical_step_granularity": "known_taxol_pathway_anchor",
            "biochemical_step_granularity_confidence": "high",
            "granularity_assignment_mode": "taxol_curated_pathway",
            "granularity_evidence_summary": "Taxol known pathway curated substrate/product pair",
            "composite_rule_flag": "false",
            "smarts_direction": "left_to_right",
            "molecular_direction": "substrate_to_product",
            "direction_evidence_type": "curated_taxol_T1_direction",
            "source_direction": "forward",
            "direction_handling": "source_columns_as_curated_main_pair",
            "reversible_group_id": "",
            "direction_variant": "forward",
            "normalized_direction": "substrate_to_product",
            "direction_qc_status": "direction_qc_ok",
            "direction_qc_note": "Curated Taxol pathway main pair is oriented substrate-to-product.",
            "reverse_transform_available": "false",
            "reverse_rule_ids": "",
            "reverse_rule_relation": "not_reversible_curated_pathway_step",
            "cofactor_or_donor_class": join_values(profile["left"]),
            "external_participant_roles": ";".join(
                [f"{x}:left_required_external_participant" for x in profile["left"]]
                + [f"{x}:right_external_product" for x in profile["right"]]
            ),
            "participant_role_confidence": profile["confidence"],
            "reaction_representation_scope": "source_curated_node_transformation",
            "transferred_group": "",
            "leaving_group_class": "",
            "source_ec_numbers": ec,
            "candidate_ec_numbers": ec,
            "template_ec_candidates": ec,
            "database_ec_candidates": ec,
            "ec_prior_candidates": ec,
            "full_ec_numbers": ec if "-" not in ec else "",
            "partial_ec_numbers": ec if "-" in ec else "",
            "supported_partial_ec_numbers": ec if "-" in ec else "",
            "prior_ec_numbers": ec,
            "ec_status": "taxol_curated_ec",
            "ranked_ec_numbers": ec,
            "top_ec_number": ec,
            "top_ec_confidence": "1.0",
            "top_ec_assignment_mode": "taxol_known_pathway_csv",
            "top_ec_evidence_types": "curated_taxol_pathway",
            "top_ec_sources": "taxol_pathway.csv",
            "top_ec_specificity": "full" if "-" not in ec else "partial",
            "top_ec_granularity": "enzyme_step",
            "top3_ec_numbers": ec,
            "top3_ec_confidences": "1.0",
            "top3_ec_assignment_modes": "taxol_known_pathway_csv",
            "top3_ec_sources": "taxol_pathway.csv",
            "top3_ec_specificities": "full" if "-" not in ec else "partial",
            "top3_ec_granularities": "enzyme_step",
            "ec_annotation_scope": "taxol_curated_step",
            "broad_ec_classes": ec.split(".")[0] if ec else "",
            "broad_ec_class_count": "1" if ec else "0",
            "ec_candidate_count": "1" if ec else "0",
            "ec_conflict_flag": "false",
            "ec_conflict_level": "none",
            "ec_reaction_type_consistency": "curated_pathway_not_recomputed",
            "ec_reaction_type_expected_classes": "",
            "ec_reaction_type_observed_classes": ec.split(".")[0] if ec else "",
            "ec_reaction_type_top_class": ec.split(".")[0] if ec else "",
            "ec_reaction_type_consistency_note": "Taxol curated augmentation carries pathway EC evidence; reaction-type EC consistency is not recomputed in this auxiliary augmentation.",
            "ec_reaction_type_consistency_mode": "curated_taxol_auxiliary_augmentation",
            "ec_directionality_scope": "curated_direction_specific_taxol_step",
            "ec_directionality_warning": "",
            "reverse_ec_inheritance_policy": "do_not_reuse_for_reverse_edge;use_explicit_reverse_directional_rule_if_present",
            "strict_ec_annotation_use": "true",
            "ec_evidence_summary_json": json.dumps({"source": "taxol_pathway.csv", "enzyme": enzyme, "ec": ec}, sort_keys=True),
            "primary_candidate_families": enzyme,
            "secondary_candidate_families": "",
            "family_assignment_mode": "taxol_pathway_enzyme_name",
            "family_evidence": "Taxol pathway enzyme column",
            "family_annotation_available": "true",
            "family_annotation_confidence": "high",
            "family_annotation_scope": "taxol_curated_step",
            "family_evidence_sources": "taxol_pathway.csv",
            "evidence_layer_best": "T1_Bio_Core",
            "evidence_layers_all": "T1_Bio_Core",
            "template_sources": "TaxolKnownPathway_Curated",
            "source_reaction_ids": enzyme,
            "curated_taxol_anchor": "true",
            "curated_pathway_name": "TaxolKnownPathway",
            "curated_pathway_step_ids": f"taxol_pathway_{idx+1}",
            "rhea_ids": "",
            "kegg_ids": "",
            "metanetx_ids": "",
            "final_rule_confidence": "1.0",
            "strict_core_use": "true",
            "expanded_use": "true",
            "exploratory_use": "true",
            "template_count": "1",
            "source_record_count": "1",
            "example_reaction_smiles": main_pair,
            "example_substrate_smiles": substrate,
            "example_product_smiles": product,
            "reaction_delta_fingerprint": "",
            "notes": "Taxol known pathway rule generated from curated main pair; external participants may be represented in augmentation columns.",
            "coverage_release_component": "taxol_known_pathway_augmented",
            "rule_form": "taxol_main_pair_exact",
            "network_application_mode": "single_substrate_main_pair_anchor",
            "multi_reactant_included": "false",
            "can_apply_single_substrate_network": "true",
            "can_apply_projected_main_pair_network": "true",
            "rdkit_compile_status": "not_recompiled",
            "n_reactants_rdkit": "1",
            "n_products_rdkit": "1",
            "required_external_participants_left": "",
            "required_external_participants_right": "",
            "external_participant_basis": "",
            "external_participant_confidence": "",
            "full_reaction_semantic": "",
            "main_pair_reaction_smiles": main_pair,
            "projected_main_pair_reaction_smarts": main_pair,
            "source_taxol_pathway_file": str(path),
            "augmentation_source": "taxol_known_pathway_csv",
        }

        main = {c: "" for c in all_cols}
        main.update(common)
        main["smarts_rule_id"] = f"TAXOL_MAINPAIR_{idx + 1:06d}"
        main["rule_id"] = f"TAXOL_{enzyme}_MAINPAIR"
        main["taxol_rule_variant"] = "main_pair_no_external_participants"
        main_rows.append(main)

        ext = {c: "" for c in all_cols}
        ext.update(common)
        ext["smarts_rule_id"] = f"TAXOL_EXTERNAL_{idx + 1:06d}"
        ext["rule_id"] = f"TAXOL_{enzyme}_WITH_EXTERNAL_PARTICIPANTS"
        ext["template_scope"] = "taxol_known_pathway_with_inferred_external_participants"
        ext["rule_form"] = "taxol_projected_main_pair_with_external_participants"
        ext["network_application_mode"] = "single_substrate_projected_with_required_external_participants"
        ext["taxol_rule_variant"] = "main_pair_with_inferred_external_participants"
        ext["taxol_external_participant_inference"] = "true"
        ext["required_external_participants_left"] = join_values(profile["left"])
        ext["required_external_participants_right"] = join_values(profile["right"])
        ext["external_participant_basis"] = profile["basis"]
        ext["external_participant_confidence"] = profile["confidence"]
        ext["full_reaction_semantic"] = semantic
        ext["reaction_representation_scope"] = "source_curated_node_transformation_with_inferred_external_participants"
        ext["exact_abstraction_qc_status"] = "taxol_curated_main_pair_with_inferred_external_participants"
        ext["notes"] = (
            "Taxol known pathway projected main-pair rule with inferred external participants; "
            "participants are semantic annotations, not source-provided stoichiometric records."
        )
        ext_rows.append(ext)

        mapping_rows.append(
            {
                "enzyme": enzyme,
                "ec": ec,
                "reaction_type": reaction_type,
                "substrate_smiles": substrate,
                "product_smiles": product,
                "main_pair_reaction_smiles": main_pair,
                "required_external_participants_left": join_values(profile["left"]),
                "required_external_participants_right": join_values(profile["right"]),
                "external_participant_basis": profile["basis"],
                "external_participant_confidence": profile["confidence"],
                "full_reaction_semantic": semantic,
                "inference_scope": "semantic_external_participant_annotation",
            }
        )

    return pd.DataFrame(main_rows, columns=all_cols), pd.DataFrame(ext_rows, columns=all_cols), pd.DataFrame(mapping_rows)


def _add_original_aug_cols(chunk: pd.DataFrame) -> pd.DataFrame:
    chunk = chunk.copy()
    for col in AUGMENT_COLUMNS:
        if col not in chunk.columns:
            chunk[col] = ""
    n_react = pd.to_numeric(chunk.get("n_reactants_rdkit", ""), errors="coerce").fillna(0).astype(int)
    is_multi = n_react > 1
    chunk["rule_form"] = is_multi.map(lambda x: "full_multi_reactant_smarts" if x else "single_reactant_smarts")
    chunk["network_application_mode"] = is_multi.map(
        lambda x: "requires_multi_reactant_or_projected_engine" if x else "single_substrate_replay"
    )
    chunk["multi_reactant_included"] = is_multi.map(lambda x: "true" if x else "false")
    chunk["can_apply_single_substrate_network"] = is_multi.map(lambda x: "false" if x else "true")
    projected = chunk.get("example_reaction_smiles", pd.Series([""] * len(chunk), index=chunk.index)).fillna("").astype(str)
    chunk["projected_main_pair_reaction_smarts"] = projected.where(projected.str.contains(">>", regex=False), "")
    chunk["can_apply_projected_main_pair_network"] = is_multi & chunk["projected_main_pair_reaction_smarts"].astype(bool)
    chunk["can_apply_projected_main_pair_network"] = chunk["can_apply_projected_main_pair_network"].map(_norm_bool)
    chunk["augmentation_source"] = "original_release"
    return chunk


def _projected_network_chunk(chunk: pd.DataFrame) -> pd.DataFrame:
    projected = chunk.copy()
    n_react = pd.to_numeric(projected["n_reactants_rdkit"], errors="coerce").fillna(0).astype(int)
    projected_smarts = projected["projected_main_pair_reaction_smarts"].fillna("").astype(str)
    is_multi_projected = (n_react > 1) & projected_smarts.str.contains(">>", regex=False)
    if is_multi_projected.any():
        projected.loc[is_multi_projected, "source_full_multisubstrate_reaction_smarts"] = projected.loc[
            is_multi_projected, "reaction_smarts"
        ]
        projected.loc[is_multi_projected, "reaction_smarts"] = projected.loc[
            is_multi_projected, "projected_main_pair_reaction_smarts"
        ]
        projected.loc[is_multi_projected, "reaction_smarts_hash"] = projected.loc[
            is_multi_projected, "reaction_smarts"
        ].map(lambda x: sha256_text(x)[:24])
        projected.loc[is_multi_projected, "smarts_rule_id"] = projected.loc[
            is_multi_projected, "smarts_rule_id"
        ].astype(str) + "_PROJECTED"
        projected.loc[is_multi_projected, "rule_id"] = projected.loc[is_multi_projected, "rule_id"].astype(
            str
        ) + "_PROJECTED"
        projected.loc[is_multi_projected, "rule_form"] = "projected_main_pair_from_multi_reactant_rule"
        projected.loc[
            is_multi_projected, "network_application_mode"
        ] = "single_substrate_projected_from_multi_reactant_rule"
        projected.loc[is_multi_projected, "can_apply_single_substrate_network"] = "true"
        projected.loc[is_multi_projected, "projected_from_multisubstrate_rule"] = "true"
        projected.loc[
            is_multi_projected, "projection_qc_status"
        ] = "example_main_pair_projection_not_generalized"
        projected.loc[
            is_multi_projected, "projection_note"
        ] = "Projected reaction_smarts was replaced by the source example substrate>>product pair; use as exact/projected edge evidence, not a generalized multi-reactant SMARTS."
    projected.loc[~is_multi_projected, "projected_from_multisubstrate_rule"] = projected.loc[
        ~is_multi_projected, "projected_from_multisubstrate_rule"
    ].replace("", "false")
    return projected


def augment_multisubstrate_taxol_release(
    *,
    input_library: str | Path,
    rule_qc: str | Path,
    taxol_pathway: str | Path,
    output_dir: str | Path,
    chunk_size: int = 50000,
) -> dict[str, Any]:
    output = ensure_dir(output_dir)
    input_library = Path(input_library)
    rule_qc = Path(rule_qc)

    qc = _read_qc_index(rule_qc).rename(
        columns={
            "status": "rdkit_compile_status",
            "n_reactants": "n_reactants_rdkit",
            "n_products": "n_products_rdkit",
        }
    )

    header = pd.read_csv(input_library, sep="\t", nrows=0)
    base_columns = list(header.columns)
    all_columns = list(dict.fromkeys(base_columns + AUGMENT_COLUMNS))
    main_taxol, external_taxol, participant_map = _taxol_rows(taxol_pathway, all_columns)

    full_out = output / "reaction_smarts_library.coverage_max_multisubstrate_taxol_augmented.tsv"
    multi_out = output / "reaction_smarts_library.full_multisubstrate_rules.tsv"
    projected_out = output / "reaction_smarts_library.network_projected_with_taxol_external.tsv"
    main_taxol_out = output / "taxol_known_pathway.rules.main_pair.tsv"
    external_taxol_out = output / "taxol_known_pathway.rules.with_external_participants.tsv"
    participant_out = output / "taxol_external_participant_mapping.tsv"
    summary_out = output / "taxol_multisubstrate_augmentation_summary.json"

    for path in [full_out, multi_out, projected_out]:
        if path.exists():
            path.unlink()

    total_rows = 0
    multi_rows = 0
    projected_rows = 0
    first_full = True
    first_multi = True
    first_projected = True

    qc_cols = ["smarts_rule_id", "rdkit_compile_status", "n_reactants_rdkit", "n_products_rdkit"]
    for chunk in pd.read_csv(input_library, sep="\t", dtype=str, chunksize=chunk_size):
        chunk = chunk.fillna("")
        # Prefer QC status from the actual network compile pass.
        chunk = chunk.drop(columns=[c for c in qc_cols[1:] if c in chunk.columns], errors="ignore")
        chunk = chunk.merge(qc[qc_cols], on="smarts_rule_id", how="left")
        chunk[["rdkit_compile_status", "n_reactants_rdkit", "n_products_rdkit"]] = chunk[
            ["rdkit_compile_status", "n_reactants_rdkit", "n_products_rdkit"]
        ].fillna("")
        chunk = _add_original_aug_cols(chunk)
        chunk = chunk.reindex(columns=all_columns, fill_value="")

        n_react = pd.to_numeric(chunk["n_reactants_rdkit"], errors="coerce").fillna(0).astype(int)
        is_multi = n_react > 1
        is_projectable = chunk["can_apply_single_substrate_network"].eq("true") | chunk[
            "can_apply_projected_main_pair_network"
        ].eq("true")

        chunk.to_csv(full_out, sep="\t", index=False, mode="w" if first_full else "a", header=first_full)
        first_full = False

        if is_multi.any():
            chunk.loc[is_multi].to_csv(
                multi_out, sep="\t", index=False, mode="w" if first_multi else "a", header=first_multi
            )
            first_multi = False
            multi_rows += int(is_multi.sum())

        if is_projectable.any():
            projected_chunk = _projected_network_chunk(chunk.loc[is_projectable])
            projected_chunk.to_csv(
                projected_out,
                sep="\t",
                index=False,
                mode="w" if first_projected else "a",
                header=first_projected,
            )
            first_projected = False
            projected_rows += int(is_projectable.sum())

        total_rows += int(len(chunk))

    # Append Taxol main-pair and participant-augmented variants to full and projected releases.
    taxol_combined = pd.concat([main_taxol, external_taxol], ignore_index=True).reindex(columns=all_columns, fill_value="")
    taxol_combined.to_csv(full_out, sep="\t", index=False, mode="a", header=False)
    taxol_combined.to_csv(projected_out, sep="\t", index=False, mode="a", header=False)

    main_taxol.to_csv(main_taxol_out, sep="\t", index=False)
    external_taxol.to_csv(external_taxol_out, sep="\t", index=False)
    participant_map.to_csv(participant_out, sep="\t", index=False)

    summary = {
        "mode": "taxol_multisubstrate_release_augmentation",
        "input_library": str(input_library),
        "rule_qc": str(rule_qc),
        "taxol_pathway": str(taxol_pathway),
        "original_rows": total_rows,
        "original_full_multisubstrate_rows": multi_rows,
        "original_projectable_rows": projected_rows,
        "original_multisubstrate_projected_rows": int(multi_rows),
        "taxol_main_pair_rows_added": int(len(main_taxol)),
        "taxol_external_participant_rows_added": int(len(external_taxol)),
        "final_augmented_library_rows": int(total_rows + len(taxol_combined)),
        "final_projected_library_rows": int(projected_rows + len(taxol_combined)),
        "outputs": {
            "augmented_library": str(full_out),
            "full_multisubstrate_rules": str(multi_out),
            "network_projected_with_taxol_external": str(projected_out),
            "taxol_main_pair_rules": str(main_taxol_out),
            "taxol_external_participant_rules": str(external_taxol_out),
            "taxol_external_participant_mapping": str(participant_out),
            "summary": str(summary_out),
        },
        "notes": [
            "Original multi-reactant SMARTS are retained as full_multi_reactant_smarts and require a multi-reactant/projected network engine.",
            "Taxol external participants are inferred from enzyme/EC annotations because taxol_pathway.csv does not provide stoichiometric co-substrates.",
            "Taxol main-pair rules are also preserved without inferred external participants.",
        ],
    }
    write_json(summary_out, summary)
    return summary
