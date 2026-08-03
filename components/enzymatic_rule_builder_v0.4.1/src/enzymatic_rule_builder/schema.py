from __future__ import annotations

NORMALIZED_COLUMNS = [
    "record_id", "source_database", "evidence_layer", "source_reaction_id", "source_file", "parser_name",
    "curated_taxol_anchor", "curated_pathway_name", "curated_pathway_step_id",
    "reaction_smiles", "substrate_smiles", "product_smiles", "reaction_smarts", "reaction_equation",
    "direction", "is_reversible", "ec_numbers", "template_ec_candidates", "database_ec_candidates", "ec_prior_candidates",
    "rhea_ids", "kegg_ids", "metanetx_ids", "reaction_type_source", "reaction_subtype_source",
    "cofactor_or_donor_class", "enzyme_name", "protein_ids", "protein_sequence", "source_evidence_text",
    "raw_row_index", "source_row_hash",
]

MAIN_PAIR_COLUMNS = NORMALIZED_COLUMNS + [
    "main_substrate_smiles", "main_product_smiles", "canonical_reaction_smiles",
    "canonical_substrate_smiles", "canonical_product_smiles", "removed_participants",
    "external_participant_roles", "participant_role_confidence", "reaction_representation_scope",
    "transferred_group", "leaving_group_class",
    "donor_class", "acceptor_atom_class", "transferred_group_class",
    "main_pair_projection_method", "main_pair_projection_note",
    "main_pair_method", "main_pair_confidence", "direction_handling", "reversible_group_id", "direction_variant",
    "normalized_direction", "direction_qc_status", "direction_qc_note",
    "reaction_delta_fingerprint", "reaction_delta_json",
    "abstracted_from_exact_reaction", "derived_from_exact_anchor", "rxnmapper_confidence",
    "rdchiral_extraction_status", "abstracted_smarts_applies_to_original_pair",
    "exact_abstraction_qc_status", "benchmark_exclusion_flag",
]

TEMPLATE_COLUMNS = MAIN_PAIR_COLUMNS + [
    "template_id", "template_hash", "template_origin", "template_scope", "template_generalization",
    "source_reaction_smarts", "reverse_template_hash",
    "predictive_rule_use", "anchor_edge_use", "template_extraction_status", "template_qc_status", "template_qc_note",
    "consensus_group_id", "consensus_generation_mode", "consensus_evidence_rows",
    "consensus_source_database_count", "consensus_evidence_layer_support", "consensus_qc_status",
    "consensus_representative_rule_ids", "consensus_supporting_reaction_types",
]

RULE_COLUMNS = [
    "rule_id", "template_hash", "template_scope", "predictive_rule_use", "anchor_edge_use",
    "template_qc_status", "template_qc_note",
    "abstracted_from_exact_reaction", "derived_from_exact_anchor", "rxnmapper_confidence",
    "rdchiral_extraction_status", "abstracted_smarts_applies_to_original_pair",
    "exact_abstraction_qc_status", "benchmark_exclusion_flag",
    "reaction_smarts", "example_reaction_smiles", "example_substrate_smiles", "example_product_smiles",
    "reaction_delta_fingerprint", "reaction_type", "reaction_subtype", "reaction_type_assignment_mode",
    "rule_application_unit", "biochemical_step_granularity", "biochemical_step_granularity_confidence",
    "granularity_assignment_mode", "granularity_evidence_summary", "composite_rule_flag",
    "reaction_center_count", "independent_reaction_center_count", "functional_group_change_count",
    "main_functional_group_changes",
    "smarts_direction", "molecular_direction", "direction_evidence_type", "source_direction",
    "direction_handling", "reversible_group_id", "direction_variant", "normalized_direction",
    "direction_qc_status", "direction_qc_note", "reverse_template_hash",
    "reverse_transform_available", "reverse_rule_ids", "reverse_rule_relation", "cofactor_or_donor_class",
    "external_participant_roles", "participant_role_confidence", "reaction_representation_scope",
    "transferred_group", "leaving_group_class",
    "donor_class", "acceptor_atom_class", "transferred_group_class",
    "main_pair_projection_method", "main_pair_projection_note",
    "source_ec_numbers", "candidate_ec_numbers", "template_ec_candidates", "database_ec_candidates", "ec_prior_candidates",
    "full_ec_numbers", "partial_ec_numbers", "supported_partial_ec_numbers", "prior_ec_numbers", "ec_status", "ranked_ec_numbers",
    "top_ec_number", "top_ec_confidence", "top_ec_assignment_mode", "top_ec_evidence_types", "top_ec_sources",
    "top_ec_specificity", "top_ec_granularity",
    "top3_ec_numbers", "top3_ec_confidences", "top3_ec_assignment_modes", "top3_ec_sources",
    "top3_ec_specificities", "top3_ec_granularities", "ec_annotation_scope",
    "broad_ec_classes", "broad_ec_class_count", "ec_candidate_count", "ec_conflict_flag", "ec_conflict_level",
    "ec_reaction_type_consistency", "ec_reaction_type_expected_classes", "ec_reaction_type_observed_classes",
    "ec_reaction_type_top_class", "ec_reaction_type_consistency_note", "ec_reaction_type_consistency_mode",
    "ec_directionality_scope", "ec_directionality_warning", "reverse_ec_inheritance_policy",
    "strict_ec_annotation_use", "ec_evidence_summary_json",
    "primary_candidate_families", "secondary_candidate_families", "family_assignment_mode", "family_evidence",
    "family_annotation_available", "family_annotation_confidence", "family_annotation_scope", "family_evidence_sources",
    "evidence_layer_best", "evidence_layers_all", "template_sources", "source_reaction_ids", "source_record_ids",
    "curated_taxol_anchor", "curated_pathway_name", "curated_pathway_step_ids",
    "rhea_ids", "kegg_ids", "metanetx_ids", "template_count", "source_record_count",
    "consensus_group_id", "consensus_generation_mode", "consensus_evidence_rows",
    "consensus_source_database_count", "consensus_evidence_layer_support", "consensus_qc_status",
    "consensus_representative_rule_ids", "consensus_supporting_reaction_types",
    "source_score", "qc_score", "family_score", "reaction_type_score", "direction_score", "scope_score", "final_rule_confidence",
    "strict_core_use", "expanded_use", "exploratory_use", "exclusive_release_tier", "leakage_risk", "notes",
]

ANCHOR_EDGE_COLUMNS = [
    "anchor_edge_id", "source_database", "source_reaction_id", "enzyme_name", "substrate_smiles", "product_smiles",
    "canonical_reaction_smiles", "ec_numbers", "evidence_layer", "direction", "source_reaction_reversibility", "source_file",
    "curated_taxol_anchor", "curated_pathway_name", "curated_pathway_step_id",
]

FAMILY_EVIDENCE_COLUMNS = [
    "evidence_id", "match_type", "match_value", "source_database", "source_reaction_id", "rhea_id", "kegg_id",
    "metanetx_id", "ec_number", "protein_id", "primary_family", "secondary_families", "family", "family_role",
    "evidence_type", "confidence", "evidence_source", "provenance",
]

LAYER_ORDER = {
    "T1_Bio_Core": 1,
    "T2_Bio_Extended": 2,
    "T3_Chem_like": 3,
    "Unknown": 9,
}

EC_CLASS_LABELS = {
    "1": "EC1_oxidoreductase",
    "2": "EC2_transferase",
    "3": "EC3_hydrolase",
    "4": "EC4_lyase",
    "5": "EC5_isomerase",
    "6": "EC6_ligase",
    "7": "EC7_translocase",
}

SMARTS_LIBRARY_COLUMNS = [
    "smarts_rule_id", "rule_id", "reaction_smarts", "reaction_smarts_hash",
    "smarts_library_tier", "exclusive_release_tier", "template_hash", "template_scope", "predictive_rule_use", "template_qc_status",
    "abstracted_from_exact_reaction", "derived_from_exact_anchor", "rxnmapper_confidence",
    "rdchiral_extraction_status", "abstracted_smarts_applies_to_original_pair",
    "exact_abstraction_qc_status", "benchmark_exclusion_flag",
    "reaction_type", "reaction_subtype", "reaction_type_assignment_mode",
    "rule_application_unit", "biochemical_step_granularity", "biochemical_step_granularity_confidence",
    "granularity_assignment_mode", "granularity_evidence_summary", "composite_rule_flag",
    "reaction_center_count", "independent_reaction_center_count", "functional_group_change_count",
    "main_functional_group_changes",
    "smarts_direction", "molecular_direction", "direction_evidence_type", "source_direction",
    "direction_handling", "reversible_group_id", "direction_variant", "normalized_direction",
    "direction_qc_status", "direction_qc_note", "reverse_template_hash",
    "reverse_transform_available", "reverse_rule_ids", "reverse_rule_relation", "cofactor_or_donor_class",
    "external_participant_roles", "participant_role_confidence", "reaction_representation_scope",
    "transferred_group", "leaving_group_class",
    "donor_class", "acceptor_atom_class", "transferred_group_class",
    "main_pair_projection_method", "main_pair_projection_note",
    "source_ec_numbers", "candidate_ec_numbers", "template_ec_candidates", "database_ec_candidates", "ec_prior_candidates",
    "full_ec_numbers", "partial_ec_numbers", "supported_partial_ec_numbers", "prior_ec_numbers", "ec_status", "ranked_ec_numbers",
    "top_ec_number", "top_ec_confidence", "top_ec_assignment_mode", "top_ec_evidence_types", "top_ec_sources",
    "top_ec_specificity", "top_ec_granularity",
    "top3_ec_numbers", "top3_ec_confidences", "top3_ec_assignment_modes", "top3_ec_sources",
    "top3_ec_specificities", "top3_ec_granularities", "ec_annotation_scope",
    "broad_ec_classes", "broad_ec_class_count", "ec_candidate_count", "ec_conflict_flag", "ec_conflict_level",
    "ec_reaction_type_consistency", "ec_reaction_type_expected_classes", "ec_reaction_type_observed_classes",
    "ec_reaction_type_top_class", "ec_reaction_type_consistency_note", "ec_reaction_type_consistency_mode",
    "ec_directionality_scope", "ec_directionality_warning", "reverse_ec_inheritance_policy",
    "strict_ec_annotation_use", "ec_evidence_summary_json",
    "primary_candidate_families", "secondary_candidate_families",
    "family_assignment_mode", "family_evidence",
    "family_annotation_available", "family_annotation_confidence", "family_annotation_scope", "family_evidence_sources",
    "evidence_layer_best", "evidence_layers_all", "template_sources",
    "curated_taxol_anchor", "curated_pathway_name", "curated_pathway_step_ids",
    "source_reaction_ids", "rhea_ids", "kegg_ids", "metanetx_ids",
    "consensus_group_id", "consensus_generation_mode", "consensus_evidence_rows",
    "consensus_source_database_count", "consensus_evidence_layer_support", "consensus_qc_status",
    "consensus_representative_rule_ids", "consensus_supporting_reaction_types",
    "final_rule_confidence", "strict_core_use", "expanded_use", "exploratory_use",
    "template_count", "source_record_count",
    "example_reaction_smiles", "example_substrate_smiles", "example_product_smiles",
    "reaction_delta_fingerprint", "notes",
]


PARTICIPANT_REGISTRY_COLUMNS = [
    "participant_smiles", "participant_hash", "heavy_atom_count", "occurrence_count",
    "reactant_count", "product_count", "source_database_count", "source_databases",
    "source_reaction_ids", "role_class", "role_assignment_mode", "role_confidence",
    "registry_class", "transferred_group", "leaving_group_class", "provenance"
]
