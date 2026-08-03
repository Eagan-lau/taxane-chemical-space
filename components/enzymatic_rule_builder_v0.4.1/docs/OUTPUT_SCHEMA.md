# Output schema overview

## Final predictive release: `04_release/reaction_smarts_library.*.tsv`

These files are the final reaction SMARTS rule libraries used by downstream network construction.

Inclusion criteria:

```text
template_scope = generalized_template
predictive_rule_use = true
reaction_smarts is non-empty
```

Tier files:

- `reaction_smarts_library.all.tsv`: all generalized predictive SMARTS rules.
- `reaction_smarts_library.T1_core.tsv`: high-confidence biochemical rules for strict network construction; mutually exclusive with T2/T3 releases.
- `reaction_smarts_library.T2_extended.tsv`: rules that did not qualify for T1 but qualify for biochemical extended use.
- `reaction_smarts_library.T3_exploratory.tsv`: rules that did not qualify for T1/T2 but qualify for exploratory use.
- `reaction_smarts_library.T1_only.tsv`, `T2_only.tsv`, and `T3_only.tsv`: explicit mutually exclusive contribution aliases.
- `reaction_smarts_library.T2_cumulative.tsv`: T1 + T2 rules for sensitivity analysis.
- `reaction_smarts_library.T3_cumulative.tsv`: T1 + T2 + T3 rules for sensitivity analysis.
- `reaction_smarts_library.core.tsv`, `expanded.tsv`, and `exploratory.tsv`: compatibility aliases of the T1/T2/T3 files.

`reaction_smarts_rules.*.tsv` files are written as aliases with the same content for older scripts.

Important columns:

```text
rule_id
smarts_library_tier
exclusive_release_tier
template_hash
reaction_smarts
abstracted_from_exact_reaction
derived_from_exact_anchor
rxnmapper_confidence
rdchiral_extraction_status
abstracted_smarts_applies_to_original_pair
exact_abstraction_qc_status
benchmark_exclusion_flag
reaction_type
reaction_subtype
rule_application_unit
biochemical_step_granularity
biochemical_step_granularity_confidence
composite_rule_flag
functional_group_change_count
main_functional_group_changes
smarts_direction
molecular_direction
direction_evidence_type
source_direction
direction_handling
direction_variant
reversible_group_id
normalized_direction
direction_qc_status
direction_qc_note
reverse_template_hash
reverse_transform_available
reverse_rule_ids
reverse_rule_relation
cofactor_or_donor_class
external_participant_roles
participant_role_confidence
reaction_representation_scope
transferred_group
leaving_group_class
candidate_ec_numbers
broad_ec_classes
primary_candidate_families
secondary_candidate_families
family_annotation_available
family_annotation_confidence
family_annotation_scope
curated_taxol_anchor
curated_pathway_name
curated_pathway_step_ids
evidence_layer_best
evidence_layers_all
template_sources
source_reaction_ids
rhea_ids
kegg_ids
metanetx_ids
consensus_group_id
consensus_generation_mode
consensus_evidence_rows
consensus_source_database_count
consensus_evidence_layer_support
consensus_qc_status
consensus_representative_rule_ids
consensus_supporting_reaction_types
final_rule_confidence
strict_core_use
expanded_use
exploratory_use
```

## Audit table: `03_rules/general_transformation_rules.audit.tsv`

This file contains all normalized rule-like records, including:

- generalized reaction SMARTS templates;
- exact substrate-product anchors;
- records that are useful for provenance but not direct SMARTS-based prediction.

Do not use this file directly for pairwise network construction.

## Exact-reaction abstraction report: `01_sources/exact_reaction_abstraction.tsv`

This file is written when `--abstract-exact-reactions` is enabled. It records RXNMapper/RDChiral conversion attempts for selected external exact reactions.

Key columns:

```text
record_id
source_database
source_reaction_id
evidence_layer
canonical_reaction_smiles
mapped_reaction_smiles
rxnmapper_confidence
extraction_status
reaction_smarts
template_qc_status
template_qc_note
```

Only rows with successful mapped extraction, SMARTS QC, and replay validation against the original substrate-product pair are promoted to generalized predictive SMARTS. Failed rows remain exact anchors or provenance records.

## Exact anchors: `04_release/curated_exact_anchor_edges.tsv`

This file stores exact curated pathway reactions such as known Taxol-pathway edges. These records can be used for:

- known-pathway backbone anchoring;
- positive controls;
- recall benchmarking;
- visualization of known reactions.

They are not generalized SMARTS rules.

## Anchor-derived SMARTS: `04_release/reaction_smarts_library.anchor_derived.tsv`

This optional file is written when `--generalize-exact-anchors` is enabled. It contains exact-anchor substrate-product pairs that were successfully generalized with RXNMapper/RDChiral and replay-validated against the original pair.

Companion files:

- `anchor_generalization.report.tsv`: one row per attempted exact anchor, including failures.
- `anchor_generalization.templates.raw.tsv`: raw promoted anchor-derived templates before deduplication.
- `anchor_generalization.templates.deduplicated.tsv`: deduplicated anchor-derived templates.
- `anchor_generalization.rules.audit.tsv`: rule annotations for anchor-derived SMARTS.
- `anchor_generalization.validation.json`: validation summary.
- `anchor_generalization.summary.json`: count summary and runtime settings.

This library is separate from `reaction_smarts_library.T1_core.tsv` and should be treated as exploratory unless manually reviewed or supported by independent evidence.

## Validation files

`reaction_smarts_rules.validation.json` confirms whether the SMARTS release contains only generalized predictive reaction SMARTS records. Empty releases are valid for anchor-only builds but insufficient for downstream network construction.


## EC confidence fields in SMARTS releases

The files under `04_release/reaction_smarts_library.*.tsv` include a de-redundant EC ranking for every SMARTS rule. The key fields are:

| Field | Meaning |
|---|---|
| `top_ec_number` | Highest-confidence EC after evidence aggregation and redundancy suppression. |
| `top_ec_confidence` | Confidence score for `top_ec_number`, based on evidence type and source layer; EC granularity is reported separately and is not a confidence penalty. |
| `top3_ec_numbers` | Up to three de-redundant EC candidates, ordered by confidence. |
| `top3_ec_confidences` | Semicolon-aligned confidence values for `top3_ec_numbers`. |
| `top_ec_assignment_mode` | Strongest assignment mode: `source_direct`, `template_direct`, `database_cross_reference`, `prior_only`, or `missing`. |
| `top_ec_evidence_types` | Evidence types supporting the top EC. |
| `top_ec_sources` | Source database names contributing evidence for the top EC. |
| `full_ec_numbers` | Full four-level ECs among supported candidates. |
| `partial_ec_numbers` | Partial ECs retained for broad-class annotation. |
| `candidate_ec_numbers` | Supported, de-redundant EC candidates from source/template/database evidence; EC priors are excluded. |
| `ec_prior_candidates` | Low-confidence EC priors retained separately. |
| `ec_conflict_level` | `none`, `low`, `medium`, or `high`, based on multiple EC candidates, close top-two scores, and broad EC-class conflicts. |
| `ec_reaction_type_consistency` | QC status comparing supported EC class with the emitted directional reaction type: `consistent`, `mixed_supported_classes`, `mixed_top_inconsistent`, `inconsistent`, `not_assessed`, or `missing_ec_not_assessed`. |
| `ec_reaction_type_expected_classes` | Conservative EC classes expected for the directional reaction type, used only for QC. |
| `ec_reaction_type_observed_classes` | EC classes observed among supported EC candidates. |
| `ec_reaction_type_consistency_note` | Human-readable explanation for the consistency status. |
| `ec_directionality_scope` | Scope of the EC annotation; production SMARTS use `direction_specific_rule`. |
| `ec_directionality_warning` | Warning for reversible, reverse-corrected, or unknown-direction source evidence. |
| `reverse_ec_inheritance_policy` | Explicit instruction that reverse edges must not reuse this EC annotation unless an explicit reverse directional rule exists. |
| `strict_ec_annotation_use` | Boolean flag for whether the EC assignment is suitable for strict enzyme annotation. |

`candidate_ec_numbers` should be treated as a ranked candidate set, while `top_ec_number` is the preferred EC annotation for downstream network construction. `strict_ec_annotation_use=false` means the EC should not be used as a high-confidence enzyme-function claim.

## Family assignment modes

`family_assignment_mode` can be:

- `external_evidence`: specific family supplied by a provenance-tracked family evidence table.
- `ec_supported_family_unassigned`: no specific family is claimed, but high-confidence source/template-supported EC evidence supports enzymatic use; EC granularity is reported in the EC fields.
- `none`: no family evidence and no strict EC fallback.


## Direction fields

`reaction_smarts` in release files is always intended to be applied left-to-right. `direction_handling` records how source direction was normalized: `kept_forward`, `reversed_from_source`, `split_reversible_forward`, `split_reversible_reverse`, or an unknown-direction left-to-right fallback. `reverse_rule_ids` links explicit paired rules when a reversible source was split; it does not imply that reverse chemistry was inferred automatically.

`direction_qc_status` is used for scoring and tiering. Known substrate-to-product, corrected reverse, and split-reversible rows can enter strict tiers. Unknown-direction rows are retained as exploratory evidence unless independently supported.

## v0.4.1 additions

The release uses T1/T2/T3 evidence tiers only:

- `T1_Bio_Core`: direct biochemical or curated pathway evidence, including known Taxol-pathway anchors tagged with `curated_taxol_anchor`.
- `T2_Bio_Extended`: biochemical extended evidence requiring more cautious interpretation.
- `T3_Chem_like`: natural-product-like or chemical-like exploratory evidence.

The SMARTS release also includes data-driven consensus fields:

- `consensus_group_id`
- `consensus_generation_mode`
- `consensus_evidence_rows`
- `consensus_source_database_count`
- `consensus_evidence_layer_support`
- `consensus_qc_status`
- `consensus_representative_rule_ids`
- `consensus_supporting_reaction_types`

Candidate and promoted consensus reports are written under `03_rules/`.

## v0.4.0 additions

The SMARTS release now includes rule-application and biochemical-granularity fields:

- `rule_application_unit`: always `single_smarts_application` for predictive SMARTS rows.
- `biochemical_step_granularity`: one of `likely_single_enzyme_step`, `possible_composite_step`, `known_composite_step`, or `uncertain`.
- `biochemical_step_granularity_confidence`: confidence of the granularity annotation, not rule validity.
- `granularity_assignment_mode`: source of the granularity call, such as `single_functional_group_delta`, `multi_functional_group_delta`, or `source_metadata_known_composite`.
- `composite_rule_flag`: explanatory flag for rules that may represent composite/net transformations.
- `functional_group_change_count` and `main_functional_group_changes`: coarse functional-change annotations used for interpreting a single SMARTS application.

Family annotation is separated from rule validity:

- `family_annotation_available`
- `family_annotation_confidence`
- `family_annotation_scope`

External participant handling is exposed as annotation:

- `external_participant_roles`
- `participant_role_confidence`
- `reaction_representation_scope`
- `transferred_group`
- `leaving_group_class`

These fields are carried into `reaction_smarts_library.*.tsv` so downstream network edges can inherit them as annotations.
