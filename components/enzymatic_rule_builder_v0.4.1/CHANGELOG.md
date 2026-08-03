# Changelog

## v0.4.1

- Replaced the separate Taxol T0 tier with `T1_Bio_Core` curated anchors tagged by `curated_taxol_anchor`, `curated_pathway_name`, and pathway-step provenance.
- Added direction QC fields and scoring so source-defined, reversed-to-substrate-product, reversible-split, and unknown-direction rows are distinguishable in releases.
- Expanded role-aware main-pair projection beyond acyl-CoA to recurring transfer donors such as nucleotide-sugar/phosphate carriers and SAM-like methyl donors, while preserving external-participant annotation.
- Added whole-library data-driven consensus extraction. QC-passing SMARTS are clustered by reaction type, functional-group delta, transfer role, donor/acceptor class, structural delta and direction QC, then promoted as auditable representative consensus rules.
- Added canonical T1/T2/T3 release files: `reaction_smarts_library.T1_core.tsv`, `reaction_smarts_library.T2_extended.tsv`, and `reaction_smarts_library.T3_exploratory.tsv`. Legacy `core/expanded/exploratory` names remain compatibility aliases.
- Updated raw-build metadata and CLI defaults for v0.4.1-native builds from raw external databases without using old `network_ready` intermediates.

## v0.4.0

- Added a separate `anchor_generalization` workflow. `--generalize-exact-anchors` attempts RXNMapper/RDChiral generalization of exact anchors across selected evidence layers and writes `reaction_smarts_library.anchor_derived.tsv` plus reports, without mixing anchor-derived SMARTS into the production core release.
- Defined the downstream network unit as `single_smarts_application`; every SMARTS-release row now carries `rule_application_unit`.
- Added biochemical step-granularity annotation fields: `biochemical_step_granularity`, confidence, assignment mode, evidence summary, `composite_rule_flag`, reaction-center/functional-change counts, and main functional-group changes.
- Kept composite database SMARTS instead of discarding them, while labeling likely single-enzyme, possible composite, known composite, and uncertain rule granularity.
- Decoupled enzyme-family evidence from rule validity and final rule-confidence scoring; family annotation is retained as auxiliary downstream genome-mining support.
- Added `family_annotation_available`, `family_annotation_confidence`, and `family_annotation_scope` fields.
- Added a data-derived participant-role registry workflow. `build-rules` can derive a provenance-tracked recurring participant registry when no cofactor registry is supplied, and `build-participant-registry` exports TSV/YAML registries for review.
- Treated cofactor/donor handling as external-participant role annotation rather than silent deletion; SMARTS releases now retain participant roles, confidence, representation scope, transferred group, and leaving-group fields.
- Production default: when curated Taxol exact anchors are present, exact-to-SMARTS abstraction is attempted by default; successful anchor-derived SMARTS still require replay validation before release and remain excluded from external-only recall benchmark.
- `build-from-raw` auto-includes Taxol anchors when `--taxol-pathway` is supplied and uses T1/T2 as default exact-abstraction layers.
- Continued SMARTS-only release validation: exact anchors remain separated in `curated_exact_anchor_edges.tsv`; network construction should consume `reaction_smarts_library.core.tsv`.

## v0.3.7 and earlier

- Added directional molecular transformation handling and removed `is_reversible` from SMARTS releases.
- Added exact-derived SMARTS replay validation.
- Added top1/top3 EC ranking with separate confidence and granularity.
- Added SMARTS-only release files and anchor separation.
