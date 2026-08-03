# Enzymatic Rule Builder v0.4.1

`enzymatic_rule_builder` builds a **database-derived, evidence-stratified, directional reaction SMARTS rule library** for molecular transformation network construction.

The package’s responsibility is deliberately narrow: it builds a **general reaction SMARTS rule library**. It does **not** connect the 600+ taxane nodes. A downstream network builder should apply `reaction_smarts_library.T1_core.tsv` to standardized molecules; a direct computational edge is one `single_smarts_application` that transforms node A into node B.

## Core design principles

1. **Final predictive releases are SMARTS-only.**
   `reaction_smarts_library.*.tsv` contains only generalized reaction SMARTS with `predictive_rule_use = true`; exact substrate-product anchors are separated.

2. **“One step” means one SMARTS application.**
   The release carries `rule_application_unit = single_smarts_application`. This does not claim that every rule is a single elementary enzyme event.

3. **Biochemical granularity is annotated, not assumed.**
   Rules are labeled as `likely_single_enzyme_step`, `possible_composite_step`, `known_composite_step`, or `uncertain` using source metadata, reaction-type text and functional-group delta evidence. Composite rules can remain in the core library, but downstream edge interpretation can distinguish them.

4. **Curated Taxol reactions are T1 curated anchors, not a separate tier.**
   Known Taxol pathway reactions are merged into `T1_Bio_Core` and tagged with `curated_taxol_anchor = true`. Exact-to-SMARTS abstraction is attempted when requested; Taxol-derived SMARTS must pass replay validation and are excluded from external-only Taxol recall benchmarks.

5. **EC confidence and EC granularity are separate.**
   Partial/prefix EC annotations can be high-confidence when supported by high-quality sources. The rule table reports `top_ec_confidence` separately from `top_ec_granularity`.

6. **EC annotations are direction-specific and reaction-type checked.**
   EC labels are not reused blindly for reverse edges. The build reports `ec_directionality_scope`, `reverse_ec_inheritance_policy`, and `ec_reaction_type_consistency` so source EC evidence can be interpreted against the emitted substrate-to-product reaction.

7. **Family evidence is auxiliary.**
   Candidate enzyme-family annotations come only from external evidence tables and are not hard-coded. Family annotation is reported separately and no longer controls whether a SMARTS rule is valid.

8. **Cofactor/donor handling is participant-role annotation.**
   External participants are not silently deleted. The build records `external_participant_roles`, `participant_role_confidence`, `reaction_representation_scope`, `transferred_group`, and `leaving_group_class`. A participant registry can be supplied or derived from database reaction participants for review.

9. **Consensus SMARTS are extracted across the whole database collection.**
   v0.4.1 clusters all QC-passing generalized SMARTS by reaction type, functional-group delta, transfer role, donor/acceptor class, structural delta and direction QC. Promoted consensus rules are data-driven representative SMARTS with full provenance; no enzyme-family mappings are hard-coded.

10. **Direction is quality-controlled.**
   Source direction is normalized to substrate-to-product where possible. Reversed and reversible sources are marked explicitly, while unknown-direction rows are downgraded to exploratory evidence.

11. **T1/T2/T3 release files are mutually exclusive by priority.**
   Contribution tiers are assigned as `T1_only` first, then `T2_only`, then `T3_only`. Cumulative sensitivity files are still written separately as `T2_cumulative` and `T3_cumulative`.

## Installation

```bash
cd enzymatic_rule_builder_v0.4.1
python -m pip install -e '.[dev,chem]'
python -m pytest -q
```

## Main build

```bash
enzymatic-rules build-rules \
  --manifest configs/source_manifest.example.yaml \
  --output-root ./rule_build \
  --data-driven-consensus \
  --require-smarts-rules
```

If Taxol exact anchors are present in the manifest, they are loaded as T1 curated anchors. If RXNMapper/RDChiral are unavailable, the build reports the optional dependency status and keeps Taxol reactions as exact anchors only.

## Build from raw external database root

```bash
enzymatic-rules build-from-raw \
  --external-db-root /path/to/external_databases \
  --taxol-pathway /path/to/taxol_pathway.csv \
  --output-root ./rule_build_full \
  --abstract-exact-reactions \
  --data-driven-consensus \
  --require-smarts-rules
```

`--include-taxol-anchors` is auto-enabled when `--taxol-pathway` is supplied.

## Build a data-derived participant registry

```bash
enzymatic-rules build-participant-registry \
  --manifest configs/source_manifest.example.yaml \
  --output-dir ./participant_registry
```

This produces a TSV/YAML recurring-participant registry that can be reviewed and supplied back to `build-rules` with `--cofactor-yaml`.

## Key outputs

```text
03_rules/general_transformation_rules.audit.tsv
03_rules/data_driven_consensus_candidate_groups.tsv
03_rules/data_driven_consensus_rules.tsv
03_rules/data_driven_consensus_summary.json
04_release/reaction_smarts_library.all.tsv
04_release/reaction_smarts_library.T1_core.tsv
04_release/reaction_smarts_library.T2_extended.tsv
04_release/reaction_smarts_library.T3_exploratory.tsv
04_release/reaction_smarts_library.T1_only.tsv
04_release/reaction_smarts_library.T2_only.tsv
04_release/reaction_smarts_library.T3_only.tsv
04_release/reaction_smarts_library.T2_cumulative.tsv
04_release/reaction_smarts_library.T3_cumulative.tsv
04_release/curated_exact_anchor_edges.tsv
05_benchmark/known_taxol_pathway_recall.external_generalized_only.tsv
```

Use `04_release/reaction_smarts_library.T1_core.tsv` as the default strict downstream network-construction input. `reaction_smarts_library.core.tsv` is still written as a compatibility alias with the same content.
