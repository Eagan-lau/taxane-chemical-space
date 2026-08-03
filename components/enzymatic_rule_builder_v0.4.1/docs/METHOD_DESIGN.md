# Method design

The package builds a **database-derived reaction SMARTS library** for one-step enzyme-executable transformations.

## Conceptual separation

The workflow separates four concepts:

1. **Source reactions/templates**: raw database records, including exact reactions and existing templates.
2. **Audit rules**: unified records with provenance, evidence, scores, and use flags.
3. **Predictive SMARTS rules**: generalized reaction SMARTS records that can be applied to molecule nodes.
4. **Exact anchors**: curated exact substrate-product reactions, such as Taxol known-pathway reactions.

Only concept 3 is used to build a predicted network. Exact anchors can be visualized or used for benchmark/known-pathway backbone, but they are not generalized prediction rules.

## Why final release is SMARTS-only

A downstream network builder must evaluate whether two nodes can be connected by applying a reaction rule to the source molecule and checking whether the target molecule is produced. This requires a generalized reaction SMARTS. Exact substrate-product pairs cannot be applied to unseen molecules and therefore cannot be part of the predictive SMARTS library.

## Recommended workflow

```text
enzymatic_rule_builder
  → discover raw Rhea/RetroRules/BioNavi/KEGG/MetaNetX files
  → generate a v0.4.1-native raw manifest
  → optionally generalize exact reactions with RXNMapper/RDChiral
  → project main substrate-product pairs with participant-role annotation
  → extract data-driven consensus SMARTS across all databases
  → evidence-stratified SMARTS-only rule library
  → exact anchors retained separately

taxane_network_builder
  → apply reaction_smarts_library.T1_core.tsv to 600+ taxane nodes
```

v0.4.1 can perform exact-reaction abstraction internally with RXNMapper/RDChiral when `--abstract-exact-reactions` is enabled. This is intended for T1/T2 external exact reactions that have not yet been generalized upstream. `build-from-raw` includes curated Taxol anchors as `T1_Bio_Core` by default when a `--taxol-pathway` file is supplied. Taxol-derived SMARTS are allowed into the production T1 library only after replay validation, but are excluded from external-only Taxol recall benchmarks to avoid leakage.

## Anchor generalization

`--generalize-exact-anchors` runs a separate exact-anchor generalization workflow. It attempts to convert exact anchor substrate-product pairs into generalized reaction SMARTS with RXNMapper/RDChiral, requires SMARTS parsing and replay validation against the original pair, and writes an independent exploratory release:

```text
04_release/reaction_smarts_library.anchor_derived.tsv
04_release/anchor_generalization.report.tsv
04_release/anchor_generalization.summary.json
```

Anchor-derived SMARTS are not mixed into `reaction_smarts_library.T1_core.tsv` by the separate anchor-generalization command. For clean separation, run anchor generalization with `--no-abstract-exact-reactions`; this keeps the production core library based on source generalized SMARTS, while `reaction_smarts_library.anchor_derived.tsv` can be used for exploratory gap filling or sensitivity analysis.

## Confidence tiers

- **T1 core**: high-confidence generalized SMARTS rules for the main predictive network.
- **T2 extended**: rules that did not qualify for T1 but qualify for biochemical extended use.
- **T3 exploratory**: rules that did not qualify for T1/T2 but qualify for exploratory use.
- **exact anchors**: curated exact reactions; not predictive SMARTS rules.

The T1/T2/T3 release files are mutually exclusive by priority for contribution analysis. Cumulative sensitivity releases are written separately as `T2_cumulative` and `T3_cumulative`.

## Data-driven consensus extraction

v0.4.1 adds whole-library consensus extraction. The algorithm starts only from QC-passing generalized SMARTS rows, then clusters them by reaction type, functional-group change, transferred-group class, donor class, acceptor atom class, reaction delta fingerprint and direction QC. A cluster can be promoted only when it meets minimum evidence-row, template-support and source-database thresholds. The promoted consensus rule stores the representative SMARTS and all consensus provenance columns; enzyme-family labels are not hard-coded.

## Data-driven constraint

The package avoids hard-coded reaction-type-to-family mappings. Candidate enzyme family annotations are imported from external family evidence tables and preserved with provenance.

If no external enzyme-family evidence is available, a rule can still enter the strict core when its SMARTS, source/evidence tier, QC/replay status, reaction annotation, and EC evidence are strong. EC specificity and granularity are reported explicitly. Family assignment can remain `not_assigned_at_rule_stage`; this does not invalidate the rule.

## Evidence-layer constraint

Evidence layers are inferred from both source identity and computability:

- direct biochemical SMARTS with EC/Rhea/KEGG/MetaNetX support can be T1;
- annotation-only KEGG/MetaNetX records are retained as support evidence but do not automatically become T1 predictive templates;
- BioNavi-NP BioChem exact reactions are T2 until generalized and validated;
- USPTO/NPL-like records are T3 chemical-like exploration.

## v0.4.1 conceptual update

Taxol known-pathway reactions are no longer modeled as a separate evidence tier. They are `T1_Bio_Core` curated anchors with explicit flags, so the release has a simple T1/T2/T3 hierarchy while preserving Taxol provenance and benchmark leakage control.

Source direction is normalized to substrate-to-product where possible. Reversed source rows and reversible source rows remain usable but are annotated; unknown-direction rows are downgraded to exploratory use.

## v0.4.0 conceptual update

A direct computational transformation is defined as a single reaction SMARTS application, not necessarily a guaranteed single-enzyme elementary event. Some database templates may encode net or composite transformations. These are retained when source/QC evidence is strong, but they are annotated with `biochemical_step_granularity` so downstream edge interpretation can distinguish likely single-enzyme steps from possible or known composite transformations.

Family evidence is auxiliary. The core rule library is gated by SMARTS availability, source quality, QC/replay validation, reaction/EC evidence, and evidence tier. Candidate enzyme-family fields support later genome mining but do not determine rule validity.

Cofactor/donor processing is treated as participant-role annotation. The node-transformation SMARTS may omit external donors/cofactors, but the omitted participants and their inferred roles are reported. Registries should be data-derived and provenance-tracked.
