# Reproducible computational workflow

This document records the executable order of operations for the taxane
reaction-grammar study. It complements `STUDY_PROTOCOL.md`, which defines the
scientific estimands and interpretation boundaries.

## 0. Database-to-SMARTS rule-library construction

The exact upstream implementation is bundled in
`components/enzymatic_rule_builder_v0.4.1/`. Install it together with the
optional atom-mapping dependencies:

```bash
python -m pip install -e \
  'components/enzymatic_rule_builder_v0.4.1[chem,mapper,dev]'
```

Prepare the external database snapshots according to the source manifest and
their original provider terms, then run:

```bash
enzymatic-rules build-from-raw \
  --external-db-root /path/to/external_databases \
  --taxol-pathway /path/to/taxol_pathway.csv \
  --output-root ./rule_build_full \
  --abstract-exact-reactions \
  --data-driven-consensus \
  --require-smarts-rules
```

This stage writes mutually exclusive `T1_only`, `T2_only`, and `T3_only`
reaction-SMARTS releases, exact-anchor audits, direction and EC quality-control
fields, and source-provenance summaries. Verify the input and output hashes
against the companion data record before using the generated rules in the
downstream workflow. Raw third-party records are not part of the software
release.

## 1. Environment and input snapshots

Install the package in an RDKit-enabled Python 3.10 environment:

```bash
python -m pip install -e '.[chem,test]'
```

Record software versions and verify all configured input hashes before any
scientific computation:

```bash
taxane-grammar-space record-environment \
  --study-config configs/study.yaml \
  --output-dir "${STUDY_ROOT}/00_environment_and_inputs"
```

## 2. Primary T1 grammar

Prepare the mutually exclusive T1 release partition and the reviewed
taxane-domain consensus rules:

```bash
taxane-grammar-space prepare-rules \
  --rules "${RULES_T1}" \
  --release-tier T1 \
  --representatives-per-group 3 \
  --max-reactant-atoms 48 \
  --require-single-center \
  --output-dir "${STUDY_ROOT}/01_T1_grammar"

taxane-grammar-space prepare-domain-rules \
  --rules "${TAXANE_DOMAIN_RULES}" \
  --output-dir "${STUDY_ROOT}/01_taxane_domain_grammar"
```

Apply the evidence and structural gates once to the complete T1
representatives and once to the subset activated by G0:

```bash
taxane-grammar-space select-grammar \
  --activated "${STUDY_ROOT}/01_T1_grammar/generative_grammar.T1_primary.tsv" \
  --representatives-per-group 3 \
  --output-dir "${STUDY_ROOT}/02_global_selection"

taxane-grammar-space screen-grammar \
  --grammar "${STUDY_ROOT}/01_T1_grammar/generative_grammar.T1_primary.tsv" \
  --nodes "${KNOWN_TAXANES}" \
  --output-dir "${STUDY_ROOT}/02_G0_screen"

taxane-grammar-space select-grammar \
  --activated "${STUDY_ROOT}/02_G0_screen/generative_grammar.T1_G0_activated.tsv" \
  --representatives-per-group 3 \
  --output-dir "${STUDY_ROOT}/02_G0_selection"
```

Combine global and G0-compatible context representatives, then add the
separately traceable reviewed domain rules:

```bash
taxane-grammar-space assemble-open-grammar \
  --global-selected "${STUDY_ROOT}/02_global_selection/taxane_activated_grammar.primary.tsv" \
  --g0-selected "${STUDY_ROOT}/02_G0_selection/taxane_activated_grammar.primary.tsv" \
  --representatives-per-group 4 \
  --output-dir "${STUDY_ROOT}/02_open_grammar"

taxane-grammar-space augment-domain-grammar \
  --open-grammar "${STUDY_ROOT}/02_open_grammar/taxane_open_grammar.tsv" \
  --domain-grammar "${STUDY_ROOT}/01_taxane_domain_grammar/taxane_domain_grammar.primary.tsv" \
  --output-dir "${STUDY_ROOT}/02_primary_grammar"
```

The frozen study grammar contains 65 external productions and 9 reviewed
domain productions. Exact SMARTS duplicates are removed after combination.

## 3. Primary G0-G3 generation

Generate the direction-specific derivation graph:

```bash
taxane-grammar-space generate-space \
  --grammar "${STUDY_ROOT}/02_primary_grammar/taxane_reaction_grammar.primary.tsv" \
  --nodes "${KNOWN_TAXANES}" \
  --max-generation 3 \
  --max-products-per-parent-rule 256 \
  --min-source-atom-retention 0.65 \
  --max-abs-formal-charge 2 \
  --output-dir "${STUDY_ROOT}/03_primary_G0_G3"
```

Use `--resume` after an interruption. Completed generations and committed
within-generation parent checkpoints are retained.

Validate all relational and generation invariants:

```bash
taxane-grammar-space validate-space \
  --database "${STUDY_ROOT}/03_primary_G0_G3/taxane_reaction_grammar_space.sqlite" \
  --output-dir "${STUDY_ROOT}/04_validation"
```

## 4. Evidence-layer sensitivity

Repeat rule preparation, selection, G0 screening, and open-grammar assembly
independently for `RULES_T2` and `RULES_T3`. Generate only G1 for each
partition. Do not pool these structures with primary counts.

```bash
taxane-grammar-space compare-g1-sensitivity \
  --primary-nodes "${STUDY_ROOT}/03_primary_G0_G3/chemical_space_nodes.tsv" \
  --primary-summary "${STUDY_ROOT}/03_primary_G0_G3/chemical_space_build_summary.json" \
  --t2-nodes "${T2_ROOT}/05_G1_space/chemical_space_nodes.tsv" \
  --t2-summary "${T2_ROOT}/05_G1_space/chemical_space_build_summary.json" \
  --t3-nodes "${T3_ROOT}/05_G1_space/chemical_space_nodes.tsv" \
  --t3-summary "${T3_ROOT}/05_G1_space/chemical_space_build_summary.json" \
  --output-dir "${STUDY_ROOT}/05_G1_evidence_sensitivity"
```

## 5. Leakage control and domain replay

Run the external-only benchmark with exact-pair leakage exclusion:

```bash
taxane-grammar-space benchmark \
  --activated-grammar "${EXTERNAL_BENCHMARK_GRAMMAR}" \
  --taxol-pathway "${TAXOL_PATHWAY}" \
  --nodes "${KNOWN_TAXANES}" \
  --decoys-per-positive 20 \
  --exclude-taxol-derived \
  --output-dir "${STUDY_ROOT}/05_external_leakage_control"
```

Run the same command with the reviewed domain grammar and
`--no-exclude-taxol-derived` for internal replay calibration. The latter is
not an independent pathway-prediction result.

## 6. Structure-based analysis

```bash
taxane-grammar-space analyze-space \
  --nodes "${STUDY_ROOT}/03_primary_G0_G3/chemical_space_nodes.tsv" \
  --events "${STUDY_ROOT}/03_primary_G0_G3/derivation_events.tsv" \
  --application-audit "${STUDY_ROOT}/03_primary_G0_G3/rule_application_audit.tsv" \
  --rejections "${STUDY_ROOT}/03_primary_G0_G3/rejection_events.tsv" \
  --projection-max-nodes-per-generation 20000 \
  --similarity-max-nodes-per-generation 50000 \
  --random-seed 1729 \
  --output-dir "${STUDY_ROOT}/06_analysis"
```

Run the same analysis command for the independent T2 and T3 G1 spaces.

## 7. Figures and manuscript

```bash
taxane-grammar-space render-study-figures \
  --analysis-dir "${STUDY_ROOT}/06_analysis" \
  --provenance-dir "${STUDY_ROOT}/01_rule_library_provenance" \
  --sensitivity-dir "${STUDY_ROOT}/05_G1_evidence_sensitivity" \
  --external-benchmark-dir "${STUDY_ROOT}/05_external_leakage_control" \
  --domain-benchmark-dir "${STUDY_ROOT}/05_domain_replay" \
  --t2-analysis-dir "${T2_ROOT}/06_analysis" \
  --t3-analysis-dir "${T3_ROOT}/06_analysis" \
  --output-dir "${STUDY_ROOT}/07_figures"

taxane-grammar-space render-manuscript \
  --template "${DATA_RELEASE}/manuscript/Main_Manuscript_with_Figures.md" \
  --analysis-dir "${STUDY_ROOT}/06_analysis" \
  --output-dir "${STUDY_ROOT}/08_manuscript"
```

Copy the frozen Methods, references, figure legends, and supplementary
information into the same manuscript release directory.

## 8. Release audit

The release is complete only after:

1. G0-G3 generation and all relational checks pass;
2. every main and supplementary panel has a source table;
3. PDF, SVG, and 400-dpi PNG are present for every figure;
4. no final manuscript placeholder remains;
5. input and output hashes are recorded; and
6. all unit and integration tests pass.

After assembling the versioned release directory, run:

```bash
taxane-grammar-space audit-release \
  --code-dir "${CODE_ROOT}" \
  --provenance-dir "${STUDY_ROOT}/00_provenance" \
  --space-dir "${STUDY_ROOT}/03_primary_G0_G3" \
  --validation-dir "${STUDY_ROOT}/04_validation" \
  --analysis-dir "${STUDY_ROOT}/06_analysis" \
  --figures-dir "${STUDY_ROOT}/07_figures" \
  --manuscript-dir "${STUDY_ROOT}/08_manuscript" \
  --output-dir "${STUDY_ROOT}/09_release_audit"
```

The command writes a machine-readable check table, file manifest, and SHA-256
checksum list.
