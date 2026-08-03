# Taxane Reaction Grammar Space

Reproducible source code for the study **A provenance-resolved enzymatic
reaction grammar delineates the generative chemical space of taxanes**.

The software constructs an evidence-stratified enzymatic reaction grammar,
applies the final T1 grammar iteratively to known taxanes, records every
directional parent-rule-product event, validates the generated chemical space,
and reproduces the reported analyses and figures.

The exact upstream database-to-SMARTS implementation used in the study is
included as `components/enzymatic_rule_builder_v0.4.1/`. Raw third-party
database files are not redistributed.

## Scientific scope

- G0 contains 648 standardized known taxanes.
- G1 contains 15,801 unique structures first reached in one accepted grammar
  event.
- G2 contains 223,823 unique structures first reached in two events.
- G3 contains 2,362,766 exploratory structures first reached in three events.
- Full stereochemistry-aware InChIKey is the primary molecular identity.
- G1 and G2 define the primary near-seed chemical space; G3 is an exploratory,
  right-censored frontier.

Generated structures are hypotheses of reaction-grammar accessibility. They
are not asserted to be experimentally observed metabolites or validated
biosynthetic pathways.

EC prediction, genomic annotation, protein prioritization, docking, TDCN, and
scaffold-group analyses are outside this repository's scientific scope.

## Repository layout

```text
src/taxane_reaction_grammar_space/  Core scientific implementation
components/enzymatic_rule_builder_v0.4.1/  Raw-reaction rule builder
tests/                              Unit and integration tests
configs/study.example.yaml          Portable study configuration
docs/                               Protocol and reproducibility documentation
scripts/                            Released-data reconstruction utilities
publication_workflows/              Frozen final-figure and release workflows
```

Large input tables and generated G0-G3 data are distributed in the companion
data release. The data DOI will be added to the tagged software release before
manuscript submission.

## Installation

The frozen analysis used Python 3.10.20 and RDKit 2026.03.2. A portable Conda
specification is provided in `environment.yml`.

```bash
conda env create -f environment.yml
conda activate taxane-grammar-space
python -m pip install -e .
python -m pip install -e 'components/enzymatic_rule_builder_v0.4.1[chem,mapper,dev]'
```

For an existing RDKit-enabled environment:

```bash
python -m pip install -e '.[chem,test]'
```

## Verification

```bash
pytest -q
pytest -q components/enzymatic_rule_builder_v0.4.1/tests
taxane-grammar-space --version
enzymatic-rules --version
```

The manuscript-release helper is an independently packaged utility. Validate
it from its own directory to avoid importing a different installed version:

```bash
cd publication_workflows/release_tools/article_release_builder_v0.3.0
PYTHONPATH=src pytest -q
```

The complete command order and interpretation boundaries are documented in
`docs/REPRODUCIBLE_WORKFLOW.md` and `docs/STUDY_PROTOCOL.md`.

## Released-data reconstruction

After downloading and extracting the companion data record:

```bash
bash restore_release_archives.sh
python scripts/rebuild_sqlite_from_release.py \
  --space-dir 03_primary_G0_G3 \
  --output 03_primary_G0_G3/taxane_reaction_grammar_space.sqlite
```

Every public release includes SHA-256 manifests. Verify those files before
running analyses or comparing numerical results.

## Reproducing the scientific workflow

Copy the portable configuration and update only the input/output paths:

```bash
cp configs/study.example.yaml configs/study.yaml
```

Do not change the frozen scientific parameters when reproducing the published
counts. Execute the stages in `docs/REPRODUCIBLE_WORKFLOW.md` in order.

## License and citation

The source code is released under the BSD 3-Clause License; see `LICENSE`.
Author-generated documentation is distributed under CC BY 4.0. Third-party
software and data remain subject to their respective licenses.

Citation metadata are provided in `CITATION.cff`. Cite both the archived
software release and the companion dataset when reusing the workflow.
