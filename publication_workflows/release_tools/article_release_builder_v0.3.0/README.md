# Taxane article release builder v0.3.0

This package assembles the frozen taxane reaction-grammar analysis into a
submission-oriented, auditable article release. It does not rerun rule
construction, reaction-SMARTS extraction, G0-G3 enumeration, fingerprint or
descriptor calculations, convergence analysis, or bridge analysis.

Version 0.3.0:

- rewrites the manuscript around four result-led scientific arguments;
- uses Figure 2 V12, in which G0-G2 are node resolved and G3 is retained as
  a complete density layer;
- discloses the bounded G2 display-coordinate relaxation as a visualization
  operation that does not alter molecular identities, parents, or edges;
- places every main-figure caption immediately below its figure;
- emits an independent figure-legends document;
- consolidates 15 supplementary figures into eight thematic figures;
- consolidates 22 supplementary tables into 11 logical workbook sheets;
- preserves every original supplementary figure and table in the source-data
  archive;
- exports complete G0-G3 molecular CSV files, including SMILES;
- validates headline counts, figure audits, captions, cross-references,
  workbook topology, molecular CSVs, and the frozen hard-linked input snapshot.

## Run

```bash
build-taxane-article-release \
  --primary-release /path/to/G0_G3_primary_release \
  --editorial-release /path/to/editorial/release \
  --output /path/to/new/article_release
```

The output directory must be empty or absent. Use the same frozen primary and
editorial releases to reproduce the assembled package.
