# Compressed Research Release

The transport release preserves the complete tabular scientific record while
avoiding slow, redundant transfer of multi-gigabyte uncompressed files.

## Scientific interpretation boundary

- G0 is the standardized known-taxane seed collection.
- G1 and G2 define the primary near-seed transformation space.
- G3 is retained as an exploratory frontier and is not interpreted as an
  observed-metabolite catalogue. Its downstream bridge status is
  right-censored because no G4 expansion was performed.
- Structural comparisons use molecular fingerprints, physicochemical
  descriptors, functional-state changes, atom-level edits, grammar usage,
  route convergence, and latent bridges. Categorical core-framework
  assignments are not used in the analysis.

## Directly accessible content

- final English manuscript, Methods, figure legends, supplement, and references;
- all main and supplementary figures in PDF, SVG, and 400-dpi PNG;
- validation, provenance, test, and release-audit reports;
- compact primary-grammar summaries and executable rules;
- complete reproducible source code.

## Archived content

Large rule, molecular-state, derivation-event, rejection, path, sensitivity,
and panel-source tables are stored in `archives/*.tar.zst`. Run:

```bash
bash restore_release_archives.sh
```

from the release root to restore every archived directory in place. The
command requires `zstd` and `tar`.

## SQLite reconstruction

The original SQLite file is not duplicated in the transport package because
its five scientific tables are exported losslessly in the G0-G3 archive.
After restoring the archives, rebuild the indexed database with:

```bash
python code/scripts/rebuild_sqlite_from_release.py \
  --space-dir 03_primary_G0_G3 \
  --output 03_primary_G0_G3/taxane_reaction_grammar_space.sqlite
```

`database_table_schema.sql` and `database_index_schema.sql` preserve the
original schema. The rebuild script verifies node and event counts against the
recorded build summary.

## Integrity

`09_release_audit` contains scientific consistency checks and SHA-256 hashes
for the uncompressed source release. `transport_manifest.tsv` and
`transport_checksums.sha256` describe the files actually stored in the
compressed transport release. `transport_exclusions.tsv` records the omitted
redundant SQLite binary and its original SHA-256 hash.
