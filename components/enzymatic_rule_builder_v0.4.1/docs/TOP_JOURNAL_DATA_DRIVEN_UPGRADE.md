# v0.4.1 Top-Journal Data-Driven Rule Build

This package builds a database-derived, auditable enzymatic reaction SMARTS
library for downstream molecular transformation networks.

## Core Policy

1. Evidence tiers are T1/T2/T3 only.
   - T1_Bio_Core: curated Taxol anchors, Rhea/RetroRules/KEGG/MetaNetX-supported biochemical rules.
   - T2_Bio_Extended: biochemical exact reactions and weaker but biological evidence such as BioNavi BioChem.
   - T3_Chem_like: exploratory chemical/NP-like evidence such as USPTO_NPL.

2. Curated Taxol pathway reactions are not a separate T0 tier.
   - They are T1_Bio_Core rows with `curated_taxol_anchor=true`.
   - Exact anchor identity is preserved with `anchor_edge_use=true`.
   - If generalized into SMARTS, the rule carries `benchmark_exclusion_flag=curated_taxol_anchor_derived_exclude_from_external_recall`.

3. Direction is explicit.
   - Released SMARTS are always left-to-right as substrate -> product.
   - Reverse source reactions are corrected.
   - Reversible source reactions are split into explicit directional members.
   - Unknown-direction rules are marked `direction_qc_unknown_exploratory_only` and cannot enter T1/T2 strict releases.

4. Multi-participant reactions are projected by role.
   - Full reaction context is retained in annotation fields.
   - Main-pair projection removes external donors/carriers such as acyl-CoA, nucleotide-sugars, ATP-like phosphate donors and SAM-like methyl donors.
   - Projection fields include `donor_class`, `acceptor_atom_class`, `transferred_group_class`, `leaving_group_class`, and `main_pair_projection_method`.

5. Consensus SMARTS are data-driven.
   - All QC-passing generalized SMARTS are clustered across databases by reaction type, functional-group delta, transferred group, donor/acceptor class, structural delta and direction QC.
   - Every cluster is written to `03_rules/data_driven_consensus_candidate_groups.tsv`.
   - Clusters passing evidence thresholds are promoted to `03_rules/data_driven_consensus_rules.tsv` and included in final SMARTS releases.
   - This flow does not hard-code enzyme-family mappings.

6. T1/T2/T3 releases are mutually exclusive for contribution analysis.
   - First assign `T1_only` when a rule qualifies for strict core use.
   - If not T1, assign `T2_only` when it qualifies for expanded biochemical use.
   - If not T1/T2, assign `T3_only` when it qualifies for exploratory use.
   - Cumulative sensitivity files are written separately as `T2_cumulative` and `T3_cumulative`.

7. EC annotations are evidence-based, reaction-type checked, and direction-specific.
   - `ec_reaction_type_consistency` records whether the supported EC class is compatible with the emitted directional reaction type.
   - Incompatible or top-inconsistent EC annotations cannot be used as strict EC claims.
   - `reverse_ec_inheritance_policy` records that reverse edges must use explicit reverse directional rule evidence rather than blindly reusing the forward EC.

## Main Outputs

- `04_release/reaction_smarts_library.T1_core.tsv`
- `04_release/reaction_smarts_library.T2_extended.tsv`
- `04_release/reaction_smarts_library.T3_exploratory.tsv`
- `04_release/reaction_smarts_library.T1_only.tsv`
- `04_release/reaction_smarts_library.T2_only.tsv`
- `04_release/reaction_smarts_library.T3_only.tsv`
- `04_release/reaction_smarts_library.T2_cumulative.tsv`
- `04_release/reaction_smarts_library.T3_cumulative.tsv`
- `04_release/reaction_smarts_library.all.tsv`
- `03_rules/data_driven_consensus_candidate_groups.tsv`
- `03_rules/data_driven_consensus_rules.tsv`
- `04_release/curated_exact_anchor_edges.tsv`
- `build_summary.json`

## Recommended Build

```bash
enzymatic-rules build-from-raw \
  --external-db-root /path/to/external_databases \
  --taxol-pathway /path/to/taxol_pathway.csv \
  --output-root ./rule_build_full \
  --abstract-exact-reactions \
  --data-driven-consensus \
  --require-smarts-rules
```

The primary downstream network input for strict analysis is:

```text
04_release/reaction_smarts_library.T1_core.tsv
```

Use `T2_extended` for evidence-expanded analysis and `T3_exploratory` only for
hypothesis generation.
