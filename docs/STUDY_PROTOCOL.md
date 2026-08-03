# Frozen Study Protocol

## Central question

How does a database-derived, evidence-stratified enzymatic reaction grammar
expand the one-, two-, and three-step chemical space accessible from known
taxanes?

## Scope

The study includes reaction grammar construction, grammar validation,
generation of G1-G3 candidate intermediates, chemical-space analysis, and
reaction-path analysis. Molecular fingerprints, physicochemical descriptors,
functional-state vectors, atom-level reaction edits, and directed derivation
topology provide complementary representations of the accessible space.

## State definitions

- G0: standardized known taxanes.
- G1: unique valid products derived from G0 by one accepted grammar event.
- G2: unique valid products first observed from G1.
- G3: unique valid products first observed from G2; interpreted as an
  exploratory frontier rather than a metabolite catalogue.
- Known recovery: a generated structure that matches a G0 structure.
- First-observation convergence: at least two distinct parent structures
  reach a generated structure during the generation in which it is first
  observed. Later rediscovery events are retained only as audit data.
- Latent bridge: an unobserved intermediate on a directed route that leaves
  one known G0 structure, traverses up to three generated layers, and returns
  to a distinct G0 structure.
- Convergent intermediate: a structure satisfying first-observation
  convergence. Rule and semantic-group multiplicity are retained as support
  annotations but do not define structural convergence.
- Structural path: a generation-increasing sequence of distinct source-target
  structure edges. All accepted arrivals to a target at its first-observed
  generation are eligible, including arrivals recorded after the target's
  insertion event.
- Terminal frontier structure: a G3 structure at the stopping horizon.
  Because G3 was not used as a parent layer, its downstream reconnection
  status is right-censored.

## Molecular identity

Full stereochemistry-aware InChIKey is the primary identity for generated
states. Connectivity InChIKey14 is retained only as an auxiliary field for
comparison with the known taxane table. Stereoisomers are never collapsed in
the primary generation counts.

## Primary grammar

The complete T1 evidence library is preserved unchanged. The executable
primary grammar is derived from T1 rules that:

- are marked predictive and strict-core;
- pass template QC and, where applicable, exact-reaction replay;
- encode one reactant and one principal product;
- have a single reaction center supported either by source annotation or by
  recomputation from mapped reaction SMARTS;
- have direction QC;
- are not composite;
- satisfy an assigned transformation semantic and either direct single-step
  evidence, validated consensus evidence, or independent biological-database
  support;
- pass structural growth, product-pattern, mapping-retention, and generic-atom
  gates;
- are compressed into semantic groups defined by reaction delta, atom/bond
  edit signature, direction, transferred-group class, donor class, and
  acceptor class.

Up to three representative SMARTS are retained per semantic group to preserve
substrate-context alternatives without counting database context variants as
independent grammar productions. Domain-seeded representatives activated by
G0 taxanes are united with global representatives so that G1 chemistry is not
lost through representative choice and new grammar productions can activate
dynamically at G2 or G3.

## Expansion and sensitivity layers

- Primary analysis: evidence-gated T1, with G1-G2 defining the near-seed
  chemical-space neighborhood and G3 defining the exploratory frontier;
  reviewed taxane-domain consensus rules remain separately traceable.
- Sensitivity analysis: mutually exclusive T2 and T3 partitions are each
  evaluated independently at G1 and are never pooled with primary counts.

## Product-level QC

Every derivation at every generation must pass:

- RDKit sanitization and valence checks;
- a single-principal-product requirement;
- observed-versus-expected elemental delta agreement;
- source-atom retention;
- configured element and charge constraints;
- global deduplication against all earlier generations;
- immediate reverse-cycle annotation.

No product cap may silently truncate accepted scientific counts. If a
parent-rule application exceeds the configured enumeration guard, the entire
application is rejected as `excessive_enumeration`, retained in the audit, and
reported as a completeness limitation.

## Leakage-controlled validation

Taxol curated anchors are removed globally before known-pathway recovery is
measured. For each benchmark reaction, external rules carrying the exact same
substrate-product connectivity pair are removed as a second leakage gate.
Recovery is summarized by reaction-center fingerprint group, and full
stereochemical recovery is reported separately from InChIKey14 connectivity
recovery. Chemically matched decoy targets from the known taxane collection
provide the negative-control rate.

## Primary outputs

1. Grammar source, QC, and semantic compression tables.
2. G0-G3 nodes, accepted derivation events, and directed edges.
3. Generation-wise attrition and novelty tables.
4. Known recovery, decoy, ablation, and sensitivity benchmarks.
5. Fingerprint and physicochemical chemical-space projections.
6. Functional-state transitions and reaction-center edit distributions.
7. Parent-defined convergence, structure-edge path multiplicity, semantic and
   rule-event audit layers, and latent bridge analyses.
8. Publication figures, source tables, manuscript, and supplementary methods.

Fingerprint projection and nearest-G0 similarity summaries use deterministic,
generation-stratified sampling when a generation exceeds the configured cap.
The population size, sample size, sampling fraction, policy, and random seed
are recorded in source tables.
