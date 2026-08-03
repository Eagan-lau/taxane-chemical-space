# Materials and Methods

## Study design

We asked how a provenance-preserving enzymatic reaction grammar expands the
chemical space accessible from known taxanes over one, two, and three
successive transformations. The analysis was deliberately restricted to
reaction representation, rule abstraction, structure generation, and
chemical-space topology. Molecular fingerprints, physicochemical descriptors,
functional-state vectors, reaction edits, and directed derivation topology
were prespecified as complementary representations of chemical-space
expansion.
The primary interpretation was assigned to the first two generated layers
(G1 and G2). The third layer (G3) was prespecified as an exploratory frontier
and was not interpreted as a catalogue of naturally occurring metabolites.

## Reaction data sources

Biochemical and reaction-chemistry records were assembled from RetroRules
[1], Rhea [2], MetaNetX/MNXref [3], KEGG [4], and the BioChem and
natural-product-like patent subsets distributed with BioNavi-NP [5]. A
curated table of 30 reported Taxol-pathway transformations, informed by
established pathway biochemistry and recent network reconstitution studies
[16,17], was retained as a domain calibration set. Database snapshots were processed independently
before cross-source deduplication, so every retained rule preserved its source
database, source reaction identifiers, evidence tier, abstraction route, and
quality-control history. Source-support counts were treated as overlapping:
a rule supported by more than one database contributed to each corresponding
source count.

Source records were not treated as interchangeable evidence. The BioNavi-NP
USPTO_NPL subset supplied reaction-chemistry support only and did not, in
isolation, establish enzymatic plausibility. KEGG and MetaNetX records in the
available snapshots lacked structure-complete reaction representations for
direct SMARTS abstraction; they were retained for reaction identity,
cross-reference, and annotation provenance but contributed no released rule
rows directly. Source-specific contributions are reported explicitly rather
than inferred from the number of parsed records.

The compiled source collection contained 630,280 normalized reaction records.
Main-substrate and main-product selection produced 630,280 directional
molecular-pair records, which yielded 631,386 raw templates or exact-reaction
anchors and 548,838 deduplicated template or anchor records. Exact reactions
with sufficient structural information were submitted to atom mapping and
reaction-center abstraction. RXNMapper [6] was used for atom correspondence,
and stereochemistry-aware SMARTS extraction and replay followed RDChiral
principles [7]. The release contained 353,524 generalized reaction SMARTS
partitioned into mutually exclusive T1 (192,907), T2 (148,679), and T3
(11,938) release sets. T1 was used for the primary analysis; T2 and T3 were
evaluated separately as evidence-layer sensitivity analyses.

## Participant-role normalization and reaction direction

Reaction records were standardized into substrate-to-product orientation
before abstraction. For reactions containing cofactors, donors, leaving
groups, salts, or other external participants, a role-aware projection
selected the principal transformed substrate and principal product. Omitted
participants remained in rule metadata and were not silently deleted from the
provenance record. Directional fields recorded the source orientation, any
reversal applied during normalization, and the evidence supporting the final
left-to-right transformation.

The executable grammar accepted one principal reactant template and one
principal product template. Multi-reactant source reactions could therefore
contribute a main-pair rule when donor or cofactor participation was
represented in the annotation and the structural edit on the principal
substrate was replayable. Rules with unresolved direction, internally
inconsistent atom mapping, or composite reaction centers were excluded from
the primary grammar and retained in the audit.

## Rule abstraction and evidence stratification

Mapped reaction SMARTS were normalized and deduplicated by canonical SMARTS
hash. For each rule, the reaction center was inferred from mapped atom and
bond changes and independently recomputed from the SMARTS when source
annotations were unavailable. Rule records included the expected elemental
delta, changed-atom signature, mapped-atom retention, reactant and product
context size, reaction-center count, product count, direction status, source
databases, replay outcome, and transformation semantic.

The release tiers were exclusive rather than cumulative. T1 represented the
highest-evidence release partition, T2 an extended biochemical partition, and
T3 an exploratory partition. Tier membership was inherited from the
database-build evidence and quality-control workflow and was not reassigned
post hoc on the basis of taxane activity. Within each tier, a rule was
eligible for semantic compression only if it encoded a readable
single-reactant/single-product SMARTS, passed template quality control, was not
flagged as composite, contained a recoverable edit signature, used no more
than 48 reactant heavy atoms, and had one verified reaction center.

Eligible rules were grouped by a semantic key combining elemental delta,
mapped atom and bond edits, direction, transferred-group class, donor class,
acceptor class, and reaction semantic. Context alternatives were retained by
selecting up to three representatives per T1 semantic group, two per T2 group,
and one per T3 group. This compression reduced database-specific context
variants while retaining local alternatives around the reacting site.

## Construction of the primary executable grammar

The complete T1 release remained unchanged as the archival high-evidence
library. A smaller executable grammar was constructed for prospective
generation. Candidate T1 rules had to satisfy the following gates: assigned
transformation semantics; one inferred reaction center; mapped-atom retention
of at least 0.70; no generic atom in the product pattern; no excessive product
pattern growth; no excessive expected heavy-atom gain; and evidence consistent
with a plausible single biochemical edit. Support for the last criterion was
provided by source annotation, a validated consensus rule, independent
biological-database support, or a recomputed single center with assigned
semantics and no explicit multi-step flag.

Global representatives were combined with representatives activated by the
known taxane collection. This domain-seeded open-grammar design prevented a
globally selected context variant from suppressing a chemically equivalent
variant that matched taxanes, while allowing rules not active at G0 to become
active at later generations. Nine reviewed taxane-domain consensus rules were
added as separately traceable calibration rules. Exact SMARTS duplicates were
removed after combination. The resulting primary grammar comprised 74 rules
in 61 semantic groups: 65 externally derived rules and 9 reviewed
taxane-domain rules. The distinction between the 353,524-rule release and the
74-rule executable grammar was retained throughout all reports.
Rules not matched by any accessible parent remained in the grammar provenance
but contributed neither nodes nor derivation events. This included
sulfur-dependent productions in a chemical space containing no accepted
sulfur-bearing structure; such rules were treated as reachability-negative
within the reported depth.

## Known taxane seed space and molecular identity

The seed collection contained 648 valid taxane structures (G0). Structures
were parsed and sanitized with RDKit [8], canonicalized as isomeric SMILES, and
assigned molecular formulae, exact masses, physicochemical descriptors, and
InChIKeys [9]. Full stereochemistry-aware InChIKey was the primary identity
used for global deduplication. The first 14 InChIKey characters were retained
only as a connectivity-level auxiliary identifier for recovery analyses.
Stereoisomers were not collapsed in generation counts.

Before generation, every executable rule was compiled and screened against
G0 with an RDKit pattern-fingerprint prefilter followed by exact substructure
matching. The screen recorded both parent-rule matches and distinct matching
sites and was checked for prefilter false negatives.

## Iterative chemical-space generation

The 74-rule primary grammar was applied in the forward, substrate-to-product
direction. G1 comprised products generated directly from G0; G2 comprised
structures first observed from G1; and G3 comprised structures first observed
from G2. A structure observed previously remained assigned to its earliest
generation. Every accepted transformation was stored as a directional
parent-rule-product event, including events that returned to G0 or to a
structure first observed in an earlier generation.

For each parent-rule pair, all unique substructure sites were enumerated.
RDKit reaction execution was limited by a guard of 256 products per
parent-rule application. This guard did not truncate accepted products: when
an application exceeded the guard, the complete application was rejected as
`excessive_enumeration` and retained in the audit. Each product tuple then
underwent sanitization, valence checking, single-fragment checking, allowed
element checking, absolute formal-charge checking (maximum 2), source-heavy-
atom retention checking (minimum 0.65), identity-change checking, and
observed-versus-expected elemental-delta checking. Products failing any gate
were written to a rejection table with the parent, rule, generation, raw
product index, and rejection reason.

Accepted products were canonicalized and globally deduplicated by full
InChIKey. Each event stored source and target identities, grammar rule,
semantic group, reaction type, evidence layer, expected and observed formula
deltas, source-atom retention, changed-source-atom count, Morgan-fingerprint
similarity to the parent, G0 full-identity and connectivity matches, and an
immediate reverse-cycle flag. The complete state and event history was stored
in SQLite and exported to tab-separated tables.

## Integrity validation

Generation outputs were checked after completion for primary-key uniqueness,
full-InChIKey uniqueness, valid endpoint references, legal generation
transitions, agreement between event and node generation assignments,
consistency of target-is-new flags, absence of undeclared product truncation,
and agreement between generation summaries and database counts. Validation
was performed before chemical-space analysis. Any failed relational or
completeness check invalidated the corresponding release.

## Leakage-controlled pathway evaluation

Two complementary evaluations were prespecified. First, an external-only
specificity control removed all rules derived from the curated Taxol pathway.
For each of the 30 curated transformations, any remaining external rule
carrying the same substrate-product connectivity pair was also removed. Rule
recovery was evaluated separately at full stereochemical identity and
connectivity identity. Six hundred chemically matched alternative targets
drawn from G0 served as negative controls. Wilson 95% confidence intervals
were calculated for recovery and decoy-match rates.

Second, a domain-informed replay calibration retained the reviewed
taxane-domain rules. Because those rules were informed by the pathway, this
analysis measured internal representational adequacy, not independent pathway
prediction. Fisher's exact test compared curated-pair recovery with decoy
matching.

## Evidence-layer sensitivity analysis

T2 and T3 were never pooled with the primary T1 generation. Each exclusive
release partition was independently subjected to the same single-center
compression, activation, G1 generation, product-level quality control, and
validation. Full-InChIKey structure sets were compared by intersection,
union, and Jaccard similarity. Membership tables retained the exact evidence
layer or combination of layers producing each G1 structure.

## Chemical-space and physicochemical analyses

Generation-wise counts included first-observed structures, cumulative
structures, directional derivation events, active rules, raw products,
rejected products, known-space reconnections, and immediate reverse cycles.
Molecular descriptors comprised exact mass, heavy-atom count, calculated
logP, topological polar surface area, hydrogen-bond donors and acceptors,
rotatable bonds, ring count, fraction sp3 carbon, and formal charge.
Generation-wise shifts were expressed relative to the G0 mean in units of the
G0 standard deviation.

Morgan fingerprints (radius 2, 2,048 bits) [10] were calculated with RDKit.
Maximum Tanimoto similarity to G0 was calculated for every structure in
generations containing no more than 50,000 structures and for a deterministic
within-generation random sample of 50,000 structures in larger generations.
The population size, sample size, sampling fraction, random seed (1729), and
sampling policy were written to the output. A two-dimensional descriptive
projection was obtained with scikit-learn [13] by TruncatedSVD of a
generation-stratified sample of at most 20,000 fingerprints per generation.
The projection was used only for visualization and not for hypothesis
testing.

## Multiscale reaction-space visualization

The complete Figure 2 display was assembled from frozen molecular-state,
derivation, and precomputed layout records; no reaction rule was reapplied and
no structure, descriptor, benchmark, or network edge was recalculated. All G0,
G1, and G2 structures were rendered as individual nodes. Ordinary G0 and G1
nodes used the same marker diameter, G0 was emphasized by colour and drawing
order, and G2 marker diameter was one-half that of G1. Established G0 edges
and accepted G0-to-G1 derivations were drawn explicitly. G3 was represented by
a component-normalized density field and a deterministic density-derived
stipple. Stipple marks encode density and do not represent additional
structures.

To reduce display overplotting without altering topology, a shared softened
empirical-cumulative-distribution coordinate transform was applied to all
plotted network layers. Existing G2 nodes assigned to the largest established
G0 component then underwent a bounded deterministic local relaxation toward
nearby underoccupied display regions. This operation changed display
coordinates only. Node identifiers, molecular structures, generation labels,
parent links, and explicit G0/G1 edges were invariant, and pre-relaxation
coordinates and per-node displacements were retained in the released source
table. Layout coordinates organize topology and local display occupancy; they
are not molecular-descriptor or chemical-distance axes.

## Functional-state and reaction-edit analyses

Each structure was enumerated for ten predefined, structure-based functional
states: free hydroxyl, ester, carboxylic acid or carboxylate, ketone or
aldehyde, ether, amide, epoxide, alkene, aromatic atoms, and phosphate.
These SMARTS-defined states were intentionally non-mutually exclusive; for
example, an ester oxygen can also satisfy the generic ether pattern. Counts
therefore describe reproducible structural features rather than a partition
of atoms into unique functional-group assignments.
For every directed derivation event, the target-minus-source state vector was
computed. Identical vectors were assigned a stable transition code and
summarized by generation, rule, semantic group, source count, and target
count.

Reaction locality was summarized by the number of source atoms whose element,
formal charge, degree, aromaticity, or mapped-neighbor set changed. Observed
molecular-formula deltas were aggregated independently of rule annotations.
These analyses were derived directly from atom-level molecular structures.

## Grammar-use concentration

Accepted events were counted by exact rule and by semantic group for each
generation and across all generations. Concentration was quantified by
Shannon entropy, normalized Shannon entropy, the effective rule number
(exponential Shannon entropy), the Herfindahl-Hirschman index, the Gini
coefficient, and cumulative event fractions contributed by the most used 1,
5, and 10 rules. These metrics separated nominal grammar breadth from the
effective number of productions controlling generated-space expansion.

## Convergence, route multiplicity, and latent bridges

A generated structure was classified as convergent only when at least two
distinct parent structures reached it during the generation in which it was
first observed. Events that rediscovered the same structure in later
generations were retained in separate all-generation audit fields but did not
define first-observation convergence. Distinct grammar-rule and semantic-group
counts were retained as support-diversity annotations but did not define route
convergence, thereby preventing redundant templates from inflating the result.
Immediate route multiplicity was measured primarily by first-observation
parent count.
Multi-generation structural-path counts were computed by dynamic programming
over generation-increasing source-target edges after deduplicating identical
structure pairs. Every accepted event recorded in the target's first-observed
generation was eligible; the insertion-only
`target_is_new` flag was not used as a path filter. Parallel counts over
source-target-semantic-group edges and all parent-rule-product events were
retained as semantic and raw audit layers, respectively.

Latent bridges were defined on a conservative directed acyclic view containing
only generation-increasing edges to newly observed structures plus direct
edges from generated structures back to G0. Same-generation links and reverse
cycles remained in the event audit but were excluded from bridge inflation.
For every generated structure, bit-set propagation recorded all upstream G0
ancestors and downstream G0 recoveries. A latent bridge required at least one
directed pair of distinct G0 structures connected through that intermediate.
Counts were summarized by generation, intermediate, and directed G0 pair.
Because G3 structures were not used as parents in an additional expansion
round, their downstream bridge status was right-censored and was not
interpreted as an observed absence of bridging.

## Statistical analysis and visualization

All analyses were deterministic given the input snapshots and random seed.
Counts are reported exactly. Wilson score intervals were used for binomial
rates, and two-sided Fisher exact tests were evaluated with SciPy [14] for
curated-versus-decoy recovery comparisons. No null-hypothesis test was
applied to descriptive chemical-space projections. Tabulation used pandas
[11], and numerical array operations used NumPy [15]. Figures were generated
with Matplotlib [12] and exported as editable SVG, Type 42 font-embedded PDF,
and 400-dpi PNG. Every panel was accompanied by a tab-separated source-data
table and provenance manifest.

The final analysis environment used Python 3.10.20, RDKit 2026.03.2, pandas
2.3.3, NumPy 2.2.6, SciPy 1.15.3, scikit-learn 1.7.2, and Matplotlib 3.10.9.
Exact package versions, command-line parameters, input hashes, and output
hashes were recorded with the release.
