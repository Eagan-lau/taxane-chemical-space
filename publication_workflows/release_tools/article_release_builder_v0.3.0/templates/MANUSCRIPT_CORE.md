# A provenance-resolved enzymatic reaction grammar delineates the generative chemical space of taxanes

**Article type:** Research Article

**Target journal:** Journal of Cheminformatics

**Authors:** AUTHOR_INPUT_REQUIRED

**Affiliations:** AUTHOR_INPUT_REQUIRED

**Corresponding author:** AUTHOR_INPUT_REQUIRED

## Abstract

### Background

Reaction databases document individual transformations but do not directly
define the chemical neighborhoods reachable from a natural-product family.
Turning heterogeneous records into a generative model requires explicit
control of participant roles, direction, atom correspondence, reaction-center
context, evidence provenance, and product validity. The central challenge is
therefore not to accumulate the largest possible template collection, but to
derive a compact, auditable grammar that can be transferred beyond the
reactions from which it was learned.

### Results

We normalized 630,280 reaction records and obtained 353,524 directional
reaction SMARTS in mutually exclusive T1, T2, and T3 evidence partitions.
Semantic compression and executability filters produced grammars of 74, 17,
and 1,568 productions, respectively. T1 combined the highest evidence level
with the broadest accepted one-step coverage and was selected prospectively
for multigeneration enumeration. Application of T1 to 648 known taxanes (G0)
yielded 15,801 one-rule products (G1), 223,823 two-rule products (G2), and
2,362,766 exploratory three-rule products (G3). The resulting space comprised
2,603,038 full-InChIKey-distinct structures connected by 7,928,209 accepted
directional derivation events. Despite this expansion, 88.7-89.4% of events
at each generation altered one source atom under a structure-derived locality
definition, and a small set of oxygenation and acyl-state productions
dominated grammar use. First-observation convergence increased from 1.1% in
G1 to 82.5% in G2 and 97.8% in G3. A conservative acyclic analysis identified
227 generated bridge candidates supporting 309 ordered connections between
known taxanes.

### Conclusions

A provenance-resolved reaction grammar transforms a static reaction archive
into an explicit generative model of taxane chemical space. The accessible
space is combinatorially large but chemically low-dimensional and
topologically convergent. The released structures, directional derivations,
quality-control outcomes, and bridge hypotheses provide an auditable basis
for targeted metabolite searching and reaction testing while preserving the
distinction between grammar accessibility, enzymatic feasibility, and natural
occurrence.

**Keywords:** taxane; reaction SMARTS; reaction grammar; chemical-space
enumeration; natural products; molecular network; latent intermediates;
cheminformatics

## Introduction

Natural-product diversity is produced by a limited repertoire of chemical
operations acting repeatedly on structurally elaborate molecular frameworks.
Oxidation, acyl transfer, hydrolysis, glycosylation, methyl transfer, and
related local edits can generate families containing hundreds of observed
molecules and a much larger set of plausible but unsampled intermediates.
Conventional molecular catalogues describe the observed members, whereas
reaction databases describe isolated precedents. Neither representation alone
answers a generative question: what chemical space becomes accessible when
transferable biochemical edits are iterated around a defined molecular family?

Reaction SMARTS offer a machine-executable representation of local molecular
rewrites, and atom-mapped reaction corpora have enabled their extraction for
retrosynthesis, pathway design, and metabolic-space exploration [1, 6, 7].
However, the apparent breadth of a rule collection can be inflated by
database-specific context, duplicated reactions, ambiguous participant roles,
reversed direction, incomplete stereochemistry, and templates that fail when
reapplied. A rule set can therefore be large yet chemically narrow, or
apparently productive because a few permissive templates dominate its output.
Prospective exploration requires a sharper distinction between an archival
rule resource, which preserves breadth and provenance, and an executable
grammar, which contains the context-aware productions admitted to generation.

Evidence quality is equally consequential. Biochemically curated reactions,
integrated cross-reference records, and reaction-chemistry corpora do not
provide interchangeable support for enzymatic transferability. Family-specific
curation can repair biologically important blind spots, but it also creates a
potential leakage path if recovery of the same pathway is presented as
independent prediction. A defensible workflow must therefore stratify evidence
before enumeration, retain source provenance through abstraction, and report
domain calibration separately from external generalization.

Taxanes provide a demanding system in which to test this framework [16, 17].
Their dense stereochemistry, polyoxygenated cores, repeated acylation states,
and large substituents make molecular identity and edit direction
consequential. The known taxane library contains many locally related
structures, yet the intermediate chemical space linking those observations
remains sparsely sampled. This combination makes taxanes suitable for asking
whether a compact reaction grammar can generate a large chemical neighborhood
without losing a transparent connection to local chemistry.

Here, we constructed a provenance-resolved, evidence-stratified reaction
resource and compiled it into executable T1, T2, and T3 grammars. We selected
T1 for primary exploration using evidence level, product acceptance, and
one-step coverage rather than grammar size alone. Starting from 648 known
taxanes, we enumerated one-, two-, and three-rule neighborhoods while retaining
every accepted parent-rule-product event and every rejected application. We
then quantified expansion, edit locality, functional-state changes, grammar
concentration, convergence, reverse cycles, known-space reconnection, and
latent bridges. The analysis does not assign enzymes or claim that generated
structures occur in nature; it maps the chemical hypotheses implied by an
auditable reaction grammar.

## Results

### A provenance-resolved build separates archival breadth from executable chemical grammar

The integrated collection comprised 630,280 normalized reaction records from
RetroRules, Rhea, MetaNetX, KEGG, BioNavi-NP BioChem, BioNavi-NP USPTO_NPL,
and a curated Taxol-pathway set (Fig. 1A; Table S1A). Source identity was
retained because these resources contributed distinct forms of evidence.
Biochemical curation, integrated reaction identity, rule-level support, and
reaction-chemistry precedent were therefore recorded separately rather than
collapsed into an undifferentiated source count.

Role-aware principal-pair projection, directional normalization,
deduplication, atom mapping, reaction-center abstraction, replay testing, and
quality control yielded 353,524 generalized directional reaction SMARTS
(Fig. 1B; Table S1B). The release was partitioned into 192,907 T1, 148,679 T2,
and 11,938 T3 rows. These mutually exclusive partitions preserve archival
breadth and row-level provenance. Subsequent semantic compression and
executability gates produced much smaller grammars containing 74 T1, 17 T2,
and 1,568 T3 productions. The difference between release size and grammar size
is substantive: the former measures traceable reaction-rule evidence, whereas
the latter measures the nonredundant productions admitted to prospective
generation.

The grammar abstraction is illustrated by six complete Taxol-pathway
substrate-product pairs (Fig. 1C). Each molecular transformation is reduced to
a directional local production that retains the reaction center and enough
neighboring context to be chemically constrained while remaining transferable
to another compatible substrate. The examples span C-H hydroxylation,
O-acetyl transfer, alcohol oxidation, O-deacetylation, O-benzoyl transfer, and
N-benzoyl transfer. Thus, a grammar production is neither a whole-reaction
lookup nor an unconditioned functional-group label; it is an executable,
context-aware molecular rewrite.

Independent G0 application differentiated the three evidence layers. T1, T2,
and T3 generated 15,801, 3,069, and 2,973 unique G1 structures, respectively,
with accepted raw-product fractions of 51.3%, 51.0%, and 41.4% (Fig. 1D;
Table S1C). The T1-T2 G1 Jaccard similarity was 0.094, whereas comparisons
involving T3 were approximately 0.002 (Fig. 1E; Table S1D). Lower tiers
therefore accessed different neighborhoods rather than simply reproducing
weaker copies of T1. T1 was selected for primary enumeration because it
combined the highest evidence tier with the broadest accepted one-step space;
T2 and T3 were retained as sensitivity and exploratory comparators (Fig. S2).

Two pathway evaluations defined the scope of that choice (Table S1E). After
removal of all Taxol-derived rules and any external rule carrying the exact
benchmark connectivity pair, the strict external-only grammar recovered none
of 30 curated pathway reactions and matched none of 600 decoys. In contrast,
the domain-informed grammar replayed 23 of 30 curated pairs at connectivity
level and 15 at full stereochemical identity, while matching 4 of 600 decoys
(odds ratio 489.6; Fisher exact P = 5.47 x 10^-32). The first result exposes
a genuine external-coverage limitation; the second establishes internal
representational adequacy after domain calibration. It is not interpreted as
independent pathway discovery.

{{FIGURE_1}}

### A compact T1 grammar unfolds a 2.6-million-structure neighborhood around known taxanes

Applying the 74-production T1 grammar to 648 known taxanes produced 15,801
unique G1 structures through 18,504 accepted directional events. Iteration
over the first-observed frontier yielded 223,823 G2 structures through 490,855
events and 2,362,766 G3 structures through 7,418,850 events (Table S3). The
combined space comprised 2,603,038 full-InChIKey-distinct structures. Full
stereochemical identity was used for global deduplication; the first 14
InChIKey characters were used only for explicitly designated
connectivity-level recovery analyses.

Figure 2 resolves all 648 G0, 15,801 G1, and 223,823 G2 structures as
individual nodes, together with 1,670 established G0 edge records and 15,991
accepted G0-to-G1 derivation records. The 2,362,766 G3 structures are retained
at their complete numerical scale as a component-normalized descendant-density
background. G2 is displayed as smaller neutral-gray nodes so that the
two-rule neighborhood remains countable without obscuring G0/G1 topology.
A bounded deterministic relaxation redistributes overplotted G2 display
coordinates within the largest G0 component; it changes neither molecular
identity, generation membership, parent linkage, nor any explicit edge.
Representative G0 and accepted G1 structures link the abstract network to
specific molecular edits, and the logarithmic generation inset exposes the
scale separation between the observed seed collection and the exploratory
frontier.

Expansion was large but not uniformly efficient. Unique new structures per
raw product declined from 0.438 at G1 to 0.215 at G2 and 0.142 at G3
(Table S3). The decline did not result from silent truncation. Parent-rule
applications that exceeded the prespecified enumeration guard were rejected
as complete applications, and all sanitization, elemental-delta, retention,
charge, and fragmentation failures were retained in the product audit
(Fig. S3; Tables S4 and S5). Increasing rediscovery and rejection therefore
accompanied the expanding frontier.

These results motivate a depth-dependent interpretation. G1 is the immediate
reaction-grammar neighborhood of known taxanes, and G2 captures combinations
of two transferable edits while remaining linked to a resolved parent
history. G3 reveals the scale and topology of continued expansion, but it is
more exposed to propagated rule promiscuity and was not expanded further as a
parent generation. It is therefore reported as exploratory rather than as an
inventory of plausible natural metabolites.

{{FIGURE_2}}

### A small vocabulary of local edits governs multigeneration expansion

The scale of the generated space did not arise from chemically opaque
replacement events. Two explicit G0-G3 trajectories show that each
multigeneration path can be decomposed into released, directional
parent-rule-product records (Fig. 3A). Every product retains its immediate
parent, transferred production, generation of first observation, and
source-product similarity. The trajectories provide molecular examples of
how repeated local oxygenation and acyl-state edits propagate through the
space.

An independent atom-neighborhood comparison reached the same conclusion at
population scale. Events with one changed source atom accounted for 88.7% of
G1, 89.2% of G2, and 89.4% of G3 derivations (Fig. 3B; Table S8B). Events
with more than three changed source atoms decreased from 9.9% at G1 to 8.3%
at G3. This locality metric does not assert that every projected edge is a
single enzymatic step; it establishes that accepted products generally differ
from their parents through spatially restricted molecular edits rather than
wholesale structural replacement.

Structure-derived functional-state transitions and elemental deltas supplied
an orthogonal chemical audit (Fig. 3C,D; Tables S7 and S8B). Free-hydroxyl
gain was the dominant transition, followed by reciprocal hydroxyl/ester
changes, carbonyl formation, and oxygenated acyl-state changes. Addition of
one oxygen atom was the leading elemental signature, accompanied by recurrent
carbon-hydrogen-oxygen gains and losses consistent with acyl transfer and
hydrolysis-like edits. Because these summaries were recomputed from source
and product structures, not copied from rule labels, their agreement with the
encoded grammar semantics supports the chemical coherence of the generated
space (Figs. S4 and S5).

Nominal grammar breadth nevertheless overstated effective generative breadth.
Only 14 of the 74 T1 productions generated accepted events across G0-G3. The
most active oxygenation production accounted for 60.5% of G1 events, and the
five most active productions accounted for 93.8%. The effective number of
rules increased only from 3.89 at G1 to 4.52 at G3 (Table S8A). The
2.6-million-structure space was therefore generated by a low-dimensional edit
vocabulary operating over a high-dimensional substrate population.

{{FIGURE_3}}

### Convergence transforms combinatorial expansion into a reticulate space with latent bridges

The derivation graph became strongly many-to-one with depth. A structure was
classified as convergent only when at least two distinct parents reached it
during its generation of first observation. By this definition, 170 of 15,801
G1 structures (1.1%), 184,570 of 223,823 G2 structures (82.5%), and
2,310,226 of 2,362,766 G3 structures (97.8%) were convergent (Fig. 4F;
Table S9). The sharp transition from a predominantly divergent first
generation to extensive later-generation convergence explains why accepted
event counts grew faster than unique-structure yield.

Returns to known space and immediate reverse cycles displayed distinct depth
profiles (Fig. 4G; Fig. S6). Connectivity-level G0 reconnections decreased
from 3.92% of accepted G1 events to 0.793% at G2 and 0.029% at G3. Immediate
reverse-cycle fractions increased from 0.865% to 1.32% and 1.66%,
respectively. The frontier therefore moved rapidly away from the observed
library even as local reversible motifs remained recurrent. Reverse cycles
were retained in the complete event audit but excluded from conservative
bridge inflation.

To identify generated structures with interpretable topological roles, we
constructed a directed acyclic view containing generation-increasing edges to
new structures and direct returns to G0. A latent bridge required a generated
intermediate to connect at least one ordered pair of distinct known taxanes.
This analysis identified 208 G1 and 19 G2 bridge candidates, for a total of
227 intermediates supporting 309 ordered G0 pairs (Fig. 4A-E; Tables S10 and
S11). Four molecule-resolved routes demonstrate how a generated intermediate
can reconnect known structures while retaining direction, rule provenance,
path multiplicity, and similarity evidence. The complete candidate gallery
and directed reconnection graph are provided in Figs. S7 and S8.

Bridge status is a prioritization criterion, not a biosynthetic claim. It
identifies structures occupying positions where a small number of
grammar-consistent edits reconnect observed chemistry. G3 bridge status is
right-censored because no fourth parent generation was enumerated. The 227
candidates are therefore most appropriately viewed as testable targets for
metabolite searching or reaction assays.

{{FIGURE_4}}

## Discussion

This study establishes a generative representation of taxane chemistry that
is distinct from both a reaction database and a molecular catalogue. The
353,524-row SMARTS release preserves the breadth, provenance, and
quality-control history of the source evidence. The 74-production T1 grammar
is the compact executable model used to generate chemical hypotheses.
Maintaining both objects resolves a common interpretive problem: release size
does not measure the number of independent chemical operations, and grammar
size does not capture the evidence breadth from which those operations were
derived.

The principal result is that a compact, high-evidence grammar opens a taxane
space more than three orders of magnitude larger than the known seed
collection within three applications. Yet the generated space is not an
unstructured combinatorial cloud. Approximately nine in ten accepted events
at every generation are localized to one changed source atom, oxygenation and
acyl-state changes dominate the structure-derived chemistry, and only a small
subset of productions controls most accepted events. Taxane chemical-space
expansion is therefore broad in molecular outcome but low-dimensional in its
operative edit vocabulary.

The evidence-stratified design is essential to this interpretation. T1, T2,
and T3 produce chemically distinct one-step neighborhoods, so pooling them
would obscure the relation between evidence quality and accessible space.
Moreover, the leakage-controlled pathway evaluation shows that external
reaction precedent alone did not recover the curated Taxol pairs under the
strict exclusion regime. Domain-informed replay restored substantial
coverage, but that recovery measures calibration fidelity rather than
independent discovery. Reporting these two results together makes the role of
family knowledge explicit and prevents circular validation.

The network topology adds a second layer of organization that cannot be
inferred from node counts. Extensive G2 and G3 convergence shows that many
different edit histories collapse onto the same molecular states. This
reticulation both moderates unique-space growth and supplies path redundancy.
At the same time, the declining frequency of returns to G0 indicates that the
exploratory frontier increasingly departs from known taxane chemistry. The
227 conservative bridge candidates occupy the informative intersection:
they are unobserved structures with traceable, short directional connections
between observed molecules.

Several limitations bound the biological interpretation. Principal-pair
projection omits explicit donor and cofactor molecules from executable
SMARTS, so donor availability, thermodynamics, compartmentation, enzyme
expression, and substrate specificity are not modeled. Product-level quality
control establishes structural consistency but not reaction kinetics.
Template application can produce incompletely specified new stereocenters,
although full InChIKey identity prevents their silent collapse with distinct
stereoisomers. Broad productions may overestimate accessible chemistry when
iterated, which is why G3 is explicitly exploratory. Finally, grammar distance
is not guaranteed to equal the number of enzymes or experimentally separable
steps.

These constraints also define productive next experiments. G1 and G2
structures can be prioritized using bridge position, path support, distance
to G0, local edit class, and analytical detectability. Candidate
intermediates can then be sought by targeted mass spectrometry or tested with
defined enzyme systems. More generally, the same provenance-resolved design
can be transferred to other natural-product families, provided that
family-specific calibration and external generalization are evaluated
separately.

## Conclusions

A provenance-resolved reaction grammar converts heterogeneous reaction
records into an auditable generative map of taxane chemistry. Seventy-four
high-evidence productions expanded 648 known taxanes into 2,603,038 unique
structures through 7,928,209 directional events. This space is simultaneously
large, locally edited, strongly convergent, and dominated by a small operative
rule vocabulary. By releasing complete molecular states, directional
derivations, rejection audits, and conservative bridge hypotheses, the study
provides a falsifiable chemical-space resource while avoiding claims that
grammar-accessible structures are observed metabolites or confirmed
biosynthetic intermediates.

{{METHODS}}

## Declarations

### Ethics approval and consent to participate

Not applicable to the computational analysis; AUTHOR_CONFIRMATION_REQUIRED.

### Consent for publication

AUTHOR_INPUT_REQUIRED.

### Availability of data and materials

The reproducible release contains the evidence-stratified rule summaries, the
primary executable grammar, complete G0-G3 molecular-state and directional
derivation tables, product rejection and application audits, benchmark
outputs, figure source tables, and analysis code. G0-G3 molecular states,
including SMILES, are supplied as generation-specific CSV files. Permanent
repository identifiers and accessions are AUTHOR_INPUT_REQUIRED.
Redistribution of source database records remains subject to the licenses and
access terms of the corresponding resources.

### Competing interests

AUTHOR_INPUT_REQUIRED.

### Funding

AUTHOR_INPUT_REQUIRED.

### Authors' contributions

AUTHOR_INPUT_REQUIRED.

### Acknowledgements

AUTHOR_INPUT_REQUIRED.

{{REFERENCES}}
