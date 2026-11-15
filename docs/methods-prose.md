> **Archived pre-audit draft, not release documentation.** Its cohort, endpoint and polymyxin claims are superseded by manuscript v10, [methods.md](methods.md), [results.md](results.md) and [final-model-release.md](final-model-release.md). Do not use the numbers below for the final 178-record model.

# Methods

## Overview

MuriVerse comprises a curated rodent in-vivo pharmacokinetic half-life dataset, a
SMILES-verified peptide-only training set, and a fixed-hyperparameter ensemble half-life
regressor built on a fork of PeptiVerse (Zhang et al., 2026).

The upstream PeptiVerse half-life task, which MuriVerse is defined against, is a
species-agnostic, route-agnostic, phase-agnostic regression over 375 combined records, of
which 130 are supplied as amino-acid sequence only and 245 as SMILES strings (Zhang et
al., 2026). Those are PeptiVerse's two molecular input representations. No
three-dimensional structure is used at any point, in PeptiVerse or in MuriVerse; the
archive of three-dimensional structures distributed with the source database was
evaluated and rejected, for the reasons given under Source database below.

MuriVerse narrows that label definition to **mouse and rat in-vivo elimination half-life
after a parenteral dose, measured by a named bioanalytical method on a plasma, serum or
whole-blood sample from a normal animal, for a single peptide chain of at most 50
residues whose SMILES has been independently decoded back to its residue sequence**.
Every record's dose route, assay method and SMILES completeness is established from the
primary source; records for which none of these can be established are not used.

That label definition produces a 320-record cohort. **The shipped model is fitted on the
213 intravenous records of it, not on all 320.** Only an intravenous dose has no
absorption phase, so an intravenous terminal half-life is an elimination half-life
throughout, whereas an extravascular one reports the absorption rate whenever absorption
is the slower step, the condition known as flip-flop kinetics. A route-pooled target
therefore mixes two physical quantities, and the mixing is measurable: every cross-route
transfer across the four species-by-route combinations is offset by 0.45 to 0.65 log₁₀ units
with the sign flip-flop kinetics predicts, while every cross-species transfer stays
inside ±0.21. The model is accordingly restricted by route and pooled over species, and
not the reverse. The 320-record cohort remains the input to development and is what
Figure 2 and the cohort sections below describe.

All reported quantities are derived from the curated cohort by a documented, repeatable
analysis pipeline; no data file is hand-edited at any stage. Two model variants are
retained: the route-restricted model, which is the shipped one, and a route-pooled model
kept for the cross-route transfer diagnostics described below, which cannot be performed
with a model that has never seen an extravascular record.

## Data provenance and cohort assembly

The data are organised as a single nested chain, each set a strict subset of the one
before it:

| Set | Records | Defined by |
|---|---|---|
| Source database | 4,500 | every shipped record of PEPlife 2.0 |
| Eligible set | 1,159 | mouse or rat, dosed in vivo, of peptide scale |
| Strict peptide set | 467 | the dosed molecule writable as a peptide SMILES and that SMILES verified |
| Analysis set | 320 | parenteral route, complete SMILES, assay named by a source |
| Intravenous set | 213 | the analysis-set rows dosed intravenously; the set the shipped model is fitted on |

Every count reported anywhere in this work refers to one of those five sets, and the
sections below define them in order. We draw the three middle stages of the chain as one
funnel, from the eligible set to the analysis set (Figure 2a).

### Source database

The initial record pool was PEPlife 2.0 (Chaudhary et al., 2024), retrieved through the
database's public interface, giving 4,500 records at one row per database entry. Mouse
and rat cohorts were defined by the database's own organism annotation: the mouse cohort
is the union of its Mouse and Mice organism queries, and the rat cohort the union of its
Rat and Rats queries. Where the archived organism response was unavailable, membership
was assigned by pattern-matching the rodent name in the free-text sample description;
cell-line matrices that name a rodent, such as 3T3, RAW 264, B16, Neuro-2a, Hepa, CT26,
4T1 and MC38, were excluded explicitly. The database's organism annotation and the
free-text match agreed exactly on every record where both could be verified.

Two additional half-life resources were incorporated, measured, and then removed. The
source database's archive of three-dimensional structures was rejected because its
identifier-to-structure mapping contradicts the sequence annotation on 171 of 394
checkable records, its coordinates do not encode stereochemistry, and 399 of its 1,837
files are experimental structures of parent proteins rather than of the dosed peptide. PepTherDia 2.0 and ThPDB
2.0 were rejected because both contain essentially no rodent half-lives: PepTherDia
none, and ThPDB five of 3,380, themselves duplicates of about two unique values.

### Full-text retrieval and curation

For every record's PubMed identifier, the corresponding full-text identifier was
resolved through NCBI's E-utilities and the full text retrieved; for sources outside the
open-access subset, the source database's own mirrored article was used instead. Each
identifier was resolved to title, journal, year, digital object identifier and
publication type. Retrieved full text was rendered as body-only plain text with the
abstract explicitly fenced, so that extraction cannot draw values from it, and with the
Methods, Results, Discussion and figure and table captions preserved. A per-record and
per-study provenance table was compiled for each species.

Study-level evidence used downstream — group size, dose, route, formulation, sample
matrix, bioanalytical method, half-life phase and numeric pharmacokinetic values — was
curated manually from the rendered full text into one structured record per source,
following a written extraction protocol. The protocol enforces body-only sourcing, with
abstracts fenced and inadmissible as evidence; a short verbatim quote of at most 25
words and a section tag for every non-null value; and no inference from outside
knowledge, with unreported values recorded explicitly as not reported. Each curated
observation is keyed by the source database's record identifier, so that arm-level
evidence is joined to the record table only on an explicit identifier and never on fuzzy
name or sequence matching. Forty-eight of the 184 source identifiers behind the strict
set carry a curated extraction.

### Label definition and phase filtering

The rules in this section define which reported half-life counts as the quantity
MuriVerse models. They are applied by both of the set builders below, which share one
implementation of them rather than each carrying its own.

The provenance table was narrowed to records satisfying four conditions: a live animal
was dosed; the readout matrix was plasma, serum, blood or whole animal; the source
database reports a non-null half-life; and the half-life phase is an elimination phase
rather than absorption, distribution or activity. The absorption-phase exclusion removes
the two papers that report a subcutaneous absorption half-life in addition to the
elimination half-life of the same animals, the latter being retained. The alpha or beta
subscript is read before the word "elimination", because several papers label both
compartments "elimination" and distinguish them only by subscript. A unit of inverse
hours is treated as a transcription error for hours: every such record is an absorption
half-life and is removed by the elimination-phase gate before any unit conversion is
applied. Duplicate database entries that transcribe the same measurement under differing
whitespace are collapsed to the first copy.

### Eligible set

We draw the whole path from the eligible set to the analysis set as one stacked funnel,
1,159 records down to 320, split into mouse and rat at every stage (Figure 2a). The
eligible set is that funnel's first bar.

The eligible set is built upstream of the figure, directly from the 4,500
source-database records, by gates the funnel does not draw: the sample text names a mouse
or rat, retaining 2,032 records; the record carries the database's in-vivo flag (1,465);
the half-life is positive, as taking its logarithm requires (1,442); a generous
pre-filter of at most 80 sequence letters (1,224); a live dosed animal rather than an
isolated matrix (1,189); not a cultured-cell matrix (1,186); and deduplication of the
database's own repeated transcriptions of a single measurement (1,159). The resulting
**1,159 eligible rodent records** are mouse or rat, dosed in vivo, and of peptide scale.
The 80-letter bound here is deliberately loose; the true 50-residue limit is applied on
parsed residues inside the drawn funnel, at the stage that removes 33 records, which is
why Figure 2a can summarise its first bar as at most 50 residues while a later bar still
enforces it.

The curated full-text verdicts and evidence columns are read from the two per-species
tables, `mouse/in_vivo_half_life.csv` and `rat/in_vivo_half_life.csv`, and transferred
onto the eligible table by record identifier, so the eligible table is an information
superset of them. Those tables identify the animal from the database's organism
annotation rather than from the sample text, so they supply curated columns for the
records they share and are not a stage of this chain.

The stages we draw, under the labels the panel prints, with the records surviving each
stage and the number it removes (Figure 2a):

| Stage | Records | Removed |
|---|---|---|
| Eligible records | 1,159 | — |
| Curated review pass | 1,152 | 7 |
| Normal animal and matrix | 1,103 | 49 |
| No PEG/Fc/albumin/polymer | 806 | 297 |
| Single peptide chain | 804 | 2 |
| Residue-level sequence | 669 | 135 |
| ≤ 50 residues | 636 | 33 |
| No other partner | 497 | 139 |
| Elimination half-life | 468 | 29 |
| Verified SMILES | 467 | 1 |
| Parenteral route | 383 | 84 |
| Whole molecule | 326 | 57 |
| Named assay | 320 | 6 |

The first nine removing stages are the strict structural gates, which take the eligible
set to the 467 kept records and remove 692 between them. The last three are the cohort
requirements described under The analysis set below, which are applied by a later
stage of the pipeline and remove a further 147. The two groups are drawn as one funnel
because leaving the second out would show 467 kept records becoming 320 analysed ones
with no step accounting for the gap. The panel closes with a final bar restating the
320-record analysis set.

Stages that remove nothing are bookkeeping rather than selection and are not drawn: the
builder records one row for retaining observations and another for curating the dose
route, and neither is a filter. Every one of the 692 records removed by a named
structural gate is recorded with the gate and the verbatim text that triggered it, and
the 147 removed by the three cohort requirements are recorded separately with the
identifiers of the records each one drops.

## The strict peptide set

The strict set is the subset of eligible rodent records for which the dosed molecule can
be written as a peptide SMILES and that SMILES independently verified. Of the 1,159
eligible rodent records, 467 were kept, comprising 111 mouse and 356 rat, under eleven
sequential gates, nine of which remove records and appear as bars in Figure 2a; every one
of the 692 excluded records is recorded with the gate that removed it and the verbatim
text that triggered it.

All eleven gates are listed below with the reasoning behind each; the nine that remove
records are the structural rows of the funnel tabulated above. In order, with records
removed in parentheses, they were: (1) not curated out (7);
(2) physiological circulating half-life, meaning the record describes a normal animal
with a plasma, serum or blood readout, excluding organ-ablated, knockout, transgenic,
hypertensive, diabetic and tumour-bearing animals and tissue, lysate and lymph samples
(49); (3) no engineered partner class, excluding polyethylene glycol, polymer, Fc,
albumin and depot or particle conjugates (297), because the clearance of such a molecule
is dominated by the partner rather than the pharmacophore and no single peptide
SMILES is the whole dosed molecule; (4) single peptide chain, excluding repeat and
polymer notation (2); (5) a faithful residue sequence that parses residue by residue,
excluding sequences reconstructed from prose or holding staple placeholders (135); (6) a
peptide-scale cap of 50 residues on the verified sequence (33), which defines the
dataset's scope; (7) no other attached partner (139); (8) elimination half-life (29);
(9) structure verified by independent decode (1); (10) observations retained (0); and
(11) dose route curated from the sources (0).

### Sequence verification

The source database ships two sequence annotations. The original marks a D residue by
lower-casing its letter and delimits a non-natural residue with hyphens; the derivative
is a flattened form in which case has been destroyed. We read the original, so that
per-residue chirality is representable. A hyphen-delimited sequence is parsed
conservatively: a chunk is accepted only when it is a known residue code, a single
letter, or an entirely upper-case run. Otherwise the record is refused, because a
four-letter chunk such as pIPhe is a single residue, not five residues each of whose
letters happens to be a valid single-letter code. A leading N-prefix is refused
outright, because splitting a name such as N-MeVal-EtPro on hyphens yields a bare N that
would be read as asparagine. Where the flattened annotation disagrees with an
upper-cased reading of the original, the record is refused.

The strict set admits a controlled vocabulary of non-standard residues whose structure
is a plain side-chain substitution: Aib, Orn, Nle, Nva, Abu, Dab, Dpr, homoarginine,
Cha, Nal and Tys. These residues are both built into SMILES and decoded back from
it. Forty-one kept records contain at least one such residue. The protein language model
used for sequence embedding has no token for any of them, so those positions enter the
sequence embedding as the unknown-residue token; which positions those are, and how many
there are, is retained on every row. Backbone-altering residues — sarcosine, any
N-methylated residue, hydroxyproline, azetidine carboxylic acid, octahydroindole
carboxylic acid, pyroglutamate and isovaline — are refused with their reason named.

### Structure building and verification

Structures were built for modified peptides from residue templates rather than from a
standard sequence-to-molecule converter, which accepts only the 20 standard L-amino
acids and always emits the L enantiomer. Residues, head-to-tail macrocyclization,
positioned disulfide bridges, per-residue chirality, and terminal modifications
including acetylation, amidation, methyl esterification and pyroglutamate are encoded
directly in the SMILES string rather than carried as separate features.

Every clause of the source database's free-text chemical, N-terminal and C-terminal
modification annotations was classified into one of three verdicts. Inert clauses add no
atoms, as for a legend, a chirality note, an isotopic substitution or a
sub-stoichiometric radiotracer. Drawable clauses are real fragments with a resolvable
attachment site, which are built and bonded on; examples include acyl lipidation,
γ-glutamyl and 8-amino-3,6-dioxaoctanoyl spacers, C16-diacid and S3-octanoyl ghrelin.
Blocking clauses are real but non-drawable modifications such as a polydisperse
kilodalton-scale polyethylene glycol, a glycan, a protein partner, a proprietary dye, a
drug payload, or a structure given only by reference to the paper. Any blocking clause
excludes the record from the strict set.

Verification runs in two independent steps on every kept record. First, the generated
SMILES is read back to a residue sequence using the same residue templates in
reverse. Residue identity, order, topology — head-to-tail macrocycle up to rotation, and
positioned disulfides — and configuration are all checked. For each alpha carbon the
stereochemical descriptor is read from the SMILES string and mapped to L or D; natural
L-series residues are S except cysteine, which is R, a fact established empirically from
the cheminformatics toolkit's output rather than assumed. A constitution-only check
would pass an all-L SMILES for a mixed-chirality peptide, which is precisely the error
this verification exists to exclude. Second, the full molecule's exact mass minus the
bare peptide's exact mass must equal the sum of the attached fragments' expected masses
to within 0.01 Da; this confirms that each conjugate went on once, at one site, losing
exactly the one hydrogen the bond costs. A record failing either check is excluded
rather than corrected.

A completeness tier is reported alongside each SMILES. Complete means the
constitution is fully represented, with chirality possibly unspecified. Partly resolved
means one of the peptide's own axes is not pinned down, as with an unpositioned
cyclization or an L/D series the source leaves open. Modification-not-drawn means a
modification clause names a real attached fragment that the builder cannot represent.
The strict set keeps all three tiers, but the analysis set requires complete
structures only.

Three records show a conflict between the source database's summary chirality annotation
and the stereochemistry encoded in sequence case. The model and Figure 2 use the
sequence-verified chirality; the database's reported value and the evidence used for
each row are both retained.

### Dose-route curation

The source database's dose-route field is blank on 453 of the 1,159 eligible rodent
records, or 39%, and on 238 of the 467 strict records, so reading it alone leaves the
route unreported for those rows even where the source states it plainly. We re-read the
primary sources — full text where reachable, patent text, the abstract, and NLM MeSH
indexing — and recorded the verdict for every source together with the verbatim quote
and its origin. That ledger is reconciled with each record to yield a route
classification and a four-level evidence grade.

On the 320 analysis-set records the evidence grades are: primary arm, where the source names
the route for that very arm, either a single-route study or a half-life matched to its
own entry in the source's pharmacokinetic table (62 records); primary source, where the
source names the routes it used, all parenteral, and the record is attributed at source
level (92 records); database route column, where the route is as the source database
transcribed it with no source re-read (166 records); and none, where no primary source
names a route for the arm (0 records, these being dropped by the parenteral-route
requirement imposed on the analysis set). Thirty of the 84 sources behind the
originally blank records could not be resolved, most because the full text is paywalled
and neither the abstract nor the MeSH indexing names a route; each is listed in the
ledger with its reason.

### Assay-method curation

The free-text assay annotation was bucketed into a closed vocabulary with an explicit
unknown level: radioimmunoassay, ELISA, mass spectrometry, high-resolution mass
spectrometry, AlphaLISA, Western blot, scintillation counting, and bioassay, the last
being a plasma concentration read off a cell-based reporter or an inhibition titre
rather than measured physically. For records whose assay text names a separation rather
than a detector, such as thin-layer chromatography; or the readout rather than the
method, such as a serum decay experiment; or nothing at all, the method was curated from
the primary source. Unknown is left standing where both the assay annotation and the
sources are silent, and the analysis set drops those records: what quantity decayed
is not a detail of a half-life, it is the half-life. Of the 33 originally unknown records,
27 were resolved from the sources or the assay text and 6 were removed.

### The analysis set

The 467 kept records were taken to the 320 analysis-set records by three further cohort
requirements, each tied to the mechanistic definition of the label. These are the last
three removing bars of the Figure 2a funnel, and the 320 records they leave are the analysis
set that panel's final bar reports. First, the dose route must be parenteral —
intravenous, subcutaneous, intraperitoneal or intramuscular — because an enteral or
mucosal half-life is absorption-limited and under flip-flop kinetics its terminal slope
is the absorption rate rather than elimination (84 records removed). Second, the SMILES must be complete, because a partly-resolved or
modification-not-drawn SMILES is not the molecule the label belongs to, and the
molecular encoder cannot distinguish the two (57 removed). Third, the assay method must
be named by a source (6 removed; these six had a median half-life of 23 h against the
cohort's 0.427 h and held three of the four longest labels in the set).

The resulting 320 rows comprise 74 mouse and 246 rat records, drawn from 129 source
identifiers, covering 166 distinct sequence-embedding inputs in 90 combined-species homology
clusters. Half-lives span 0.0028 to 62.0 h with a median of 0.427 h; sequences run from
3 to 50 residues with a median of 29; and chirality is 291 all-L, 24 mixed and 5 all-D.

Repeated treatment arms and independent reports of the same molecule are retained.
Candidate repeats sharing structure, species and label are flagged for evidence-based
review, but no record is removed on the basis of sequence and half-life alone. Half-life
values keep the source's qualifier: 309 exact, 5 approximate, 5 a range reported as the
median of its endpoints, and 1 an upper bound. Unusually long values are flagged for
audit, and that flag is not a model feature.

Half-life labels are expressed as log₁₀ hours, which is the quantity every reported
metric is computed on.

### The intravenous set

The analysis set is the input to model development, not the set the shipped model is
fitted on. One further restriction defines that set: the dose route must be intravenous.
This leaves the **213-record intravenous set**, comprising 179 rat and 34 mouse records,
and drops the 73 subcutaneous, 22 intraperitoneal and 12 intramuscular rows of the
analysis set.

The restriction is not a data-quality filter but a definition of the quantity being
modelled, and it is the one curation step driven by pharmacokinetics rather than by
evidence or structure. Only an intravenous dose has no absorption phase, so an
intravenous terminal half-life is an elimination half-life throughout; after an
extravascular dose the terminal phase reports the absorption rate instead whenever
absorption is the slower step. A route-pooled label is therefore a mixture of two
physical quantities. The measurements that establish this, and the alternative of
restricting by species instead, are reported in Results rather than assumed here.

The 213 records carry 68 homology clusters and 89 source identifiers, with half-lives from
0.0028 to 59.0 h and a median of 0.242 h — a shorter median than the analysis set's
0.427 h, which is the composition effect the restriction removes. Of them, 176 fall in
the training split and 37 in the validation split, both inherited from the cluster-
preserving split defined below rather than redrawn. Every constant the model is compared
against, and every quantity described as coming from the training split, is taken from
those 176 records; taking any of them from the 257-record route-pooled training split would let
the comparison borrow extravascular labels the model was never shown.

A route-pooled fit over all 320 analysis-set records is retained alongside the shipped one,
for the transfer diagnostic named in the Overview. It supplies no reported performance
number.

### Homology clustering and train/validation split

Sequences were clustered by profile-based all-against-all search (HMMER 3.3 or newer;
Eddy, 2011) over the 320 analysis-set sequences, with an edge retained when both the
per-domain independent E-value was at or below 10⁻³ and the aligned region covered at
least 80% of both sequences; clusters are the connected components of the resulting
graph. The 10⁻³ threshold reproduces the conventional homology cutoff, and the 0.8
coverage requirement reproduces the bidirectional coverage imposed by the clustering
previously used by PeptiVerse. Profile-based search was chosen over the faster
alternative because the cohort is small and its sequences short, and we observed the
profile method linking peptide families — short, divergent and often non-standard — that
the alternative missed. Clustering was performed on the combined mouse-and-rat set rather
than within species, so that a peptide family present in both animals cannot cross the
split.

An 80/20 train/validation split was then produced for consumers of the training-data
artifact. The split preserves whole clusters, so that no cluster spans both sides; it is
balanced within species; and it uses a fixed seed (1986) for reproducibility. The split is a
static artifact and is not used for any reported performance metric, as described under
Evaluation protocol below.

## Features

The feature vector is written out in six separately measured blocks.

**Sequence embedding (1,280 dimensions).** ESM-2, the 650-million-parameter variant (Lin
et al., 2023), frozen, with token embeddings mean-pooled over the residue sequence and
special and padding tokens excluded. Encoding is performed once over the 166 distinct
sequences of the 320-record analysis set and memoised per unique sequence; the
intravenous fit of 213 records draws its rows from that same encoding rather than from a
separate pass. Non-standard residues in the strict vocabulary enter as the
unknown-residue token. Nothing is truncated: the peptide-scale gate holds the cohort to
50 residues, well inside the model's maximum context of 1,022, and every record carries a
complete structure, so neither encoder has a coverage gap. The superseded route-pooled
tier, which still admitted Fc and albumin fusions and chains of protein length, did
contain three sequences past that context — an antibody heavy chain and two albumin
fusions — which were embedded from their first 1,022 residues only; the
engineered-partner and peptide-scale gates remove every one of them before the strict
cohort is formed.

**SMILES embedding (768 dimensions).** PeptideCLM (Feller et al., 2024), frozen,
mean-pooled over the verified peptide SMILES. On the strict cohort the string is the
whole dosed molecule for every record: all 320 are of complete structural coverage and
none carries an attached partner, because the engineered-partner and
other-attached-partner gates remove every record bearing a polyethylene glycol, an Fc, an
albumin, a polymer, a depot, a dye, a drug or a protein partner. Every record is
therefore encoded and none takes a missing SMILES block.

The convention that the string is the bare peptide, with a partner's identity and mass
carried by the terminal-class block instead, is inherited from the route-pooled tier,
where partnered records were still in the cohort and the division of labour was arrived
at by measurement. Drawing a partner at its true size swamped the peptide, since a
20 kDa polyethylene glycol is roughly 1,360 heavy atoms against a peptide of about 300,
and the model pools length-weighted windows, so the embedding became mostly generic
polyether. Dropping the partner but keeping the record instead collapsed polyethylene
glycol, Fc and bare peptide onto an identical SMILES. The strict tier settles the same
problem by exclusion rather than by representation, so the convention binds on none of
the records it ships.

**Molecule annotation (9 columns).** One-hot encodings of the engineering class, which is
unmodified, backbone-modified, label or imaging, lipidated, acetylated, amidated or
free; of the topology, which is linear, head-to-tail cyclic, disulfide-bridged or
bicyclic; and of the chirality, which is all-L, mixed or all-D. Each is read from the
source database's free text rather than from the half-life.

**Study context (33 columns).** One-hot encodings of the dose route (intravenous,
subcutaneous, intraperitoneal or intramuscular), the species (mouse or rat), the sample
matrix (plasma, serum, blood, tissue or unknown), the assay method as curated, the strain,
the sex, the formulation class and the dose unit, together with the reported dose on a
log scale. These are study-design variables, not molecular
properties; the model framing is explicitly what half-life a given molecule will show in
a given animal by a given route and assay, with all design variables available at
study-planning time.

The species encoding is stated separately because it has not always been present. An
earlier version of this block listed the strain but not the species, so the animal reached
the model only through the strain categories; on the large minority of records whose strain
is unknown, and on every molecule from outside the cohort — where a caller knows the animal
long before it knows the vendor strain — a mouse and a rat were the same input and the model
returned a single number for both. Two indicator columns now carry the animal directly. In a
single-species cohort they never vary and are dropped automatically, and under the
intravenous restriction they remain informative, the restricted cohort being 34 mouse and
179 rat records.

**Terminal classes and partner masses (13 columns).** One-hot encodings of the
N-terminal and C-terminal modification classes and of the heaviest-partner class, plus
two continuous masses for the polyethylene glycol and the partner. Only the terminal
classes vary on the strict cohort: at the N-terminus 297 records are free, 12 acetylated,
4 label or imaging, 3 lipidated, 3 other and 1 amidated, and at the C-terminus 218 are
free, 91 amidated, 6 other, 4 a linker or small conjugate and 1 acetylated. The partner
columns are constant, with no partner recorded and both masses zero on all 320 records,
because the engineered-partner gate excludes everything that would populate them; they
are retained so that the block's layout matches the pooled tier on which it was
defined.

**Physicochemical descriptors (16 columns).** Computed from the verified SMILES:
molecular weight, Crippen octanol-water partition coefficient, Crippen molar
refractivity, topological polar surface area, Labute approximate surface area, hydrogen
bond donor and acceptor counts, rotatable bonds, rings, aromatic rings, fraction of
sp³-hybridised carbon, heavy-atom count, amide bond count, counts of basic groups
(primary or secondary amines not adjacent to a carbonyl, plus guanidine) and acidic
groups (carboxylic acids), and the signed difference of the two as a net charge at pH
7.4. Ionisation state at pH 7.4 is inferred from substructure matches rather than from a
formal-charge annotation, which is zero for every record because the builder emits
neutral SMILES. This block exists because the annotation covariates are bucketed
class labels, so two different peptides of the same class receive an identical vector,
and because both learned representations score below the annotation block alone on the
headline grouping. The 16 interpretable descriptors are the quantities on which renal
filtration, protease exposure and plasma-protein binding depend; being computed from the
verified SMILES, they are defined for molecules with no faithful single-letter
sequence, which is the case Figure 4 turns on.

Missing values are left missing and no block is imputed. The linear and kernel members
use a median imputer fitted on training rows only, with the medians travelling with the
fitted object; the boosted-tree member receives the raw matrix and learns a default split
direction per feature. Standardisation ignores missing values during fitting and
preserves them in its output. Encoder outputs are memoised per unique sequence and per
unique SMILES, so a rebuild that does not change the encoder identity reads from cache
and completes in seconds.

## Model

The MuriVerse half-life regressor is an equal-weight mean of three members, with every
hyperparameter fixed in advance so that nothing is tuned on the records it is scored on:

* **Ridge regression**, with the penalty set to 10.0 times the larger of 1 and the
  feature width divided by 80, so that the per-column penalty is constant at the
  reference width of 80 and scales with dimensionality. It is preceded by a median
  imputer and standardisation.
* **Radial-basis-function support-vector regression**, with regularisation constant 1.0
  and tube width 0.2, using the same imputer and standardisation as ridge.
* **Gradient-boosted decision trees** (XGBoost), with a squared-error objective, learning
  rate 0.05, row subsampling 0.9, column subsampling 0.5, minimum child weight 5,
  maximum depth 2, and 150 rounds. This member sees the raw matrix with no imputation and
  routes missing values natively.

The three members fail differently: a heavily regularised linear model, a smooth kernel
model, and a shallow tree ensemble. On a few hundred records from about 100 publications,
a single boosted-tree model overfits the publication it was trained on; averaging the
three was the most stable combination under held-out family and publication
cross-validation. No model-selection step runs on the data: all three hyperparameter sets
are fixed constants.

The shipped model's feature set is molecule annotation, study context, terminal classes
and physicochemical descriptors, giving 71 columns on the strict cohort. Neither learned
embedding is used by the shipped model.

This is **not** the configuration with the highest Spearman correlation under the
headline grouping. On the intravenous cohort that is physicochemical descriptors plus
study context plus terminal classes, at ρ = 0.587 against 0.574 for the shipped set. The
shipped configuration is the best-scoring one that retains the molecule-annotation
block, which the analysis requires. Applicability is a hard constraint on a model that will be
applied to molecules the cohort does not contain, and the annotation block is the only one
definable for any molecule that can be described in words. Several physicochemical columns
cannot be computed at all for a peptide with no faithful single-letter sequence — six of
the ten residues of a polymyxin are L-2,4-diaminobutyric acid, which has no single-letter
code — and without the annotation block nothing holds such a molecule inside the trained
range: the unconstrained winner returns polymyxin half-lives of about 10⁻¹²⁰ h against an
observed 1 to 2 h. The constraint costs 0.013 in ρ, well inside the seed-to-seed spread.
The ablation still scores and reports every configuration including the unconstrained
winner, so its price is visible rather than hidden.

Under the route restriction, the route columns of the molecule-annotation and
study-context blocks are constant and therefore carry no information; standardisation
leaves a zero-variance column at zero. The column count is unchanged; the information
content of four of those columns is not.

## Evaluation protocol

The saved 80/20 split is not used for any reported number; it is a training-data artifact
only.

**Row restriction.** The cohort is restricted to intravenous records before anything
derives folds, groups or strata from the metadata, so that every downstream object
describes the restricted cohort and none can be left describing the full one. A
route-restricted run writes to its own namespace, leaving the route-pooled artifacts in
place for comparison. The figure pipeline selects between the two and asserts that the
out-of-fold rows and the cohort rows are the same set, so that a mismatch fails loudly
rather than producing panels describing two different cohorts.

**Grouped cross-validation.** Every record receives an out-of-fold prediction under
five-fold grouped cross-validation, with the fold assignment reshuffled for each of three
random seeds; each seed is a complete out-of-fold pass. Every reported metric is the mean
of the three per-seed scores. We plot the scatter and curves from a single seed for
legibility (Figure 3).

**Grouping keys.** Three groupings are scored: held-out homology clusters; held-out
source identifiers, using patent identifiers where no article identifier exists; and
held-out connected groups sharing neither a homology cluster nor a source identifier.
The third is the strict transfer estimate and the headline grouping. It is built by a
union-find over the row index of each key, returning the connected components, so that a
record whose cluster is linked to any row sharing its source — and so on transitively —
ends up in one component. On the 320-record strict cohort this produces 75 connected
groups; on the 213-record intravenous cohort the shipped model is fitted on, 55.

**Metrics.** Pearson r, Spearman ρ, root-mean-square error, mean absolute error and
coefficient of determination are evaluated on log₁₀ half-life. The area under the
receiver operating characteristic curve against the pre-stated 6-hour screen is retained
as a secondary ranking summary rather than a Figure 3 performance panel, because its
dichotomisation does not measure the accuracy of a continuous half-life prediction. Each
metric is accompanied by a 95% group bootstrap over the headline groups, using 500
resamples, resampling whole groups with replacement, under a fixed seed for
reproducibility. In every resample each seed's out-of-fold predictions are scored
separately and the three scores averaged, so the point estimate and the interval
summarise the same statistic.

**Accuracy on the absolute scale.** We report the cumulative share of held-out records
falling within a given fold of the observed half-life, over both species and by species,
together with each stratum's median absolute fold error and signed bias (Figure 3d). This
replaces the decile-binned calibration curve of earlier builds, which pooled over routes whose offsets
point in opposite directions, so that a well-centred aggregate curve could be assembled
from routes that were each systematically wrong. With the cohort restricted to one route
that failure mode is gone, and the panel splits by species instead.

A leave-group-out residual-coverage panel is also constructed but is not in the layout.
Where it is used, the 80% or 95% radius for each connected group is the corresponding
quantile of absolute out-of-fold residuals in the other groups, coverage is averaged
equally over groups, and its 95% interval is a group bootstrap. These empirical bands
describe uncertainty in this retrospective cohort; they are not prospective
pharmacokinetic prediction intervals, nor evidence that any dosing regimen is safe.

**Error stratification.** We report the median and interquartile range of the
single-seed absolute fold error, defined as ten raised to the absolute difference between
log₁₀ observed and log₁₀ predicted, by species, bioanalytical method, source contribution
size, observed-half-life tertile, peptide-length tertile and engineering class
(Figure 3f).
Administration route is no longer a stratum: every record in the shipped cohort is
intravenous, so the section would hold a single level carrying all 213 records. Source
contribution is the count of cohort records sharing a source identifier, grouped as the
largest source, other sources with at least four records, sources with two or three
records, and singleton sources. Strata holding fewer than five records are dropped rather
than drawn, because a handful of records has no stable interquartile range and invites
the eye to weigh it against strata a hundred times its size; under the route restriction
this removes the bioassay stratum, with one record, and the lipidated class, with three.
These are descriptive post-hoc slices of one shared held-out model; the analysis does not
fit or claim performance for separate subgroup models.

**Pre-stated targets.** The model is scored against Pearson r above 0.70 and Spearman ρ
above 0.50, with an area under the curve above 0.80 at the 6-hour threshold retained as a
secondary screening analysis rather than a primary target for a regression model. On the
intravenous cohort, ρ = 0.574 (met), area under the curve 0.818 (met, having been missed
at 0.785 on the route-pooled cohort) and r = 0.568 (not met). The classification metric
rests on only 17 positives of 213, so it is coarse.

**Block ablation.** The ablation prices every named configuration — embedding-only,
annotation-only, context-only, every pair, and the compact 71-column shipped set — under
every grouping, so that which blocks ship is a measurement rather than an assertion.
Twenty-nine configurations are scored on the intravenous cohort. The ablation is also
what prices the applicability constraint above, since it reports the unconstrained winner
alongside the shipped set.

**Nested selection.** The optimism of choosing the feature configuration on the same
folds that score it was estimated directly: at every outer fold, the configuration is
picked by three-fold inner grouped cross-validation on the training rows, and the chosen
configuration is scored on the held-out outer fold. The candidate set is every
non-diagnostic configuration of at most 100 columns; the embedding configurations are
left out of the inner selection but kept in the ablation, since every embedding
configuration scores below the compact set under the headline grouping and searching
2,000-column configurations inside each fold would cost hours. On the intravenous cohort
this gives r = 0.456 and ρ = 0.460 against a headline ρ of 0.574, so roughly 0.11 of the
headline is selection optimism — a larger gap than the route-pooled cohort showed, and
one that should be quoted alongside the headline rather than instead of it.

**Refit for application.** After measurement, the shipped configuration and the two
structural configurations — the physicochemical-descriptor-best and the SMILES-embedding-best —
are refit on the scored records with the same fixed hyperparameters and saved. The shipped
models are therefore refit on the 213 intravenous records. These refit models are used only
for the held-out polymyxin application and for the route and species transfer
diagnostics; no reported metric is computed on them.

## Held-out application: colistin and polymyxin B

The shipped and structural-embedding refit models were applied to four congeners —
polymyxin B1 and B2, and colistin A and B — whose structures were taken from PubChem
isomeric SMILES records. This is a different provenance from the reconstructions used
in training, and is labelled as such. Neither drug appears in the strict cohort as a
named record; the only mention in the eligible set, searched across its curated free-text
columns, is colistin as a comparator in
the free-text dosing notes of the ZY4 record, whose half-life is ZY4's own. No peptide in
the training set has detectable homology to any of the congeners.

Each congener was predicted under two input-encoding variants crossing the N-terminal
acyl wording, giving 4 shipped-model predictions per drug and species and 16 in total. The
dose route is not varied: under the intravenous shipped model the route columns are
constant in training, so a subcutaneous input is not a different point in a learned
dimension but an extrapolation along an axis of zero training variance, and the number
returned would be an artefact of the scaler. Figure 4c measures what such a read would be
worth.

The C-terminus is annotated as unmodified only. An amidated variant was run in earlier
builds on the reasoning that the macrolactam bond is itself an amide; that reasoning was
incorrect and the variant has been removed. A C-terminal amide is –C(=O)NH₂, whereas the
ring is Thr10–C(=O)–NH–Dab4, an isopeptide bond to the Dab4 side-chain amine. Both are
amides, but only one is a C-terminal amide, and colistin and polymyxin B carry the ring,
not the terminal amide. Because the isomeric SMILES already encodes the macrocycle and is
identical under either setting, the variant changed only the one-hot terminal class,
asserting a modification the molecules do not have; it moved the prediction 1.7-fold,
further than the drug, the congener or the species, and so dominated every contrast in the
analysis. The cyclisation is carried by the topology block, and free is used in its sense
of no terminal modification recorded.

The remaining encoding variants test sensitivity to how a molecule is written down and are
not independent measurements: the N-terminal wording changes nothing at all, its two rows
being bit-identical. The 16 predictions therefore carry 8 distinct values across the four
congeners, and we draw one mark per molecule, being the median of its two encoding variants
with a bar spanning the congeners (Figure 4b). Species is a model input, so the mouse and
rat predictions for a given molecule and encoding differ; strain and sex are set to
unknown, and study dose is left missing, which the median imputer fills for the ridge and
support-vector members while the boosted-tree member routes it natively. Six of the ten residues of each congener are
L-2,4-diaminobutyric acid, which has no single-letter code; this no longer affects the
shipped model's prediction, since neither saved configuration reads the sequence block,
but it does affect the homology diagnostic, which is why four surrogate spellings are
retained there and reported as a spread.

### Comparison with observed half-lives

Observed half-lives were taken from fifteen published studies. Five of them —
Landersdorfer 2018, Jiao 2023, Yapa 2013, and the two Lin 2017 papers — were added from a
Europe PMC search for population pharmacokinetics or humanized dosing of colistin or
polymyxin B in mice or rats; the remaining ten are a convenience sample curated during
the project. No search string, date, hit count, screening log, protocol or second
reviewer was recorded, and the full text of Lin 2018 could not be retrieved. The
selection is therefore not a systematic review, and the sample may be biased toward
studies that report a half-life or a population model.

Values are plasma or serum terminal half-lives. The colistimethate prodrug arms of He
2013 and Yapa 2013, the liposomal and tissue values of He 2013, Jiao 2023's 35.8 h value
from a moribund-mouse study in which clearance had collapsed, and the kidney
concentrations of Wang 2023 are all excluded. For five studies reporting population
pharmacokinetic parameters rather than a half-life — van der Meijden 2023, Landersdorfer
2018, both Lin 2017 papers and Yapa 2013 — we derived a half-life as the slow eigenvalue
of the published compartmental model and labelled it as derived. The central-volume
half-life from the same models is reported for comparison but kept out of the main
comparison. For van der Meijden 2023 the half-life was derived at each of eight dose
levels from a one-compartment model with first-order absorption and dose-dependent
clearance. The derivation reproduces Landersdorfer 2018's own published value of 4.1 h.

Mean residence times, namely 1.58 h from Li 2003 and 4.07 h from Wang 2023, are reported
in Figure 4d as asterisks but are not compared with the model and not averaged into the
observed combined value, because a mean residence time is not a half-life. For Wang 2023
the mean residence time is drawn alongside that study's 2.53 h terminal half-life.

Per-study methodology was extracted from the published full text of all fifteen studies,
covering the half-life type; the analysis method, which is noncompartmental, one
compartment, population pharmacokinetic, derived slow-phase, or central-volume; the
animal state; the matrix basis of unbound against total drug; and the sampling window in
hours. The sampling window for Matar 2020 is read from that paper's own Figure 10 at
about 6 h, because its text does not state the sampling schedule; every other row carries
a value read directly from the methods text. Figure 4d renders these fields as the two
text columns beside each row.

Observed combined values per drug and species were computed as the animal-weighted
geometric mean across dose groups, with the animal-weighted standard deviation on the
log₁₀ scale across dose groups. Components of a single dose group were averaged first, so
that animals are counted once. Sensitivity versions were computed with every study
counted once, with reported values only, with intravenous doses only, and with each study
left out in turn. The animal-weighted choice was declared in advance and allows a study
with many dose groups to dominate: van der Meijden 2023 carries 67.5% of the polymyxin B
mouse weight. The same weighting gives a value per source of observation. Sources are
grouped into Jiao 2023; the Tam lab, being the University of Houston group of V. H. Tam,
comprising Abdelraouf 2012, He 2013 and Manchandani 2016; the Monash group, being the
Monash University group of R. L. Nation, J. Li and colleagues, comprising Dudhani 2010,
Cheah 2015, Landersdorfer 2018, Sivanesan 2017, Yapa 2013 and both Lin 2017 papers, with
Li 2003 by the same investigators from the University of South Australia; and other.

### Baselines

The model is compared with three constants read from the training split: the median
half-life of the same species, of both species combined, and of the lipidated training
peptides.

**The constants are read from the same records the applied model was fitted on**, which
under the intravenous shipped model is the intravenous training split of 176 records: mouse
0.485 h, rat 0.194 h, both species 0.218 h, and the three lipidated intravenous training
peptides 1.75 h. A route-pooled baseline would let the constant borrow the extravascular
labels — longer by the 0.45 to 0.65 log₁₀ flip-flop offset — and so hold labels the model
itself was not permitted to use. In mouse that choice alone moves the constant from
0.485 h to 0.894 h, onto the observed polymyxin values, which is enough to decide whether
the constant or the model is the better predictor. The analysis verifies that the
reported baseline split size and medians equal those of the intravenous training split.

The SMILES-embedding variant predicts 0.11 to 1.17 h for these congeners and the
top-scoring compact configuration that drops the molecule-annotation block (physicochemical
descriptors plus study context plus terminal classes) predicts 0.5 to 1.7 h — physically
plausible but weaker than the shipped set on held-out ranking, which is why neither is
drawn. A previous build of this analysis reported both configurations returning polymyxin
half-lives of 10⁻¹¹⁶ to 10⁻¹³⁷ h; that was a block-order mismatch between training (which
concatenated blocks in configuration-list order) and the inference assembler (which
reordered them by layout), applied to any configuration whose own list disagreed with the
layout. The shipped configuration happens to list its blocks in layout order and was
unaffected.

### Why the transfer diagnostic uses the route-pooled model

The route and species transfer diagnostic applies the **route-pooled** models, against
the shipped default. Its purpose is to test externally whether a model mis-transfers
across route. A model fitted on intravenous records alone cannot answer that: it has no
extravascular level to predict under, so an extravascular prediction from it is an
extrapolation along an axis of zero variance, and the route effect it appeared to report
would be an artefact of the scaler. The instrument for the diagnostic has to be the model
being diagnosed. This is the same relation that Figure 3's panels a and b bear to its
panels c onward, and it is why the transfer results are unaffected by the route
restriction.

### Dose sensitivity

A dose-versus-prediction sweep is run over both animals, since species is a model input,
and only ratios to the dose-missing prediction are reported. Over 5 to 40 mg/kg the
shipped model's prediction moves by 0.6%, a log-log slope of about zero, against 2 to
4.7-fold observed in the three dose-ranging studies, a log-log slope of 0.26 to 0.74. **Only one of the 213
intravenous training records carries a stated dose**, and three of the 320 do, so the
model's dose response is an extrapolation from a single observation rather than a learned
effect. The dose panel of earlier builds is no longer drawn, because a response resting
on one training record should not be given a panel; the sweep is still reported in
tabular form.

## Figures

All main-text figures and their captions were generated from the cohort by a figure
pipeline that separates data preparation from drawing. Supplementary Figure 1 draws the
full preprocessing pipeline that feeds Figure 2 — the species composition of the three
source databases (PEPlife2, ThPDB2 and PepTherDia2), the build funnel from the 4,500
PEPlife2 records down to the 467 structurally verified rodent in-vivo records, and the
three reconciliation gates (route, structure, assay) that take that set to the 320
modelling rows.

In Figure 2's funnel, panel a, each species's shade lightens toward the eligible end and
deepens toward the analysis set, so that the colour deepens in step with the cohort
shrinking; every other Figure 2 panel draws the two species in those same solid tones.
Supplementary Figure 1 reuses the Figure 2 mouse-blue and rat-green pair in the one
panel split by species (mouse / rat across stages); the gate panels draw kept-by-gate in
maroon and dropped-by-gate in grey.

The dimensional reduction in Figure 2h uses uniform manifold approximation and projection
with cosine distance, 15 neighbours, a minimum distance of 0.3 and a fixed seed (1986). The
smoothed Gaussian-kernel-weighted mean of per-sequence median log₁₀ half-life is drawn in
light grey, with surface opacity rising with height, and the sequence points are coloured
by species, with grey for sequences dosed in both animals.

## Data availability

The curated cohort, the curation ledgers, the filtering counts, the exclusion list, the
held-out predictions behind every reported score, the figure source tables and the refit
production models are all available in the project repository. Both the shipped
intravenous model and the route-pooled model retained for the transfer diagnostics are
included. Pre-computed embedding datasets are regenerated on build and are not
distributed.

## Statistical reporting and limitations

All reported metrics are means over three random seeds. Intervals are 95% group
bootstraps with 500 resamples over the 55 connected cluster-and-source groups of the
intravenous cohort the shipped model is fitted on; the route-pooled 320-record cohort has
75 such groups. Pre-stated targets are reported as pass or fail. No multiple-testing
correction is applied: the three metrics describe three facets of one prediction rather
than three independent hypothesis tests, and the ablation prices every configuration
separately and is not interpreted as a hypothesis test.

The headline evaluation is limited by between-source measurement disagreement. Taking
every sequence-and-species combination that more than one source reports, reducing each
source to its mean log₁₀ half-life for that combination, and forming every unordered pair
of sources, the intravenous cohort holds 33 such pairs, agreeing at r = 0.75 with a
median disagreement of 0.35 log units against a label spread of 0.82. Two laboratories
measuring the same peptide in the same animal differ by about a factor of two, which
bounds what any predictor can score against a single held-out reported label. An earlier
version of this analysis quoted 66 pairs agreeing at r = 0.57, from a pairing rule that
could not be reproduced; the rule above is the one now executed and verified. About three
quarters of the residual variance of the held-out predictions is a constant offset per
source; if each source-level offset were known, Pearson r would be 0.93.

Selection optimism is material: nested selection of the feature configuration inside each
training fold gives ρ = 0.460 against the headline 0.574, so roughly 0.11 of the headline
is the configuration having been chosen on these data.

The route restriction improves every metric on a like-for-like comparison, raising ρ from
0.514 to 0.574 on the same 213 records with the same shipped configuration, but the gain in
ρ alone is not statistically resolvable: a paired bootstrap over the 55 connected groups
gives +0.063 with an interval of −0.005 to +0.157. The restriction is justified by the
pharmacokinetic rationale of flip-flop kinetics and by the 0.45 to 0.65 log₁₀ cross-route
bias, not by that correlation.

The dose, strain and formulation covariates are absent or unparsed on most records — only
one of the 213 intravenous training records carries a dose — and the retrieval ledger
marks 85 of the 120 article identifiers in the cohort as not retrieved, while nine patent
identifiers account for a further 91 records. Species is reached only through strain in
the shipped feature set, and the polymyxin application shows a true species gap of
0.43 to 0.61 log₁₀ that the model predicts as approximately zero, so Figure 3's small
cross-species bias should be read as evidence that the model does not distinguish the
species, not that the species do not differ. The polymyxin held-out application is a
small test of four molecules from outside the training chemistry, not a general estimate
of performance on novel chemistry.
