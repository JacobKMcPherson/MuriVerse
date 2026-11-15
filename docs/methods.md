# Methods

## Final post-audit model

The final distributed model is trained on **178 source-verified intravenous records**
(150 rat, 28 mouse; 76 sources). The 213-record model is an archived **pre-audit**
development benchmark, not the release default. The endpoint is source-reported
circulating pharmacokinetic half-life after intravenous dosing, not uniformly
verified terminal elimination half-life. See [the release record](final-model-release.md).

The fixed 71-feature configuration has grouped cross-validated Spearman rho 0.445
(95% CI 0.153-0.563); nested feature-configuration selection gives 0.297
(0.096-0.450). Both use five folds, three seeds and 50 connected cluster-source groups.
Table 1 reports these post-audit estimates. Figure 3 and Table 2 document historical
pre-audit development; Figure 4 and Tables 5-7 use the post-audit models.

## Overview

MuriVerse is a curated rodent in-vivo pharmacokinetic (PK) half-life dataset,
a structure-verified peptide-only training set, and a fixed-hyperparameter
ensemble half-life regressor built on top of a fork of PeptiVerse
(Zhang et al., 2026). PeptiVerse's half-life task is a species-agnostic,
route-agnostic, phase-agnostic regression on 375 combined records (130
sequence-only, 245 SMILES; Zhang et al., 2026). MuriVerse narrows that
label definition to **mouse and rat source-reported circulating pharmacokinetic half-life after a
parenteral dose, measured by a named bioanalytical method on a plasma,
serum or whole-blood sample from a normal animal, for a single peptide
chain of at most 50 residues whose SMILES has been independently decoded
back to its residue sequence**. Every record's dose route, assay method
and structure tier is established from the primary source; records in
which none of these can be established are not used.

That label definition produces a 320-record development cohort. **The final model
is fitted on 178 post-audit intravenous records; the original 213-record fit is pre-audit.**
An intravenous dose avoids an absorption phase, but does not establish that a reported
half-life is a terminal elimination estimate. An extravascular terminal slope
one reports the absorption rate whenever absorption is the slower step
(flip-flop kinetics). A route-pooled target therefore mixes two physical
quantities, and the mixing is measurable: every cross-route transfer in
the four species x route combinations is offset by 0.45-0.65 log10 with the sign
flip-flop predicts, while every cross-species transfer stays inside
+/-0.21 (`evaluate_species_route_cells.py`, `table36`). The model is
accordingly restricted by route and pooled over species, and not the
reverse. The 320-record cohort remains the input to development and is
what Figure 2 and the cohort sections below describe. All pipeline
outputs — the cohort tables, the model, every figure and every table —
are reproducible from scripts in the repository; nothing in the data
directory is hand-edited.

The code and the committed data files are available at
<https://github.com/JacobKMcPherson/MuriVerse>; the commands that
rebuild the historical development analyses downstream of the curated cohort are
`build_strict_peptide_set.py`, `build_halflife_training_set.py --tier
strict`, `evaluate_halflife_model.py --cohort rodent --tier strict
--route iv --nested` and `figures/code/run_all.R`. Dropping `--route iv`
reproduces the route-pooled model, which is kept on disk for comparison
and is what `evaluate_polymyxin_species_route.py` deliberately applies.
`python MuriVerse/audit_figure_numbers.py` re-derives every verifiable
number in the Figure 2, 3 and 4 captions from the data and exits non-zero
if any has drifted.

## Data provenance and cohort assembly

The data form a single nested chain, each set a strict subset of the one
before it:

| Set | Records | Defined by |
|---|---|---|
| Source database | 4,500 | every shipped record of PEPlife 2.0 |
| Eligible set | 1,159 | mouse or rat, dosed in vivo, of peptide scale |
| Strict peptide set | 467 | dosed molecule writable as a peptide SMILES, and that SMILES verified |
| Analysis set | 320 | parenteral route, complete SMILES, assay named by a source |
| Intravenous set | 213 | the analysis-set rows dosed intravenously; the set the shipped model is fitted on |

Every count reported anywhere in this work refers to one of those five
sets. We draw the three middle stages as one funnel, from the eligible set
to the analysis set (Figure 2a).

### Source database

The initial record pool was PEPlife 2.0 (Chaudhary et al., 2024), pulled
via the database's REST API by `build_peplife2.py` to produce
`data/peplife2/peplife2_records.csv` (4,500 records; one row per PEPlife2
entry). Mouse and rat cohorts were defined by PEPlife2's own organism
annotation (`fetch_peplife2_api.py --species mouse|rat`): the mouse
cohort is the union of the `Mouse` and `Mice` organism queries and the
rat cohort is the union of the `Rat` and `Rats` queries. Where the
archived organism response was unavailable, membership was assigned by
pattern-matching the rodent name in the `Test Sample` free-text column;
cell-line matrices that name a rodent were excluded by an explicit
cell-line regex (e.g. 3T3, RAW 264, B16, Neuro-2a, Hepa, CT26, 4T1,
MC38). PEPlife2's own `Organism` annotation and the free-text match
agreed exactly on every record where both could be verified.

Two additional half-life resources were wired into the pipeline, measured,
and removed: PEPlife2's PDB structure archive (because
its ID-to-structure mapping contradicts the sequence column on 171 of 394
checkable records, its coordinates do not encode stereochemistry, and 399
of its 1,837 files are experimental PDB entries of parent proteins
rather than the dosed peptide) and PepTherDia 2.0 and ThPDB 2.0 (because
both contain essentially no rodent half-lives: PepTherDia zero, ThPDB
five of 3,380, duplicates of ~2 unique values). The measurements that
decided their removal are in `docs/data-provenance-quality.md`.

### Full-text retrieval and curation

For every record's PMID, `build_rodent_cohort.py` resolved the PMC
full-text identifier via NCBI's E-utilities and retrieved the BioC
full text; `fetch_peplife2_papers.py` additionally mirrored PEPlife2's
own PDF per publication for sources outside the PMC open-access subset.
`fetch_pubmed_metadata.py` resolved each PMID to title, journal, year,
DOI, PMCID and publication type. Retrieved full text was rendered as
body-only plain text (`render_fulltext.py` for PMC BioC;
`render_pdf_text.py` for publisher PDFs) with the abstract explicitly
fenced so that extraction cannot source values from it, and the body
preserved from the Methods, Results, Discussion and figure and table
captions. A per-record and per-study provenance table was written by
`build_study_dossier.py` as
`<species>/records_provenance.csv` and `<species>/studies.csv`.

Study-level evidence used downstream — group size, dose, route,
formulation, sample matrix, bioanalytical method, half-life phase, and
numeric PK values — was curated manually from the rendered full text
into one JSON file per source following `EXTRACTION_PROTOCOL.md`. The
protocol enforces body-only sourcing (abstracts are fenced and may not
be quoted as evidence), a short verbatim quote (≤25 words) and section
tag for every non-null value, and no inference from outside knowledge:
unreported values are recorded as `"not reported"`. Each curated
observation is keyed by the PEPlife2 record identifier so that arm-level
evidence (`pk_observations.csv`) is joined to the record table only on
an explicit identifier, not on fuzzy name or sequence matching. Forty-eight of
the 184 source identifiers behind the strict set carry an extraction JSON.

### Label definition and phase filtering

The rules in this section define which reported half-life counts as the
quantity MuriVerse models. Both set builders below apply them through one
shared implementation rather than each carrying its own.

`subset_in_vivo_half_life.py` narrowed `records_provenance.csv` to
records satisfying four conditions: (i) a live animal was dosed
(`animal_dosed`); (ii) the readout matrix was plasma, serum, blood or
whole animal; (iii) PEPlife2 reports a non-null half-life
(`half_life_h`); (iv) the half-life phase is an elimination phase, not
absorption, distribution or activity. The absorption-phase exclusion
removes the two papers that report a subcutaneous absorption half-life
in addition to the elimination half-life of the same animals; the latter
is retained. The alpha/beta subscript is read before the word
"elimination" because several papers label both compartments
"elimination" and distinguish them only by subscript. Hours-1 is treated
as a transcription typo for hours (every "Hours-1" row is an absorption
half-life and is removed by the elimination gate before any unit
conversion). Duplicate PEPlife2 entries that transcribe the same
measurement under different whitespace are collapsed to the first copy.

### Eligible set

`build_rodent_eligible_set.py` applies the rules above directly to the
4,500 source-database records and writes
`rodent/eligible_in_vivo_half_life.csv` with the per-stage counts in
`rodent/eligible_filtering_counts.csv`. Its gates, which Figure 2a does
not draw, are: the sample text names a mouse or rat (2,032 records
retained); the record carries PEPlife2's in-vivo flag (1,465); the
half-life is positive, as taking its logarithm requires (1,442); a
generous pre-filter of at most 80 sequence letters (1,224); a live dosed
animal rather than an isolated matrix (1,189); not a cultured-cell matrix
(1,186); and deduplication of PEPlife2's own repeated transcriptions of a
single measurement (1,159). The resulting **1,159 eligible rodent
records** are mouse or rat, dosed in vivo, and of peptide scale, and are
the first bar of Figure 2a. The 80-letter bound is deliberately loose;
the true 50-residue limit is applied to parsed residues inside the drawn
funnel, at the stage that removes 33 records.

The curated full-text verdicts and evidence columns are read from the
per-species tables `mouse/in_vivo_half_life.csv` and
`rat/in_vivo_half_life.csv` and transferred onto the eligible table by
record identifier, so the eligible table is an information superset of
them. Those tables identify the animal from PEPlife2's organism
annotation rather than from the sample text, so they supply curated
columns for the records they share and are not a stage of this chain.

## The strict peptide set

The strict set, produced by `build_strict_peptide_set.py`, is the subset
of eligible rodent records for which the dosed molecule can be written
as a peptide and the SMILES independently verified. Of 1,159 eligible
rodent records, 467 were kept (111 mouse, 356 rat) under eleven
sequential gates recorded in `rodent/strict/filtering_counts.csv`;
every one of the 692 excluded records is listed in `exclusions.csv`
with the gate that removed it and the verbatim text that triggered it.

Gates were (in order; records removed in parentheses): (1) not curated
out (7); (2) physiological circulating half-life — the record describes
a normal animal with a plasma, serum or blood readout, excluding
organ-ablated, knockout, transgenic, hypertensive, diabetic and
tumour-bearing animals and tissue, lysate and lymph samples (49); (3) no
engineered partner class — PEG, polymer, Fc, albumin or depot/particle
(297), because the clearance of such a molecule is dominated by the
partner rather than the pharmacophore and no single peptide SMILES is
the whole dosed molecule; (4) single peptide chain, excluding repeat
and polymer notation such as `(VPGAG)120` (2); (5) a faithful residue
sequence that parses residue by residue, excluding sequences
reconstructed from prose or holding staple placeholders (135); (6)
peptide-scale cap at 50 residues on the verified sequence (33), which
defines the dataset's scope; (7) no other attached partner (139); (8)
elimination half-life (29); (9) structure verified by independent decode
(1); (10) observations retained (0); (11) dose route curated from the
sources (0).

### Sequence verification

PEPlife2 ships two sequence columns. `Sequence` is the original, which
marks a D residue by lower-casing its letter (`rPKPfQwFwLL`) and
delimits a non-natural residue with hyphens (`H-Aib-QGTFTSDYSKYL`).
`sequence` is a flattened derivative in which case has been destroyed.
We read `Sequence`, so per-residue chirality is representable. A
hyphen-delimited sequence is parsed conservatively: a chunk is accepted
only when it is a known residue code, a single letter, or an entirely
upper-case run; otherwise the record is refused, because a four-letter
chunk such as `pIPhe` is a single residue, not five residues each of
whose letter happens to be a valid single-letter code. An `N-` prefix is
refused outright, because splitting `N-MeVal-EtPro` on hyphens yields a
bare `N` that would be read as asparagine. Where the flattened column
disagrees with an upper-cased reading of the original, the record is
refused.

The strict set admits a controlled vocabulary of non-standard residues
whose structure is a plain side-chain substitution: Aib, Orn, Nle, Nva,
Abu, Dab, Dpr, hArg, Cha, Nal and Tys. These residues are both built
into SMILES and decoded back from it. Forty-one kept records contain
at least one such residue. ESM-2 has no token for any of them, so those
positions enter the sequence embedding as `X`;
`nonstandard_positions` and `n_nonstandard` are retained on every row.
Backbone-altering residues (Sar, N-methyl-anything, Hyp, Aze, Oic, pGlu,
Iva) are refused with their reason named.

### Structure building and verification

`peptide_smiles.py` builds a SMILES for a modified peptide from residue
templates (not `Chem.MolFromSequence`, which accepts only the 20
standard L-amino acids and always emits the L enantiomer). Residues,
head-to-tail macrocyclization, positioned disulfide bridges, chirality
(per residue) and terminal modifications (acetylation, amidation,
methylester, pyroglutamate) are encoded directly in the SMILES string
rather than carried as separate features. `peptide_mods.py` classifies
every clause of PEPlife2's free-text `Chem-MOD`, `N-ter MOD` and
`C-ter MOD` columns into one of three verdicts — `inert` (no atoms
added, e.g. a legend, chirality note, isotopic substitution or
sub-stoichiometric radiotracer), `drawable` (a real fragment with a
resolvable attachment site, built and bonded on; e.g. acyl lipidation,
γ-glutamyl and 8-amino-3,6-dioxaoctanoyl spacers, C16-diacid, S3
octanoyl ghrelin), or `blocking` (a real but non-drawable modification
such as a polydisperse kDa-scale PEG, a glycan, a protein partner, a
proprietary dye, a drug payload, or a structure given only as "see
paper"). Any blocking clause excludes the record from the strict set.

Verification runs in two independent steps on every kept record. First,
`reconstruct_residues` reads the generated SMILES back to a residue
sequence using the same residue templates in reverse. Residue identity,
order, topology (head-to-tail macrocycle up to rotation; positioned
disulfides) and configuration are checked: for each alpha carbon the CIP
descriptor is read from the string and mapped to L or D (natural L-series
residues are S except cysteine, which is R, established empirically from
RDKit output rather than assumed). A constitution-only check would pass
an all-L SMILES for a mixed-chirality peptide, which is precisely the
error this verification exists to exclude. Second, the full molecule's
exact mass minus the bare peptide's exact mass must equal the sum of the
attached fragments' expected masses to 0.01 Da; this checks that each
conjugate went on once, at one site, losing exactly the one hydrogen the
bond costs. A record that fails either check is excluded rather than
corrected.

`structure_tier` is reported alongside the SMILES: `complete` (the
constitution is fully represented, with chirality possibly unspecified);
`partly resolved` (one of the peptide's own axes — e.g. an
unpositioned cyclization or an L/D series the source leaves open — is
not pinned down); `modification not drawn` (a `Chem-MOD` clause names a
real attached fragment that the builder cannot represent). The strict
set keeps all three tiers, but the modelling table requires
`structure_tier == complete`.

Three records show a conflict between PEPlife2's summary `Chirality`
column and the stereochemistry encoded in sequence case; the model and
Figure 2 use the sequence-verified chirality, and PEPlife2's reported
value is retained in `chirality_reported`. `stereo_source` records the
evidence used for each row.

### Dose-route curation

PEPlife2's `Dose Route` field is blank on 453 of the 1,159 eligible
rodent records (39%), and on 238 of the 467 strict records, so reading
it alone leaves `route_class == "not reported"` for those rows even
where the source states the route plainly. We re-read the primary
sources — full text where it is reachable, patent text, the abstract
and NLM MeSH indexing — and
recorded the verdict for every source in
`data/curation/route_curation.csv` with the verbatim quote and where it
came from. `halflife_condition_curation.py` reconciles that ledger with
each record and writes `route_class` plus a four-level
`route_evidence_level`. On the 320 modelling records the counts are:
`primary_arm` (the source names the route for this very arm — a
single-route study, or a half-life matched to its own entry in the
source's PK table; 62 records); `primary_source` (the source names the
routes it used, all parenteral; the record is attributed at source
level; 92); `database_route_column` (the route PEPlife2 transcribed,
with no source re-read; 166); `none` (no primary source names a route
for the arm; 0 — these records are dropped by the parenteral-route
requirement imposed on the modelling table, see "The modelling table"
below). Thirty of the 84 sources behind the originally-blank records
could not be resolved, most because the full text is paywalled and
neither the abstract nor the MeSH indexing names a route; each is
listed in the ledger with its reason.

### Assay-method curation

`study_conditions.py` buckets the free-text `Assay` column into a closed
vocabulary with an explicit `unknown` level: radioimmunoassay, ELISA,
mass spectrometry, high-resolution mass spectrometry, AlphaLISA, Western
blot, scintillation counting, and bioassay (plasma concentration read
off a cell-based reporter or an inhibition titre rather than measured
physically). For records whose `Assay` text names a separation rather
than a detector (e.g. "TLC"), the readout rather than the method (e.g.
"Serum decay experiment"), or nothing at all, the method was curated
from the primary source into `data/curation/assay_curation.csv`.
`unknown` is left standing where the `Assay` column and the sources are
both silent, and the modelling table drops those records: what quantity
decayed is not a detail of a half-life, it IS the half-life. Of the
33 originally-unknown records, 27 were resolved from the sources or the
`Assay` text and 6 were removed.

### The modelling table

`build_halflife_training_set.py --tier strict` takes the 467 kept records
to the 320 modelling records (74 mouse, 246 rat; 129 source identifiers; 166
distinct ESM-2 input sequences in 90 combined-species HMMER homology clusters;
half-life 0.0028–62.0 h, median 0.427 h; 3–50 residues, median 29; 291
all-L, 24 mixed-chirality, 5 all-D) by three further cohort requirements,
each tied to the mechanistic definition of the label: (i) the dose route
must be parenteral — iv, sc, ip or im — because an enteral or mucosal
half-life is absorption-limited and under flip-flop kinetics its terminal
slope is the absorption rate rather than elimination (84 records
removed); (ii) `structure_tier == complete`, because a partly-resolved or
modification-not-drawn SMILES is not the molecule the label belongs to
and the molecular encoder cannot distinguish the two (57); (iii) the
assay method is named by a source (6; these six had a median half-life
of 23 h against the cohort's 0.427 h and held three of the four longest
labels in the set).

Repeated treatment arms and independent reports of the same molecule are
retained; candidate repeats with the same structure, species and label
are flagged in `repeat_observation_candidates.csv` for evidence-based
review, but no record is removed on the basis of sequence and half-life
alone. Half-life values keep the source's qualifier: 309 exact,
5 approximate (e.g. `~4`), 5 range (median of endpoints) and 1 upper
bound (`<5`); the qualifier is recorded in `half_life_qualifier`.
Unusually long values are flagged in `half_life_extreme_qc` for audit;
the QC flag is not a model feature.

Half-life labels are expressed as log10(hours), which is the quantity
every reported metric is computed on.

### Homology clustering and train/validation split

Sequences were clustered by `phmmer` (HMMER 3.3+; Eddy, 2011) all-against-all
over the 320 analysis-set sequences, with an edge retained when both the
per-domain i-Evalue was at or below 1e-3 and the aligned region covered
at least 80% of both sequences; clusters are the connected components
of the resulting graph. The 1e-3 i-Evalue reproduces the conventional
homology cutoff, and 0.8 coverage reproduces the bidirectional coverage
that `mmseqs cluster -c 0.8 --cov-mode 0` (the clustering previously
used by PeptiVerse) imposes. HMMER was chosen over MMseqs2 because the
cohort is small and short, and we observed phmmer linking peptide
families — short, divergent and often non-standard — that
MMseqs2 missed. Clustering was performed on the combined mouse-and-rat
set rather than within species, so that a peptide family present in
both animals cannot cross the split.

An 80/20 train/validation split was then produced for consumers of the
training-data artifact. The split preserves whole clusters (no cluster
spans both splits), is balanced within species, and uses a fixed seed
(1986) for reproducibility. The split is a static artifact; it defines no
cross-validation fold and no cross-validated score (the evaluation scores all 213
intravenous records out of fold, including the 37 in the validation split; see
Section "Evaluation protocol"). It supplies only the baseline constants of the
polymyxin application (from its post-audit records) and the label-shift panel of Figure 3i.

## Features

The feature vector is written out in six separately-measured blocks
whose column counts, computations and layout are recorded in
`feature_names.json`:

**Sequence embedding (`sequence`; 1,280 dimensions).** ESM-2 (650M-parameter
`facebook/esm2_t33_650M_UR50D`; Lin et al., 2023), frozen; mean-pooled token
embeddings over the residue sequence excluding special tokens and
padding. Encoding runs over the 166 distinct sequences of the 320-record
analysis set (the intravenous-only fit of 213 records is a row subset of
it, not a separate encoding pass), and the vectors are memoised per
unique sequence. Non-standard residues in the strict vocabulary enter as
`X`. No sequence is truncated: the peptide-scale gate caps the cohort at
50 residues (`MAX_RESIDUES` in `build_strict_peptide_set.py`), so the
longest input is 50 residues against ESM-2's 1,022-residue context, and
every record carries a `complete` SMILES, so there is no coverage gap on
either encoder. The superseded route-pooled tier, which admitted Fc and albumin fusions and chains of protein length, did contain longer sequences (the embedding cache still holds four, of 1,130 to 3,611 residues, embedded from their first 1,022 residues); the engineered-partner and peptide-scale gates remove all of them before the strict cohort is formed. The eligible set's longest sequence is 80 residues.

**SMILES embedding (`smiles`; 768 dimensions).** PeptideCLM
(`aaronfeller/PeptideCLM-23M-all`; Feller et al., 2024), frozen;
mean-pooled over the verified peptide SMILES, with the SMILES-SPE
tokenizer vocabulary at `tokenizer/`. On the strict cohort the string is
the whole dosed molecule: `smiles_coverage` is `complete` and
`partner_class` is `none` on all 320 records, because the
engineered-partner and other-attached-partner gates remove every record
carrying a PEG, Fc, albumin, polymer, depot, dye, drug or protein
partner (436 records between them; see `filtering_counts.csv`). Every
record is therefore encoded, and no row takes the NaN SMILES block.

The convention that the string is the bare peptide, with a partner's
identity and mass carried by `termini` instead, is inherited from the
route-pooled tier, where partnered records were still in the cohort and
the division of labour was arrived at by measurement: drawing a partner
at its true size swamped the peptide (a 20 kDa PEG is ~1,360 heavy atoms
against a ~300-atom peptide, and PeptideCLM pools length-weighted
windows, so the embedding became mostly generic polyether), whereas
dropping the partner but keeping the record collapsed PEG, Fc and bare
peptide onto an identical SMILES. The strict tier resolves the same
problem by exclusion rather than by representation, so the convention
binds on no record it actually ships.

**Molecule annotation (`molecule_annotation`; 9 columns).** One-hot
`hl_engineering` (4 columns: unmodified / none stated, backbone modification
only, lipidated, radio/fluorescent label only), `topology` (2: cyclic, linear)
and `chirality` (3: d, l, mix), each bucketed from PEPlife2 free text rather
than from the half-life. A `smiles_coverage` indicator is built but is
`complete` on every strict record, so `_drop_constant` removes it; terminal
modifications are in `termini`, not here.

**Study context (`study_context`; 33 columns).** One-hot `route_model` (4:
im, ip, iv, sc), `species` (2: mouse, rat), `sample_matrix` (4: blood, plasma,
serum, unknown), `assay_method` (6: ELISA / immunoassay, HPLC, LC-MS, bioassay,
radioimmunoassay, radiometric / imaging), `strain` (8, including unknown),
`sex` (3), `formulation_class` (3: depot, nondepot, not reported) and
`pk_dose_unit` (2: mg/kg, not reported), plus the one continuous column
`log10_dose_mg_kg`. Age and body-weight columns are built but empty on this
cohort and are dropped. On the 213 intravenous records the three formulation
columns are constant (`not reported` on all of them) and the dose is reported
on one record, so the dose column is constant too. These are study-design
variables, not molecular properties; the model framing is explicitly "what
half-life will this molecule show in this animal by this route and assay",
with all design variables available at study-planning time.

`species` is the block's load-bearing column for that framing and is
stated separately because it has not always been present. Until the
species fix the block listed `strain` but not `species`, so the animal
reached the model only through the strain one-hots; on the large
minority of records with `strain=unknown`, and on every molecule from
outside the cohort (where a caller knows the animal but not the vendor
strain), mouse and rat were the same input vector and the model returned
one number for both. Two one-hot columns, `species=mouse` and
`species=rat`, now carry it directly. Under a single-species cohort they
are constant and are dropped automatically by `_drop_constant`; under
the intravenous restriction they remain the only study-design axis with
variance besides matrix, assay, strain and sex.

**Terminal classes (`termini`; 13 columns).** One-hot N-terminal class (6:
acetylation, amidation, free, label / imaging, lipidation, other / unparsed)
and C-terminal class (5: acetylation, amidation, free, linker / small
conjugate, other / unparsed), plus the binary flags `has_n_mod` and
`has_c_mod`. On the 320 strict records the N-terminus is 297 free, 12
acetylated, 4 label or imaging, 3 lipidated, 3 other and 1 amidated, and the
C-terminus 218 free, 91 amidated, 6 other, 4 linker or small conjugate and 1
acetylated. The block is DESIGNED to carry a heaviest-partner class and two
log-scaled masses (`peg_kda`, `partner_kda`), but `partner_class` is `none` and
both masses are 0.0 on all 320 records, because the engineered-partner gate
excludes everything that would populate them, so `_drop_constant` removes them
and they are **not** among the 13 columns.

**Physicochemical descriptors (`physchem`; 16 columns).** Computed from
the verified SMILES with RDKit (version in `requirements.txt`):
molecular weight, Crippen clogP, Crippen molar refractivity, topological
polar surface area, Labute approximate surface area, number of hydrogen
bond donors and acceptors, rotatable bonds, rings, aromatic rings,
fraction of sp3 carbon, heavy-atom count, amide bond count, counts of
basic (SMARTS: primary/secondary amine not adjacent to a carbonyl, plus
guanidine) and acidic (SMARTS: carboxylic acid) ionizable groups, and
their difference as a signed `net_charge_ph74`. Ionisation state at
pH 7.4 is inferred from substructure matches rather than from
`formal_charge`, which is zero for every record (the builder emits
neutral SMILES). The block exists because annotation covariates are
bucketed class labels, so two different peptides of the same class
receive the identical vector, and the two learned representations (ESM-2
and PeptideCLM) score below the annotation block alone on the headline
grouping. The 16 interpretable descriptors are the quantities on which
renal filtration, protease exposure and plasma-protein binding depend;
being computed from the verified SMILES means they are defined for
molecules with no faithful single-letter sequence, which is the case
Figure 4 turns on.

Missing values are left as NaN; no block is imputed. Scikit-learn
members use a median imputer (fit on training rows only; the medians
travel with the fitted object), and XGBoost receives the raw matrix and
learns a default split direction per feature. `StandardScaler` ignores
NaN during fitting and preserves it in output.

Encoder outputs are memoised per unique sequence (`esm2_embedding_cache.npz`)
and per unique SMILES (`peptideclm_embedding_cache.npz`), so a rebuild
that does not change the encoder identity reads from cache and takes
seconds.

Features are written as HuggingFace `datasets.DatasetDict` saved to disk
at `../results/half_life_muriverse/`
under two names: `halflife_rodent_strict_with_embeddings` (the shipped
block set) and `halflife_rodent_strict_ablation_with_embeddings` (every
block, used by `evaluate_halflife_model.py` so that any block can be
ablated).

## Model

The MuriVerse half-life regressor (`halflife_model.HalfLifeEnsemble`)
is an equal-weight mean of three members, with every hyperparameter
fixed in advance so that nothing is tuned on the records it is scored
on:

* **Ridge regression** (`sklearn.linear_model.Ridge`) with
  `alpha = 10.0 × max(1, d/80)` for a feature matrix of width `d`, so the
  per-column penalty is constant at the reference width of 80 and scales
  with dimensionality. Preceded by a `SimpleImputer(strategy="median")`
  and a `StandardScaler`.
* **RBF-kernel support-vector regression** (`sklearn.svm.SVR`) with
  `C = 1.0` and `epsilon = 0.2`, same imputer and scaler as ridge.
* **XGBoost** (`xgb.train`; version in `requirements.txt`) with the
  `reg:squarederror` objective, `eta = 0.05`, `subsample = 0.9`,
  `colsample_bytree = 0.5`, `min_child_weight = 5`, `max_depth = 2` and
  150 rounds. XGBoost sees the raw matrix (no imputation) and routes
  missing values natively.

The three members fail differently: a heavily regularised linear model,
a smooth kernel model, and a shallow tree ensemble. On a few hundred
records from ~100 publications, a single boosted-tree model overfits
the publication it was trained on; averaging the three was the most
stable combination under held-out family and publication cross-
validation. No model-selection step runs on the data: all three
hyperparameter sets are literal constants in `halflife_model.py`.

The shipped model's feature set is `annotation + context + termini +
physchem` (71 columns on the strict cohort). No ESM-2 or PeptideCLM
embedding is used by the shipped model.

It is **not** the configuration with the highest Spearman ρ under the
cluster+PMID grouping. On the intravenous cohort that is `physchem +
study context + termini` at ρ = 0.587, against 0.574 for the shipped
set. The shipped configuration is the best-scoring one that satisfies
`REQUIRED_SHIPPED_BLOCKS` in `evaluate_halflife_model.py`, which demands
the `molecule_annotation` and `study_context` blocks. Study context is
required because it is the only block that carries the animal; molecule
annotation is required so that a prediction always rests on the engineering
class, topology and chirality a curator recorded rather than on bulk
physicochemistry alone. Both configurations return physically plausible
half-lives for the polymyxins (post-audit refits: shipped 0.57–0.92 h,
unconstrained winner 0.45–0.79 h over the 16 rows; pre-audit refits 0.98–1.47 and 0.48–0.86 h), so the constraint rests on definability, not on
a numerical pathology. An earlier build reported the unconstrained winner
returning about 10⁻¹²⁰ h; that was a block-order mismatch between training
and inference (see the comment above `REQUIRED_SHIPPED_BLOCKS`) and is fixed.
The constraint costs 0.013 in ρ, well inside the seed-to-seed standard
deviation of 0.048 (per-seed ρ 0.533, 0.561, 0.627), and the ablation still
scores and reports every configuration including the unconstrained winner, so
its price is visible rather than hidden.

Both learned encoders are published pretrained checkpoints used frozen and not
fine-tuned here: ESM-2 650M (`facebook/esm2_t33_650M_UR50D`, revision
`08e4846e537177426273712802403f7ba8261b6c`) and PeptideCLM-23M-all
(`aaronfeller/PeptideCLM-23M-all`, revision
`a0847d8231d236645a2c4f629590118716c6fdda`). The weight-file SHA-256 values, the
embedding-cache hashes and a recomputation of cached vectors with the pinned
weights are in
`results/revision_analysis/snapshot/embedding_provenance.json`
(`analysis/embedding_provenance.py`).

Under the route restriction 11 of the 71 columns are constant on the 213
records and carry no information, and `StandardScaler` leaves a zero-variance
column at zero: the four route columns, the three formulation columns, the dose
column (one reported value), `strain=C57BL/6J`, and the N-terminal lipidation
and other / unparsed classes. The species columns are not constant — the
intravenous cohort is 34 mouse and 179 rat. The column count is unchanged; the
effective width is 60. The per-column listing, with the saved scaler and
imputer values checked against the matrix, is
`results/revision_analysis/snapshot/feature_manifest.csv`
(`analysis/feature_manifest.py`).

## Evaluation protocol

`evaluate_halflife_model.py` is where the model is measured. The saved
80/20 split defines no cross-validation fold and no performance score: all 213
intravenous records are held out once per seed, the 176 training-split and the
37 validation-split records alike (`analysis/split_reconciliation.py` replays
the folds and reproduces the saved out-of-fold predictions to 2e-16). It is a
training-data artifact, used only for the label-shift panel of Figure 3i and for
the baseline constants of the polymyxin application, which are medians of the
146 post-audit intravenous rows of the training split (27 mouse, 119 rat; the pre-audit
analysis used the 176 training-split rows).

**Row restriction.** `--route iv` restricts the cohort to intravenous
records before anything derives folds, groups or strata from the
metadata, so every downstream object describes the restricted cohort and
none can be left describing the full one. A route-restricted run writes
to its own namespace — `oof_predictions_iv.csv`, `table31..._iv.csv`,
`table31b..._iv.csv` and `rodent/strict/model_iv/` — leaving the
route-pooled artifacts in place for comparison. The figure pipeline
selects between them with `HL_ROUTE` in
`figures/code/05_data_halflife_model.R` (default `"iv"`), and
`05_data_halflife_model.R` asserts that the out-of-fold rows and the
cohort rows are the same set, so a mismatch fails loudly rather than
producing panels that describe two different cohorts.

**Grouped cross-validation.** Every record receives an out-of-fold
prediction under five-fold `GroupKFold` with the fold assignment
reshuffled for each of three random seeds (0, 1, 2); each seed is a
complete out-of-fold pass. Every reported metric is the mean of the
three per-seed scores; the plotted scatter and curves in Figure 3 are
drawn from seed 0 only, for legibility.

**Grouping keys.** Three groupings are scored:

* `cluster` — held-out HMMER homology clusters;
* `PMID` — held-out source identifiers (patent identifiers where no
  PMID exists);
* `cluster_PMID` — held-out connected groups sharing neither a homology
  cluster nor a source identifier. This is the strict transfer estimate
  and the headline grouping.

The `cluster_PMID` key is built by `union_groups`, which runs a
union-find over the row index of each key and returns the connected
components, so a record whose cluster is linked to any row sharing its
source — and so on transitively — ends up in one component. On the
320-record strict cohort this produces 75 connected groups; on the
213-record intravenous cohort the shipped model is fitted on, 55.

**Metrics.** Pearson r, Spearman ρ, RMSE, MAE and R² are evaluated on
log10 half-life. AUC against the pre-stated ≥ 6 h screen is retained as
a secondary ranking summary in the result tables, not as a Figure 3
performance panel; its dichotomisation does not measure the accuracy of
a continuous half-life prediction. Each metric is accompanied by a 95%
group bootstrap over the `cluster_PMID` groups (500 resamples; whole
groups resampled with replacement; identical seed for reproducibility).
In every resample each seed's out-of-fold predictions are scored
separately and the three scores averaged, so the point estimate and the
interval summarise the same statistic (`bootstrap_ci` in
`evaluate_halflife_model.py`; `analysis/metrics_iv.py` reproduces it from
the saved predictions). Predictions are never averaged across seeds before
scoring.

**Accuracy on the absolute scale.** We report the cumulative share of
held-out records falling within a given fold of the observed half-life,
over both species and by species, together with each stratum's median absolute fold
error and signed bias (Figure 3d). This replaces the decile-binned calibration curve of
earlier builds, which pooled over routes whose offsets point in opposite
directions, so that a well-centred aggregate curve could be assembled from
routes that were each systematically wrong. With the cohort restricted to one
route that failure mode is gone, and the panel splits by species instead.
A leave-group-out residual-coverage panel is still constructed by `fig3.R`
but is not in the layout; where it is used, the 80% or 95% radius for each
connected cluster/source group is the corresponding quantile of absolute
out-of-fold residuals in the other groups, coverage is averaged equally over
groups, and its 95% interval is a group bootstrap. These empirical bands
describe uncertainty in this retrospective cohort; they are not prospective
pharmacokinetic prediction intervals or evidence that a dosing regimen is
safe.

**Error stratification.** We report the median and interquartile
range of the seed-0 absolute fold error, `10^|log10(observed) −
log10(predicted)|`, by species, bioanalytical method, source contribution
size, observed-half-life tertile, peptide-length tertile and engineering
class (Figure 3f). Administration route is no longer a stratum: every record in the
shipped cohort is intravenous, so the section would hold one level carrying
all 213 records. Source contribution is the count of cohort records with the
same source identifier, grouped as the largest source, other sources with at
least four records, sources with two or three records, and singleton sources.
Strata holding fewer than five records are dropped rather than drawn, because
a handful of records has no stable interquartile range and invites the eye to
weigh it against strata a hundred times its size; under the route restriction
this removes `bioassay` (n = 1) and the lipidated class (n = 3). These are
descriptive post-hoc slices of one shared held-out model; the analysis does
not fit or claim performance for separate subgroup models.

**Pre-stated targets.** The model is scored against Pearson r > 0.70
and Spearman ρ > 0.50, with AUC(≥ 6 h) > 0.80 retained as a secondary
screening analysis rather than a primary target for a regression model.
On the intravenous cohort ρ = 0.574 (met), AUC = 0.818 (met, having been
missed at 0.785 on the route-pooled cohort) and r = 0.568 (not met). The
AUC rests on only 17 positives of 213, so it is coarse. See the Figure 3
caption for the continuous-outcome results.

**Block ablation.** The ablation prices every named configuration
(embedding-only, annotation-only, context-only, every pair, and a
compact 71-column `annotation + context + termini + physchem` set) under
every grouping, so that which blocks ship is a measurement rather than
an assertion. 29 configurations are scored on the intravenous cohort.
Results are written to
`figures/tabs/table31_muriverse_halflife_feature_ablation_strict_iv.csv`
(route-restricted) or `..._strict.csv` (route-pooled). The ablation is
also what prices the applicability constraint above: it reports the
unconstrained winner alongside the shipped set.

**Nested selection.** `--nested` additionally estimates the optimism of
choosing the feature configuration on the same folds that score it: at
every outer fold, the configuration is picked by three-fold inner
grouped cross-validation on the training rows and the chosen
configuration is scored on the held-out outer fold. The candidate set is
every non-diagnostic configuration of at most 100 columns (the embedding
configurations are left out of the inner selection but kept in the
ablation; every embedding configuration scores below the compact set
in the ablation under cluster_PMID, and searching 2,000-column
configurations inside each fold would cost hours). The result is
written to `table31b_muriverse_halflife_nested_selection_iv.csv` and is
the historical pre-audit nested-selection estimate. On
the 213-record pre-audit cohort it gives r = 0.456 and ρ = 0.460 (95% CI 0.267–0.548)
against the pre-audit fixed-configuration ρ of 0.574, so roughly 0.11 of that estimate is selection
optimism — a larger gap than the route-pooled cohort showed, and one that
should be quoted alongside the headline rather than instead of it. It nests
only the choice of feature configuration; the homology clustering, the
connected-group construction, the route restriction and the retention of the
molecule-annotation block were decided with the full cohort in view, so
residual optimism may remain.

**Refit for application.** After measurement, the shipped and the two
structural configurations (physicochemical-descriptor–best and
embedding-best) are refit on the scored records with the same fixed
hyperparameters and saved as `.joblib` files to `rodent/strict/model_iv/`
for a route-restricted run, or `rodent/strict/model/` for a route-pooled
one. These historical models are refit on the 213 intravenous-labelled records;
these are the pre-audit models, retained for reproducibility, and they include 37
records the primary-source audit later found discordant. After the audit the shipped
configuration and the two structural configurations were refit with the same fixed
hyperparameters on the 178 source-verified records and saved to
`rodent/strict/model_iv_postaudit/` (`shipped`, `structural`, `structural_embedding`,
`training_membership.csv`; post-audit, distributed model;
`MuriVerse/revision_analysis/refit_postaudit_model.py`).
`predict_polymyxins.py` loads the post-audit models by default (`--model-suffix _iv` selects
the pre-audit models in `model_iv/`), and its baseline constants are read from the post-audit
records of the training split. The pre-audit route-pooled model is used by
`evaluate_polymyxin_species_route.py` (which loads
`model/`, deliberately — see below); no reported metric is computed on
them.

## Held-out application: colistin and polymyxin B

`predict_polymyxins.py` applies the shipped and structural-embedding
refit models to four congeners (polymyxin B1 and B2, colistin A and B)
whose structures were taken from PubChem's isomeric SMILES, which is a
different provenance from the RDKit reconstructions used in training and
is labelled as such. Neither drug appears in the strict cohort as a
named record (`table24_polymyxin_leakage_check.csv`); the only mention
in the eligible set is colistin as a comparator in the free-text
dosing notes of the ZY4 record, whose half-life is ZY4's. No peptide in
the training set has detectable homology to any of the congeners
(`table25_polymyxin_surrogate_homology.csv`).

Each congener was predicted under two input-encoding variants that
cross the N-terminal acyl wording, giving 4 shipped-model predictions per
drug and species and 16 in total. The
dose route is no longer varied: under the intravenous shipped model the
route columns are constant in training, so a subcutaneous input is not a
different point in a learned dimension but an extrapolation along an axis
with zero training variance, and the number returned would be an artefact
of the scaler. Figure 4c measures what such a read would be worth.

The C-terminus is annotated as unmodified only. An `Amidation` variant was
run in earlier builds on the reasoning that the macrolactam bond is itself
an amide; that reasoning was incorrect and the variant has been removed. A
C-terminal amide is –C(=O)NH₂, whereas the ring is Thr10–C(=O)–NH–Dab4, an
isopeptide bond to the Dab4 side-chain amine. Both are amides, but only one
is a *C-terminal* amide, and colistin and polymyxin B carry the ring, not
the terminal amide. Because the isomeric SMILES already encodes the
macrocycle and is identical under either setting, the variant changed only
the one-hot terminal class, asserting a modification the molecules do not
have; it moved the prediction 1.7-fold, further than the drug, the congener
or the species, and so dominated every contrast in the analysis. The
cyclisation is carried by the `topology` block, and `Free` is used in its
sense of no terminal modification recorded.

The remaining encoding variants test sensitivity to how a molecule is
written down and are not independent measurements: the N-terminal wording
changes nothing at all, its two rows being bit-identical, so the
16 predictions carry 8 distinct values across the four congeners and we draw
one mark per molecule, the median of its two encoding variants with a bar
spanning the congeners (Figure 4b). Species is a model input, so the mouse
and rat predictions for a given molecule and encoding differ; strain and sex
are set to `unknown`, and study dose is left missing, which the median
imputer fills for the ridge and support-vector members while the
boosted-tree member routes it natively. Six of the ten residues of each congener are
L-2,4-diaminobutyric acid, which has no single-letter code; this no
longer affects the shipped model's prediction (neither saved
configuration reads the ESM-2 sequence block) but it does affect the
HMMER homology diagnostic, which is why four surrogate spellings are
retained there and reported as a spread.

### Comparison with observed half-lives

Observed half-lives were taken from eighteen published studies curated in
`data/polymyxin_studies.csv` (flag `in_figure4`; the methods and
inclusion rules are set out in `figures/figs/fig4_methods_note.md`).
Five of these studies — Landersdorfer 2018, Jiao 2023, Yapa 2013, and
the two Lin 2017 papers — were added from a Europe PMC search for
population PK or humanized dosing of colistin or polymyxin B in mice or
rats; the remaining ten are a convenience sample curated during the
project. Three studies were added in revision (2026-10-17) to give colistin in
mouse an intravenous reference: Lin 2017 (colistin; intravenous arm, 2.64 mg
base/kg; its intratracheal arm is not used), Jin 2009 (5 mg/kg IV bolus; the paper
reports no half-life, so one was derived from the plotted plasma means) and Liao
2023 (a rat study; mouse experiments in that paper were safety and survival).
The primary comparison uses intravenous dose groups only; the all-route comparison
is secondary (Supplementary Figure 3). No search string, date, hit count, screening log, protocol or
second reviewer is recorded in the repository, and the full text of Lin
2018 could not be retrieved. The selection is therefore not a
systematic review and the sample may be biased toward studies that
report a half-life or population model.

Values are plasma or serum terminal half-lives. The CMS prodrug arms of
He 2013 and Yapa 2013, liposomal and tissue values of He 2013, Jiao
2023's 35.8 h value from a moribund-mouse study (clearance had
collapsed), and the kidney concentrations of Wang 2023 are not used.
For five studies that report population PK parameters rather than a
half-life (van der Meijden 2023, Landersdorfer 2018, both Lin 2017
papers and Yapa 2013), we derived a half-life as the slow-eigenvalue
of the published compartmental model and labelled it `derived`;
`ln2 × V_c/CL` from the same models (the central-volume half-life) is
reported for comparison but kept out of the main comparison. For van
der Meijden 2023 the half-life was derived as `ln2 × V_c/F ÷ CL/F` at
each of eight dose levels from a one-compartment model with first-order
absorption and dose-dependent clearance. The derivation reproduces
Landersdorfer 2018's own published 4.1 h value.

Mean residence times (Li 2003, 1.58 h; Wang 2023, 4.07 h) are reported
in Figure 4d as asterisks labelled "MRT" but are not compared with the
model and not averaged into the observed combined value: a MRT is not a
half-life. For Wang 2023 the MRT is drawn alongside the study's 2.53 h
terminal half-life; both MRTs render in the "other"-source grey, with
no study-row highlight.

Per-study methodology was extracted from the published full text of all
eighteen studies and is recorded in
`MuriVerse/data/polymyxin_study_groups.csv` (half-life type and the
analysis — NCA, one-compartment, population PK, derived slow-phase, or
`ln2 × V/CL` — animal state, matrix basis of unbound vs total drug, and
the sampling window in hours). The sampling window for Matar 2020 is
read from its Figure 10 (~6 h) because the paper's text does not state
the sampling schedule; every other row carries a value read directly
from the methods text. Figure 4d renders these fields as the two text
columns beside each row ("half-life · analysis" and "animals · matrix ·
sampled").

Observed combined values per drug and species were computed as the
animal-weighted geometric mean across dose groups, with the
animal-weighted standard deviation on log10 scale across dose groups
(`table49_polymyxin_observed_by_drug_species.csv`). Components of a
single dose group (B1, B2, ...) were averaged first so animals are
counted once. Sensitivity versions were computed with every study
counted once, with reported values only, with intravenous doses only,
and with each study left out in turn
(`table52`, `table55_polymyxin_pooled_sensitivity.csv`). The animal-
weighted choice was declared in advance and lets a study with many dose
groups dominate (van der Meijden 2023 carries 68% of the polymyxin B
mouse weight). The same weighting gives a value per source of
observation (`table52_polymyxin_observed_by_source.csv`); sources are
grouped into Jiao 2023, the Tam lab (University of Houston group of
V. H. Tam: Abdelraouf 2012, He 2013, Manchandani 2016), the Monash
group (Monash University group of R. L. Nation, J. Li and colleagues:
Dudhani 2010, Cheah 2015, Landersdorfer 2018, Sivanesan 2017, Yapa 2013
and both Lin 2017 papers; Li 2003 by the same investigators from the
University of South Australia) and `other`
(`MuriVerse/data/polymyxin_study_groups.csv`).

### Baselines

The model is compared with three constants read from the training split:
the median half-life of the same species, of both species combined, and of
the lipidated training peptides. These are reported in
`table27_polymyxin_fold_error_vs_baseline.csv`, `table26` and panel 4b.

**The constants are read from the records the applied model was fitted
on.** The applied models are the post-audit ones, so the constants are the
medians of the 146 post-audit intravenous records of the training split (27
mouse, 119 rat; corrected labels and species): mouse 0.547 h, rat 0.196 h,
both species 0.245 h, and the single lipidated intravenous training peptide
4.3 h. Recomputed on all 178 post-audit records they are 0.614, 0.245, 0.259
and 4.3 h, and no comparison changes. A route-pooled baseline would let the
constant borrow the extravascular labels — longer by the 0.45-0.65 log10
flip-flop offset — and so hold labels the model itself was not permitted to
use. The pre-audit run (`--model-suffix _iv`) uses the 176-record
intravenous training split (mouse 0.485 h, rat 0.194 h, both 0.218 h,
lipidated 1.75 h). `audit_figure_numbers.py` checks that `table26`'s
baseline split size and medians equal those of the post-audit intravenous
training split.

The post-audit PeptideCLM-embedding refit predicts 0.08–0.63 h for these
congeners and the post-audit refit of the top-scoring compact configuration
(physchem + study context + termini) predicts 0.45–0.79 h — both physically
plausible; on the 178 post-audit records their grouped cross-validated ρ is
0.292 and 0.456 against 0.445 for the shipped set (pre-audit: PeptideCLM
ρ = 0.247, unconstrained compact ρ = 0.587, shipped 0.574); the
embedding predictions are reported in
`table51_polymyxin_embedding_extrapolation.csv` and neither variant is
drawn. A previous build of this analysis reported both configurations
returning polymyxin half-lives of 10⁻¹¹⁶ to 10⁻¹³⁷ h; that was a
block-order mismatch between training (which concatenated blocks in
configuration-list order) and `predict_polymyxins.assemble` (which
reordered them by layout). The fix is in `_blocks_of` in
`evaluate_halflife_model.py` and `assemble` in `predict_polymyxins.py`;
the shipped configuration happens to list its blocks in layout order and
was unaffected.

### Why one script uses the route-pooled model

`evaluate_polymyxin_species_route.py` sets `MODEL_TAG = ""` and applies
the **route-pooled** models, against the shipped default. Its job is to
test externally whether a model mis-transfers across route. A model
fitted on intravenous records alone cannot answer that: it has no
subcutaneous level to predict under, so an extravascular prediction from
it is an extrapolation along an axis with zero variance and the route
effect it appeared to report would be an artefact of the scaler. The
instrument for the diagnostic has to be the model being diagnosed. This
is the same relation Figure 3's panels a-b have to its panels c onward,
and it is why `table56` is unaffected by the route restriction.

### Dose sensitivity

A dose-vs-prediction sweep
(`table54_polymyxin_prediction_vs_dose.csv`) is run over both animals,
since species is a model input, and only ratios to the dose-missing
prediction are reported. Over 5–40
mg/kg the post-audit shipped model's prediction moves by 0.3–0.4% (log–log slope about
zero), against 2–4.7-fold observed in the three dose-ranging studies
(log–log slope 0.26–0.74). **Only one of the 213 intravenous training
records carries a stated dose** (three of the 320 do), so the model's
dose response is an extrapolation from a single observation rather than
a learned effect. The dose panel of earlier builds is no longer drawn,
because a response resting on one training record should not be given a
panel; the sweep is still written to its table.

## Figures, tables and reproducibility

All main-text figures and their captions were generated by R scripts
under `figures/code/`, orchestrated by `figures/code/run_all.R`. Each
figure separates data preparation from drawing: a `NN_data_*.R` script
loads and summarises, and a `figN.R` script draws. We draw Figure 2 from
`04b_data_muriverse_strict.R`, Figure 3 from `05_data_halflife_model.R`
and Figure 4 from `06_data_polymyxin.R`.
Supplementary Figure 1 draws the full preprocessing pipeline that
feeds Figure 2 — the species composition of the three source databases
(PEPlife2, ThPDB2, PepTherDia2), the build funnel from 4,500 PEPlife2
records down to the 467 structurally verified rodent in-vivo records,
and the three reconciliation gates (route, structure, assay) that take
that set to the 320 modelling rows. It is drawn from
`supp_fig1.R`, which uses its own data loader at
`02_data_preprocessing.R` and does not source `fig2.R`. The mouse-blue
/ rat-green pair of Figure 2 (`#4169E1` / `#2E8B57`) is reused in the
one Supplementary Figure 1 panel split by species (f, mouse / rat
across stages); the gate panels (g, h, i) use the shared maroon
`COL_ANALYSIS` for kept-by-gate and grey for dropped-by-gate.
In Figure 2's funnel (panel a) each species's shade lightens toward
the eligible end and deepens toward the analysis set, so the colour
deepens in step with the cohort shrinking; every other Figure 2 panel
draws the two species in those same solid tones.
The UMAP in Figure 2h uses `uwot` with cosine distance,
15 neighbours, minimum distance 0.3 and seed 1986; the smoothed
Gaussian-kernel-weighted mean of per-sequence median log10 half-life
is drawn in light grey (surface opacity rises with height), and the
sequence points are coloured by species (grey for sequences dosed in
both animals).

Four commands rebuild everything downstream of the curated cohort:

```bash
python MuriVerse/build_strict_peptide_set.py
python MuriVerse/build_halflife_training_set.py --tier strict
python MuriVerse/evaluate_halflife_model.py --cohort rodent --tier strict \
    --jobs 16 --nested
Rscript figures/code/run_all.R
```

Encoder weights (`facebook/esm2_t33_650M_UR50D`;
`aaronfeller/PeptideCLM-23M-all`) are pulled from the Hugging Face
model cache. PeptideCLM additionally needs the SPE tokenizer vocabulary
at `tokenizer/`. The pipeline requires `torch`, `transformers`,
`xgboost`, `scikit-learn`, `joblib`, `rdkit`, `datasets` and `pandas`
(Python ≥ 3.11; the exact versions used are in `requirements.txt`), the
`phmmer` binary from HMMER 3.3 or newer, and R with `tidyverse`,
`ggpubr`, `uwot` and `ggnewscale`. The exact session information under which the
figures in the paper were produced is in
`figures/tabs/sessionInfo.txt`.

## Data and code availability

The MuriVerse repository (code, curation ledgers, filtering counts,
exclusions, OOF predictions and figure tables) is available at
<https://github.com/JacobKMcPherson/MuriVerse>. The curated 320-record
modelling table is at
`MuriVerse/rodent/strict/halflife_meta_with_split.csv`; the complete
467-record kept set with provenance at
`MuriVerse/rodent/strict/in_vivo_half_life.csv`; the two curation
ledgers at `MuriVerse/data/curation/`; the arm-level evidence at
`MuriVerse/rodent/pk_observations.csv`; and the refit production models
at `MuriVerse/rodent/strict/model_iv/` (pre-audit shipped, intravenous),
`MuriVerse/rodent/strict/model_iv_postaudit/` (post-audit, refit on the 178
source-verified intravenous records) and
`MuriVerse/rodent/strict/model/` (route-pooled, retained for the
transfer diagnostics); the held-out predictions behind every reported
score at `oof_predictions_iv.csv` and `oof_predictions.csv` beside them.
The HuggingFace Arrow datasets of
pre-computed embeddings are intermediate build products and are not checked in;
the evaluation uses the published embedding caches.

## Statistical reporting and limitations

All reported metrics are means over three random seeds; intervals are
95% group bootstraps (500 resamples), scoring each seed separately before
averaging. The final 178-record post-audit cohort has 50 connected
`cluster_PMID` groups; the historical 213-record pre-audit cohort has 55,
and the route-pooled 320-record cohort has 75.
Pre-stated targets are reported pass/fail in the Figure 3 caption. No
multiple-testing correction is applied: the three metrics (Pearson r,
Spearman ρ, AUC) describe three facets of one prediction, not three
independent hypothesis tests, and the ablation prices every
configuration separately and is not interpreted as a hypothesis test.

The headline evaluation is limited by between-source measurement
disagreement. Taking every sequence x species combination that more than
one source reports, reducing each source to its mean log10 half-life for
that combination and forming every unordered pair of sources, the
intravenous cohort holds 33 such pairs, agreeing at r = 0.75 with a
median disagreement of 0.35 log units against a label spread of 0.82.
Two laboratories measuring the same peptide in the same animal differ by
about a factor of two, which bounds what any predictor can score against
a single held-out reported label. (Earlier versions of this document
quoted "66 pairs, r = 0.57" from a pairing rule no script in the
repository generated and no reconstruction reproduced; the rule above is
executed by `audit_figure_numbers.py`.) About three quarters of the
residual variance of the held-out predictions is a constant offset per
source; if each source-level offset were known, Pearson r would be 0.93.

Selection optimism is material: nested selection of the feature
configuration inside each training fold gives ρ = 0.460 against the
headline 0.574, so roughly 0.11 of the headline is the configuration
having been chosen on these data.

The route restriction improves every metric on a like-for-like
comparison (ρ 0.514 → 0.574 on the same 213 records with the same shipped
configuration), but the ρ gain alone is not statistically resolvable: a
paired bootstrap over the 55 connected cluster/source groups gives
+0.063 [-0.005, +0.157]. The restriction is justified by the pharmacokinetic
rationale of flip-flop kinetics and by the 0.45–0.65 log₁₀ cross-route
bias, not by that correlation.

The dose, strain and formulation covariates are absent or unparsed on
most records — only one of the 213 intravenous training records carries a
dose — and the retrieval ledger marks 85 of the 120 article PMIDs in
the cohort as not retrieved; nine patent identifiers account for a
further 91 records. Species in the model is reached only through strain,
and the polymyxin application shows a true species gap of 0.31–0.55 log₁₀
that the model predicts as approximately zero, so Figure 3's small
cross-species bias should be read as evidence that the model does not
distinguish the species, not that the species do not differ. The
polymyxin held-out application is a small test of four molecules from
outside the training chemistry, not a general estimate of performance on
novel chemistry.
