> **Archived pre-audit draft, not release documentation.** Its cohort, endpoint and polymyxin claims are superseded by manuscript v10, [methods.md](methods.md), [results.md](results.md) and [final-model-release.md](final-model-release.md). Do not use the numbers below for the final 178-record model.

# Results

We describe the training set (Figure 2, Supplementary Figure 1), the
cross-validated evaluation (Figure 3) and the held-out application to colistin and
polymyxin B (Figure 4).

The sets reported below form a single nested chain, each a strict subset of the one before
it: 4,500 source-database records, 1,159 eligible rodent records, 467 strict peptide
records, and the 320-record analysis set. Every count here refers to one of those four.

## A structure-verified rodent half-life cohort of 320 records

The modelling unit is one source observation. Applying the MuriVerse label definition to
the PEPlife 2.0 pool leaves **320 records — 74 mouse and 246 rat — covering 166 distinct
sequence-embedding inputs, 129 sources and 90 combined-species homology clusters**
(Figure 2b).
Every record is a single peptide chain of 3 to 50 residues, parenterally dosed, measured by
a bioanalytical method that a source names, and carrying a peptide SMILES for the whole dosed
molecule that has been independently decoded back to its residue sequence. Reported
elimination half-lives span 0.0028 to 62.0 h with a median of 0.427 h (Figure 2c).

The selection funnel runs from 1,159 eligible records to 320 (Figure 2a). The largest
single losses are the exclusion of polyethylene glycol, Fc, albumin and polymer
conjugates, removing 297; the requirement for a residue-level sequence, removing 135; and
the removal of records with another dosed partner, removing 139. Three cohort
requirements then take the 467 records that clear the structural gates down to 320: the
dose route must be intravenous, subcutaneous, intraperitoneal or intramuscular, removing
84; the SMILES must be complete, removing 57; and the assay must be named, removing 6.

Species overlap is small but real: of the 166 unique input sequences, 12 occur in
both animals, and 7 of the 90 homology clusters span both species. We count 246
records, 128 sequences and 62 clusters for rat against 74, 50 and 35 for mouse
(Figure 2b). Mouse and
rat are therefore combined for model development, and whole clusters are kept on one side
of any split.

Peptide length has a median of 29 residues (Supplementary Figure 1j). The route
distribution is dominated by intravenous dosing — 213 intravenous, 73 subcutaneous, 22
intraperitoneal and 12 intramuscular, with **no record reading "not reported"** (Figure
2e). Filtering also collapses the engineering vocabulary: 4 of the 9 classes present in
the eligible set survive into the analysis set (Figure 2d), where they are unmodified or
none stated, with 189 records; backbone modification only, with 101; radio or fluorescent
label only, with 21; and lipidated, with 9. The 176-record intravenous training split the
shipped model is fitted on holds 114, 50, 9 and 3 of them respectively. A dimensional
reduction of the mean-pooled sequence embeddings over the 166 sequences shows that long-
and short-half-life peptides occupy distinguishable but overlapping regions of sequence
space (Figure 2h).

### Audit and provenance of the cohort

We provide the panels that let a reader audit each record's origin in Supplementary
Figure 3.
Sources are concentrated but not dominated by one study: the ten largest sources
contribute 126 records, or 39% of the set, and the largest single source, a patent,
contributes 61 records, or 19% (Supplementary Figure 1a). The homology clusters used as
the cross-validation unit are correspondingly fine-grained, with 90 clusters of which 57
are singletons and the largest holds 40 records (Supplementary Figure 1b).

Chemically the set is more diverse than the "natural peptide" label suggests. 113 of 320
records carry an annotated terminal modification, overwhelmingly C-terminal amidation
with 95 records and N-terminal acetylation with 12 (Supplementary Figure 1c). 67 records
are cyclic against 253 linear, and chirality is 291 all-L, 24 mixed L and D, and 5 all-D
(Supplementary Figure 1d); 268 of 320 are composed only of natural L-residues
(Supplementary Figure 1e).

The whole-molecule requirement is drawn over the 467 structurally kept records rather
than the 320 analysed ones, because every analysed record is complete by construction:
406 of the 467 carry a SMILES for the whole dosed molecule and 61 do not
(Supplementary Figure 1f). Those 61 fail for reasons the SMILES builder records — an
unconfirmed substitution in 23 cases, an unrecognised modification in 16, unspecified
cyclisation in 8, a conflict between the annotation and the sequence in 7, a non-standard
residue that could not be drawn in 3, a backbone surrogate in 1, and other causes in 3
(Supplementary Figure 1g).

Two panels bear directly on the modelling result. Half-life differs systematically by
assay method, with medians spanning 0.16 to 1.85 h across bioassay, liquid
chromatography-mass spectrometry, ELISA or immunoassay, high-performance liquid
chromatography, radiometric or imaging methods, and radioimmunoassay (Supplementary
Figure 3h) — assay is a study-context feature, not a molecular one. And route provenance
is uneven: **154 of 320 routes read back to a primary source, with the weaker
database-column provenance concentrated in the longest-half-life quintiles**
(Supplementary Figure 1i).

## The shipped model is restricted to intravenous records (Figure 3a, b)

The cohort of 320 records is the *input* to model development, not the set the shipped
model is fitted on. **The shipped model is fitted on the 213 intravenous records**, being
179 rat and 34 mouse. The restriction is a result in its own right, and Figure 3's first
two panels are the measurement behind it.

The physical argument is that a route-pooled half-life is not one quantity. After an
intravenous dose the terminal half-life is an elimination half-life. After an
extravascular dose it is one only if absorption is the faster step; when absorption is
slower, the terminal phase reports the absorption rate instead. Pooling the routes
therefore gives the regressor a target that mixes two different physical quantities.

The data agrees, and it agrees asymmetrically. Cutting the cohort into four
species-by-route combinations and refitting the same ensemble on four different training pools
per combination separates the two axes that are confounded in any marginal comparison, rat being
73% intravenous against mouse's 46%. We report the four combinations in Publication Table 1.

To separate the two axes, we partitioned the 320 records into four species-by-route
combinations — rat IV with 179 records, rat EV with 67, mouse IV with 34, and mouse EV with 40 —
and refitted the same ensemble on four different training pools per combination:
route/species-matched, which is the same species and the same route as the test combination and
carries the evaluator key `same cell`; cross-species, which is the other animal by the
same route; cross-route, which is the same animal by the other route; and pooled, which
is every record outside the test combination. Each pool was fitted under grouped five-fold
cross-validation that held out whole connected groups of homology cluster and source
identifier, so neither family nor publication reached the scoring model, and the
protocol was repeated over five seeds because the smaller combinations are thin. All four pools
were scored on exactly the same held-out records of each combination, so the four colours in each
row of Figure 3a are directly comparable; the mean ± 1 SD over seeds is the point and
bar, with the per-seed values as small dots. We report ranking (Spearman ρ,
left), bias (mean(predicted − observed) in log₁₀ h, middle) and evidence (−log₁₀ of a
one-sided t p-value for positive ρ, right), with sub-chance correlations clipped to zero
rather than registered as significant; because the held-out records of a grouped pass are not
independent, the evidence facet is a guide to scale rather than an exact test.

**Only rat IV was modellable on its own**, at route/species-matched ρ = 0.62. The three smaller combinations
each sat below ρ = 0.25 when trained only on their matching records, and the mouse-EV-only
pool was worse than chance at ρ = −0.19, with five per-seed values of −0.09, −0.42,
−0.36, −0.19 and +0.13 that span almost the whole plausible range for 40 records.
**Pooling every record outside the test combination raised three of the four combinations and cost rat IV almost
nothing**: mouse IV rose from 0.23 to 0.41, mouse EV from −0.19 to 0.18 and rat EV from
0.19 to 0.39, while rat IV gave up only 0.62 to 0.57. *Transfer behaviour was
combination-specific.* For mouse IV the species axis carried, at cross-species ρ = 0.43 against
cross-route 0.14, so a sufficient volume of rat intravenous data substituted for the
thin mouse intravenous one. For rat EV and mouse EV the route axis carried, at 0.44
against −0.06 for rat EV, and 0.27 against 0.18 for mouse EV; extravascular combinations were
therefore better served by extravascular data from the other animal than by intravenous
data from the same one. These rankings justify combining species, and the bias measurement
of the middle facet, which we redraw against the two transfer axes, is what then forces
restricting by route (Figure 3b).

**Every cross-route transfer is offset by 0.45 to 0.65 log₁₀ — a factor of 3 to 4 — and
every cross-species transfer stays inside ±0.21.** The signs are the ones flip-flop
kinetics predicts: models trained on intravenous records under-predict extravascular
labels, which are absorption-limited and therefore longer, and models trained on
extravascular records over-predict intravenous ones. That this reproduces independently
in mouse and in rat, at near-identical magnitude and with opposite signs either way,
means it cannot be an artefact of species composition.

The within-molecule version of the same result is Figure 3h: across the 32
sequence-and-species combinations measured by both routes, so that the molecule and the
animal are fixed and only the route varies, subcutaneous exceeds intravenous in 24 of 32
pairs, or 75%, at a paired Wilcoxon p = 0.001. The paired shift is modest, at a median
+0.068 log₁₀ or a 1.17-fold ratio, and much smaller than the marginal gap between the
route medians, which is 0.24 h intravenous against 0.96 h subcutaneous, or roughly
4-fold. Most of that marginal difference is composition: the two routes were used on
different peptides.

There is a second, more prosaic reason. Route is a free lever in a route-pooled cohort:
scoring every record by its own route's median half-life, with no molecule involved at
all, reaches Spearman ρ = 0.42 against the route-pooled label. Part of any route-pooled
headline is therefore the model recovering a study-design variable.

**What the restriction is worth, measured like for like.** Publication Table 2 scores the
same 213 intravenous records, under the same shipped configuration and the same grouped
protocol, once as a route-pooled fit and once as an intravenous fit. Spearman rho rises
from 0.507 to 0.574, the area under the curve at 6 h from 0.785 to 0.818, mouse rho from
0.427 to 0.567, and the root-mean-square error falls from 0.699 to 0.675 log10 hours.

The gain in ρ is **not on its own statistically resolvable**: a paired bootstrap over the
55 connected cluster-and-source groups puts it at +0.066 with an interval of −0.016 to
+0.138. The case for restricting by route is therefore the physics, the bias measurement
above, and the removal of the route lever — not the headline correlation. A route-pooled
model scores partly by recovering which route a record came from, and its predictions
have no single route at which to be read.

An earlier version of this analysis reported the gain as ρ rising from 0.514 to 0.587
with a bootstrap interval excluding zero. That comparison was against the configuration
of physicochemical descriptors plus study context plus terminal classes, which is the
top-scoring configuration in the ablation but is **not** the shipped one. It is recorded
here because the distinction matters and the superseded numbers should not reappear.

## Cross-validated half-life prediction on the intravenous cohort (Figure 3c–g)

The model is an equal-weight ensemble of ridge regression, a radial-basis-function
support-vector regressor and a shallow gradient-boosted tree model, with every
hyperparameter fixed in advance, over **71 compact columns**: molecule annotation with 9,
study context with 33, including the mouse and rat indicator pair, terminal classes with
13 and physicochemical descriptors with 16.
No sequence or SMILES embedding is used. Under the route restriction the route columns
are constant in training and carry no information, and standardisation leaves a
zero-variance column at zero. The species columns are not constant: the intravenous
cohort is 34 mouse and 179 rat records.

Evaluation holds out whole groups under five-fold cross-validation, reshuffling the fold
assignment for each of three seeds; each seed is a complete out-of-fold pass and the
headline numbers are the mean of the three per-seed scores. The grouping is the **union
of homology cluster and source identifier, which yields 55 connected groups over the 213
records**, so that a held-out peptide shares neither a family nor a publication with the
records it is predicted from. Intervals are 95% group bootstraps over those groups.

### Performance over both species and by species

We report the held-out metrics over both species and by species in Publication
Table 3.

Against the three pre-stated targets, the strict evaluation gives Spearman ρ = 0.574
with an interval of 0.322 to 0.655, so the target above 0.50 is **met**; an area under
the curve for half-life at or above 6 h of 0.818, so the target above 0.80 is **met**;
and Pearson r = 0.568 with an interval of 0.381 to 0.651, so the target above 0.70 is
**not met**. The classification target, missed at 0.785 on the route-pooled cohort, is
met on the intravenous one; the correlation target is not.

Selection optimism is larger than it was on the route-pooled cohort. Choosing the feature
configuration inside each training fold rather than on the folds that score it gives
r = 0.456 and ρ = 0.460, so roughly **0.11 of the headline ρ is selection optimism** and
the honest nested estimate is ρ of about 0.46. This should be quoted alongside the
headline, not instead of it: the shipped configuration is fixed in advance of any single
fold, but it was chosen on these data.

The usable result remains a ranking, not a calibrated absolute prediction. We obtain
mouse ρ = 0.56, with a seed range of 0.50 to 0.68 over 34 records, and rat
ρ = 0.60, with a range of 0.55 to 0.65 over 179 records, both as three-seed means
(Figure 3c). Mouse is the volatile stratum at 34 records: it reads above rat at seed 0,
0.68 against 0.65, and below it in the other two seeds, which is why the three-seed mean
is the number quoted. The
per-species gap is much narrower than it was on the route-pooled cohort, where the mouse
stratum reached only ρ = 0.378. Note that this route-pooled figure is computed over all 320
records rather than over the 213 compared in the table above, where the like-for-like
mouse pair is 0.427 against 0.567; the comparison here is a statement about composition
rather than a performance delta. The reason is visible in that composition: mouse was 54%
extravascular, so the mouse stratum of a route-pooled model was mostly records whose label
was a different physical quantity.

Ordering and absolute accuracy come apart in mouse, which is the distinction a
correlation alone hides. We find that **57% of held-out records fall within
3-fold of the observed half-life**, being 41% in mouse and 60% in rat, with a median
absolute fold error of 4.1-fold in mouse against 2.2-fold in rat and a signed bias of
−0.09 and −0.04 log₁₀ respectively (Figure 3d). The model is close to centred in both animals; the
mouse penalty is spread, not offset. Read as a screen at the pre-declared 6 h target
(Figure 3e) the area under the curve is 0.82, ranging 0.77 to 0.85 across seeds, but only
17 of the 213 records are positive, being 8 mouse and 9 rat, so the curve is coarse and
the interval wide.

Error varies by about three-fold across descriptive strata of the same held-out
predictions (Figure 3f). By assay, radiometric and imaging records are worst at 5.8-fold
over 13 records and radioimmunoassay best at 1.7-fold over 90. By engineering class,
radio- or fluorescent-label-only peptides reach 4.9-fold over 17 records against 2.0-fold
for unmodified ones over 135. The largest single source is the best-predicted stratum in
the panel, at 1.4-fold over 54 records, and singleton sources the worst, at 4.3-fold over
53 — the source-offset effect seen from the error side. A single global error or
prediction interval is therefore insufficient for pharmacokinetic decision-making.

### Assay method is the largest single input; molecular descriptors contribute jointly

The feature ablation under the strict grouping (Figure 3g) scores 29 configurations, of
which Publication Table 4 ranks the twelve the readings below turn on.

Three readings follow, unchanged in direction from the route-pooled analysis and sharper
in magnitude. First, **study context carries most of the signal**: the 33 context columns
alone reach ρ = 0.360 and adding the physicochemical block reaches ρ = 0.521, whereas the
9 molecule-annotation columns alone reach only ρ = 0.104. Second, **the large pretrained
embeddings do worse, not better**: every configuration containing the 1,280-dimensional
sequence embedding or the 768-dimensional SMILES embedding scores below the compact set
under held-out clusters and sources, while rising when only sources are held out — the
embeddings recognise families and publications, and the strict grouping removes exactly
that. Third, 71 hand-built columns beat 2,048 learned dimensions on this cohort, which is
a statement about the cohort's size, not about the embeddings.

### The shipped configuration is not the top-scoring one

The configuration of physicochemical descriptors plus study context plus terminal classes
outscores the shipped set by 0.013 in ρ (0.587 against 0.574) but omits the
molecule-annotation block. The shipped configuration keeps the annotation block because
several physicochemical columns cannot be computed at all for a peptide with no faithful
single-letter sequence, and without it nothing holds such a molecule inside the trained
range: the unconstrained winner returns polymyxin half-lives of about 10⁻¹²⁰ h against an
observed 1–2 h. The 0.013 cost is well inside the seed-to-seed spread. The ablation still
scores and reports every configuration — including the unconstrained winner — so the
constraint's price is visible rather than hidden.

### The ceiling is between-source disagreement

Performance degrades as the holdout gets stricter (Publication Table 5).

The last line of the table is the diagnostic. With each source's mean removed from both the prediction
and the label, the model still orders peptides inside a single source at r = 0.433 and
ρ = 0.462 — appreciably less than its headline score. The difference between the strict
evaluations and this one **is** the source-level offset: a constant per publication that the
model implicitly recovers when it can see the source's other records and cannot know for a
genuinely new source. Directly, **about three quarters of the residual variance is a
constant offset per source**: supply that offset and r rises from 0.57 to 0.93.

The literature sets the same bound from the other side. Taking every
sequence-and-species combination that more than one source reports, reducing each source
to its mean log₁₀ half-life for that combination, and forming every unordered pair of
sources, the intravenous cohort contains **33 such cross-source pairs, agreeing at
r = 0.75 with a median disagreement of 0.35 log units against a label spread of 0.82**,
the last being the standard deviation of log₁₀ half-life. Two laboratories measuring the
same peptide in the same animal differ by about a factor of two, which caps what any
predictor can score against a single held-out reported label. An earlier version quoted
66 pairs agreeing at r = 0.57 with a disagreement of 0.39, from a pairing rule that could
not be reproduced; the rule above is the one now executed and verified.

We make the same point from the data side (Figure 3i). The saved 80/20 split is
cluster-preserving, so whole homology clusters — and the label ranges that travel with
them — land on one side, and the training labels, with 176 records at a mean of −0.53,
and the validation labels, with 37 records at a mean of −0.57, differ in density. This is
why the headline evaluation reshuffles folds across three seeds rather than relying on
that single draw, and no performance number in this analysis comes from it.

The practical consequence is that the unmet correlation target is unlikely to be closed
by more features of the kind already present. Route is now established for all 320 records,
154 of them from a primary source, and establishing it did not raise r — what it did was
identify which records belong to one physical quantity. What remains missing is the rest of
the source-level design. Dose, strain and formulation are absent or unparsed on most
records; **only one of the 213 intravenous training records carries a dose at all**.
Recovering them requires the source documents themselves, of which the retrieval ledger
marks 85 of the 120 article identifiers in this cohort as not retrieved, with 9 patent
identifiers accounting for a further 91 records.

## Held-out application to colistin and polymyxin B (Figure 4)

### Setup and leakage

Nothing is retrained for the case study. The prediction script loads the
**post-audit** models and applies them to four
molecules — polymyxin B1, polymyxin B2, colistin A and colistin B — each built from its own PubChem
isomeric SMILES. The earlier pre-audit models remain available with
the pre-audit option; the manuscript's case study uses the post-audit models.

The test is genuinely held out. The leakage check reports
**zero** records whose name is a polymyxin in the strict training set and zero in
the modelled cohort; the only mention anywhere in the eligible set is
colistin appearing as a comparator in the free-text dosing notes of the ZY4
record, whose half-life is ZY4's. No training peptide has detectable homology
to any congener: every surrogate query in
the homology table reports no sequence at or below the
cohort's clustering E-value of 10⁻³ (best E-value 1.2–3.0 across 236–256 sequences
reported at E ≤ 1000), and the query is run against all 320 analysis-set
sequences, a superset of the 178 post-audit (and 213 pre-audit) intravenous records the models are
fitted on, which is the conservative choice.

Each congener is predicted under two input variants crossing the N-terminal
acyl wording (literal or normalised), giving 16 shipped-model predictions across
the two species. **Only intravenous inputs are predicted**: under the route
restriction the route columns are constant in training, so a subcutaneous input
would be an extrapolation along an axis with zero training variance rather than a
different point in a learned dimension. **The C-terminus is annotated as
unmodified only.** An Amidation variant was previously run alongside it, on the
reasoning that the macrolactam bond is itself an amide; that reasoning was wrong
and the variant has been removed. A C-terminal amide is –C(=O)NH₂, whereas the
ring is Thr10–C(=O)–NH–Dab4, an isopeptide bond to a side-chain amine: both are
amides, but only one is a *C-terminal* amide, and these molecules have the ring.
Annotating amidation therefore asserted a modification the molecules do not carry
while leaving the SMILES — which already encodes the macrocycle — untouched, and
because it moved the prediction further than the drug, the congener or the
species did, it dominated every real contrast in the analysis. The cyclisation is
carried by topology, and Free is used in its sense of no terminal
modification recorded. The N-terminal variants test sensitivity to how a molecule
is written down and are **not** independent measurements — the wording changes
nothing at all, the two rows being bit-identical — so the 16 predictions carry 8
distinct values. Species is a model input, so the mouse and rat predictions for a
molecule differ; strain and sex are set to unknown and study dose is left missing.

### Against intravenous observations the model is long in mouse and short in rat

The comparison is intravenous only, because the models are fitted on intravenous records and predict
intravenous inputs. The comparison against observations from every route is secondary (Supplementary
Figure 3; the _pooled tables). Across the 16 shipped predictions the span is **0.57–0.92 h**; per drug
and species the median is colistin 0.92 h in mouse and 0.78 h in rat, polymyxin B 0.74 h and 0.62 h
(Figure 4b, Table 6). Within a species the four congeners spread 1.38-fold,
and that spread is not symmetric: colistin A and B agree to within 1.02-fold, while polymyxin B1 and B2
differ by 1.20-fold although the two pairs differ by the same single methylene in the acyl chain.

No median over all 16 is quoted, because the predictions split by drug (polymyxin B 0.57–0.80 h, colistin
0.77–0.92 h) and a pooled median would hide the species and drug structure the comparison is about; the
comparison is made per drug and species throughout.

The animal-weighted geometric means of the intravenous dose-group half-lives are 0.45 h for colistin in
mouse (two studies, 37 animals), 1.59 h for colistin in rat (six studies, 36
animals), 0.58 h for polymyxin B in mouse (two studies, 56 animals) and 1.90 h for
polymyxin B in rat (three studies, 33 animals). The model reads **long in both mouse
cells** — 2.05-fold for colistin and 1.26-fold for polymyxin B — and **short in both rat
cells**, 2.05-fold for colistin and 3.06-fold for polymyxin B (median per-congener fold errors; the
range across the four rows of a cell is 2.04–2.07, 1.15–1.38, 2.03–2.06 and 2.78–3.34). The mouse
agreement reported before the comparison was made intravenous came from subcutaneous half-lives, which the
model was not fitted on. The pattern follows from the species term being far too small: intravenous
half-lives are observed 3.6-fold (0.55 log₁₀) longer in rat than in mouse for colistin and 3.3-fold
(0.51 log₁₀) for polymyxin B, while the post-audit model places rat about 15% (0.072 log₁₀) *shorter*
than mouse. The column is present and its sign is wrong.

The model does not reproduce the direction of the between-drug difference. It
predicts colistin at 1.24× polymyxin B in both species, where the intravenous
observations give polymyxin B the longer, by 1.31× in mouse and 1.19× in rat. Those observed differences sit
inside the spread between dose groups, so no ranking claim is supported either way;
the 1.20-fold spread between polymyxin B1 and B2, two congeners the model should
treat as near-identical, is as large as the between-drug margin it gives.

### Against constants that need no model

We price the predictions against constants taken from the training
split alone (Figure 4b). **The baseline must come from the data the applied model was
allowed to use**, which under the route restriction is the intravenous training
split: the 146 post-audit intravenous records of the saved training split (27 mouse and 119 rat records,
with the audit's corrected labels and species). The constants are the median of the same species
(mouse 0.547 h, rat 0.196 h), the median of both species (0.245 h), and the one lipidated intravenous
training peptide (4.3 h). The applied model was fitted on all 178 post-audit records, so the constants
come from a subset of the rows it saw; recomputed on all 178 they are 0.614, 0.245, 0.259 and 4.3 h, and
the comparison is unchanged. A route-pooled baseline would instead let the constant borrow the
extravascular labels — longer by the 0.45–0.65 log₁₀ flip-flop offset established in Figure 3a–b — and
the comparison would then be against labels the model itself was not permitted to see.

**A constant was closer than the model in three of the four intravenous cells: the mouse-median constant
in both mouse cells (1.22-fold from the colistin reference against the model's 2.05, 1.07-fold from the
polymyxin B reference against 1.26) and the lipidated constant in rat polymyxin B (2.27-fold against 3.06).
In rat colistin the model was closer than every constant** (2.05-fold against 2.70 for the lipidated
constant, 2.91 for the mouse median and 8.13 for the rat median). Against its own species median the model is closer
in 8 of the 16 shipped predictions (every rat row, where that constant is 8.1–9.7-fold off), and in none
of the mouse rows. The lipidated constant rests on a single training record and is to be read with great
caution. Against observations from every route (Supplementary Figure 3) the picture is more favourable to
the model, but that reading rests on subcutaneous and intraperitoneal half-lives the model was not fitted
on.

### The verdict depends on which observed half-life is the reference

The intravenous references are few. Polymyxin B in mouse combines He 2013
(0.27 h; the model is 2.71-fold long against it) and Lin 2017 aerosolized
polymyxin B (1.03 h, derived by us from its population model; 1.42-fold short);
the cell is 1.26-fold long with both. Colistin in mouse combines Lin 2017 colistin (0.62 h, infected mice,
reported) and Jin 2009 (0.29 h, healthy mice; **the paper reports no half-life, so
we derived it by log-linear regression of the last three plasma means plotted in
its Figure 2, read by eye; 0.26–0.33 h under alternative readings**): the model is
1.48-fold long against Lin alone, 3.14-fold against Jin alone and 2.05-fold
against both. Colistin in rat spans six studies, from Liao 2023 (0.21 h; a rat study whose
3 mg/kg dose was infused over 15 minutes and whose half-life is shorter than the
infusion) to the slow phase of dry-powder colistin in Lin 2017 (10.3 h, a weakly
determined second compartment); the fold error is 2.05 with all six, 2.47 without
Liao 2023 and 1.66–2.47 across every single deletion. Polymyxin B in rat is the
most stable cell, 2.31–3.36-fold across deletions.

Across every route, laid side by side (Supplementary Figure 3), the spread is wider still:

| Drug, species | Reference | Observed t½ | Fold error | n animals (studies) |
|---|---|---|---|---|
| Polymyxin B, mouse | Jiao 2023 | 3.78 h | 5.13× **under** | 6 (1) |
| Polymyxin B, mouse | Tam lab, reported | 0.27 h | 2.71× over | 24 (1) |
| Polymyxin B, mouse | Monash group, derived | 1.61 h | 2.18× **under** | 47 (2) |
| Polymyxin B, mouse | other studies, derived | 1.20 h | 1.64× **under** | 160 (1) |
| Polymyxin B, rat | Tam lab, reported | 2.07 h | 3.31× **under** | 25 (2) |
| Polymyxin B, rat | Monash group, reported | 1.45 h | 2.32× **under** | 8 (1) |
| Colistin, mouse | Monash group, reported | 1.13 h | 1.23× **under** | 205 (3) |
| Colistin, mouse | Monash group, derived | 0.29 h | 3.14× over | 16 (1) |
| Colistin, rat | Monash group, reported | 1.36 h | 1.75× **under** | 13 (2) |
| Colistin, rat | Monash group, derived | 4.07 h | 5.24× **under** | 12 (2) |
| Colistin, rat | other studies, reported | 1.09 h | 1.40× **under** | 17 (3) |

Against Jiao 2023 (subcutaneous) the model is 5.1-fold **short**; against the Tam-lab
mouse value it is 2.7-fold **long**. A single prediction of about 0.7 h cannot be
right about both, and the model has no infection or sampling-design
input with which to choose. That a prediction lying between the longest and the shortest of the
published values agrees with neither extreme is a statement about those values, not evidence that the
model has captured the pharmacokinetics.

The underlying disagreement is methodological rather than simply disputed.
Reading the studies side by side separates the short Tam-lab values from Jiao
2023 along four axes — observation window and estimator (designs ending at
6–7.5 h with four to six sampling times, fitted per animal or to destructively
sampled averages, against a population fit to every sample of every mouse),
health and route (healthy and intravenous against infected and subcutaneous,
with clearance reported about two-fold lower in the lung model), species, and
the choice of reported versus derived quantity. The effect is not a laboratory
signature: within the Monash group, noncompartmental half-lives are short
(Sivanesan 2017, 1.3–1.6 h in rats; Li 2003, 1.24 h for colistin) while the
slow phases of that group's own population compartmental models are longer for
three of four (Landersdorfer 2018, 4.1 h; Yapa 2013, 3.0 h; dry-powder colistin
Lin 2017, 10.3 h), and the central-volume half-life ln2 × Vc/CL of those same
models is short again (0.46, 0.49, 1.46 and 0.60 h). Which of these is "the"
half-life depends on whether the terminal phase was identified, which depends
on how long and how often the animals were sampled. No study measured both
designs in the same animals, so the discrepancy is not resolved here.

### The combined reference is itself fragile

The pooled-sensitivity table shows that the all-route fold errors
are a property of the weighting, not of the data. For polymyxin B in mouse, van
der Meijden 2023 carries 67.5% of the animal weight (eight dose groups of one
study) and values we derived carry 87.3%. Re-pooling moves the fold error:

| Drug, species | Animal-weighted | Each study once | Reported values only | IV doses only |
|---|---|---|---|---|
| Polymyxin B, mouse | 1.13 h → 1.55× | 1.39 h → 1.91× | 0.46 h → 1.60× | 0.58 h → 1.26× |
| Polymyxin B, rat | 1.90 h → 3.06× | 1.72 h → 2.77× | 1.90 h → 3.06× | 1.90 h → 3.06× |
| Colistin, mouse | 1.02 h → 1.11× | 0.72 h → 1.27× | 1.13 h → 1.23× | 0.45 h → 2.05× |
| Colistin, rat | 1.70 h → 2.19× | 1.63 h → 2.10× | 1.20 h → 1.54× | 1.59 h → 2.05× |

Leaving any single study out moves the animal-weighted polymyxin B mouse value only
between 0.99 and 1.33 h, so it is stable to single deletions. It is **not** stable to a
change in what counts as an observation: the polymyxin B mouse fold error stays between 1.26 and 1.91
across the four weightings, but the sign flips (the 0.74 h prediction lies between the 0.46 h
reported-only reference and the 1.13–1.39 h pooled references), and for colistin in mouse the fold error
runs from 1.11 (animal-weighted) to 2.05 (intravenous doses only). The headline range is the
intravenous-only choice, which is route-matched to the model.

### Dose does not rescue the prediction, and the model cannot represent it

In the observations, half-life rises with dose within studies that varied it,
in the direction expected of saturable clearance. For colistin SC, Dudhani 2010
reports the half-life rising from 0.30 h at 5 mg/kg to 1.39 h at 40 mg/kg
(4.7-fold over an 8-fold dose range), and Cheah 2015 reports it rising from
1.21 h at 10 mg/kg to 2.63 h at 40 mg/kg. For polymyxin B SC, van der Meijden
2023 reports it rising from 0.64 h at 0.5 mg/kg to 2.26 h at 64 mg/kg. The
log–log slopes across the three studies span 0.26–0.74. Dose does not,
however, explain the disagreement between sources: Jiao 2023 and van der
Meijden 2023 both dose mice subcutaneously at 8–16 mg/kg, yet their reported
half-lives differ by 2.1–3.2-fold; He 2013 and Manchandani 2016 both dose at
3 mg/kg IV, in different species, yet their reported half-lives differ by
9-fold.

In the model, **only one of the 213 pre-audit intravenous training records carries a
dose** (the 178 post-audit records are a subset), so any dose response is an extrapolation from a single observation and
is not interpretable. Supplying a dose barely moves the prediction at all:
from 1.000-fold at 0.1 mg/kg to 0.991–0.993-fold at 64 mg/kg, a total span of under 1%
across a 640-fold dose range, in the same direction for both drugs and both
animals. Over the
5–40 mg/kg range it moves by under 0.5% against 2–4.7-fold in the
observations — a log–log slope of about zero against 0.26–0.74. The dose panel of
earlier builds is no longer drawn, because a response resting on one training
record should not be given a panel.

### What the external route and species test shows (Figure 4c, d)

Figure 3a–b made two claims from inside the training cohort: that cross-route
transfer is biased by 0.45–0.65 log₁₀, and that cross-species transfer stays
inside ±0.21. Figure 4c tests both against the polymyxin data, from outside the
cohort. We deliberately apply the **route-pooled** model, because a model fitted on
intravenous records alone has no extravascular level to be tested across; the
instrument that justifies the restriction has to be the unrestricted model.

Route is at least directional. The model predicts a constant +0.19 log₁₀ for the
intravenous-to-extravascular shift against observed effects that vary by case —
+0.22 for colistin in rat, +0.46 for colistin in mouse and +0.70 for polymyxin B in
mouse — close to the smallest and 2.4- and 3.7-fold smaller than the other two,
because a single one-hot route column can produce only one offset. The scatter
around the diagonal is consistent with the 0.45–0.65 log₁₀ cross-route offset
reported internally.

Species it cannot represent. The literature shows rat half-lives 0.31–0.55 log₁₀
(2.1–3.6-fold) above mouse for the same drug and route, and the shipped feature
set predicts approximately zero for every one of them. Species reaches this model
only through strain, and the strain levels are completely interleaved across the
two animals, so the full span across all strain levels (0.300 log₁₀) is
narrower than the 0.51 and 0.55 log₁₀ intravenous species gaps the literature shows
for polymyxin B and colistin. This qualifies Figure 3's species result: the small
cross-species offset reported there is not evidence that species barely matters,
but that the model does not distinguish species from the strain encoding alone.

### The embedding and physchem-only variants, applied here but not drawn

The configuration that reads the 768-dimensional PeptideCLM embedding of the
SMILES alongside the physicochemical and study-context blocks (817 columns in
total) was also applied to the four congeners. Its post-audit refit predicts 0.08–0.63 h — physically
plausible but short of the observed 0.45–1.90 h, and its grouped cross-validation Spearman ρ on the
178-record post-audit cohort is 0.292 against 0.445 for the selected compact configuration (0.247 and
0.574 in the 213-record pre-audit analysis). That correlation describes development-cohort performance,
not ranking among the four congeners. The top-scoring compact configuration that drops the
molecule-annotation block, physchem + study context + termini, predicts 0.45–0.79 h on the same
molecules. Neither configuration is drawn in Figure 4 as a headline, and neither is shipped.

A previous build of this analysis reported both configurations returning polymyxin
half-lives of 10⁻¹¹⁶ to 10⁻¹³⁷ h. That was a bug in the saved per-role block list, which
was written in layout order while training had concatenated the slices in
configuration-list order. Inference therefore assembled the polymyxin feature vector in a
different column order than training for any configuration whose own list disagreed with
the layout; the scaler statistics for a sparse one-hot column (mean 0.005, scale 0.068)
were then applied to polymyxin's molecular weight (1203), producing z-scores in the
thousands and ridge predictions of −400 log₁₀ h that averaged with the other members to
about −131. The shipped configuration happens to list its blocks in layout order and was
unaffected (predictions on polymyxin were correct throughout); the fix is noted here for
the record.

## Summary

MuriVerse assembles 320 rodent half-life records in which the dose route, the
bioanalytical assay and the whole-molecule structure are each established from the
primary source. **Only 213 of them are used to fit the shipped model**: the cohort splits
into four species-by-route combinations, and while pooling across species costs at most
±0.21 log₁₀ of bias, every cross-route transfer is offset by 0.45 to 0.65 log₁₀ with the
sign flip-flop kinetics predicts. A route-pooled half-life label is a mixture of an
elimination rate and an absorption rate, so the model is restricted by route and pooled
over species rather than the reverse.

On the intravenous cohort, a 71-column fixed-hyperparameter ensemble reaches Spearman
ρ = 0.574 and an area under the curve of 0.818 at the 6-hour threshold under the
strictest available holdout, meeting the pre-stated ranking and classification targets
and missing the correlation target; nested selection puts the honest ρ nearer 0.46. The
ablation shows that most of the usable signal is study context rather than molecular
structure, and that 1,280- and 768-dimensional pretrained embeddings lose to 71
hand-built columns once families and publications are held out. The within-source
analysis locates the ceiling in a source-level offset the model cannot know for a new
publication: about three quarters of the residual variance is such an offset, and two
laboratories measuring the same peptide in the same animal agree at only r = 0.75.

Applied without retraining to colistin and polymyxin B — four molecules with no
name, sequence, SMILES or detectable homology in the training set — the post-audit model
predicts 0.57–0.92 h against intravenous animal-weighted observed values of
0.45–1.90 h: long in mouse (2.05-fold for colistin, 1.26-fold for polymyxin B) and short in rat
(2.05-fold and 3.06-fold). Species is a model input, but the difference it applies is far too small
(rat about 15% shorter than mouse, against 3.3–3.6-fold longer observed), which is what the pattern
reflects. A constant drawn from the same training data is closer than the model in three of the four
cells — the mouse-median in both mouse cells and the lipidated-class median, which rests on one record,
in rat polymyxin B — and the model is closer than every constant only in rat colistin; it beats the
species-median constant in 8 of 16 predictions (every rat row), not all 16. It gets the direction of the
colistin/polymyxin B difference wrong in both species, and whether it looks accurate at all depends on
which published half-life is taken as the reference, a choice the literature does not settle. The
colistin mouse reference rests on two studies, one of whose half-lives we derived from a plotted
profile.

This is a two-molecule held-out test from outside the training chemistry, not an estimate
of general performance. The cross-validated figures of Figure 3 should not be read as
expected performance on novel chemistry, and the absolute predictions of either analysis
should not be read as half-lives.
