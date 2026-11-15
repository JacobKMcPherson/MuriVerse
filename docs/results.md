# Results

## Final post-audit performance

The primary model is the **178-record post-audit model** (150 rat, 28 mouse;
76 sources; 50 connected evaluation groups). Table 1 reports fixed-configuration
Spearman rho **0.445 (95% CI 0.153-0.563)** and nested feature-configuration selection
rho **0.297 (0.096-0.450)**. Fixed-configuration Pearson r is 0.485, RMSE 0.703
log10 h and MAE 0.550 log10 h. These are grouped cross-validation estimates, not
an independent test of the final all-record refit.

The endpoint is **source-reported circulating pharmacokinetic half-life after
intravenous dosing**, not uniformly verified terminal elimination half-life.
The post-audit polymyxin application predicts colistin/polymyxin B medians of
0.92/0.74 h in mouse and 0.78/0.62 h in rat. Fold errors are 2.05/1.26 in mouse
and 2.05/3.06 in rat; a constant baseline is closer in three of four cells.

The development-cohort results below, including Figure 3, Table 2, the 0.574/0.460
estimates and 213-record sensitivity analyses, are **historical pre-audit results**.
They must not be used as the final model's headline performance. The polymyxin
section, Figure 4 and Tables 5-7 use the post-audit models.
See [the release record](final-model-release.md) for artifacts and checksums.

Every number below is read from the tables in [`figures/tabs/`](../figures/tabs)
and the figures in [`figures/figs/`](../figures/figs), all of which are
regenerated from the cohort by `Rscript figures/code/run_all.R`. We describe
the training set (Figure 2, Supplementary Figure 1), the cross-validated
evaluation (Figure 3) and the held-out application to colistin and polymyxin B
(Figure 4).

## A structure-verified rodent half-life cohort of 320 records

The modelling unit is one source observation in
`MuriVerse/rodent/strict/halflife_meta_with_split.csv`. Applying the MuriVerse
label definition to the PEPlife 2.0 pool leaves **320 records — 74 mouse and
246 rat — covering 166 distinct ESM-2 input sequences, 129 sources and 90
combined-species HMMER homology clusters** (Figure 2b). Every record is a single peptide
chain of 3–50 residues, parenterally dosed, measured by a bioanalytical method
that a source names, and carrying a peptide SMILES for the whole dosed molecule
that has been independently decoded back to its residue sequence. Reported
elimination half-lives span 0.0028–62.0 h with a median of 0.427 h
(Figure 2c).

The selection funnel runs from 1,159 eligible records to 320 (Figure 2a). The
largest single losses are the exclusion of PEG, Fc, albumin and polymer
conjugates (−297), the requirement for a residue-level sequence (−135) and the
removal of records with another dosed partner (−139). Three cohort
requirements then take the 467 records that clear the structural gates down to
320: the dose route must be intravenous, subcutaneous, intraperitoneal or
intramuscular (−84), the structure tier must be `complete` (−57), and the assay
must be named (−6).

Species overlap is small but real: of the 166 unique input sequences, 12
occur in both animals, and 7 of the 90 homology clusters span both species
(Figure 2b: 246 records / 128 sequences / 62 clusters for rat
against 74 / 50 / 35 for mouse). Mouse and rat are therefore combined for model
development, and whole clusters are kept on one side of any split.

Peptide length has a median of 29 residues (Supplementary Figure 1j). The
route distribution is dominated by intravenous dosing — 213 intravenous, 73
subcutaneous, 22 intraperitoneal and 12 intramuscular, with **no record
reading "not reported"** (Figure 2e). Filtering also collapses the
engineering vocabulary: 4 of the 9 classes present in the eligible set
survive into the analysis set (Figure 2d), where they are
`unmodified / none stated` (189), `backbone modification only` (101),
`radio/fluorescent label only` (21) and `lipidated` (9); the 176-record
pre-audit intravenous split used for historical baseline constants holds 114, 50, 9
and 3 of them. A UMAP of mean-pooled ESM-2 embeddings over the 166 sequences
shows that long- and short-half-life peptides occupy distinguishable but
overlapping regions of sequence space (Figure 2h).

### The cohort is concentrated in a few sources but chemically varied

We provide the panels that let a reader audit each record's origin in
Supplementary Figure 1. Sources are concentrated but not dominated by one study: the ten
largest sources contribute 126 records (39% of the set) and the largest single
source, a patent, contributes 61 records (19%) (Supplementary Figure 1a). The
homology clusters used as the cross-validation unit are correspondingly
fine-grained — 90 clusters, of which 57 are singletons and the largest holds 40
records (Supplementary Figure 1b).

Chemically the set is more diverse than the "natural peptide" label suggests.
113 of 320 records carry an annotated terminal modification, overwhelmingly
C-terminal amidation (95 records) and N-terminal acetylation (12)
(Supplementary Figure 1c). 67 records are cyclic against 253 linear, and
chirality is 291 all-L, 24 mixed L/D and 5 all-D (Supplementary Figure 1d);
268 of 320 are composed only of natural L-residues (Supplementary Figure 1e).

The whole-molecule requirement is drawn over the 467 structurally kept records
rather than the 320 analysed ones, because every analysed record is `complete`
by construction: 406 of the 467 carry a SMILES for the whole dosed molecule and
61 do not (Supplementary Figure 1f). Those 61 fail for reasons the SMILES
builder records — an unconfirmed substitution (23), an unrecognised
modification (16), unspecified cyclisation (8), a conflict between the
annotation and the sequence (7), a non-standard residue that could not be drawn
(3), a backbone surrogate (1) and other causes (3) (Supplementary Figure 1g).

Two panels bear directly on the modelling result. Half-life differs
systematically by assay method, with medians spanning 0.16–1.85 h across
bioassay, LC-MS, ELISA/immunoassay, HPLC, radiometric/imaging and
radioimmunoassay (Supplementary Figure 1h) — assay is a study-context feature,
not a molecular one. And route provenance is uneven: **154 of 320 routes read
back to a primary source, with the weaker `database column` provenance
concentrated in the longest-half-life quintiles** (Supplementary Figure 1i).

## The shipped model is restricted to intravenous records (Figure 3a, b)

The cohort of 320 records is the *input* to model development, not the set the
shipped model is fitted on. **The shipped model is fitted on the 213
intravenous records** (179 rat, 34 mouse). The restriction is a result in its
own right, and we report the measurement behind it in the first two panels of
Figure 3.

The physical argument is that a route-pooled `half_life_h` is not one quantity.
After an intravenous dose the terminal half-life is an elimination half-life.
After an extravascular dose it is one only if absorption is the faster step;
when absorption is slower, the terminal phase reports the absorption rate
instead. Pooling the routes therefore gives the regressor a target that mixes
two different physical quantities.

The data agrees, and it agrees asymmetrically. Cutting the cohort into four
species × route combinations and refitting the same ensemble on four different
training pools per combination (`MuriVerse/evaluate_species_route_cells.py`;
`table36`) separates the two axes that are confounded in any marginal
comparison — rat is 73% intravenous against mouse's 46%:

| Test combination | n | Matched ρ | Pooled ρ | Cross-species ρ | Cross-route ρ | Cross-species bias | Cross-route bias |
|---|---|---|---|---|---|---|---|
| Rat IV | 179 | 0.62 | 0.57 | 0.24 | 0.20 | +0.21 | +0.65 |
| Rat EV | 67 | 0.19 | 0.39 | −0.06 | 0.44 | +0.07 | −0.45 |
| Mouse IV | 34 | 0.23 | 0.41 | 0.43 | 0.14 | −0.21 | +0.45 |
| Mouse EV | 40 | −0.19 | 0.18 | 0.18 | 0.27 | −0.18 | −0.65 |

To separate the two axes, we partitioned the 320 records into four
species × route combinations (rat IV, 179 records; rat EV, 67; mouse IV, 34;
mouse EV, 40) and refitted the same ensemble on four different training
pools per combination: *route/species-matched* (same species and same route as
the test combination; the evaluator key is `same cell`), *cross-species* (other
animal, same route), *cross-route* (same animal, other route) and
*pooled* (every record outside the test combination). Each pool was
fitted under grouped five-fold cross-validation that held out whole
connected groups of homology cluster and source identifier, so neither
family nor publication reached the scoring model, and the protocol was
repeated over five seeds because the smaller combinations are thin. All four
pools were scored on exactly the same held-out records of each combination, so
the four colours in each row of Figure 3a are directly comparable; the
mean ± 1 SD over seeds is the point and bar, with the per-seed values
as small dots. We report ranking (Spearman ρ, left), bias
(mean(predicted − observed) in log₁₀ h, middle) and evidence (−log₁₀ of
a one-sided t p-value for positive ρ, right), with sub-chance
correlations clipped to zero rather than registered as significant
(Figure 3a); because the held-out records of a grouped pass are not independent, the
evidence facet is a guide to scale rather than an exact test.

**Only rat IV was modellable on its own** (route/species-matched
ρ = 0.62). The
three smaller combinations each sat under ρ = 0.25 trained only on their
matching records, and the mouse-EV-only pool was *worse than chance* at
ρ = −0.19, with per-seed values of −0.09, −0.42, −0.36, −0.19 and +0.13
that span almost the whole plausible range for 40 records. **Pooling
every record outside the test combination raised three of the four combinations and cost rat IV almost
nothing**: mouse IV rose 0.23 → 0.41, mouse EV −0.19 → 0.18 and rat EV
0.19 → 0.39, while rat IV gave up only 0.62 → 0.57. *Transfer
behaviour was combination-specific.* For mouse IV the species axis carried
(cross-species ρ = 0.43 against cross-route 0.14), so a sufficient
volume of rat intravenous data substituted for the thin mouse
intravenous one. For rat EV and mouse EV the route axis carried (0.44
against −0.06, and 0.27 against 0.18); extravascular combinations were
therefore better served by extravascular data from the other animal
than by intravenous data from the same one. These rankings justify
combining species, and the bias measurement of the middle facet, which we
redraw against the two transfer axes, is what then forces restricting by
route (Figure 3b).

**Every cross-route transfer is offset by 0.45–0.65 log₁₀ — a factor of
3 to 4 — and every cross-species transfer stays inside ±0.21.** The signs are the
ones flip-flop kinetics predicts: IV-trained models under-predict extravascular
labels, which are absorption-limited and therefore longer, and EV-trained
models over-predict intravenous ones. That this reproduces independently in
mouse and in rat, at near-identical magnitude and with opposite signs either
way, means it cannot be an artefact of species composition.

We also measured the within-molecule version of the same result (Figure 3h):
across the 32
sequence × species combinations measured by both routes, so that the molecule
and the animal are fixed and only the route varies, SC exceeds IV in 24 of 32
pairs (75%, paired Wilcoxon p = 0.001). The paired shift is modest — median
+0.068 log₁₀, a 1.17-fold ratio — and much smaller than the marginal gap
between the route medians (IV 0.24 h against SC 0.96 h, ≈4-fold), because most
of that marginal difference is composition: the two routes were used on
different peptides.

There is a second, more prosaic reason. Route is a free lever in a route-pooled
cohort: scoring every record by its own route's median half-life, with no
molecule involved at all, reaches Spearman ρ = 0.42 against the route-pooled label.
Part of any route-pooled headline is therefore the model recovering a
study-design variable.

**What the restriction is worth, measured like for like.** On the same 213
intravenous records, the same shipped configuration and the same grouped protocol:

| | Route-pooled fit | IV-restricted fit |
|---|---|---|
| Spearman ρ | 0.507 | **0.574** |
| AUC (≥ 6 h) | 0.785 | **0.818** |
| RMSE (log₁₀ h) | 0.699 | **0.675** |
| Mouse ρ | 0.427 | **0.567** |

The ρ gain is **not on its own statistically resolvable**: a paired bootstrap
over the 55 connected cluster/source groups puts it at +0.066 [−0.016, +0.138].
The case for restricting by route is therefore the physics, the bias
measurement above, and the removal of the route lever — not the headline
correlation. A route-pooled model scores partly by recovering which route a
record came from, and its predictions have no single route at which to be read.

An earlier build of this analysis reported the gain as ρ 0.514 → 0.587 with a
bootstrap interval excluding zero. That comparison was against
`physchem + study context + termini`, the top-scoring configuration in the
ablation, which is **not** the shipped one; it is recorded here because the
distinction matters and the superseded numbers should not reappear.

## Cross-validated half-life prediction on the intravenous cohort (Figure 3c–g)

The model (`MuriVerse/halflife_model.py`) is an equal-weight ensemble of ridge
regression, an RBF support-vector regressor and a shallow XGBoost model, with
every hyperparameter fixed in advance, over **71 compact columns**: molecule
annotation (9), study context (33, including the mouse/rat indicator pair),
terminal classes (13) and RDKit
physicochemical descriptors (16). No sequence or structure embedding is used.
Under the route restriction the route columns are constant in training and
carry no information; `StandardScaler` leaves a zero-variance column at zero.
The species columns are not constant: the intravenous cohort is 34 mouse and
179 rat.

`MuriVerse/evaluate_halflife_model.py --route iv` holds out whole groups under
5-fold cross-validation, reshuffling the fold assignment for each of three
seeds; each seed is a complete out-of-fold pass and the headline numbers are
the mean of the three per-seed scores. The grouping is the **union of homology
cluster and source identifier, which yields 55 connected groups over the 213
records**, so a held-out peptide shares neither a family nor a publication with
the records it is predicted from. Intervals are 95% group bootstraps over those
groups.

### Ranking is moderate in both species; absolute accuracy is poor, especially in mouse

| Stratum | n | Spearman ρ | Pearson r | AUC (≥ 4 h) | R² | RMSE (log₁₀ h) | MAE (log₁₀ h) |
|---|---|---|---|---|---|---|---|
| Both species | 213 | 0.574 | 0.568 | 0.829 | 0.315 | 0.675 | 0.517 |
| Rat | 179 | 0.600 | 0.572 | 0.828 | 0.315 | 0.646 | 0.490 |
| Mouse | 34 | 0.562 | 0.500 | 0.787 | 0.247 | 0.811 | 0.658 |

Source: `table10_muriverse_halflife_branch_metrics.csv`.

Against the three pre-stated targets, the strict evaluation gives Spearman
ρ = 0.574 [0.322, 0.655] (target > 0.50: **met**), AUC for half-life ≥ 6 h =
0.818 (target > 0.80: **met**) and Pearson r = 0.568 [0.381, 0.651]
(target > 0.70: **not met**) (`table31c`). The AUC target, missed at 0.785 on
the route-pooled cohort, is met on the intravenous one; the correlation target
is not.

Selection optimism is larger than it was on the route-pooled cohort. Choosing the
feature configuration inside each training fold rather than on the folds that
score it gives r = 0.456 and ρ = 0.460 (95% CI 0.267–0.548; `table31b`), so roughly **0.11 of the
headline ρ is selection optimism** and the honest nested estimate is ρ ≈ 0.46.
This should be quoted alongside the headline, not instead of it: the shipped
configuration is fixed in advance of any single fold, but it was chosen on
these data.

The usable result remains a ranking, not a calibrated absolute prediction.
We obtain mouse ρ = 0.56 (seed range 0.50–0.68, n = 34) and rat ρ = 0.60
(0.55–0.65, n = 179) as three-seed means (Figure 3c). Mouse is the volatile
stratum at 34 records: it reads *above* rat at seed 0 (0.68 against 0.65) and
below it in the other two seeds, which is why the three-seed mean is the number
quoted. The per-species gap is much narrower
than it was on the route-pooled cohort, where mouse reached only ρ = 0.378 (over all
320 records; the like-for-like mouse pair on the same 213 records is 0.427 against
0.567), and the reason is visible in the composition: mouse was 54% extravascular, so the
mouse stratum of a route-pooled model was mostly records whose label was a
different physical quantity.

Ordering and absolute accuracy come apart in mouse, which is the distinction a
correlation alone hides. We find that **57% of held-out records fall
within 3-fold of the observed half-life** (mouse 41%, rat 60%), with a median
absolute fold error of 4.1× in mouse against 2.2× in rat and a signed bias of
−0.09 and −0.04 log₁₀ (Figure 3d). The model is close to centred in both animals; the mouse
penalty is spread, not offset. Read as a screen at the pre-declared 6 h target
(Figure 3e) the AUC is 0.82 (0.77–0.85 across seeds), but only 17 of the 213
records are positive (mouse 8, rat 9), so the curve is coarse and the interval
wide.

Error varies by about three-fold across descriptive strata of the same held-out
predictions (Figure 3f). By assay, radiometric/imaging records are worst at
5.8-fold (n = 13) and radioimmunoassay best at 1.7-fold (n = 90); by
engineering, radio/fluorescent-label-only peptides reach 4.9-fold (n = 17)
against 2.0-fold for unmodified ones (n = 135). The largest single source is
the best-predicted stratum in the panel (1.4-fold, n = 54) and singleton
sources the worst (4.3-fold, n = 53) — the source-offset effect seen from the
error side. A single global error or prediction interval is therefore
insufficient for pharmacokinetic decision-making.

### Assay method is the largest single input; molecular descriptors contribute jointly

The feature ablation under the strict grouping (Figure 3g;
`table31_muriverse_halflife_feature_ablation_strict_iv.csv`) scores 29
configurations:

| Feature block | Dims | ρ | r | AUC |
|---|---|---|---|---|
| physchem + study context + termini (top-scoring, **not shipped**) | 62 | 0.587 | 0.582 | 0.826 |
| annotation + context + termini + physchem (**shipped**) | 71 | 0.574 | 0.568 | 0.818 |
| physchem + study context | 49 | 0.522 | 0.512 | 0.782 |
| study context + termini | 46 | 0.456 | 0.460 | 0.772 |
| molecule annotation + study context | 42 | 0.462 | 0.431 | 0.743 |
| study context only | 33 | 0.347 | 0.323 | 0.648 |
| SMILES + physchem + study context | 817 | 0.247 | 0.278 | 0.772 |
| sequence + SMILES | 2,048 | 0.183 | 0.224 | 0.573 |
| SMILES (PeptideCLM) only | 768 | 0.177 | 0.206 | 0.712 |
| sequence (ESM-2) only | 1,280 | 0.145 | 0.176 | 0.523 |
| molecule annotation only | 9 | 0.104 | 0.166 | 0.516 |
| baseline (train mean) | 1 | −0.290 | −0.302 | 0.267 |

Three readings follow, unchanged in direction from the route-pooled analysis
and sharper in magnitude. First, **study context carries most of the signal**:
the 33 context columns alone reach ρ = 0.347 and adding the physicochemical
block reaches ρ = 0.522, whereas the 9 molecule-annotation columns alone reach
only ρ = 0.104. Second, **the large pretrained embeddings do worse, not
better**: every configuration containing the 1,280-dimensional ESM-2 or
768-dimensional PeptideCLM embedding scores below the compact set under
held-out clusters and sources while rising when only sources are held out — the
embeddings recognise families and publications, and the strict grouping removes
exactly that. Third, 71 hand-built columns beat 2,048 learned dimensions on
this cohort, which is a statement about the cohort's size, not about the
embeddings.

### The shipped configuration is not the top-scoring one

`physchem + study context + termini` outscores the shipped set by 0.013 in ρ
(0.587 against 0.574) but omits the molecule-annotation block. The shipped
configuration keeps the annotation block, enforced by `REQUIRED_SHIPPED_BLOCKS`
in `evaluate_halflife_model.py`, for a reason that rests on definability: study
context is the only block that carries the animal, and the annotation block is
what makes a prediction rest on the engineering class, topology and chirality a
curator recorded rather than on bulk physicochemistry alone. Both
configurations return physically plausible half-lives for the polymyxins
(post-audit refits: shipped 0.57–0.92 h, unconstrained winner 0.45–0.79 h over the 16 rows; the pre-audit refits gave 0.98–1.47 and 0.48–0.86 h). An
earlier build reported the unconstrained winner returning about 10⁻¹²⁰ h; that
was a block-order mismatch between training and inference, now fixed, and it is
not a reason for the constraint. The 0.013 cost is well inside the seed-to-seed
standard deviation of 0.048, and the ablation still scores and reports every
configuration — including the unconstrained winner — so the constraint's price
is visible rather than hidden.

### The ceiling is between-source disagreement

Performance degrades as the holdout gets stricter
(`table31c_muriverse_halflife_stringency.csv`):

| Evaluation | n | r | ρ | 
|---|---|---|---|
| Sources held out only | 213 | 0.614 | 0.633 |
| Clusters held out only | 213 | 0.587 | 0.585 |
| Clusters + sources held out (headline) | 213 | 0.568 | 0.574 |
| Configuration chosen inside folds | 213 | 0.456 | 0.460 |
| Within source (source mean removed) | 160 | 0.433 | 0.462 |

The last line of the table is the diagnostic. With each source's mean removed from both the
prediction and the label, the model still orders peptides inside a single
source at r = 0.433, ρ = 0.462 — appreciably less than its headline score. The
difference between the strict evaluations and this one **is** the source-level offset:
a constant per publication that the model implicitly recovers when it can see
the source's other records and cannot know for a genuinely new source. Directly,
**about three quarters of the residual variance is a constant offset per
source**: supply that offset and r rises from 0.57 to 0.93.

The literature sets the same bound from the other side. Taking every
sequence × species combination that more than one source reports, reducing each
source to its mean log₁₀ half-life for that combination and forming every
unordered pair of sources, the intravenous cohort contains **33 such
cross-source pairs, agreeing at r = 0.75 with a median disagreement of 0.35
log units against a label spread of 0.82** (SD of log₁₀ half-life). Two
laboratories measuring the same peptide in the same animal differ by about a
factor of two, which caps what any predictor can score against a single
held-out reported label. (Earlier builds quoted "66 pairs, r = 0.57,
disagreement 0.39" from a pairing rule no script in the repository generated
and no reconstruction reproduced; the rule above is executed by
`MuriVerse/audit_figure_numbers.py` and the result is reproducible.)

We make the same point from the data side (Figure 3i). The saved 80/20 split is
cluster-preserving, so whole homology clusters — and the label ranges that
travel with them — land on one side, and the training (n = 176, mean −0.53) and
validation (n = 37, mean −0.57) label densities differ. This is why the
headline evaluation reshuffles folds across three seeds rather than relying on
that single draw, and no performance number in this analysis comes from it.

The practical consequence is that the unmet Pearson target is unlikely to be
closed by more features of the kind already present. Route is now established
for all 320 records, 154 of them from a primary source, and establishing it did
not raise r — what it did was identify which records belong to one physical
quantity. What remains missing is the rest of the source-level design — dose,
strain and formulation are absent or unparsed on most records; **only one of the
213 intravenous training records carries a dose at all** — and recovering it
requires the source documents themselves, of which the retrieval ledger marks
85 of the 120 article PMIDs in this cohort as not retrieved, with 9 patent
source IDs accounting for a further 91 records.

## Held-out application to colistin and polymyxin B (Figure 4)

### The polymyxins were held out of training, and the case study is not independent validation

Nothing is retrained for the case study. `MuriVerse/predict_polymyxins.py` loads the
**post-audit** models (`MuriVerse/rodent/strict/model_iv_postaudit/`: the shipped configuration and the
two structural configurations, each fitted with the same fixed hyperparameters on the 178
source-verified intravenous records, see the clean-cohort refit below) and applies them to four
molecules — polymyxin B1, polymyxin B2, colistin A and colistin B — each built from its own PubChem
isomeric SMILES. The earlier pre-audit models (`model_iv/`, fitted on the 213 intravenous-labelled
records, which include the 37 records the primary-source audit found discordant) remain available with
`--model-suffix _iv`; the manuscript's case study uses the post-audit models.

The test is genuinely held out. `table24_polymyxin_leakage_check.csv` reports
**zero** records whose name is a polymyxin in the strict training set and zero in
the modelled cohort; the only mention anywhere in the eligible set is
colistin appearing as a comparator in the free-text dosing notes of the ZY4
record, whose half-life is ZY4's. No training peptide has detectable homology
to any congener: every surrogate query in
`table25_polymyxin_surrogate_homology.csv` reports no sequence at or below the
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
unmodified only.** An `Amidation` variant was previously run alongside it, on the
reasoning that the macrolactam bond is itself an amide; that reasoning was wrong
and the variant has been removed. A C-terminal amide is –C(=O)NH₂, whereas the
ring is Thr10–C(=O)–NH–Dab4, an isopeptide bond to a side-chain amine: both are
amides, but only one is a *C-terminal* amide, and these molecules have the ring.
Annotating amidation therefore asserted a modification the molecules do not carry
while leaving the SMILES — which already encodes the macrocycle — untouched, and
because it moved the prediction further than the drug, the congener or the
species did, it dominated every real contrast in the analysis. The cyclisation is
carried by `topology`, and `Free` is used in its sense of no terminal
modification recorded. The N-terminal variants test sensitivity to how a molecule
is written down and are **not** independent measurements — the wording changes
nothing at all, the two rows being bit-identical — so the 16 predictions carry 8
distinct values. Species is a model input, so the mouse and rat predictions for a
molecule differ; strain and sex are set to unknown and study dose is left missing.

### Against intravenous observations the model is long in mouse and short in rat

The comparison is intravenous only, because the models are fitted on intravenous records and predict
intravenous inputs. The comparison against observations from every route is secondary (Supplementary
Figure 3; the `_pooled` tables). Across the 16 shipped predictions the span is **0.57–0.92 h**; per drug
and species the median is colistin 0.92 h in mouse and 0.78 h in rat, polymyxin B 0.74 h and 0.62 h
(Figure 4b, Table 6). Within a species the four congeners spread 1.38-fold,
and that spread is not symmetric: colistin A and B agree to within 1.02-fold, while polymyxin B1 and B2
differ by 1.20-fold although the two pairs differ by the same single methylene in the acyl chain.

No median over all 16 is quoted, because the predictions split by drug (polymyxin B 0.57–0.80 h, colistin
0.77–0.92 h) and a pooled median would hide the species and drug structure the comparison is about; the
comparison is made per drug and species throughout.

The animal-weighted geometric means of the intravenous dose-group half-lives
(`table49_polymyxin_observed_by_drug_species_iv.csv`) are 0.45 h for colistin in
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
observations give polymyxin B the longer, by 1.31× in mouse and 1.19× in rat
(`table29_polymyxin_drug_discrimination.csv`). Those observed differences sit
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
constant, 2.91 for the mouse median and 8.13 for the rat median;
`table27_polymyxin_fold_error_vs_baseline_iv.csv`). Against its own species median the model is closer
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
most stable cell, 2.31–3.36-fold across deletions
(`results/revision_analysis/polymyxin_iv_leave_one_study_out.csv`).

Across every route, laid side by side (Supplementary Figure 3;
`table53_polymyxin_fold_error_by_source.csv`), the spread is wider still:

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

`table55_polymyxin_pooled_sensitivity.csv` shows that the all-route fold errors
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
animals (`table54_polymyxin_prediction_vs_dose.csv`). Over the
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
instrument that justifies the restriction has to be the unrestricted model
(`MuriVerse/evaluate_polymyxin_species_route.py`).

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
total) was also applied to the four congeners
(`table51_polymyxin_embedding_extrapolation.csv`). Its post-audit refit predicts 0.08–0.63 h — physically
plausible but short of the observed 0.45–1.90 h, and its grouped cross-validation Spearman ρ on the
178-record post-audit cohort is 0.292 against 0.445 for the selected compact configuration (0.247 and
0.574 in the 213-record pre-audit analysis). That correlation describes development-cohort performance,
not ranking among the four congeners. The top-scoring compact configuration that drops the
molecule-annotation block, `physchem + study context + termini`, predicts 0.45–0.79 h on the same
molecules. Neither configuration is drawn in Figure 4 as a headline, and neither is shipped.

A previous build of this analysis reported both configurations returning
polymyxin half-lives of 10⁻¹¹⁶ to 10⁻¹³⁷ h. That was a bug in
`predict_polymyxins.assemble` and `_blocks_of` in `evaluate_halflife_model.py`:
the saved per-role block list was in layout order while training had
concatenated the slices in configuration-list order, so inference assembled the
polymyxin feature vector in a different column order than training for any
configuration whose own ordering disagreed with the layout. The scaler stats
for a sparse one-hot column (mean 0.005, scale 0.068) were then applied to
polymyxin's molecular weight (1203), producing z-scores in the thousands and
ridge predictions of −400 log₁₀ h that averaged with the other two members to
about −131. The shipped configuration happens to list its blocks in layout
order and was unaffected (predictions on polymyxin were correct throughout);
the fix is noted here for the record.

## Summary

MuriVerse assembles 320 rodent half-life records in which the dose route, the
bioanalytical assay and the whole-molecule structure are each established from
the primary source. **Only 213 of them are used to fit the shipped model**: the
cohort splits into four species × route combinations, and while pooling across species
costs at most ±0.21 log₁₀ of bias, every cross-route transfer is offset by
0.45–0.65 log₁₀ with the sign flip-flop kinetics predicts. A route-pooled
half-life label is a mixture of an elimination rate and an absorption rate, so
the model is restricted by route and pooled over species rather than the
reverse.

On the intravenous cohort, a 71-column fixed-hyperparameter ensemble reaches
Spearman ρ = 0.574 and AUC(≥ 6 h) = 0.818 under the strictest available
holdout, meeting the pre-stated ranking and classification targets and missing
the correlation target; nested feature-configuration selection puts ρ at 0.460
(0.267–0.548). Assay method has the largest individual permutation importance,
while molecular descriptors contribute jointly (study context alone ρ = 0.347;
with the physicochemical descriptors 0.522), and 1,280- and 768-dimensional
frozen pretrained embeddings lose to 71 hand-built columns once families and
publications are held out. The
within-source analysis locates the ceiling in a source-level offset the model
cannot know for a new publication: about three quarters of the residual
variance is such an offset, and two laboratories measuring the same peptide in
the same animal agree at only r = 0.75.

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

This is a two-molecule held-out test from outside the training chemistry, not
an estimate of general performance. The cross-validated figures of Figure 3
should not be read as expected performance on novel chemistry, and the absolute
predictions of either analysis should not be read as half-lives.

## Revision analyses (external review, round 1)

These numbers match the revised manuscript; scripts are in
`MuriVerse/revision_analysis/` and outputs in `figures/tabs/` and
`results/revision_analysis/`.

- **Primary-source audit of the 213 intravenous-labelled records**
  (`MuriVerse/data/curation/halflife_label_verification.csv`; Supplementary
  Table 4): 170 verified in full primary-source text, 3 from the abstract only,
  3 with a confirmed value but no stated intravenous route, 9 with a database
  value not supported by the source and 28 contradicted by it (37 discordant).
  The audit confirmed an intravenous route for 197 records, including all 170.
  The 213-record analysis is retained as the original/pre-audit benchmark.
- **Fixed-configuration sensitivity analyses** (Supplementary Table 5,
  Supplementary Figure 2e,f): Spearman ρ 0.540 (0.282–0.617; n = 176) without the
  discordant records; 0.510 (n = 179) with three corrected values restored;
  0.541 (n = 173) full text or abstract; 0.474 (0.242–0.561; n = 170) full text;
  0.442 (n = 117) full text without the largest patent; 0.535 (n = 203) exact
  labels. The configuration was not reselected inside these subsets, so they do
  not remove the configuration-selection optimism.
- **Post-audit model** (`MuriVerse/rodent/strict/model_iv_postaudit/`, the
  distributed final model): 178 records (170 full text, 3 abstract, 3 corrected
  values, 2 species corrections). Fixed configuration ρ 0.445 (0.153–0.563), MAE
  0.550; nested feature-configuration selection ρ 0.297 (0.096–0.450). The
  pre-audit model's out-of-fold predictions on the same 178 records give ρ 0.579.
  The same fixed hyperparameters were also refit on the 178 records for the two structural
  configurations; grouped cross-validation ρ is 0.456 (physicochemical-descriptor-best) and 0.292
  (SMILES-embedding-best). The polymyxin case study (Figure 4, Tables 5–7) applies these post-audit
  models: median predictions 0.92 h (colistin) and 0.74 h (polymyxin B) in mouse, 0.78 and 0.62 h in rat;
  fold errors against intravenous references 2.05, 1.26, 2.05 and 3.06; a constant was closer than the
  model in three of the four cells. The pre-audit models remain available with
  `predict_polymyxins.py --model-suffix _iv`.
- **Ensemble against its members** (Supplementary Table 3): MAE 0.517
  (ensemble), 0.529 (ridge), 0.570 (SVR), 0.592 (XGBoost), 0.702 (training-fold
  species median). The ensemble had the lowest point-estimate error of the five
  estimators, but an advantage over ridge regression alone was not established
  (paired difference −0.012; −0.086 to +0.065).
- **Ranking in thirds** (`results/revision_analysis/ranking_utility.csv`): 59%
  (36–71%) of the observed top third fall in the predicted top third (33% by
  chance); 15% (10–28%) of the predicted top third come from the observed bottom
  third.
- **Permutation importance** (Supplementary Table 2): assay method 0.074
  (0.046–0.101) log₁₀ h increase in MAE, strain 0.013; no single physicochemical
  descriptor above 0.011.
- **Short peptides**: joining records with an identical canonical sequence or
  structure changed no group (55 before and after); without the ten records of
  six residues or fewer, ρ was 0.566 (0.194–0.670; n = 203).
- **Bootstrap**: each seed is scored separately in every resample and the three
  scores averaged, so the interval matches the reported three-seed mean.

## Figures and tables

| Item | File | Caption |
|---|---|---|
| Figure 2 | [`figures/figs/fig2.pdf`](../figures/figs/fig2.pdf) | [`fig2_caption.md`](../figures/figs/fig2_caption.md) |
| Figure 3 | [`figures/figs/fig3.pdf`](../figures/figs/fig3.pdf) | [`fig3_caption.md`](../figures/figs/fig3_caption.md) |
| Figure 4 | [`figures/figs/fig4.pdf`](../figures/figs/fig4.pdf) | [`fig4_caption.md`](../figures/figs/fig4_caption.md) |
| Supplementary Figure 1 | [`figures/figs/supp_fig1.pdf`](../figures/figs/supp_fig1.pdf) | [`supp_fig1_caption.md`](../figures/figs/supp_fig1_caption.md) |
| Supplementary Figure 3 (all-route comparison) | [`figures/figs/supp_fig3.pdf`](../figures/figs/supp_fig3.pdf) | [`supp_fig3_caption.md`](../figures/figs/supp_fig3_caption.md) |

Principal source tables, all in [`figures/tabs/`](../figures/tabs):
`table10_muriverse_halflife_branch_metrics.csv`,
`table10b_muriverse_halflife_ranking_metrics.csv`,
`table31_muriverse_halflife_feature_ablation_strict_iv.csv`,
`table31b_muriverse_halflife_nested_selection_iv.csv`,
`table36_muriverse_halflife_species_route_cells.csv`,
`table32_muriverse_halflife_route_flipflop_pairs.csv`,
`table31c_muriverse_halflife_stringency.csv`,
`table45_muriverse_performance_by_structure_tier.csv`,
`table24_polymyxin_leakage_check.csv`,
`table25_polymyxin_surrogate_homology.csv`,
`table27_polymyxin_fold_error_vs_baseline_iv.csv` (and the all-route `table27_...csv`),
`table29_polymyxin_drug_discrimination.csv`,
`table49_polymyxin_observed_by_drug_species_iv.csv` (and the all-route `table49_...csv`),
`table51_polymyxin_embedding_extrapolation.csv`,
`table53_polymyxin_fold_error_by_source.csv`,
`table54_polymyxin_prediction_vs_dose.csv`,
`table55_polymyxin_pooled_sensitivity.csv`.
