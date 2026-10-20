# MuriVerse strict rodent peptide half-life set

> **Historical peptide-scale build:** an earlier tier kept 1,044 cohort records and 977
> modelling rows, including 219 without verified SMILES. That cohort is no longer active;
> the current modelling input is the peptide-only, structure-verified set below.

> **Current build (peptide-only):** the eligible set contributes 1,159 candidate records.
> A record is kept only if it is a single peptide chain of at most 50 residues that parses
> residue by residue, carries no attached half-life partner (no PEG, Fc, albumin, polymer,
> particle, dye, drug or protein partner), reports an elimination half-life, and has a SMILES
> that decodes back to its sequence. Every source observation is retained, including repeated
> treatment arms and independent studies of the same molecule. Possible duplicate entries are
> flagged with their source-arm descriptors for review. It must also describe a normal
> circulating-PK setting: records in organ-ablated, knockout, transgenic, hypertensive, diabetic or
> tumour-bearing animals, and tissue, lysate or lymph samples, are excluded by the
> `physiological circulating half-life` gate. **467 records remain: 111 mouse and
> 356 rat.** The stage counts are in `filtering_counts.csv`; every excluded record, with its
> gate and the text that triggered it, is in `exclusions.csv`. Three further cohort
> requirements - a parenteral dose route, a whole-molecule SMILES and a named bioanalytical
> method - take these 467 kept records to the **320** the model is fitted on; all three are
> applied in `build_halflife_training_set.py` and described under "Why this set exists"
> below.

| stage | records | removed |
|---|---|---|
| eligible rodent records | 1159 | - |
| not curated out | 1152 | 7 |
| physiological circulating half-life (normal animal; plasma, serum or blood) | 1103 | 49 |
| no engineered partner class (PEG, polymer, Fc, albumin, depot/particle) | 806 | 297 |
| single peptide chain (no repeat or polymer notation) | 804 | 2 |
| faithful residue sequence | 669 | 135 |
| peptide scale (<= 50 aa) | 636 | 33 |
| no other attached partner | 497 | 139 |
| elimination half-life (no absorption, distribution or activity phase; no magnitude cutoff) | 468 | 29 |
| structure verified | 467 | 1 |
| observations retained; repeat candidates flagged | 467 | 0 |
| dose route curated from the sources | 467 | 0 |

Half-life values: ranges are the median of their endpoints, approximate (`~2 h`) and
one-sided (`<5 min`) values keep their stated number, and `half_life_qualifier` records the
source qualifier. Valid extreme labels are retained and carry `half_life_extreme_qc`; that
QC flag is not a model feature. "Hours-1" is a source typo for hours; those rows are
absorption half-lives and are excluded by phase.

467 mouse and rat in-vivo pharmacokinetic half-life records with a sequence-verified peptide
SMILES, each a chain of 50 residues or fewer. The peptide sequence, topology and chirality
are independently checked; reported details that are unspecified or could not be drawn are
recorded in `structure_tier` and `unresolved`. Built by
`MuriVerse/build_strict_peptide_set.py` from the pooled analysis set in
`../in_vivo_half_life.csv`; nothing here is hand-edited.

| file | contents |
|---|---|
| `in_vivo_half_life.csv` | the current 467 kept source observations, with verified structure, modification, repeat-candidate and extreme-label QC fields, and the source-curated `route_class` with its `route_evidence_level` |
| `exclusions.csv` | every one of the 692 gate-excluded source records, with the gate that removed it and the text that triggered it |
| `filtering_counts.csv` | the stage-by-stage funnel below, machine-readable |
| `repeat_observation_candidates.csv` | within-study duplicate candidates and same-structure cross-source observations, preserved for evidence review |
| `../pk_observations.csv` | arm-level full-text PK evidence keyed to PEPlife2 IDs; conditions are joined only when the ID is explicit |
| `halflife_meta_with_split.csv` | the 320 modelling rows - kept records that are parenterally dosed, whose SMILES is the whole molecule, and whose assay a source names - with HMMER cluster and train/validation split |
| `oof_predictions.csv` | one out-of-fold prediction per modelling row under cluster, PMID and combined holdouts; this is what Figure 3 plots |
| `hmmer/` | the `phmmer` all-vs-all search that produced those clusters |

## Why this set exists

The eligible set is a faithful record of what PEPlife2 reports, and that is why it is the
wrong thing to train on. Over its 1,159 records the half-life-extension class moves the
median by a factor of 93:

| `hl_engineering` | n | median t½ (h) | median length |
|---|---|---|---|
| unmodified / none stated | 430 | 0.35 | 28 |
| backbone modification only | 247 | 0.61 | 37 |
| radio/fluorescent label only | 122 | 0.86 | 26 |
| pegylated | 141 | 5.25 | 30 |
| lipidated | 60 | 6.92 | 32 |
| polymer conjugate | 74 | 9.05 | 39 |
| albumin-based | 31 | 10.0 | 35 |
| depot / particle formulation | 29 | 11.3 | 13 |
| **fc fusion** | **25** | **32.3** | **46** |

An Fc fusion's clearance is the partner's rather than the pharmacophore's, so a model fitted
on the eligible set learns the partner. (The lengths above understate the point: the eligible
set applies a generous 80-letter pre-filter, so the longest fusions are already gone.) Worse,
for those records no SMILES is true — the pipeline deliberately draws the peptide only and
passes partner identity and mass as separate features — so the molecular encoder reads a
string that is not the molecule.

The rule here is instead: **a peptide sequence and a verified peptide structure.** A record
is kept only when the sequence column is a real sequence of template-bearing residues and
the emitted SMILES independently decodes to that peptide, including its resolved topology
and chirality. The builder also checks known terminal and `Chem-MOD` details; an unspecified
feature or a reported modification that cannot be drawn is retained in the quality tier and
unresolved annotation rather than being presented as complete.

That admits acetylation and amidation, head-to-tail and disulfide macrocycles, D and mixed
series, side-chain lipidation with its γ-glutamyl and OEG spacers, and exactly-structured
covalent labels. PEG conjugates are excluded at every size, along with Fc fusions, albumin
binders, other polymer conjugates, depot and particle formulations, glycosylation, drug
payloads and proprietary dyes.

## The funnel

The funnel above is the current one. The strict set has **467 records** from **184 source
identifiers**, including 11 patent identifiers (356 rat, 111 mouse). `exclusions.csv` lists
the 692 excluded records of the eligible set, each with its gate.

Not every kept record is a modelling row. The concentration-decay phase filter and the
feature-level de-duplication that earlier builds applied still remove none, but two cohort
requirements apply in `build_halflife_training_set.py`:

| requirement | removed | why |
|---|---|---|
| parenteral dose route (iv, sc, ip, im) | 84 | an enteral or mucosal half-life is absorption-limited, and under flip-flop kinetics its terminal slope is the absorption rate rather than elimination; records whose route no source establishes go with them, because nothing says which of the two they are |
| whole-molecule structure (`structure_tier == complete`) | 57 | a `partly resolved` or `modification not drawn` SMILES is not the molecule the label belongs to, and the molecular encoder cannot tell the difference |
| named bioanalytical method (`assay_method != unknown`) | 6 | what quantity decayed is not a detail of a half-life, it IS the half-life; these six had a median label of 23 h against the cohort's 0.43 h and held three of the four longest values in the set |

That leaves **320 modelling rows** from **129 source identifiers** (246 rat, 74 mouse;
half-life 0.0028-62.0 h, median 0.427 h; 166 distinct ESM-2 input sequences in 90 HMMER
clusters; 3-50 residues, median 29; 291 all-L, 24 mixed-chirality, 5 all-D). The route
requirement is the largest of the three and is the one that most changes the cohort's
character: 213 of the 320 are intravenous, 73 subcutaneous, 22 intraperitoneal and 12
intramuscular, and none now reads "not reported".

The 166 sequence count refers to unique normalized ESM-2 inputs in the model metadata.
Uppercasing chirality-marked residues and mapping unsupported residues to `X` collapse some
distinct source strings to the same ESM-2 input, so the source export holds more.
`table34_strict_cohort.csv` reports both counts, while `table42_strict_modelling_set.csv`
reports the model-input count. Species-level counts in table 34 overlap: five source IDs and
13 source sequence strings occur in both species.

## The 50-residue cap

The peptide-scale gate is **50 residues**. It removes 33 eligible records in the current
funnel. The cap defines this dataset's peptide-scale scope; it is applied to the verified
residue sequence before structure verification.

## The other sequence column

PEPlife2 ships **two** sequence columns and they are not the same data. The pipeline had been
reading the wrong one.

`Sequence` is the original. It marks a D residue by lower-casing its letter (`rPKPfQwFwLL`)
and delimits a non-natural residue with hyphens (`H-Aib-QGTFTSDYSKYL`). `sequence` is a
flattened derivative: upper-cased and punctuation-stripped, which destroys the D/L pattern
and can obscure non-natural residue notation. The original case-bearing column is used for
the verified sequence.

Reading the original makes per-residue mixed chirality representable. The current strict set
has 410 all-L, 52 mixed-chirality and 8 all-D records. Sequence case determines chirality
where available, and the round-trip checked result is the canonical `chirality` feature.
Three Apraglutide rows disagree with PEPlife2's `Chirality` annotation; the source value is
preserved in `chirality_reported`, while the sequence-verified value is used by the model and
Figure 2. `stereo_source` records the evidence used for each row.

Delimited sequences are parsed conservatively, because case cannot be trusted to mean D
there: `pIPhe` is **one** residue whose letters p, I, P, h and e are each individually a
valid code, so expanding it would invent five residues from one. A delimited sequence is
accepted only when every chunk is a known residue code, a single letter, or an all-upper-case
run; anything else is refused. An `N-` prefix is also refused outright, since splitting
`N-MeVal-EtPro` on hyphens yields a bare `N` that would be read as asparagine.

Where the flattened column disagrees with an upper-cased reading of the original, the record
is refused: two columns of one database contradicting each other about which residues a
peptide contains is not a basis for drawing it.

### Non-standard residues

The builder templates Aib, Orn, Nle, Nva, Abu, Dab, Dpr, hArg, Cha, Nal and Tys — each a
plain side-chain substitution, so one table entry makes the residue both buildable and
decodable. 41 kept records contain at least one. **ESM-2 has no token for any of them**, so those
positions reach the sequence branch as `X`; `nonstandard_positions` and `n_nonstandard`
record exactly where, so the resulting quality gradient is measurable rather than hidden.

Residues that alter the *backbone* rather than the side chain — Sar, N-methyl-anything, Hyp,
Aze, Oic, pGlu, Iva — are still refused, and say so by name.

## Verification

Every kept record passes two independent checks, because the value of this set is that its
strings can be trusted:

* the **bare** peptide — the molecule without its side-chain conjugates — is decoded back by
  `reconstruct_residues`, an independent SMILES→sequence reader, and must return the recorded
  residues (up to rotation for a head-to-tail macrocycle). This checks residue identity,
  order, topology **and configuration**: the decoder reads each alpha carbon's CIP descriptor
  out of the string and maps it to L or D (every residue is S in the natural series except
  cysteine, which is R — established empirically, not assumed). A constitution-only check
  would pass an all-L string for a mixed-chirality peptide, which is the exact error this set
  exists to exclude.
* the **full** molecule's exact mass minus the bare peptide's must equal the sum of the
  attached fragments' masses, to 0.01 Da. This checks that each conjugate went on once, at one
  site, losing exactly the one hydrogen the bond costs.

A record failing either is excluded, not corrected.

## Two sources assessed, then removed

Both PEPlife2's PDB structure archive and PepTherDia 2.0 were wired into this pipeline,
measured, and taken back out. Neither is used now, and `docs/data-provenance-quality.md`
records the measurements that decided it. In short:

* the structure archive's ID-to-structure mapping contradicts the sequence column on 171 of
  394 checkable records, its coordinates do not encode stereochemistry at all (records
  annotated `D` come back entirely S, because PEPstrMOD built everything in the natural
  series), and 399 of its 1,837 files are experimental PDB entries of parent proteins and
  complexes rather than of the dosed peptide. A gate that excluded records for disagreeing
  with a source that unreliable was not defensible.
* PepTherDia 2.0 and ThPDB 2.0 cannot supply rodent half-life records - PepTherDia has zero
  rodent mentions anywhere, ThPDB has 5 of 3,380, duplicates of ~2 unique values - so both
  are a species mismatch for this label. While PepTherDia was briefly used as an external
  check on the builder, 19 of its 24 comparable approved drugs came out as the identical
  molecule by InChIKey (oxytocin, glucagon, exenatide, lixisenatide and angiotensin II all
  with their literature keys), and of the 5 exceptions four were omissions or errors in the
  reference. That is reassuring about the builder but it is not a reason to keep a human
  clinical dataset in a rodent pipeline, so the check is gone too.

## Judgement calls worth knowing about

**Radiotracers are treated as structurally inert.** 40 kept records carry one, flagged in
`tracer_label` (`isotope`, `radiohalogen`, `radiometal`, `radiolabel`). For ³H, ¹⁴C, ¹³C and
¹⁵N this is simply true — an isotope changes mass, not constitution. For radioiodination it is
a deliberate choice about the dosed material rather than about drawing difficulty: tracer
radioiodination runs far below one iodine per peptide, so most administered molecules carry no
iodine, and drawing it would describe the minority species. A radiometal named *without* its
chelator is excluded instead, because the chelator is a real covalent adduct of a few hundred
daltons whose absence from the text means the structure is unknown.

**No PEG, in any form.** An earlier revision admitted a discrete, countable oligomer — and
`PEG12` genuinely is an exact molecule — while excluding only the polydisperse kDa-scale
material. In this cohort that distinction bought exactly **one** record out of 134 PEGylated
ones, because 42 state no size at all and the rest are 5–46 kDa chains, frequently at an
unstated or side-chain site, so no single SMILES is true for them. One surviving record was
not worth a policy that needed a caveat to state, so PEG is now excluded outright (20
records blocked as `PEG conjugate`) and the oligomer builder is gone with it. Zero kept
records mention PEG anywhere in `Name`, `Origin`, `Chem-MOD` or either terminal column.

Note this does **not** touch the OEG / `PEG2` *spacer unit* inside a fatty-acyl recipe:
8-amino-3,6-dioxaoctanoyl is a single residue-sized linker with an exact structure, reached
only during conjugate assembly, and it never stands alone as a polymer conjugate.

**Substitution claims are recorded and reviewed.** When a source describes a residue
substitution, the builder compares the claim with the encoded sequence. Clauses that cannot
be confirmed are retained in `substitution_unconfirmed`; they are not a separate exclusion
gate. The current strict input contains 19 such records. These rows still pass sequence
parsing and structure round-trip checks; interpret the annotation as source uncertainty rather
than proof that the named substitution was built.

The earlier hand-exclusion example for PEPlife2 ID 5555 is historical: that ID is not in the
current strict set. Active exclusions and their reasons are listed in `exclusions.csv`.

**The `Name` blind spot.** One record, `Au‐AR pep‐PROTAC`, is a gold-nanoparticle conjugate
with a 26.3 h half-life — roughly 100× this set's median — whose every annotation column reads
`Free` and `unmodified`. The carrier exists only in `Name`, which no other gate reads, and
`hl_engineering` misses it too because its patterns have no term for a gold carrier. `Name` is
therefore screened against the macromolecular vocabulary. `Origin` is deliberately **not**
screened: it is provenance prose about where the sequence came from (`Derived from the tenth
domain of fibronectin type III`), which describes nothing attached to the molecule.

**13 side-chain conjugates are drawn**, at a lysine Nζ or a serine oxygen: fatty-acyl and
γ-glutamyl lipidations, a C16-diacid, and an octanoyl ghrelin at S3. Conjugate recipes with two different
spacer types are excluded rather than guessed: PEPlife2 writes
`N(eps)hexadecanedioyl-gGlu-2xOEG` distal-first and `(ADO-linker-Glu-C18)` proximal-first, and
nothing in the text distinguishes them.

## What it scores

The 320 modelling rows are the INPUT to model development. **The shipped model is fitted on
the 213 intravenous rows of them** (179 rat, 34 mouse): only an intravenous dose has no
absorption phase, so a route-pooled `half_life_h` mixes an elimination rate with an absorption
rate, and every cross-route transfer in the four species x route cells is offset by 0.45-0.65
log10 against +/-0.21 for every cross-species one. The shipped model averages ridge, RBF-SVR
and shallow XGBoost predictions over 69 annotation and physicochemical features. The underlying
Arrow dataset also contains sequence and molecular embeddings; the shipped feature set is
selected and evaluated with grouped holdouts. In cluster+PMID holdout evaluation over 55
connected groups, Pearson `r=0.568`, Spearman `ρ=0.574` and AUC(>=6 h) `0.818` (five folds,
three seeds), with `r=0.576`, `ρ=0.595` in rat against `r=0.533`, `ρ=0.567` in mouse - the
mouse branch is 34 rows and its interval is correspondingly wide.

Measured like for like on the same 213 rows with the same configuration, restricting the route
takes `ρ` from 0.507 to 0.574 and AUC from 0.785 to 0.818, though the `ρ` gain alone sits inside
its own paired bootstrap interval (+0.063 [-0.005, +0.157]); the case for the restriction is
the physics and the cross-route bias, not the correlation.

Choosing the feature configuration inside each training fold rather than on the folds that
score it gives `r=0.470` and `ρ=0.461`, so roughly 0.12 of the headline is selection optimism.

The route-pooled model is still built and kept on disk (`model/`, `oof_predictions.csv`,
`table31..._strict.csv`) because the transfer diagnostics need a model that has seen both
routes; drop `--route iv` to reproduce it.

The route-pooled 320-row scores these replace were `r=0.574`, `ρ=0.579`, AUC(>=6 h) `0.750`,
and before that a 470-row build gave `r=0.597`, `ρ=0.603`. The drop from 470 to 320 is
the expected direction. The cleaning removed three covariate levels - an unreported route, an
incomplete structure and an unnamed assay - that each marked a particular subset of sources,
and a grouped-holdout model can lean on a source-identifying level the way it leans on any
other. What is left is a cohort in which every row's route, structure and assay is
established, so the estimate has less to lean on and is a better guide to a new peptide. See the [Figure 3 caption](../../figures/figs/fig3_caption.md) and
`figures/tabs/table31c_muriverse_halflife_stringency.csv` for current results and uncertainty.

## Rebuilding

```bash
conda run -n peptiverse-local python MuriVerse/build_strict_peptide_set.py
conda run -n peptiverse-local python MuriVerse/build_halflife_training_set.py --tier strict
conda run -n peptiverse-local python MuriVerse/evaluate_halflife_model.py --tier strict
```

`--tier strict` is the default for both of the latter two. `--tier pooled` builds an optional
comparison cohort; the current final model and Figure 2 use the strict 320-row cohort. The
shipped feature vector includes a 13-dimensional `termini` block alongside molecule
annotation, study context and physicochemical features.
