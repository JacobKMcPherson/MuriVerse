# MuriVerse

> **Current training data:** The [mouse and rat training index](rodent/strict/training_data/README.md)
> describes the current **320-row**, peptide-only species-pooled model (74 mouse, 246 rat): single
> peptide chains of at most 50 residues with no attached half-life partner, an elimination
> half-life from a normal animal's circulation after a **parenteral** dose, measured by a
> **named bioanalytical method**, and a SMILES that is the **whole dosed molecule**, decoded
> back to its sequence, for every row.

Rodent in-vivo peptide half-life: a curated cohort, a sequence-verified peptide training set,
and one model. It builds on the encoders and conventions of
[PeptiVerse](https://github.com/synbioai-au/Peptiverse), which is not part of this repository.

PeptiVerse predicts peptide properties from sequence and SMILES across nine tasks. Its
half-life task is trained on a pooled, species-agnostic set with no route, matrix or PK-phase
annotation. MuriVerse narrows that to **mouse and rat in-vivo pharmacokinetics**, carries the
provenance of every value back to its paper or patent source, and trains on peptide
structures that have been decoded back to their residue sequence. Structure completeness is
tracked per row, and the model table requires the complete tier, so a verified peptide SMILES
is never mistaken for a complete representation of every reported modification. The dose route
is curated from the primary sources rather than taken from the database column alone.

Nothing in this directory is hand-edited. Every CSV is the output of a script in it, and every
figure and table under [`../figures`](../figures) is regenerated from those CSVs.

---

## The short version

```bash
# one-time setup: Python environment (see Environment) and the PeptideCLM SMILES tokenizer
bash scripts/fetch_tokenizer.sh
VENV=python   # the interpreter of the environment that has torch, transformers, xgboost, scikit-learn, rdkit and datasets

$VENV MuriVerse/build_strict_peptide_set.py                              # -> rodent/strict/    (467 records)
$VENV MuriVerse/build_halflife_training_set.py --tier strict             # -> 320 modelling rows + embeddings
$VENV MuriVerse/evaluate_halflife_model.py --cohort rodent --tier strict --jobs 16 --nested  # ~25 min -> oof_predictions.csv + ablation + saved models
Rscript figures/code/run_all.R                                           # -> every figure and table
```

Those four commands rebuild everything downstream of the curated cohort. The steps that
*build* the cohort (PEPlife2 download, full-text retrieval, curation) are slow and network-
bound, and their outputs are committed — see [Stage 1](#stage-1-cohort) if you need to rerun
them.

---

## Layout

| path | what it holds |
|---|---|
| `*.py` | every pipeline step; one script per stage, run top to bottom |
| `data/` | the raw sources: the PEPlife2 record table, its mirrored PDFs, PubMed metadata |
| `mouse/`, `rat/` | per-species **cohorts**: provenance, retrieved full text, curated extractions, and the in-vivo half-life analysis set |
| `rodent/` | the two species **pooled**: the analysis set, its SMILES audit, and the pooled modelling artifacts |
| `rodent/strict/` | the **strict set** — what the model actually trains on. Has its own [README](rodent/strict/README.md) |
| `EXTRACTION_PROTOCOL.md` | the rules a curator follows when reading a paper into an extraction JSON |

Large or re-downloadable things are gitignored: `data/peplife2/papers/`, the site export
zips, `*/pk_figures/`, and the `*_embedding_cache.npz` memoisation files.

---

## The pipeline

### Stage 1: cohort

Build the record table, select a species, fetch and render the full text. Slow, network-bound,
and its outputs are committed — you do not need to rerun this to retrain.

| script | does |
|---|---|
| `build_peplife2.py` | Step 0. Builds the PEPlife 2.0 record table from the database's own REST API → `data/peplife2/peplife2_records.csv` |
| `fetch_peplife2_api.py --species mouse\|rat` | Selects a cohort by PEPlife2's own organism annotation (`Mouse` ∪ `Mice`, `Rat` ∪ `Rats`) |
| `build_rodent_cohort.py` | Step 1. Assembles the species cohort and fetches PMC full text per PMID |
| `fetch_pubmed_metadata.py` | Resolves each PMID to title, journal, year, DOI, PMCID, publication type |
| `fetch_peplife2_papers.py` | Mirrors PEPlife2's own PDF per publication |
| `render_fulltext.py` | Step 2. Renders retrieved BioC full text as body-only plain text |
| `render_pdf_text.py` | The same, for publications absent from the PMC OA subset |
| `list_missing_pdfs.py` | Worklist of dosed studies whose PDF could not be mirrored |
| `build_study_dossier.py` | Per-record provenance + per-study dossier → `<species>/records_provenance.csv`, `studies.csv`, `STUDIES.md` |
| `subset_by_nuance.py` | Step 4. Collapses each curated extraction JSON to one row of paper-level design flags |

Curation itself is manual: a curator reads the rendered text and writes
`<species>/extractions/<PMID>.json` following `EXTRACTION_PROTOCOL.md`. 48 of the strict set's
source identifiers have an extraction JSON.

### Stage 2: analysis set

| script | does |
|---|---|
| `subset_in_vivo_half_life.py` | Narrows `records_provenance.csv` to live-dosed-animal circulating half-lives → `<species>/in_vivo_half_life.csv` |
| `build_rodent_eligible_set.py` | Selects the rodent in-vivo peptide records from the sample text → `rodent/eligible_in_vivo_half_life.csv` (**1,159 records**) and `rodent/eligible_filtering_counts.csv` |

`rodent/eligible_in_vivo_half_life.csv` is the **eligible set**: every mouse or rat record
dosed in vivo at peptide scale. It is the input to the strict gates. We draw that selection
funnel in Figure 2a.

### Stage 3: the strict set

```bash
python MuriVerse/build_strict_peptide_set.py
```

One rule: **a peptide sequence and a verified peptide structure.** A record is kept only when
the sequence column holds real template-bearing residues and the chain is **50 residues or
fewer**. The emitted peptide SMILES is independently decoded to check residue order, topology
and stereochemistry. Reported but unspecified features and modifications that could not be
drawn remain visible in `structure_tier` and `unresolved`; they are not described as complete.
Where the builder adds a conjugate, its mass delta is independently checked against the
fragment masses. Three records have a disagreement between the source's summary Chirality
field and the stereochemistry encoded in sequence case; the strict model uses the verified
sequence-derived chirality and retains the reported value in `chirality_reported`.

1,159 eligible → **467 kept**, 692 gate-excluded (each listed in `exclusions.csv` with the gate and the
verbatim text that triggered it). A `physiological circulating half-life` gate removes records
in organ-ablated, knockout, transgenic, hypertensive, diabetic or tumour-bearing animals and in
tissue, lysate or lymph samples (49 records); the engineered-partner gate is the largest (297). The full rationale, the funnel and the
residue tables are in [`rodent/strict/README.md`](rodent/strict/README.md).

#### Dose route, curated from the sources

PEPlife2 leaves `Dose Route` blank on half the cohort, so `route_class` read "not reported"
for 241 of the modelling rows even where the source states the route plainly. Those records
are resolved against the primary sources themselves - full text where it is reachable, patent
text, the abstract, and NLM's MeSH indexing - and the verdict for every source is recorded in
[`data/curation/route_curation.csv`](data/curation/route_curation.csv) with the verbatim quote
and where it came from. `halflife_condition_curation.py` reconciles that ledger with the
record and writes `route_class` plus a `route_evidence_level`:

| level | meaning | records |
|---|---|---|
| `primary_arm` | the source names the route for THIS arm - a single-route study, or a half-life matched to its own row in the source's PK table | 66 |
| `primary_source` | the source names the routes it used, all parenteral; the record is attributed at source level | 106 |
| `database_route_column` | the route PEPlife2 transcribed; no source was re-read for it | 227 |
| `none` | no primary source names a route for the PK arm | 69 |

The model table then keeps **iv, sc, ip and im only**. An oral, intranasal or rectal half-life
is an absorption-limited composite - under flip-flop kinetics the terminal slope is the
absorption rate, not elimination - so it answers a different question from the one the model
is fitted on. A record whose route no source establishes is dropped for the same reason:
nothing says which of the two it is. 30 of the 84 sources behind the blank records could not
be resolved, most because the full text is paywalled and neither the abstract nor the MeSH
indexing names a route; each is listed in the ledger with why. Supplying those PDFs and
re-running the stage is all that is needed to recover them.

#### Whole-molecule structure

Every kept record has a verified SMILES, but not every SMILES is the whole dosed molecule:
`partly resolved` leaves a reported ring or configuration unspecified, and
`modification not drawn` emits the peptide WITHOUT a clause that alters it. Both hand the
molecular encoder a string that is not the molecule the label belongs to, so the model table
requires `structure_tier == complete`. `build_strict_peptide_set.py` resolves what can be
resolved first - 18 of the 81 incomplete records were recovered by extending the builder's
vocabulary rather than by discarding them - and what remains is listed with its reason in the
kept file's `unresolved` column and drawn in Figure 2's panels n and o.

#### Assay method

`assay_method` is bucketed from the free-text `Assay` column, which for some records names a
separation rather than a detector ("TLC"), the readout rather than the method ("Serum decay
experiment"), or nothing at all. The vocabulary in `study_conditions.py` now covers the named
methods it was missing (EIA, LC/MS, LCHRMS, AlphaLISA, western blot) and carries `bioassay` as
its own level, for records whose plasma concentration was read off a cell-based reporter or an
inhibition titre rather than measured physically. What the column still cannot say is taken
from the sources and recorded in
[`data/curation/assay_curation.csv`](data/curation/assay_curation.csv). `unknown` is left
standing where the `Assay` column and the sources are both silent, and those records are then
dropped: what quantity decayed is not a detail of a half-life, it IS the half-life, so a record
whose method nothing names does not say which quantity it measured. The six that reached that
state had a median half-life of 23 h against 0.43 h for the cohort, and three of them were the
three longest labels in it, so keeping them would have put the top of the label range behind
the one covariate that cannot be checked. 33 of the original 470 rows read `unknown`; 27 were
resolved and 6 removed.

Two modules do the chemistry and are used by this stage and by `build_rodent_smiles.py`:

| module | does |
|---|---|
| `peptide_smiles.py` | Builds a SMILES for a modified peptide: residues, topology, chirality, termini; and decodes one back to a sequence |
| `peptide_mods.py` | Classifies PEPlife2 `Chem-MOD` free text and turns the drawable part into SMILES fragments |

### Stage 4: features and the model

```bash
python MuriVerse/build_halflife_training_set.py --tier strict   # default tier
python MuriVerse/evaluate_halflife_model.py --cohort rodent --tier strict
```

`build_halflife_training_set.py` builds the current **320-row** model table without dropping
repeated treatment arms or independent reports. Three cohort requirements take it there from
the 467 kept records: the dose route must be **parenteral** (iv, sc, ip or im - 84 records
removed), the SMILES must be the **whole dosed molecule** (57 removed), and the bioanalytical
method must be **named by a source** (6 removed). All three are explained below. It then
clusters the sequences with `phmmer` into **90 homology clusters**, assigns a cluster-preserving 80/20 split, and writes six separately
measured feature blocks:

```
[ ESM-2 (1,280) | PeptideCLM (768) | molecule annotation (11) | study context (36)
  | termini (13) | physicochemistry (16) ]
```

Arm-level full-text conditions and PK values are linked by PEPlife2 ID into
`rodent/pk_observations.csv` and row metadata. Label evidence tiers are audit metadata, not
model features. The outputs are `rodent/strict/halflife_meta_with_split.csv` and the Arrow
datasets under `../results/half_life_muriverse/`
(`halflife_rodent_strict_with_embeddings` and `halflife_rodent_strict_ablation_with_embeddings`).
Encoder outputs are memoised per unique sequence and per unique SMILES, so a rebuild that does
not change the model id takes seconds.

The Arrow datasets are generated by the build command and are not stored in this repository. They may be absent from a checkout until the build command above
is run; the canonical CSV metadata and fitted models are stored under `MuriVerse/rodent/strict/`.

**The model** is `halflife_model.HalfLifeEnsemble`: an equal-weight ensemble of ridge, an RBF
SVR and a shallow XGBoost, with fixed hyperparameters, on the 76 compact columns (molecule
annotation, study context, termini, physchem). The ESM-2 and PeptideCLM blocks are measured in
the ablation but score lower under held-out clusters and source identifiers, so they are not in the
shipped model.

`evaluate_halflife_model.py` is where the model is *measured*. It does not use the saved 80/20
split for reported numbers. Instead it rotates the held-out group under 5-fold grouped
cross-validation (fold assignment reshuffled for each of three seeds) and reports Pearson r,
Spearman rho and the AUC against half-life >= 6 h, each with a group-bootstrap interval, under
peptide-cluster, source-identifier, and combined cluster-plus-source holdouts:

* **`cluster`** - held-out homology clusters.
* **`PMID`** - held-out source identifiers (the field includes patent identifiers where no PMID exists).
* **`cluster_PMID`** - held-out connected groups sharing neither a homology cluster nor a
  source identifier. This is the strict transfer estimate and the headline.

Reported numbers are the mean of the per-seed scores (each seed is a full out-of-fold pass).
`--nested` additionally chooses the feature configuration inside each training fold
(`table31b`), which is the number to quote when the configuration itself was selected from the
ablation. Current headline (species-pooled strict cohort, clusters + source identifiers held out, 320
rows): Pearson 0.57, Spearman 0.58, AUC 0.75 (0.55 / 0.56 / 0.74 with nested selection);
mouse AUC 0.60, rat 0.79. Most of the remaining
error is a per-source offset the model cannot know for a new source - see the Figure 3
caption.

The headline is lower than the 470-row build it replaces (0.60 / 0.60 / 0.82), which is the
expected direction and not a regression in the model. Three of the removed covariate levels -
`route_class == "not reported"`, `structure_tier != "complete"` and `assay_method ==
"unknown"` - were each a near-constant marker of a particular subset of sources, so a
grouped-holdout model could lean on them the way it leans on any source-identifying feature.
The 320-row estimate is over a cohort where every row's route, structure and assay are
established, so there is less to lean on and the number is a better estimate of transfer to a
new peptide.

Figure 3 presents continuous-outcome performance rather than thresholded ROC or
precision–recall panels. The 4 h screening-label sensitivity analysis (56/320 positives)
remains available in `figures/tabs/table10b_muriverse_halflife_ranking_metrics.csv`, but it
does not change the regression fit or the continuous predictions used for Figure 4.

Outputs: `rodent/strict/oof_predictions.csv` (one out-of-fold prediction per record, under
each grouping) and `../figures/tabs/table31_muriverse_halflife_feature_ablation_strict.csv`.

After the ablation, it refits the selected feature set (and the best SMILES-embedding variant) on all
rows and saves them to `rodent/strict/model/` as `.joblib` files; `--models-only` refits them from the existing ablation table without re-scoring. These models are for application; all reported
performance comes from the held-out grouped folds. `predict_polymyxins.py` is their consumer.

### Tiers

`--tier strict` is the default and is what the figures use. `--tier pooled` can run the same
machinery over the two per-species `<species>/in_vivo_half_life.csv` tables instead, writing
to `rodent/` rather than `rodent/strict/`; it is a broader comparison that includes
conjugates and unresolved partners.
The current final model and its reported metrics use only the strict 320-row dataset.

### Side branches

| script | does |
|---|---|
| `build_rodent_smiles.py` | Runs the SMILES builder over **every** pooled record, not just the modelling subset → `rodent/smiles_all_records.csv`, a historical all-records structure audit |
| `build_pk_parameters.py` | Per-agent dose/route/matrix PK parameter table for the dosed cohort |
| `extract_pk_figures.py` | Crops each dosed study's concentration-time figure out of its mirrored PDF |
| `analyze_data_provenance.py` | Schema-tolerant provenance and quality audit → [`../docs/data-provenance-quality.md`](../docs/data-provenance-quality.md) |
| `predict_polymyxins.py` | Held-out application of the saved strict models to colistin and polymyxin B in mouse and rat (Figure 4) |
| `study_conditions.py` | Sample matrix, assay method, strain, sex and recovered dose route read from PEPlife2's free-text design columns (model inputs) |
| `halflife_condition_curation.py` | Dose route and assay method curated from the primary sources, against the ledgers in `data/curation/` |

---

## Figures

The main figures, supplementary figures and their tables are generated from the data
scripts; `Rscript figures/code/run_all.R` rebuilds those outputs. Data preparation is
separated from drawing — a `NN_data_*.R` script loads and summarises, a `figN.R` script
draws.

| figure | subject | data script |
|---|---|---|
| 1 | MuriVerse study overview (schematic) | not generated by code in this repository |
| 2 | The final dataset, as delivered (326 modelling rows) | `04b_data_muriverse_strict.R` |
| 3 | The half-life regressor, cross-validated, with the block ablation | `05_data_halflife_model.R` |
| 4 | Held-out application to the polymyxins, shipped vs structural | `06_data_polymyxin.R` |
| S1 | Audit panels for the strict rodent half-life set (sources, clusters, structure features, dose-route and assay provenance) | `supp_fig1.R` (panel objects inherited from `fig2.R`) |
| S2 | Performance and label reproducibility | `supp_fig2.R` (inputs in `results/supp_fig2/`) |
| S3 | Polymyxin comparison against every administration route | `fig4.R` |

Supplementary Figure 1 also draws the source-database composition of PEPlife2, ThPDB2 and
PepTherDia2. ThPDB2 and PepTherDia2 are not redistributed; their species counts are committed in
`results/database_comparison/` and are regenerated from the raw files if you place
them under `MuriVerse/data/` (see [`../DATA_SOURCES.md`](../DATA_SOURCES.md)).

`04_data_muriverse.R` still describes the pooled **mouse** cohort and writes tables 3–9 and
15–19. It is background provenance, not Figure 2's input, since Figure 2 uses the strict
rodent set. The strict selection funnel is in Figure 2a and the before/after comparison is
Figure 2d; the build history is documented in `rodent/strict/README.md`.

---

## Environment

The pipeline needs `torch`, `transformers`, `xgboost`, `scikit-learn`, `joblib`, `rdkit`, `datasets` and `pandas`, plus
the `phmmer` binary (`hmmer` package) and R with `tidyverse`, `ggpubr`, `uwot` (the UMAP in Figure 2e) and `ggnewscale` (the two fill scales on Figure 2e). The versions are pinned in [`../requirements.txt`](../requirements.txt).

Model weights (`facebook/esm2_t33_650M_UR50D`, `aaronfeller/PeptideCLM-23M-all`) are pulled
from the Hugging Face cache; PeptideCLM additionally needs the SPE tokenizer vocabulary, which `scripts/fetch_tokenizer.sh` downloads from the upstream
PeptiVerse repository into `external/PeptiVerse/` (or set `PEPTIVERSE_DIR`).
