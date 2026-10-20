# mouse

Mouse pharmacokinetic (PK) evidence for MuriVerse: the mouse half-life cohort and the
structured, evidence-quoted extractions derived from the primary full-text articles behind it.

Everything here is mouse (`Mus musculus`) data only. The source is currently
[PEPlife 2.0](https://webs.iiitd.edu.in/raghava/peplife2/) — other mouse half-life sources may
be added to this folder later, which is why files are named after their source
(`peplife2_rows.csv`) rather than after the species.

**This is one of two cohorts.** A rat cohort built by the same scripts lives in
[`../rat`](../rat/README.md), and the half-life training set is now the two POOLED, not mouse
alone — see [`../build_halflife_training_set.py`](../build_halflife_training_set.py), which
carries a `species` column, clusters over the combined mouse-and-rat set and balances the
split within each species. The full filtered mouse+rat eligible table is
[`../rodent/eligible_in_vivo_half_life.csv`](../rodent/eligible_in_vivo_half_life.csv),
built by
[`../build_rodent_eligible_set.py`](../build_rodent_eligible_set.py); it precedes
model-specific filtering. Read the rat README for the selection defect that cohort exposed (the PEPlife2 organism
query matches per publication, not per record) and for the curation findings that apply to both
cohorts, including the sequence-encoding problems fixed in `clean_sequence`.

## Why this exists

PEPlife2 records a half-life value, a `Test Sample` matrix string (e.g. `"mouse serum"`,
`"mouse plasma"`), and a `Vivo/Vitro` flag per entry, but these fields are short free-text
annotations that don't reliably distinguish, for example, a half-life measured by dosing a
live mouse (in vivo pharmacokinetics) from one measured by spiking a peptide into mouse
serum in a tube (ex vivo stability). MuriVerse needs that distinction, so this pipeline goes
back to the cited article for every mouse entry, reads the actual Methods/Results body, and
records a structured, evidence-quoted judgment of study design, strain detail, dosing, and the
reported PK values.

## Pipeline

Run from `MuriVerse/`. Steps 0-1 define the cohort, steps 2-4 curate it.

### 0. Record table from the REST API — `build_peplife2.py`

```
python build_peplife2.py
```

Pulls the complete PEPlife 2.0 record set (4,500 entries) from the database's REST endpoint
(`api/api.php?dataType=organism&dataValue=`), archives the raw JSON at
`../data/peplife2/api/records_all.json`, and writes `../data/peplife2/peplife2_records.csv`.
The API is the source of truth: it carries `dose_route`, which the site's .xlsx exports omit.
It does not return the exports' precomputed `half-life(in Sec)`, so that column is re-derived
from `Half Life` + `Units Half Life` (midpoint of a leading range, first number otherwise).

### 1. Cohort selection + full-text retrieval — `fetch_peplife2_api.py`, `build_rodent_cohort.py`

```
python fetch_peplife2_api.py --species mouse    # defines the cohort
python build_rodent_cohort.py --species mouse
```

- Cohort membership is PEPlife2's own organism annotation, queried over the REST API, not our
  text matching. PEPlife2 splits the mouse literature over two organism values, so the cohort
  is the union of `dataValue=Mouse` (369 entries) and `dataValue=Mice` (514), archived at
  `../data/peplife2/api/organism_{Mouse,Mice}.json` and merged into
  `../data/peplife2/peplife2_api_mouse.csv`. `api_query` on each row records which spelling
  returned it. If an organism query has not been archived yet, `build_rodent_cohort.py` falls
  back to matching `\bmice\b|\bmouse\b` in `Test Sample` for the uncovered entries and says so
  in `cohort_source`.
- Flags entries whose matrix is actually a mouse-derived **cell line** (3T3, RAW264, L929,
  B16, Neuro-2a, C2C12, Hepa, P388, CT26, 4T1, MC38, ...) via `matrix_is_cell_line`, since
  those are in vitro, not mouse PK, despite the "mouse" name.
- Writes the selected cohort to [peplife2_rows.csv](peplife2_rows.csv) (883 rows, one per
  PEPlife2 record; a single PMID can contribute many rows for different peptides, matrices, or
  assays).
- Resolves each distinct numeric PMID to a PubMed Central record via NCBI `elink`
  (`pubmed_pmc`), then fetches the open-access full text as BioC JSON from the PMC OA BioC
  service. A record with no PMC linkout, or a PMC record with no body sections (title/
  abstract only), is marked unavailable rather than silently skipped.
- Also fetches the PubMed abstract for provenance/QC (never used as an extraction source,
  see below).
- Records the outcome for every PMID in
  [fulltext_retrieval_status.csv](fulltext_retrieval_status.csv)
  (`pmcid`, `fulltext_retrieved`, `fulltext_status`, `fulltext_source`, `fulltext_chars`,
  `n_peplife_rows`), and stores the raw BioC document at `fulltext/<PMID>.json` for every
  PMID where retrieval succeeded.
- Supports `--skip-existing` to resume without re-fetching PMIDs already resolved.

### 2. Rendering to plain text — `render_fulltext.py`

```
python render_fulltext.py --dir mouse
```

- Converts each `fulltext/<PMID>.json` BioC document into a single plain-text file at
  `fulltext_txt/<PMID>.txt` for the extraction pass to read.
- Drops reference lists, acknowledgments/funding, author contributions, competing-interest
  statements, and supplementary blurbs (`REF`, `ACK_FUND`, `AUTH_CONT`, `COMP_INT`, `SUPPL`) —
  they add bulk and carry no PK evidence.
- Keeps figure and table captions, because concentration-time profiles and PK summary
  tables usually live there rather than in running text.
- Emits the abstract in a clearly fenced block —
  `<<<ABSTRACT - DO NOT USE AS AN EXTRACTION SOURCE>>> ... <<<END ABSTRACT>>>` — separate
  from `=== BODY ===`, so the extraction step can enforce a body-only evidence rule (see
  below). Abstracts often round off or omit detail (e.g. "half-life was significantly
  extended" instead of a number), and are excluded as a source of truth for that reason.
- Prepends a `=== PEPLIFE2 ROWS FOR THIS PMID ===` block listing the distinct
  `Sequence, Name, Half Life, Units Half Life, Test Sample, Vivo/Vitro, Assay` values PEPlife2
  recorded for that PMID. This is database metadata used only to tell the extractor which
  peptides/values to look for in the body — it is explicitly *not* an allowed extraction
  source itself.

### 3. Structured extraction — [`../EXTRACTION_PROTOCOL.md`](../EXTRACTION_PROTOCOL.md) → `extractions/<PMID>.json`

Each rendered `fulltext_txt/<PMID>.txt` is read in full (by a human or an LLM following the
protocol) and turned into one JSON file per PMID at `extractions/<PMID>.json`. The protocol's
absolute rules:

1. **Body only** — no value may be taken from the fenced abstract block; if a fact appears
   only in the abstract, it's recorded as `"not reported in body"`.
2. **No outside knowledge** — never fill a value (dose, n, half-life) from what the extractor
   already knows about the drug; only what the body states.
3. **Every non-null value needs evidence** — a verbatim quote (≤25 words) and the section it
   came from.
4. The PEPlife2-rows header is metadata for orientation only, never an extraction source.

Each `extractions/<PMID>.json` follows a fixed schema (see the protocol for the full field
list); the fields most relevant to nuance-based subsetting are:

| Field | Meaning |
|---|---|
| `study_types` | Any of `clinical`, `in vivo`, `ex vivo`, `in vitro` — all types the body supports for this paper. A single paper is very often a mix (e.g. in vivo dosing **and** an in vitro mechanism assay). |
| `study_type_evidence` | Per-type quote + section backing each entry in `study_types`. |
| `in_vivo_species` | One entry per species actually dosed as an intact live animal, with strain/sex/age when stated. A mouse-cohort paper can still dose other species; only the mouse entries are mouse PK. |
| `ex_vivo_origin` | `"animal"`, `"human"`, `"human and animal"`, or `"not reported"` — whose tissue/blood/serum was used for an ex vivo assay. |
| `ex_vivo_detail` | Free text describing the ex vivo matrix/procedure (e.g. "peptide spiked into mouse plasma, incubated at 37°C, sampled over time"). |
| `in_vitro_detail` | Free text describing any cell-culture/buffer/cell-free assay. |
| `agents[].dose`, `.route`, `.group_n` | Per-agent dosing detail — only populated for agents that were actually administered to a living animal. |
| `agents[].half_life`, `.half_life_matrix`, `.cmax`, `.tmax` | The PK read-out itself and the matrix it was measured in (e.g. `"mouse plasma"` vs `"mouse serum, ex vivo"`). |
| `concentration_time_figure` | Whether a concentration/level-vs-time (or %-remaining-vs-time) figure exists, and whether it's an in vivo PK curve or an ex vivo stability curve. |
| `bioanalytical_method` | Assay used to quantify the peptide (LC-MS/MS, ELISA, RIA, HPLC/UV, fluorescence, radiolabel counting, imaging, ...). |
| `notes` | Free-text caveats — e.g. PK reported for only one of several peptides, internal inconsistencies in the source paper, or missing Methods (e.g. Methods in supplementary information not retrieved). |

The key study-design distinction this protocol enforces (crucial for downstream subsetting):

- **in vivo** = an intact, living mouse was dosed and sampled — true animal pharmacokinetics.
- **ex vivo** = tissue/blood/serum/plasma was taken *out of* a body first, then treated/spiked
  and monitored — e.g. serum/plasma stability assays. **These are classified as ex vivo, not
  in vitro**, even though they resemble a bench assay.
- **in vitro** = cell culture, buffer, purified enzyme, or other fully cell-free/animal-free
  system.

### 4. Nuance subsetting — `../subset_by_nuance.py`

```
python subset_by_nuance.py --species mouse
```

Collapses every completed `extractions/<PMID>.json` into one paper-level row of
[nuance_subset.csv](nuance_subset.csv), flagging `is_in_vivo`/`is_ex_vivo`/`is_in_vitro`/
`is_clinical`, and categorizing `agents[].route` into `in_vivo_routes` (iv/sc/ip/im/oral/icv/
infusion_pump/not_administered/not_reported) and `agents[].half_life_matrix` +
`ex_vivo_detail` into `ex_vivo_matrices` (plasma/serum/blood/tissue_organ/urine/csf/...).
PMIDs with no completed extraction still get a row (`extraction_status` =
`rendered_not_yet_classified` / `retrieved_not_yet_rendered` / `fulltext_unavailable`), so the
file is a complete 227-PMID inventory, not just the curated subset. This is paper-level, not
record-level: a PMID's tags summarize everything reported anywhere in that paper, so they
cannot be joined back onto individual PEPlife2 rows one-to-one.

This curated, paper-level nuance data is what feeds panels (b) and (c) of `figures/figs/fig2`
(`figures/code/04_data_muriverse.R` reads `nuance_subset.csv` directly) — see
[Record-level design reclassification](#record-level-design-reclassification-full-cohort)
below for how the same question is answered across the full cohort, not just the 48 curated
papers.

### 5. Provenance assembly — `../build_study_dossier.py`

```
python fetch_pubmed_metadata.py --species mouse   # PubMed record per PMID
python fetch_peplife2_papers.py --species mouse   # PEPlife2's mirrored source PDFs
python build_study_dossier.py --species mouse
```

Answers "which study produced this half-life value?" for every record by joining the PEPlife2
fields, the recovered study design, the PubMed record behind the PMID, the curation stage, and
the curated per-agent evidence into one table per unit of analysis. This module owns the design
classification; `figures/code/04_data_muriverse.R` reads `design`/`matrix` from
[records_provenance.csv](records_provenance.csv) rather than re-deriving them.

`fetch_pubmed_metadata.py` pulls title, journal, year, DOI, PMCID, publication types and MeSH
terms for all 227 numeric PMIDs (archived efetch XML under `../data/pubmed/`).
`fetch_peplife2_papers.py` downloads the author PDF that PEPlife2 mirrors at
`peplife2/Papers/<PMID>.pdf`; 157 of 227 are on the server, which supplies readable text for 68
publications that have no PubMed Central open-access body.

Evidence now available per publication:

| Tier | Publications | Records |
|---|---:|---:|
| Curated extraction (`extractions/*.json`) | 90 | 296 |
| PMC full text retrieved, not yet curated | 10 | 54 |
| Publisher PDF only (PEPlife2 mirror) | 68 | 219 |
| Metadata only (no retrievable text) | 67 | 311 |

## File inventory

| File / folder | Produced by | Contents |
|---|---|---|
| [peplife2_rows.csv](peplife2_rows.csv) | `build_rodent_cohort.py` | Every PEPlife2 entry the organism API assigns to mouse (883 rows); includes `api_query`, `cohort_source`, `matrix_is_cell_line` and `pmid_valid` flags. |
| [records_provenance.csv](records_provenance.csv) | `build_study_dossier.py` | One row per PEPlife2 entry with the recovered `matrix`/`design`, the PubMed record, the curation stage (`evidence`), the curated per-paper evidence, and `is_repeat`/`identical_copies` duplicate flags. |
| [studies.csv](studies.csv) | `build_study_dossier.py` | One row per publication (235 identifiers): records and sequences contributed, design mix, half-life range, evidence tier, citation. |
| [sequences.csv](sequences.csv) | `build_study_dossier.py` | One row per unique sequence (407): records, contributing studies, half-life spread and a `multi_study_conflict` flag. |
| [STUDIES.md](STUDIES.md) | `build_study_dossier.py` | Readable dossier: one section per publication with its citation, contribution, design, curated evidence and a table of every record it backs. |
| [fulltext_retrieval_status.csv](fulltext_retrieval_status.csv) | `build_rodent_cohort.py` | One row per distinct numeric PMID (227) with PMC resolution + retrieval outcome. |
| `fulltext/<PMID>.json` | `build_rodent_cohort.py` | Raw PMC OA BioC document, for PMIDs where full text was retrievable (101 files). |
| `fulltext_txt/<PMID>.txt` | `render_fulltext.py` | Plain-text rendering of the BioC document, fenced abstract + body, for extraction (100 files). |
| `extractions/<PMID>.json` | manual/LLM extraction per [EXTRACTION_PROTOCOL.md](../EXTRACTION_PROTOCOL.md) | Structured, evidence-quoted PK/study-design record per PMID (90 files; every dosed PMID with body text is done). |
| [nuance_subset.csv](nuance_subset.csv) | `../subset_by_nuance.py` | One row per PMID (227): curated in-vivo/ex-vivo/in-vitro flags, routes and matrices for the 90 extracted papers, plus curation status for the rest. |
| [in_vivo_half_life.csv](in_vivo_half_life.csv) | `../subset_in_vivo_half_life.py` | Analysis subset: 513 deduplicated `in vivo PK` records with a circulating half-life measured in a live dosed mouse. 190 studies, 263 sequences, median 4.0 h. |
| [pk_figures.csv](pk_figures.csv) + `pk_figures/` | `../extract_pk_figures.py` | 124 concentration-time figures cropped out of the publisher PDFs of the dosed studies (91 of 194), with the matched caption and the crop method used. |
| [pk_parameters.csv](pk_parameters.csv) | `../build_pk_parameters.py` | One row per dosed study arm (570 arms, 194 publications): dose, route, Cmax, tmax, AUC, absorption half-life, elimination half-life, the arm's figure, and a `needs_curation` gap list. |
| [PDF_WORKLIST.md](PDF_WORKLIST.md) + [pdf_retrieval_worklist.csv](pdf_retrieval_worklist.csv) | `../list_missing_pdfs.py` | The 59 dosed publications with no local PDF (168 dosed records), grouped by retrieval route, with PubMed/DOI/PMC links. |
| [batches.json](batches.json) | manual batching | PMIDs grouped into small batches for incremental extraction work. |
| `../data/pubmed/` | `fetch_pubmed_metadata.py` | Archived efetch XML plus `pubmed_mouse.csv`. |
| `../data/peplife2/papers/` | `fetch_peplife2_papers.py` | PEPlife2-mirrored source PDFs (157 files, git-ignored) and `paper_download_status.csv`. |

As of this writing, 90 of the 100 rendered full texts have a completed `extractions/*.json`.
All 85 **dosed** publications with body text are done; the 10 remaining are non-dosed
(ex vivo / in vitro only).

## What the cohort actually rests on

- **The literature is recent and thinly spread.** 227 publications back 880 records; the
  median study contributes 2 records and the top 10 studies contribute 27% of the cohort.
  Two thirds of the records come from papers published 2016 or later, with a spike in 1999
  (82 records, mostly one antimicrobial-pseudopeptide paper).
- **Repeats are real.** 81 records are exact duplicates of another record in the same study
  (same PMID, sequence, `Test Sample` and reported value). PMID 10461747 alone contributes 64
  records that are 8 measurements repeated 8 times. Any unweighted mean or median over the raw
  rows counts those studies many times over; use `is_repeat` to collapse them.
- **Most sequences are measured once.** 386 of 407 sequences appear in a single publication;
  only 21 appear in more than one, and 13 of those disagree by more than 10x between studies
  (see `multi_study_conflict` in [sequences.csv](sequences.csv)). Exenatide
  (`HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS`) spans 13 studies and 0.14-730 h, because the
  records mix native peptide with PEGylated and depot formulations under one sequence.
- **Three sources are secondary literature.** PMIDs 35259149, 35177945 and 36997003 are a
  systematic review/meta-analysis and two reviews, contributing 52 records whose values were
  compiled from other papers rather than measured in the cited work.
- **Eight identifiers are not publications at all** (patent numbers such as
  `EP 2022074607 W` and internal identifiers such as `138-792-131-275-411`), contributing 30
  records with no retrievable source.
- **A third of in vivo records name only a fluid.** 276 of the in vivo PK records give a
  `Test Sample` such as `"Mice plasma"` with no dosing detail, so the route, dose and sampling
  scheme are only recoverable from the article itself.

## Known caveats

- A PMID with `fulltext_retrieved: False` in `fulltext_retrieval_status.csv` has no PMC
  open-access body (either no PMC record at all, or a title/abstract-only PMC deposit) and
  will never get a `fulltext_txt/` or `extractions/` file under this pipeline.
- One PEPlife2 PMID can back many rows (different peptides, matrices, or assays in the same
  paper) and one `extractions/<PMID>.json` can likewise list many `agents`, each with its own
  study-type/dose/half-life evidence — extraction is per-paper, not per-PEPlife2-row.
- Some papers' Methods live in supplementary information that isn't part of the retrieved
  PMC body; when that happens it's called out in `notes` and the corresponding detail fields
  are `"not reported in body"` rather than guessed.
- `study_types` is frequently a mix within a single paper — always filter/join on the specific
  `agents[]` entry and its own evidence, not just the paper-level `study_types` list.
- Cohort membership is PEPlife2's organism annotation. A paper that dosed mice but was
  annotated to another organism value will not appear here.

## Record-level design reclassification (full cohort)

The curated extraction pipeline above (Steps 2-4) gives ground-truth study design, but only
for the 48/227 PMIDs curated so far. For the full 880-record cohort (883 minus the 3
cell-line rows), PEPlife2's own `Vivo/Vitro` flag is the only design signal available — and it
is not reliable: it has no `ex vivo` category, so a peptide spiked into isolated mouse
plasma/serum/blood and monitored for degradation (a stability assay, not animal PK) is filed
as `in vitro`; and it is not even internally consistent (some tissue-homogenate stability
assays are flagged `in vivo`).

`figures/code/04_data_muriverse.R` recovers the actual per-record design by pattern-matching
the raw `Test Sample` free text instead of trusting the flag, writing one row per design
bucket to `figures/tabs/table9_muriverse_record_level_design.csv`:

| Design (text-recovered) | Records | With half-life | Evidence in `Test Sample` |
|---|---:|---:|---|
| in vivo PK | 571 | 558 | `in vivo`-flagged, matrix is plasma/serum/blood/whole-animal and not described as diluted |
| ex vivo (fluid stability) | 239 | 227 | `in vitro`-flagged serum/plasma/blood, e.g. `"Mouse serum + 1µM catropril"`, `"50% v/v mouse plasma"`; plus explicitly diluted or neat matrix such as `"25% mouse serum"` / `"Full mouse serum"` even when flagged `in vivo` — spiked into isolated matrix, not dosed into a live mouse |
| ex vivo (tissue homogenate) | 65 | 65 | `"Mouse liver homogenate"`, `"Mouse brain proteases"` — including entries PEPlife2 mislabeled `in vivo` (e.g. `"Mice brain homogenize"`) |
| in vitro (other) | 3 | 3 | e.g. cultured primary neurons |
| ambiguous / no flag | 2 | 2 | `Vivo/Vitro` = `N.A.` |

So of the 880 mouse records, only **571 (65%)** represent a half-life measured by dosing an
intact, living mouse; **295 records that do carry a reported half-life value are ex vivo or in
vitro** (227 fluid-stability + 65 tissue-homogenate + 3 other), not animal pharmacokinetics,
despite often sharing the same "mouse serum"/"mouse plasma" matrix label as genuine in vivo
records. The classification lives in `../build_study_dossier.py` and is checked against the
curated extractions: at the paper level it agrees with the body text for 86 of the 90 curated
publications. Where a paper has been curated, prefer its per-agent evidence
(`extractions/<PMID>.json`, surfaced in `nuance_subset.csv`) over this reclassification. The
resulting inclusion workflow is plotted in panel (b) of `figures/figs/fig2` (see
[fig2_caption.md](../../figures/figs/fig2_caption.md)).

### Live-dosed half-life subset

`../subset_in_vivo_half_life.py` builds the analysis set — a circulating half-life measured
in a live, dosed mouse, one row per distinct measurement:

```
conda run -n peptiverse-local python ../subset_in_vivo_half_life.py --species mouse
```

| Filter | Records | Dropped |
|---|---:|---:|
| cohort (`in_cohort`) | 880 | |
| a live animal was dosed (`animal_dosed == "dosed"`) | 584 | −296 |
| a half-life value is reported | 571 | −13 |
| `design == "in vivo PK"` | 558 | −13 |
| absorption-phase half-lives dropped | 534 | −24 |
| exact repeats collapsed (`is_repeat`) | **513** | −21 |

#### Historical single-species model build: 513 published records → 459 rows

This section documents the superseded mouse-only model build and is retained as its audit
trail. It is **not** the current training dataset. The final MuriVerse model uses the strict
pooled mouse-and-rat set of 320 rows in [`../rodent/strict/README.md`](../rodent/strict/README.md).

`in_vivo_half_life.csv` holds 513 records. The training set holds **459**. The gap is real and
fully accounted for, but until this update nothing printed the chain, so three different totals
for the same cohort were in circulation — 513 here, 467 from `subset_in_vivo_half_life.py`'s
summary line, and 459 from the training build.
`../build_halflife_training_set.py` now emits the whole chain, with the PEPlife2 ID of every
removed record, to `halflife_cohort_reconciliation.csv`, and asserts that it starts at the
published count and ends at the training row count:

| Filter (applied by `build_halflife_training_set.py`) | Records | Dropped |
|---|---:|---:|
| `in_vivo_half_life.csv` | 513 | |
| concentration-decay `pk_phase` only | 467 | −46 |
| a sequence exists | 463 | −4 |
| `half_life_h > 0` | 463 | −0 |
| identical features + label collapsed | **459** | −4 |

What each step removes, concretely:

* **−46, `pk_phase`.** 28 `distribution / alpha` and 18 `non-concentration` records. A
  distribution half-life and a plasma-concentration terminal half-life answer different
  questions, and a non-concentration readout (activity, excretion, radioactivity) is a
  different parameter again. `CORE_PK_PHASES` keeps `unlabelled` and `terminal / beta`. **This
  467 is the number `subset_in_vivo_half_life.py` used to print as "the analysis set"; it is
  an intermediate count and that line now says so.**
* **−4, no sequence.** IDs 4178, 4179 (37875481, IL-15 constructs) and 4520, 4521 (36075899,
  S-20-1). PEPlife2 leaves `Sequence` empty on these, so a sequence model has no input.
* **−0, non-positive label.** No mouse record fails this; it guards `log10` of the label.
* **−4, de-duplication.** IDs 1869, 2052, 2992, 2994. An identical feature vector carrying an
  identical label is one training example. Note these are **not** the same thing as
  `subset_in_vivo_half_life.py`'s `CURATED_DUPLICATES`, which are database transcription
  repeats: two of these four are *two independent papers* reporting the same value for the
  same peptide (GLP-1 at 0.033 h in 15012592 and 18085721; BNP at 0.052 h in 15587934 and
  19836183), and the other two are separate rows within 8218482 that differ only in the `Name`
  PEPlife2 gave them.

Unique sequences are **258 at both ends** of this chain: de-duplication only ever removed a row
whose sequence was already present, so the 513 → 459 reduction costs records, not peptides.

Both `animal_dosed` and `design` are needed: `design_of` buckets every tissue matrix as ex
vivo before it reads the flag, so `design` alone drops 11 dosed animals whose organs were
harvested, while `animal_dosed` alone keeps those 11 plus 2 `in vitro (other)` records —
neither is a circulating-plasma PK measurement, so both are excluded here.

**Absorption vs elimination.** PEPlife2 appends the PK parameter to the value as free text
(`"2.12 ± 0.3 (T1/2 Elimination)"`); the script parses it into a `pk_label` column and drops
the 24 absorption-phase rows, which measure entry from a subcutaneous depot into circulation
rather than persistence in it. PMIDs 28989813 and 36630826 also report the elimination
half-life of the same animals, which is kept; PMID 31281949 reports *only* SC absorption
half-lives and therefore leaves the subset entirely (hence 190 studies, not 191). PEPlife2
also spells it `"Abosrption"`, so the match is `ab[os]{2}rption`.

> **Not a rate constant.** 11 of those rows carry `Units Half Life == "Hours-1"`, which looks
> like a first-order rate constant needing `t½ = ln2/k`. It is not. The source paper
> ([PMID 36630826](../data/peplife2/papers/36630826.pdf), Table 2) mislabels *both* of its
> half-life columns `(h−1)` while its body text gives the same numbers in hours
> (`"t1/2 absorption=9.05 ± 0.9h"`) and its Methods state `t1/2 = 0.693/k` was already
> applied. The values are hours and need no conversion — they are excluded as the wrong PK
> parameter, not as the wrong units.

**Deduplication.** Repeats are collapsed on `PMID + sequence + Test Sample + Half Life +
Units`, recomputed after the other filters so the surviving copy is one that passed them. The
key ignores whitespace: PEPlife2 transcribed some tables twice with inconsistent spacing and
stray non-breaking spaces (`"23.14 ± 3.2(T1/2 Elimination)"` vs `"23.14 ± 3.2 (T1/2
Elimination)\xa0"`), which a literal key reads as distinct measurements — that accounts for 7
of the 21 repeats, including a second full copy of PMID 36630826's Table 2.

One record is corrected against the source: ID 4269 is Table 2's *absorption* value for
Ex-POEGMAopt in immunized mice (12.1 ± 1.0 h) filed under the elimination label; that arm's
elimination half-life is 105.7 h (IDs 4270/4281, both retained).

The result, [in_vivo_half_life.csv](in_vivo_half_life.csv), is **513 records / 190 studies /
263 unique sequences**; 419 report a `Dose Route` (`N.A.` on the rest). Half-lives span 0.5
min to 2440 h with a median of 4.0 h.

**The PK phase is resolved as of 2026-09-29** (was an open caveat before that): the script
now adds `pk_phase` (`unlabelled` / `terminal / beta` / `distribution / alpha` /
`non-concentration`, parsed from the same free-text label as `pk_label`, reading the
alpha/beta subscript *before* the word "elimination" because several papers label both
phases that way), plus `route_class`, `topology` and `hl_engineering`. The core analysis set
is `pk_phase in {unlabelled, terminal / beta}` — 467 of 513 records — which drops
distribution-phase values and the handful of non-concentration readouts (enzyme activity,
excretion, radioactivity).

### Curated full-text cross-check (2026-09-29)

Until this pass, `half_life_h` was PEPlife2's own transcription of the source paper's value —
re-derived from PEPlife2's free text (see step 0 above), not from a curator reading the paper.
The 90 curated papers' verbatim quotes (`extractions/*.json`) existed only as the
`curated_half_lives` evidence string, never checked against the modelling label. The script
now does that check via `add_curated_crosscheck()`, but applies a correction only where it is
unambiguous: exactly one curated agent reporting a half-life for that paper, matched to
exactly one surviving row here. Looser name matching was tried first and rejected — PEPlife2
and the curated schema both strip parenthetical detail from a construct name, so
`ELP(0FA)-GFP` and `ELP(5FA)-GFP` collide onto one key while their half-lives differ 18-fold.

| `half_life_provenance` | Records |
|---|---:|
| PEPlife2 transcription (no curated full text) | 496 |
| matches curated full text (within 15%) | 9 |
| curated full text available (ambiguous match, not auto-applied) | 7 |
| corrected to curated full text | 1 |

The one correction: PMID 28010844, Octreotide, 0.467 h → 0.628 h (the curated body text's
"37.7 min in mouse intestine (iMScope); 28 min in mouse stomach", midpointed). The original
PEPlife2 value is kept at `half_life_h_peplife2` for audit. This is a conservative first pass,
not a full reconciliation — extending it past single-agent papers needs the same route/dose/
matrix arm-matching `build_pk_parameters.py` already does for Cmax/tmax/dose, applied with the
same caution the ELP collision argues for.

## PK parameter extraction (dosed cohort)

Two scripts turn the 584 dosed records into a per-arm PK table with its source figure. They
run over the **dosed** cohort, not the `in_vivo_half_life.csv` subset, because absorption
half-lives and two-phase fits are wanted here rather than filtered out.

```
conda run -n peptiverse-local python ../extract_pk_figures.py --species mouse
conda run -n peptiverse-local python ../build_pk_parameters.py --species mouse
```

### Concentration-time figures

`../extract_pk_figures.py` crops each study's concentration-time plot out of its
PEPlife2-mirrored publisher PDF (the BioC full text has captions but no images). A figure is
found by its caption — a block that opens like a caption *and* describes a drug-concentration
readout over time, which excludes the PD time courses (blood **glucose**, tumour volume) that
sit beside it. The artwork is then the union of raster images and vector drawings above the
caption and overlapping its column. Two fallbacks matter:

- a thin union means only a rule or an axis line cleared the area threshold, so it is extended
  up to the previous text block (`crop == "artwork+gap"`);
- accepted-manuscript PDFs list every legend on one text page and put the artwork on later
  pages labelled only `"Figure 2"`, so those plate pages are indexed and matched by figure
  number (`crop == "plate p23"`).

| | Studies | |
|---|---:|---|
| figure cropped | **91** | 124 images in `pk_figures/<PMID>_<figN>.png` (37 MB, not yet tracked or ignored) |
| caption matched, crop failed | 7 | legend-only page, artwork not embedded |
| no concentration-time caption | 40 | PDF has no such figure, or it is in the SI |
| no publisher PDF | 59 | incl. 8 patents; PEPlife2 does not mirror these |

### The parameter table

`../build_pk_parameters.py` writes [pk_parameters.csv](pk_parameters.csv) — **570 arms across
all 194 dosed publications**, one row per (publication, agent, matrix, route, dose). PEPlife2
transcribed several tables twice with different spellings (`"Ex-PEGMw"` at `"1000 nmol/kg"`
and `"Ex-PEGMw(46.3KDa)"` at `"1000 nmol kg⁻¹"`), so arms are grouped on normalised keys —
parentheticals stripped from the agent name, `kg⁻¹` folded to `/kg` — and `n_records` records
how many PEPlife2 rows backed each arm.

`ka_t_half_reported` and `ke_t_half_reported` come from the phase label PEPlife2 buries in the
half-life free text, pivoted into columns, then backfilled from the curated
`half_life_absorption` / `half_life_elimination` fields. `ke_basis` says which backed ke:

| `ke_basis` | Arms | |
|---|---:|---|
| `unlabelled (assumed elimination)` | 362 | single half-life, no phase named |
| `terminal/beta` | 57 | `T1/2b`, `t1/2 second phase` |
| `elimination` | 36 | named `T1/2 Elimination` |
| `curated extraction` | 26 | phase resolved by reading the paper |
| *(none)* | 89 | no half-life for that arm |

Coverage after the curation round below:

| Field | Arms with a value | Publications | Source |
|---|---:|---:|---|
| dose | 544 / 570 | 185 / 194 | PEPlife2 `Concentration`, then curation |
| route | 485 / 570 | 159 / 194 | PEPlife2 `Dose Route`, then curation |
| elimination t½ (ke) | 474 / 570 | 188 / 194 | PEPlife2 phase label, then curation |
| Cmax | 79 / 570 | 30 / 194 | curation only |
| tmax | 73 / 570 | 31 / 194 | curation only |
| absorption t½ (ka) | 17 / 570 | 3 / 194 | PEPlife2 phase label |
| concentration-time figure | 292 / 570 | 91 / 194 | cropped from publisher PDF |

397 arms across 150 publications carry dose + route + ke; 44 arms across 22 publications
carry the full dose + route + Cmax + tmax + ke set.

**Cmax and tmax can only come from curation — PEPlife2 has no field for either.** The per-row
`needs_curation` column lists which of `dose|route|cmax|tmax|ka|ke|figure` that arm lacks, so
it doubles as the worklist.

**Absorption half-lives are close to nonexistent in this literature.** Only 3 of 194 dosed
publications report one, and every one of the 27 papers curated in the round below reported
none — including subcutaneous, intraperitoneal, oral and intranasal studies where an
absorption phase demonstrably exists. Papers overwhelmingly fit a single decay phase and
report tmax instead. Treat ka as unavailable rather than merely uncurated.

### Curation status

All dosed publications that have retrievable full text have now been curated:

| | Publications |
|---|---:|
| dosed publications | 194 |
| with BioC full text | 85 |
| **curated (`extractions/*.json`)** | **85 of 85** |
| PDF only — not yet curated | 50 |
| metadata only — nothing to curate from | 59 |

Cohort-wide that is 90 curated papers of 227 (the extra 5 are non-dosed). The ceiling is
retrieval, not curation effort: the remaining 109 dosed publications have no body text in
`fulltext_txt/`. 50 of them do have a publisher PDF, so rendering PDF text would be the next
unlock; [PDF_WORKLIST.md](PDF_WORKLIST.md) covers the 59 with nothing at all.

[EXTRACTION_PROTOCOL.md](../EXTRACTION_PROTOCOL.md) gained `half_life_phase`,
`half_life_absorption`, `half_life_elimination` and `auc`, plus an explicit warning never to
convert a reported rate constant into a half-life and not to trust an `h⁻¹` column header over
the body text.
