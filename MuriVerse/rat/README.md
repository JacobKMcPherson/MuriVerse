# Rat cohort

The rat half of MuriVerse. Built by the same scripts as [`../mouse`](../mouse/README.md) — every
step takes `--species rat` — so this file records only what differs, plus the curation findings
that are specific to these papers.

The mouse README is the reference for the pipeline itself. Read it first.

This file documents the rat source cohort and its curation. Some pooled training-set counts
below describe an earlier model build; the current final model uses the 320-row strict set in
[`../rodent/strict/README.md`](../rodent/strict/README.md), not those historical pooled rows.

## Why this exists

MuriVerse was a mouse dataset. Rat was added because PEPlife 2.0 holds more rat than mouse
records for the same question, drawn from the same database under the same filters, and the
combined set is what the training pipeline now consumes:

| | mouse | rat |
|---|---|---|
| cohort records | 883 | 1149 |
| distinct PMIDs | 235 | 256 |
| live-dosed in-vivo half-life records | 513 | 804 |
| studies behind those | 190 | 209 |
| unique sequences | 263 | 344 |
| median half-life | 4.0 h | 0.9 h |

The two are pooled — not modelled separately — by
[`../build_halflife_training_set.py`](../build_halflife_training_set.py), which carries a
`species` column, clusters over the pooled set, and balances the train/val split within each
species. Its module docstring explains why, and what it does *not* handle (species is not a
feature of either representation). 44 sequences are measured in both animals, which is the
reason clustering has to be joint.

Clustering is **HMMER**, not MMseqs2, as of this update: `phmmer` all-vs-all over the unique
sequences, an edge where the domain i-Evalue is ≤ 1e-3 and the aligned region covers ≥ 80% of
both sequences, and clusters as the connected components. `cluster_sequences` records how the
two tools compare on this cohort and where they disagree.

## Cohort selection

The rat cohort exposed a defect in the selection logic that mouse could not.
`fetch_peplife2_api.py --species rat` queries the organism values `Rat` and `Rats`; `Rat`
substring-matches `Rats`, so the second query returns a strict subset of the first and the union
is 1181 entries rather than a partition.

More importantly, the organism query matches per **publication**, not per record. 32 of those
1181 entries come from rat papers whose individual measurement was made in human, pig, dog or
NHP material ("Intravenous administration to human plasma"). For mouse the distinction is
invisible — all 883 entries returned by `Mouse`+`Mice` name a mouse in their own `Test Sample`.

`build_rodent_cohort.rodent_cohort` now intersects the API selection with the record's own
`Test Sample`. That drops exactly those 32 and yields 1149 — identical to the independent
text-matched cohort, so the two selection routes agree once the per-publication/per-record
mismatch is removed. The rejected entries are written to
`data/peplife2/peplife2_api_rat_rejected.csv` for audit.

## Full-text retrieval is the binding constraint

Only 57 of 235 numeric PMIDs (24%) are in the PMC OA subset, against 100 of 227 (44%) for mouse:

```
[done] full text retrieved for 57/235 PMIDs
[missed] 166 no PMC record; 12 not in PMC OA subset
```

The publisher-PDF fallback that recovered 61 mouse papers did not work here. PEPlife 2.0 mirrors
one PDF per PMID, but the host stopped responding partway through the rat run and then began
refusing connections outright (`ConnectTimeout` on every request). Only 15 rat PMIDs already had
a mirrored PDF from the mouse run, 4 of which are in the in-vivo set and are rendered into
`pdf_txt/`. **Re-running `fetch_peplife2_papers.py --species rat` when that host recovers is the
single highest-value outstanding step** — it would roughly triple the curatable fraction.

Net readable text: 57 PMC + 4 PDF = 61 papers, of which 50 contribute in-vivo half-life records.

## Curation status

All 50 papers with retrievable text have been read and curated to
[`../EXTRACTION_PROTOCOL.md`](../EXTRACTION_PROTOCOL.md), in `extractions/`:

```
extraction_status
fulltext_unavailable           174
extracted                       50
rendered_not_yet_classified     11
```

That is 155 of 804 in-vivo records (19%). The remaining 649 records come from 164 papers with no
retrievable full text and rest on PEPlife2's own transcription, flagged as such in
`half_life_provenance`. This is the cohort's main weakness and it is a retrieval limit, not a
curation backlog.

Where a paper reports exactly one non-absorption half-life, the value is cross-checked
automatically (`subset_in_vivo_half_life.add_curated_crosscheck`): 9 records match the curated
full text and 6 are curated but ambiguous to match. No auto-applicable disagreement was found —
every value that could be matched one-to-one was transcribed correctly.

## What curation caught

Reading the 50 papers found eight classes of defect. Each is fixed by a named, documented
per-record override in `../subset_in_vivo_half_life.py` rather than a blanket rule, because in
every case a rule general enough to catch it would also have discarded valid data.

| Defect | Where | Effect if uncorrected |
|---|---|---|
| **p-values read as half-lives** | 35703463 IDs 4577, 4580, 4583 | Table 3's `Sex` column has levels male / female / **P**, where P is the male-vs-female t-test significance. The three p-values 0.57, 0.65 and 0.22 were transcribed as elimination half-lives in hours. |
| **Wrong species** | 31818212 IDs 5301, 5302 | All five records filed as "Wistar rats plasma", but 5.58 h and 7.06 h are the **mouse** arm. Three mouse half-lives, all at the long end, would have entered the rat cohort. Reassigned via `CURATED_SPECIES`. |
| **Payload, not peptide** | 31194563, 33665501, 36047255, 34358122, 38555444 (11 records) | The measured analyte is a small molecule carried by a peptide-decorated vehicle — fasudil, cisplatin (as platinum by ICP-MS), calycosin, milrinone, tigecycline. Each paper's own free-drug control shows the label is set by the formulation. Two of these also store a polymer *name* in the `Sequence` field. |
| **Total radioactivity** | 36588717 ID 4405 | `t1/2 = 102 ± 32 h` is the half-life of total [3H] label; the same paper states the parent peptide is undetectable after 4 h with a half-life under 30 min. A ~200-fold overstatement. Reclassified via `PHASE_OVERRIDES`. |
| **Units** | 38092894 ID 4377 | Body says the distribution half-life is 11.55 **minutes**; PEPlife2 recorded 11.55 **hours**, making the distribution phase eight times longer than the elimination phase it precedes. Fixed via `UNITS_CORRECTIONS`. |
| **Unlabelled absorption phase** | 28605180 IDs 5998, 6000 | The subcutaneous 50 min and 48 min values are the *absorption* half-life of a flip-flop model whose elimination half-life is 18 min, recorded without the word "absorption" so the existing text filter could not see them. |
| **Double transcription** | 36588717 IDs 4704–4707 | One table transcribed twice under two spellings of `Test Sample` ("SD rats plasma" vs "Rats plasma"). `DUP_KEY` cannot catch this without also collapsing genuinely distinct conditions — "Immunized" vs "Naïve mice plasma" in 36630826, duodenum vs ileum in 7561632 — so it is a curated exclusion. |
| **Misspelt phase label** | 32082972 ID 5209 | `Eliminatoin Half Life` fell through to `unlabelled` while its five table-mates were `terminal / beta`. Third observed spelling of "elimination"; added to the phase regex. |

Two further classes are recorded but **not** corrected, because they are properties of the
database's representation rather than errors in a value:

* **Reviews as sources.** Four review articles contribute 44 records across both cohorts
  (35259149 alone contributes 32). Their values are re-reports of other papers and may duplicate
  primary records already in the cohort. Now flagged by `source_is_review`, added in
  `build_study_dossier.py`.
* **Patents as sources.** 139 of this cohort's 804 in-vivo records (17%) cite a patent rather
  than a paper, against 18 of 513 for mouse. `EP2277890A1` alone contributes 65 records — the
  largest single source in the rat cohort, larger than any journal article. These have no PubMed
  record and no retrievable full text, so they can never be curated through this pipeline and
  their values rest entirely on PEPlife2's transcription. A further 6 records cite a bare
  internal identifier. Now recorded in a `source_type` column:

  | `source_type` | mouse | rat |
  |---|---|---|
  | journal article (PubMed) | 485 | 657 |
  | patent | 18 | 139 |
  | unidentified (no PMID) | 12 | 6 |
* **Sequences that do not describe the dosed molecule.** This is pervasive and is the most
  important caveat for anyone training on the pooled set — see below.

## The sequence does not encode what drives the label

Curation turned up repeated cases where PEPlife2 stores a sequence that cannot explain its own
half-life. Some are outright contradictions:

* **26400022** gives peptide 1 and peptide S2 the *same* string and half-lives of 13 min and
  212 min. The difference is a stearic acid on a Lys that replaced a linker Gly.
* **34393794** gives K17F and LIT01-196 the same string and three half-lives spanning 187-fold
  (50 s, 28 min, 156 min) — a fluorocarbon chain plus a route change, neither encoded.
* **27240277** gives ENF and PEG2k-ENF the same string, 1.5 h and 16.1 h apart.
* **36232550** and **27217590** are different papers whose compounds share one string,
  `CYIQNCGKG`, with half-lives of 2.3 h and 8.5 h.

Others store a fragment of a much larger molecule: **33374407** stores 28 residues of thymosin
α1 for a construct whose half-life comes from the 600-residue PAS polypeptide fused to it;
**31156041** likewise omits a 600-residue PA polymer; **33706686** stores a 30–39mer motif for a
dosed IgG whose ~1-week half-life is FcRn-driven. And **36557850** represents a
leuprolide–oleic-acid conjugate as leuprolide *minus* its N-terminal pyroglutamate, so the
conjugate looks like a truncation of its own parent.

The practical consequence: on this cohort, half-life extension chemistry is largely invisible to
a sequence-only representation, and identical inputs carry order-of-magnitude different labels.
`hl_engineering` is the column that recovers the mechanism, and the `chem`/`fused` branches
(removed from the current builder) were where it used to enter the model.

## Sequence-alphabet defects, now fixed upstream

Three encoding problems were found and two are fixed in
`build_halflife_training_set.clean_sequence`:

* **Non-natural residues written as three-letter codes.** `H-Aib-QGTF...` upper-cased and
  stripped becomes `HAIBQGTF...` — one residue silently expanded into three, lengthening the
  sequence, inventing residues and corrupting both the ESM-2 embedding and the HMMER cluster.
  Aib appears in 52 records and Sar in 46. Hyphen-delimited codes (Aib, Nle, Orn, Dab, Sar, Cha,
  Abu, Nal, MeLeu, hArg) are now collapsed to a single `X`. A bare `AIB` inside a one-letter
  sequence is left alone.
* **Non-standard residue letters.** 218 records (16%) contain O, U, X, Z or B — ambiguity codes
  or PEPlife2's letters for non-natural residues. ESM-2 accepts all of them as rare tokens, so
  nothing crashes; it is a silent quality gradient. Now recorded in a
  `sequence_nonstandard` column rather than left invisible.
* **Stereochemistry, still lost.** Lowercase residues mark D-amino acids — **30165332**'s DV1 is
  written entirely in lowercase and the body confirms it is "a novel D-peptide antagonist of
  CXCR4". Upper-casing makes an all-D peptide indistinguishable from its all-L enantiomer, while
  protease resistance is exactly what the D-configuration buys. `Chirality` carries the
  annotation but the sequence does not.

## Pipeline

Identical to mouse, with `--species rat`:

```bash
python fetch_peplife2_api.py    --species rat   # defines the cohort
python build_rodent_cohort.py   --species rat   # + PMC full-text retrieval
python render_fulltext.py       --dir rat
python fetch_pubmed_metadata.py --species rat
python fetch_peplife2_papers.py --species rat   # publisher PDFs; host currently down
python render_pdf_text.py       --species rat
python subset_by_nuance.py      --species rat
python build_study_dossier.py   --species rat
python subset_in_vivo_half_life.py --species rat
```

Then, for the combined dataset:

```bash
python build_halflife_training_set.py --species mouse rat
```

`build_halflife_training_set.py` reads the two per-species `in_vivo_half_life.csv` tables
directly, 513 mouse and 804 rat records after the same live-dosed half-life filters.
`species` is the curated animal assignment, so two records sourced from the rat file are
correctly labeled mouse. The result is distinct from `rodent/halflife_meta_with_split.csv`,
which applies additional phase, sequence and duplicate filters for model training.

## File inventory

Same layout as `../mouse`, minus the artefacts that need text this cohort does not have:

| path | what it is |
|---|---|
| `peplife2_rows.csv` | the 1149-record cohort, with `cohort_source` and `api_query` |
| `fulltext_retrieval_status.csv` | per-PMID PMC outcome |
| `fulltext/`, `fulltext_txt/` | 57 BioC deposits and their rendered text |
| `pdf_txt/` | 4 papers rendered from a mirrored publisher PDF |
| `extractions/` | 50 curated per-paper JSON files |
| `nuance_subset.csv` | per-PMID setting/route/matrix nuance |
| `records_provenance.csv` | one row per cohort record, with design and dosing classification |
| `studies.csv`, `STUDIES.md` | per-publication dossier |
| `sequences.csv` | per-sequence summary |
| `in_vivo_half_life.csv` | the 804-record live-dosed analysis set |

The current final training set and its audit trail live in `../rodent/strict/`, not here:

| path | what it is |
|---|---|
| `../rodent/strict/halflife_meta_with_split.csv` | 320 current modelling rows with `species`, `cluster` and `split` |
| `../rodent/strict/halflife_cohort_reconciliation.csv` | current model-row reconciliation; no additional rows are removed after the strict-set build |
| `../rodent/strict/hmmer/all.fasta`, `phmmer.domtbl` | the all-vs-all search used for the current cluster assignments |

The 1,212-row pooled training table described in older revisions is no longer present in this
checkout and is not used for the reported model or figures.

No `pk_parameters.csv` or `pk_figures/` yet — `build_pk_parameters.py` and
`extract_pk_figures.py` both need the publisher PDFs.

## Known caveats

1. **81% of records are uncurated**, because their papers are not retrievable. Retrieval, not
   effort, is the limit.
2. **`route_class` is "not reported" for 389 of 804 records (48%)** — PEPlife2 leaves `Dose
   Route` blank far more often here than for mouse. Several studies whose IV and SC arms differ
   several-fold are therefore indistinguishable in the database; 28605180, 30794654, 34393794 and
   38721043 are documented instances. This caveat describes the SOURCE cohort and still holds
   of it. Downstream it is addressed rather than inherited: the strict set resolves the blank
   routes against the primary sources and records the evidence level for each
   (`MuriVerse/data/curation/route_curation.csv`), and the model table keeps only records whose
   route is parenteral and established. No modelling row reads "not reported".
3. **`pk_phase` is "unlabelled" for 647 of 804 records (80%)**. PEPlife2 carries the phase only
   as free text inside the value, and most rat records have none, so two-compartment α and β
   values cannot generally be separated.
4. **The median half-life is 0.9 h against mouse's 4.0 h.** Some of that is real allometry, but
   part is composition: the rat cohort is richer in unmodified peptides (356 of 804) and in
   short radiotracer studies. Do not read the difference as a pure species effect.
5. **Two records in this file carry `species == "mouse"`** — 31818212's mislabelled pair, kept
   here so the correction is auditable but assigned to the right animal for pooling.
