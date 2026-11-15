# Source curation ledgers

Two study-condition columns in the rodent half-life cohort read "not reported" or "unknown"
far more often than the sources warrant, because the PEPlife2 column they are derived from is
blank. These two files record what the PRIMARY sources say, so that the gap is closed from the
papers and patents rather than from the database's own prose.

Each file is keyed on the source identifier (a PMID, or a patent number), carries the verbatim
quote the verdict rests on, and names where that quote came from. Nothing here is inferred
from the database columns: `Test Sample` often states the route too, but it is a transcription,
and a transcription is what these ledgers exist to check rather than to trust.

`MuriVerse/halflife_condition_curation.py` reconciles both ledgers with the records, and its
docstring is the specification. Both are consumed by `build_strict_peptide_set.py` (route) and
`build_halflife_training_set.py` (assay).

## `route_curation.csv`

One row per source behind a record whose `route_class` read "not reported" - 84 sources
covering 241 records.

| column | meaning |
|---|---|
| `source` | PMID, or patent number for the sources that are patents |
| `evidenced_routes` | the routes the primary source evidences, `\|`-separated; empty when it evidences none |
| `evidence_level` | `primary_arm`, `primary_source` or `none` - see below |
| `evidence_source` | `fulltext`, `patent`, `abstract`, `mesh`, or a combination |
| `evidence_quote` | the sentence the verdict rests on, verbatim |
| `note` | why a source could not be resolved, or what qualifies the verdict |

`primary_arm` means the source names the route for the specific arm the record came from:
either the study used one route throughout, or the record's half-life was matched to its own
row in the source's PK table. US6656906 is the clearest case - its Table 10 reports
`t 1/2, terminal 2.02 2.00 2.46 1.86` across two SC and two IV dose groups, and all four
records match an arm exactly.

`primary_source` means the source names the routes it used and all of them are parenteral, so
the record is attributed at source level rather than arm level.

`none` means no primary source names a route for the PK arm. 30 sources are in this state.
Most are paywalled journal articles whose abstract and MeSH indexing are both silent; a few
are more interesting, and the `note` says which:

* **22286034** (15 records, the largest) - Elsevier full text unreachable.
* **1646465** - the database annotates i.p., which the sub-minute clearance the paper reports
  contradicts; no full text to settle it.
* **12083977** - a review of human clinical depot-leuprorelin PK, with no rodent dosing arm.
* **18085721** - the intravenous dose is the adenoviral vector, not the peptide.
* **WO2012138867A2** - the only "1.5 hours" in the document is a claim, not a measured rat
  half-life.
* **8100840** - the source reports anaesthetised toads, not a rodent.

A record whose route cannot be established is dropped from the model table, so supplying the
missing PDFs and re-running the stage is all that is needed to recover these.

## `assay_curation.csv`

One row per source behind a record whose `assay_method` bucketed to "unknown" and whose
`Assay` text could not yield one - 11 sources.

| column | meaning |
|---|---|
| `source` | PMID, or patent number |
| `assay_method` | the method the source states; empty when it states none |
| `evidence_source` | `fulltext` or `abstract` |
| `evidence_quote` | the sentence the verdict rests on, verbatim |
| `note` | what the database column said instead, or why it could not be resolved |

The ledger is only consulted where the `Assay` column came out "unknown". A method the column
names is read from the column, so `study_conditions.py` stays the single place that decides
what "LC/MS" means and the ledger holds only what the column could not say.

Four sources resolved. 26381279 is the instructive one: its `Assay` column says only "TLC",
which is a separation and not a detector, while the paper states that "the time dependent
amount of radiolabelled RD2 was measured in plasma and several organ homogenates by liquid
scintillation counting" - radiometric.

Unlike the route, an unresolved assay is not a reason to drop a record. The half-life is still
a half-life, and `unknown` is an honest level for it.

## `halflife_label_verification.csv`

One row per record of the 213-record intravenous-labelled cohort (reviewer 1, major 1; Supplementary
Table 4). Each label was read against the primary source: the full primary-source text (article body,
patent or open supporting information) or, where no full text could be had, the PubMed abstract,
with a verbatim quote for the value and for the route. Built by
`MuriVerse/revision_analysis/consolidate_label_verification.py` from the per-source
verdicts; every quote is re-checked against the source text (`quotes_checked` is true for all
213), and a blind second reading of 32 records agreed on all 32 label verdicts (values within 1%).

`evidence_tier` puts each record in exactly one category:

| `evidence_tier` | n | Meaning |
|---|---|---|
| `verified_body` | 170 | verified in full primary-source text: value, species, half-life phase and an intravenous route agree (167 for the arm, 3 at study level) |
| `verified_abstract` | 3 | the same, from the abstract only (full text not obtainable) |
| `value_only` | 3 | value confirmed, no intravenous route stated |
| `discrepant` | 37 | discordant with the source: database value not supported (`label_verdict` `not_located`, 9: 2 abstract only, 7 full text) or contradicted (28: different value 3, non-elimination phase 4, ambiguous arm 1, matching value for another species, arm or route, or cited rather than measured 20) |

`route_verdict` is the audit's own reading of the route for every record: `iv_confirmed_arm` 184,
`iv_confirmed_study_level` 13, `iv_not_stated` 13, `not_iv` 3. `corrected_hours` is filled only
for the three `corrected` records. `evidence_level` says whether the reading was of the full text
(`body`) or the abstract. The post-audit cohort (178 records: the 170, the 3 abstract-verified, the 3
corrected values and 2 species corrections) is listed in
`rodent/strict/model_iv_postaudit/training_membership.csv`.
