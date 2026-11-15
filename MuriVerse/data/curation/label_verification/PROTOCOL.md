# Primary-source verification of the 213 intravenous half-life labels

Reviewer 1, major 1. One ledger row per record (`verdicts/<source>.json`, a JSON array).

## Inputs

* `worklist_all.csv`: one row per record. `label_hours` is the value being checked. The
  `peplife2_*` columns are the database transcription. They tell you WHICH peptide, species,
  route and value to look for. They are never evidence.
* `MuriVerse/data/curation_cache/label_sources/<source_file>.txt`: the source text (article body
  or patent). Articles from PDFs may have scrambled two-column order or split table cells.

## Rules

1. Evidence comes only from the source text. Not from the worklist, not from memory of the drug.
2. Article abstracts are not evidence. Text before `=== BODY ===` (or inside the
   `<<<ABSTRACT ... >>>` fence) is ignored. A fact that is only in the abstract is `not_located`.
3. A patent claim or prophetic example is not a measurement. Only a reported measured value counts.
4. Every non-null verdict needs a verbatim quote of at most 40 words, copied exactly.
5. If a number cannot be read with confidence (garbled table), say so; do not guess.
6. Do not edit any file except your own verdict JSON files.

## Per-record fields

| field | values |
|---|---|
| `id` | worklist id |
| `label_verdict` | `confirmed` (source value within 15% of `label_hours`, same molecule and species), `corrected` (source gives a different value for this molecule and arm; give `corrected_hours`), `phase_not_elimination` (the matching value is a distribution, absorption, alpha-phase, MRT or activity half-life, not terminal elimination), `ambiguous_arm` (several arms or doses and nothing ties this record to one of them), `not_located` (peptide or value cannot be found in the text) |
| `source_value_text` | the value exactly as the source writes it, with unit |
| `source_value_hours` | numeric, converted to hours, or null |
| `half_life_phase_in_source` | `terminal`, `elimination`, `unspecified`, `alpha`, `absorption`, `mrt`, `other` |
| `source_quote` | verbatim, at most 40 words |
| `source_location` | e.g. `Table 2`, `Results, para 3`, `Example 5` |
| `species_in_source` | species the value was measured in |
| `species_matches` | true or false |
| `route_verdict` | `iv_confirmed_arm` (the source ties this value to an intravenous arm), `iv_confirmed_study_level` (the source used only intravenous dosing, or names iv among several routes without tying this value), `iv_not_stated`, `not_iv` (the source says another route for this value) |
| `route_quote` | verbatim, at most 40 words, or null |
| `corrected_hours` | number or null |
| `notes` | one short line, or empty |

`confirmed` needs the molecule, species, value and phase to agree. If the source lists the value
only as a range or a mean of arms, say so in `notes`.

## Evidence tiers assigned afterwards (not by the reader)

* `body_verified`: label confirmed and route iv confirmed from source body text.
* `value_only`: label confirmed, route not stated in the source.
* `abstract_only`: no body text was available; only the abstract was read (assigned separately).
* `unverifiable`: no source text could be retrieved.
