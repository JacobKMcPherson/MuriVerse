You are extracting pharmacokinetic evidence from open-access article FULL TEXT.

## Absolute rules
1. **Body only.** Each file contains a block fenced by `<<<ABSTRACT - DO NOT USE AS AN EXTRACTION SOURCE>>>` ... `<<<END ABSTRACT>>>`. You MUST NOT take any extracted value from inside that block. Use only the text after `=== BODY ===` (Methods, Results, Discussion, figure captions, table captions/contents). If a fact appears ONLY in the abstract, record it as `"not reported in body"`.
2. **No inference from outside knowledge.** If the body does not state something, write `"not reported"`. Never fill a value from what you know about the drug. Never guess a sample size, dose, or PK value.
3. **Every non-null value needs evidence.** Give a short verbatim quote (<= 25 words) plus the section it came from.
4. The `=== PEPLIFE2 ROWS FOR THIS PMID ===` header at the top is database metadata, NOT the article. Use it only to know which peptides to look for. Never use it as an evidence source.

## Task
For each assigned PMID, read `MuriVerse/<species>/fulltext_txt/<PMID>.txt` in full and write
`MuriVerse/<species>/extractions/<PMID>.json` with exactly this schema:

Some studies are not in PMC and are instead rendered from the publisher PDF into
`<species>/pdf_txt/<PMID>.txt`. Those files open with a `=== SOURCE: publisher PDF ===` header.
They follow the same rules with two extra cautions:
- Two-column reading order is reconstructed heuristically and table cells may be split across
  lines. **Verify any value that looks garbled**, and if a number cannot be read with confidence,
  record `"not reported"` rather than a guess.
- A few PDFs have no separable abstract; their header says `ABSTRACT NOT FENCED`, which means the
  abstract is still sitting at the top of the body. Identify it yourself and do not quote from it.

Set `"text_source"` at the top level to `"PMC BioC"` or `"publisher PDF"` to match which folder
the text came from.

```json
{
  "pmid": "12345678",
  "fulltext_read": true,
  "text_source": "PMC BioC",
  "study_types": ["in vivo"],
  "study_type_evidence": [{"type": "in vivo", "quote": "...", "section": "METHODS"}],
  "in_vivo_species": [{"species": "mouse", "strain_sex_age": "male C57BL/6, 8 wk", "quote": "...", "section": "METHODS"}],
  "ex_vivo_origin": null,
  "ex_vivo_detail": null,
  "in_vitro_detail": null,
  "clinical_detail": null,
  "agents": [
    {
      "name": "DualAG",
      "role": "experimental",
      "modality": "peptide",
      "sequence": "HSQGTFTSDYSKYLDSRRAQDFVQWLMNTKRNRNNIA",
      "sequence_source": "Fig. 1A sequence alignment",
      "smiles": null,
      "chemical_identifier": null,
      "group_n": "n = 8 per group",
      "group_n_quote": "...",
      "dose": "10 nmol/kg",
      "route": "subcutaneous",
      "dose_route_quote": "...",
      "half_life": "1.7 h",
      "half_life_phase": "elimination",
      "half_life_absorption": "not reported",
      "half_life_elimination": "1.7 h",
      "half_life_matrix": "mouse plasma",
      "cmax": "not reported",
      "tmax": "not reported",
      "auc": "not reported",
      "pk_quote": "...",
      "pk_section": "RESULTS"
    }
  ],
  "pk_observations": [
    {
      "peplife2_id": "1234",
      "agent_name": "DualAG",
      "species": "mouse",
      "strain": "C57BL/6",
      "sex": "female",
      "age": "8 weeks",
      "body_weight": "not reported",
      "group_n": "n = 8 per group",
      "dose": "10 nmol/kg",
      "route": "subcutaneous",
      "formulation": "saline",
      "formulation_class": "solution",
      "matrix": "plasma",
      "sampling_design": "serial sampling through 24 h",
      "bioanalytical_method": "LC-MS/MS",
      "half_life": "1.7 h",
      "half_life_phase": "unspecified",
      "cmax": "not reported",
      "tmax": "not reported",
      "auc": "not reported",
      "clearance": "not reported",
      "volume_of_distribution": "not reported",
      "mrt": "not reported",
      "bioavailability": "not reported",
      "evidence": [
        {"field": "dose, route, group_n", "quote": "...", "section": "METHODS"},
        {"field": "formulation", "quote": "...", "section": "METHODS"},
        {"field": "half_life, matrix", "quote": "...", "section": "RESULTS"}
      ]
    }
  ],
  "concentration_time_figure": {
    "available": true,
    "locator": "Fig. 2B",
    "caption_quote": "...",
    "note": "plasma concentration vs time after single s.c. dose"
  },
  "bioanalytical_method": "RP-HPLC with ESI-MS detection",
  "bioanalytical_quote": "...",
  "bioanalytical_section": "METHODS",
  "notes": "free-text caveats, e.g. PK reported only for one of three peptides"
}
```

## Field guidance
- `study_types`: array, any of `"clinical"`, `"in vivo"`, `"ex vivo"`, `"in vitro"`. A paper may have several; list all that the body supports.
  - `clinical` = human subjects/patients dosed in a trial.
  - `in vivo` = intact living animal (or human) dosed; fill `in_vivo_species` with every species dosed.
  - `ex vivo` = tissue/blood/serum/plasma taken out of a body and then treated (e.g. stability in mouse serum, isolated perfused organ). Set `ex_vivo_origin` to `"human"` or `"animal"` (name the animal in `ex_vivo_detail`).
  - `in_vitro` = cell culture, buffer, purified enzyme, cell-free system.
  - Serum/plasma stability assays are **ex vivo**, not in vitro — classify them that way and say so in `ex_vivo_detail`.
- `agents`: one object per distinct administered/tested substance, **including controls** (vehicle, saline, PBS, scrambled peptide, parent/unmodified peptide, positive comparator). Set `role` to `"experimental"`, `"control"`, `"vehicle"`, `"comparator"`, or `"reference standard"`.
  - `modality`: `"peptide"`, `"peptide conjugate"`, `"small molecule"`, `"protein/biologic"`, `"nanoparticle/formulation"`, `"vehicle"`, `"other"`.
  - `sequence`: one-letter amino acid sequence **only if the body/figure states it**. Preserve stated modifications in a parenthetical. Otherwise `"not reported"`.
  - `smiles`: only if a SMILES string is literally printed in the body; otherwise `null`. Put any CAS number, catalogue number, or IUPAC name in `chemical_identifier`.
  - `group_n`: animals (or subjects) per group for the PK/study arm, verbatim style e.g. `"n = 5"`, `"n = 3 per timepoint"`. If only a total is given, say so.
  - `half_life`/`cmax`/`tmax`: exactly as reported, with units and error term, e.g. `"1.7 ± 0.3 h"`. Use `"not reported"` when absent. If a value exists for several matrices/routes, give the most specific and note the rest in `pk_quote`.
  - `half_life_phase`: which phase the headline `half_life` is, one of `"absorption"`, `"distribution/alpha"`, `"elimination"`, `"terminal/beta"`, `"unspecified"`. Use `"unspecified"` when the paper reports a single half-life without naming a phase — do **not** infer one.
  - `half_life_absorption` / `half_life_elimination`: fill both only when the paper fits a two-phase model and reports them separately. These are **half-lives in time units**, not the rate constants `ka`/`ke`. If the paper reports only a rate constant (e.g. `"ka = 0.35 h⁻¹"`), record it verbatim in `pk_quote` and leave the half-life `"not reported"` — do not do the `t½ = ln2/k` conversion yourself.
    - Beware the reverse error: several papers print a half-life under an `h⁻¹` column header (see PMID 36630826 Table 2). Trust the parameter name and the body text over the unit string, and say which you used in `pk_quote`.
    - An absorption half-life measures entry from the injection depot, not persistence in circulation. Never put it in `half_life`.
  - `auc`: area under the concentration-time curve as reported, with units.
- `pk_observations`: one object per **species × treatment arm × PEPlife2 record**. This is
  the row-level table used for training and is separate from the molecule-level `agents`
  summary. Preserve every arm, even when the same molecule appears at several doses, routes,
  formulations or sampling plans. Set `peplife2_id` only as a join key from the header; it is
  not evidence. If a row cannot be linked to exactly one arm from the body, leave the ID
  blank and do not propagate its conditions to a training row. Keep the reported text and
  units; do not infer age from body weight, dose from dose volume, or formulation from a
  product name. Use `half_life_phase: "unspecified"` unless the body names the phase. Each
  populated field needs a `quote` and `section` entry under `evidence`. The PK fields include
  Cmax, Tmax, AUC, clearance, volume of distribution, MRT and bioavailability when reported;
  write `"not reported"` otherwise. Keep `formulation` verbatim and assign
  `formulation_class` only when the body explicitly identifies a solution/vehicle, depot,
  particle or other formulation class.
- `concentration_time_figure.available`: `true` only if a figure (or figure panel) plots concentration/level/%-remaining against time. Give the panel in `locator`. Serum-stability %-remaining vs time counts — say so in `note`.
- `bioanalytical_method`: e.g. LC-MS/MS, ELISA, RIA, RP-HPLC/UV, fluorescence, scintillation counting of radiolabel, IVIS fluorescence imaging. `"not reported"` if absent.

## Output
- Write one JSON file per PMID. Valid JSON, no markdown fences, no trailing commas.
- If a PMID's file has no usable body text, still write the JSON with `"fulltext_read": false` and a `notes` explanation.
- Do not modify any other file. Do not write summaries to disk.

Report back only: the list of PMIDs you wrote, and any PMID where the body carried no PK data at all.
