# Half-life data-quality implementation

Branch: `improve-halflife-data-quality`

## Decisions implemented

### Label evidence

The training metadata now has an explicit quality contract:

- **A** — exact PEPlife2-ID link to a primary-source PK arm, agreeing numeric half-life, exact label.
- **B** — exact PEPlife2-ID link to a primary-source PK arm, agreeing numeric half-life, but the source value is approximate, bounded, or a range.
- **C** — database transcription without exact arm-level primary-source verification.
- **D** — primary-source conflict or an arm whose half-life cannot be resolved.

Only A+B are the gold label set. C is retained for sensitivity analyses; D requires adjudication.

No narrative match or heuristic extraction can upgrade a row from C to A/B.

### Exact quantitative labels

The current manifest contains 447 exact and 23 non-exact labels
(9 approximate, 8 ranges, 5 upper bounds, 1 lower bound).

`exact_label` identifies the 447 quantitative regression labels. The 23 non-exact records
remain available as a sensitivity set and are not silently converted into equally certain
ground truth.

### Full-text arm evidence

MuriVerse already joins `pk_observations` to the model table only by exact PEPlife2 ID.
That exact-arm evidence remains authoritative over heuristic route/matrix/assay/strain/sex
recovery. The new A/B/C/D grading makes that precedence explicit: heuristic context can
improve covariates, but it cannot create a gold label.

### Structure fidelity

The former `structure_tier` is mapped to:

- `exact_whole_molecule`
- `exact_peptide_partial_modifier`
- `peptide_core_only`
- `unresolved`

`smiles_for_model` is populated only for `exact_whole_molecule`. The original `smiles`
is retained for provenance/audit. This creates an explicit model-safe structure column so an
incomplete string is never presented as though it were the administered molecule.

Under the current accounting this mapping corresponds to 406 whole molecules, 5 partly
resolved structures, and 56 peptide-core structures with an unresolved modification, across
the 467 kept records. The model table now REQUIRES the whole-molecule tier, so all 320 of its
rows are `exact_whole_molecule` and `smiles_for_model` is populated for every one of them; the
other two tiers survive only in the kept file, which is where Figure 2's panels n and o draw
them from.

### Observation identifiers

The quality metadata carries five distinct identifiers:

- `molecule_id` — canonical whole-molecule structure when available; otherwise a conservative
  annotation composite that does not claim the peptide-core SMILES is the whole molecule.
- `experiment_id` — study + species + route + matrix + assay + strain + sex + dose +
  formulation + sampling design, preferring exact-arm evidence where available.
- `study_id` — source/publication identifier.
- `arm_id` — study + exact PEPlife2 source-record identifier.
- `homology_cluster` — the existing HMMER cluster.

These keep legitimate repeated PK experiments while allowing later loss balancing by study,
molecule, or peptide family.

### Priority-source curation

`build_priority_source_curation_worklist.py` ranks sources from the **current strict model
table**, selects the largest ten, and writes their record counts, current A/B/C/D status,
exact-label count, structure-fidelity count, available local full text/extraction files, and
the remaining PEPlife2 IDs requiring arm-level curation.

The current manuscript accounting says those ten sources contribute 126/320 records.
The largest, EP2277890A1, contributes 62 kept records and 61 modelling rows, and is a
patent/non-PMID source. Its patent text HAS now been read for the dose route - it states
"each was administered to the rat intravenously (tail) or subcutaneously (back) in a dose of
0.1 mg/kg", which is recorded in `MuriVerse/data/curation/route_curation.csv` - but that
settles the route only. The label itself remains C until the patent's own PK tables are
curated arm by arm, rather than being promoted by inference.

## Commands

After the normal strict training-data build:

```bash
python MuriVerse/apply_halflife_data_quality.py
python MuriVerse/build_priority_source_curation_worklist.py
```

Outputs:

- `MuriVerse/rodent/strict/data_quality_manifest.json`
- `MuriVerse/rodent/strict/halflife_quality_metadata.jsonl`
- `MuriVerse/rodent/strict/priority_source_curation_worklist.json`

## Still to execute

The repository connection exposes the generated CSVs as Git LFS pointers and does not expose
a runnable model environment. Therefore the following must not be reported as completed until
the repository is executed in its normal conda environment:

1. primary-text curation of all remaining rows in the dynamically identified top ten sources;
2. the 447-exact-label headline rerun;
3. non-exact-label sensitivity predictions for the other 23 rows;
4. A+B gold-model evaluation once enough independent A/B groups exist;
5. study-balanced and molecule-balanced loss experiments.

No new performance number should be quoted before those runs are completed.
