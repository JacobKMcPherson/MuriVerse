# Mouse and rat training data index

The **canonical training table** is [`../halflife_meta_with_split.csv`](../halflife_meta_with_split.csv).
The model trains on its pooled mouse and rat rows so homology clusters shared across
species stay in one split. The two CSV files here are small species-specific indexes,
not separate model inputs or copies of the molecular data.

The current table has **320 rows** (74 mouse and 246 rat), from 1,159 eligible records: 467
survive the strict gates, and of those the model table keeps the parenterally dosed records
whose SMILES is the whole dosed molecule and whose bioanalytical method a source names.
Every row is a single peptide chain of at most 50 residues with no attached half-life
partner, an elimination half-life after a parenteral dose, a bioanalytical method a source
names, and a SMILES that was decoded back to its sequence and is the whole dosed molecule. A
record that fails one of the structural gates is listed with its reason in
[`../exclusions.csv`](../exclusions.csv); one that passes those but fails a cohort requirement
is listed in [`../halflife_cohort_reconciliation.csv`](../halflife_cohort_reconciliation.csv)
with the step that removed it. `manifest.json` records the current split counts
and source checksum. Repeated observations stay in the canonical table; candidates are
annotated and listed in [`../repeat_observation_candidates.csv`](../repeat_observation_candidates.csv)
for review with source-arm context.

| File | Use |
|---|---|
| [`mouse_index.csv`](mouse_index.csv) | Find mouse rows in the canonical metadata table |
| [`rat_index.csv`](rat_index.csv) | Find rat rows in the canonical metadata table |
| [`manifest.json`](manifest.json) | Counts, metadata checksum, and index schema |

`metadata_row` is the zero-based row position in the canonical CSV; `id` is the
join key for this build. Use `id` to join with the pooled Arrow datasets in
`results/half_life_muriverse/` and with
`../oof_predictions.csv`. `split` is the saved cluster-preserving train/validation
assignment (257 train, 63 validation rows in this build; mouse 60/14, rat 197/49).
This saved assignment keeps homology clusters intact; it does not hold out source identifiers.
Four source identifiers occur on both sides, so the per-split source counts in
[`table12_muriverse_halflife_split.csv`](../../../../figures/tabs/table12_muriverse_halflife_split.csv)
are not additive. Reported model performance comes from grouped cross-validation in
`evaluate_halflife_model.py`, not from that single assignment.

The source cohort is [`../in_vivo_half_life.csv`](../in_vivo_half_life.csv). Its rows are
filtered into the training table without automatic deduplication; see
[`../halflife_cohort_reconciliation.csv`](../halflife_cohort_reconciliation.csv) for the
step-by-step row counts. `structure_tier` describes how completely the peptide SMILES
captures the reported structure: `complete`, `partly resolved` or `modification not drawn`.
Every row of this table is `complete` - the model requires it - so the column is constant
here and varies only in the source cohort.
Every row must have a nonempty, sequence-verified peptide SMILES; the tier records whether
all reported modifications and structural features are represented. The generated
`../pk_observations.csv` retains arm-level full-text evidence and links it to rows only by
PEPlife2 ID; 3 of the current model rows have a full-text arm link and label match. Unknown
labels/conditions remain unknown, and label evidence tier is audit metadata, not a model input.
The model's `chirality` feature is derived from the sequence stereochemistry verified by the
SMILES round trip; PEPlife2's original summary field is retained separately as
`chirality_reported`. Three Apraglutide rows conflict between those sources and are treated as
mixed chirality in the model and Figure 2.

Regenerate these indexes after rebuilding the training table:

```bash
python MuriVerse/training_data_index.py
```

`build_halflife_training_set.py` also updates them during a pooled mouse and rat strict
build. The generator checks unique IDs and rejects a split that puts one homology cluster
on both sides.
