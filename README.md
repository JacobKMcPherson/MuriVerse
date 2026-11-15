# MuriVerse

Code and data for **MuriVerse**, a curated mouse and rat in-vivo peptide half-life dataset and
the model trained on it, applied to colistin and polymyxin B.

This repository holds everything needed to inspect the curated evidence, rebuild the modelling
set and model, and regenerate the figures and tables of the manuscript.

> Paper: _citation and DOI to be added_

## Layout

| path | contents |
|---|---|
| [`MuriVerse/`](MuriVerse) | Pipeline scripts, curated cohorts (`mouse/`, `rat/`, `rodent/`, `rodent/strict/`), source data and curation ledgers (`data/`). Start with [`MuriVerse/README.md`](MuriVerse/README.md) |
| [`figures/`](figures) | R scripts (`code/`), figures (`figs/`) and the tables behind them (`tabs/`) |
| [`results/`](results) | Small frozen inputs to the figure scripts (`supp_fig2/`, `database_comparison/`); generated model datasets are written to `results/half_life_muriverse/` |
| [`docs/`](docs) | Data provenance and quality notes |
| [`scripts/`](scripts) | `fetch_tokenizer.sh` |
| [`DATA_SOURCES.md`](DATA_SOURCES.md) | Where every input comes from, and what is and is not redistributed |

## Reproduce

```bash
# 1. environment: Python 3.11, HMMER (phmmer on PATH), R >= 4.3 (packages in requirements.txt)
pip install -r requirements.txt

# 2. the PeptideCLM tokenizer (from the upstream PeptiVerse repository, pinned and checksummed)
bash scripts/fetch_tokenizer.sh

# 3. rebuild everything downstream of the curated cohort
python MuriVerse/build_strict_peptide_set.py                          # -> MuriVerse/rodent/strict/
python MuriVerse/build_halflife_training_set.py --tier strict         # -> results/half_life_muriverse/
python MuriVerse/evaluate_halflife_model.py --cohort rodent --tier strict --route iv --jobs 16 --nested   # ~25 min; without --route iv: the route-pooled model
Rscript figures/code/run_all.R                                        # -> figures/figs, figures/tabs
```

Step 3 downloads ESM-2 and PeptideCLM weights from Hugging Face. The committed
`MuriVerse/rodent/strict/` tables and saved models let you skip it and run only the figure
scripts. [`MuriVerse/README.md`](MuriVerse/README.md) describes each stage, including the
network-bound cohort-building steps whose outputs are committed.

## Revision analyses (external review, round 1)

**Final model: 178 post-audit intravenous records**, 150 rat and 28 mouse, from 76 sources.
Fixed-configuration Spearman rho is **0.445 (95% CI 0.153-0.563)**; nested feature-configuration
selection gives **0.297 (0.096-0.450)**. These are grouped cross-validation estimates, not
independent test-set performance. The endpoint is source-reported circulating PK half-life
after intravenous dosing, not uniformly verified terminal elimination half-life.
See [the final-model release record](docs/final-model-release.md) and `SHA256SUMS-postaudit`.
Older Hugging Face weights and the demo are not the final model until the separately prepared
hosting update is published; this GitHub release is authoritative.

The revised manuscript adds a primary-source audit of the 213 intravenous-labelled records, sensitivity
analyses on source-verified subsets, a post-audit model, and comparisons the reviewers asked for. Scripts
are in `MuriVerse/revision_analysis/`, their outputs in `results/revision_analysis/`, and the frozen
inputs they read (the shipped model's out-of-fold predictions, fold assignments and feature manifest) in
`results/revision_analysis/snapshot/`.

| Manuscript item | File | Script |
|---|---|---|
| Supplementary Table 2 (permutation importance) | `figures/tabs/supp_table6_feature_permutation_importance.csv` | `feature_permutation_iv.py` |
| Supplementary Table 3 (ensemble against its members) | `figures/tabs/supp_table7_ensemble_member_comparison.csv` | `model_member_baselines.py`, `model_member_paired.py` |
| Supplementary Table 4 (label-verification ledger) | `MuriVerse/data/curation/halflife_label_verification.csv` | `consolidate_label_verification.py` |
| Supplementary Table 5 (cohort sensitivity) | `figures/tabs/supp_table8_cohort_sensitivity.csv` | `cohort_sensitivity.py`, `verification_scoring.py` |
| Post-audit model and its performance | `MuriVerse/rodent/strict/model_iv_postaudit/`, `results/revision_analysis/clean_cohort_*.csv`, `figures/tabs/table31b_muriverse_halflife_nested_selection_postaudit.csv` | `clean_cohort_refit.py`, `refit_postaudit_model.py`, `nested_postaudit.py` |
| Ranking in thirds | `results/revision_analysis/ranking_utility.csv` | `ranking_utility.py` |
| Short-peptide grouping and sensitivity | printed | `short_peptide_grouping.py`, `short_peptide_sensitivity.py` |

The file names `supp_table6/7/8` predate the manuscript numbering (Supplementary Tables 2, 3 and 5); the
older `supp_table1`-`supp_table5` files in `figures/tabs/` are the curation-gate tables.

Run order: `consolidate_label_verification.py` (needs the source texts, which are not redistributed:
re-fetch them with `fetch_curation_sources.py`; the ledger is provided as built), then
`cohort_sensitivity.py` and `verification_scoring.py`, `model_member_baselines.py` then
`model_member_paired.py`, `feature_permutation_iv.py`, `ranking_utility.py`, `metrics_iv.py`,
`clean_cohort_refit.py` then `refit_postaudit_model.py` and `nested_postaudit.py`,
`short_peptide_grouping.py`, `short_peptide_sensitivity.py`, and last `build_round1_supplementary.py`.
`MuriVerse/audit_figure_numbers.py` checks the numbers in the figures and tables against these outputs.

Two intravenous models are distributed. `MuriVerse/rodent/strict/model_iv/` is the **pre-audit** model,
fitted on all 213 intravenous-labelled records, including 37 the audit found discordant with their primary
source; it is retained unchanged for reproducibility.
`MuriVerse/rodent/strict/model_iv_postaudit/` is the **post-audit** (distributed final) model, the same
71-column configuration refit on the 178 source-verified records; the directory also holds the two structural
configurations refit on the same records. The polymyxin case study (`predict_polymyxins.py`) applies the
post-audit models by default; `--model-suffix _iv` applies the pre-audit models.

The reported values reproduce, to the reported precision, only with Python 3.11, NumPy 2.4.6 and
XGBoost 3.2.0 (`requirements.txt`); a different environment moved individual Spearman estimates by up
to 0.02.

## Not included

- The PeptiVerse code, models and training data. MuriVerse uses its tokenizer and builds on its
  encoders and conventions; cite PeptiVerse (Zhang et al., bioRxiv 2025,
  doi:10.64898/2025.12.31.697180).
- Mirrored publisher PDFs and full text, and the third-party database files; see
  [`DATA_SOURCES.md`](DATA_SOURCES.md). The curated extractions and the verbatim quotes behind
  each curation verdict are included.
- Figure 1 is a study schematic and is not generated by code.

## License

Code is released under the [MIT License](LICENSE). Data and figures created for this work (the
curation ledgers, extractions, cohort tables, model tables and figures) are released under
[CC BY 4.0](LICENSE-DATA). Records derived from third-party sources (PEPlife2, PubMed metadata)
remain under those sources' own terms; see [`DATA_SOURCES.md`](DATA_SOURCES.md).
