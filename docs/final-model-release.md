# Final Post-Audit Model

Release scope: manuscript v10 reconciliation, 2026-11-15.

## Primary Artifact

The final model is [model_iv_postaudit/shipped.joblib](../MuriVerse/rodent/strict/model_iv_postaudit/shipped.joblib),
fitted on **178 source-verified intravenous records** (150 rat, 28 mouse; 76 sources).
The [membership ledger](../MuriVerse/rodent/strict/model_iv_postaudit/training_membership.csv)
records each inclusion decision, corrected value and corrected species. Of the original
213 records, 170 are full-text verified, three abstract verified, three have corrected
values and two have corrected species. Thirty-five records are excluded: 32 of the 37
discordant records, and the three records whose value is confirmed without a stated
intravenous route. Five discordant records are corrected and retained.

The endpoint is **source-reported circulating pharmacokinetic half-life after intravenous
dosing**, in hours (modeled in log10 hours), not uniformly verified terminal elimination
half-life. Predictions are not validated exposure estimates, dose recommendations or
prospective prediction intervals.

## Evaluation

| Evaluation on 178 records | Spearman rho | 95% group-bootstrap CI |
|---|---|---|
| Fixed 71-feature configuration | 0.445 | 0.153-0.563 |
| Nested feature-configuration selection | 0.297 | 0.096-0.450 |

Five folds and three seeds jointly withhold homology clusters and publication sources
(50 connected groups). Point estimates and bootstrap intervals summarize the mean of
per-seed scores. Other development decisions are not nested. These are cross-validation
estimates, not training scores or an independent test of the all-record refit.
Fixed-configuration MAE is 0.550 log10 h; RMSE is 0.703 log10 h.

The [primary performance table](../figures/tabs/pub-tab__table1_benchmark_metrics.csv)
and [nested results](../figures/tabs/table31b_muriverse_halflife_nested_selection_postaudit.csv)
are authoritative. In the public repository, fixed-configuration predictions, fold
assignments and full metric intervals are in `results/revision_analysis/clean_cohort_*.csv`.

## Historical Analyses

`MuriVerse/rodent/strict/model_iv/` is the **pre-audit** model trained on all 213
intravenous-labelled records, including 37 subsequently found discordant. Its
0.574 fixed-configuration and 0.460 nested estimates are not final-model performance.
The corrected historical nested CI is 0.267-0.548, not 0.293-0.636.
Figure 3, Table 2, feature-importance and subset analyses are explicitly pre-audit
development evidence. The route-pooled model is retained for historical transfer diagnostics.

## Polymyxin Application

Figure 4, Supplementary Figure 3 and Tables 5-7 use post-audit models. Median predictions
for colistin and polymyxin B are 0.92 and 0.74 h in mouse, and 0.78 and 0.62 h in rat.
Relative to intravenous references, fold errors are 2.05 and 1.26 in mouse, 2.05 and
3.06 in rat. A constant baseline is closer in three of four drug-species cells.
This is exploratory application, not independent external validation.

## Reproduction And Hosting

Use the pinned model environment (Python 3.11, NumPy 2.4.6, XGBoost 3.2.0).
In the public repository run:

```bash
python MuriVerse/revision_analysis/clean_cohort_refit.py
python MuriVerse/revision_analysis/refit_postaudit_model.py
python MuriVerse/revision_analysis/nested_postaudit.py
python MuriVerse/build_publication_tables.py
python MuriVerse/audit_figure_numbers.py
```

`predict_polymyxins.py` selects the post-audit models by default; `--model-suffix _iv`
explicitly selects the pre-audit archive. The public root README gives installation
and tokenizer prerequisites. Evaluation scripts use the published embedding cache.

The release includes the post-audit models and membership ledger on GitHub and the default `postaudit` weights on Hugging Face. The unchanged 213-record model remains available as the explicitly labeled pre-audit archive. Verify `SHA256SUMS-postaudit` on GitHub and `SHA256SUMS` on Hugging Face before use. Private manuscript and review documents and source-paper full texts are excluded from the public export.
