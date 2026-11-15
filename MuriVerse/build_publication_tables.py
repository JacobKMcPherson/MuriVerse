"""Render the manuscript's publication tables from the pipeline's own result tables.

The tables in docs/results.md were prose-embedded markdown. This script emits them as
standalone publication tables under figures/tabs/ with the `pub-tab__` prefix, so a table
quoted in the text has one generated source and cannot drift from it.

Every value here is read or re-derived from a committed result table, or from the
out-of-fold predictions, in exactly the way audit_figure_numbers.py checks it. Nothing is
transcribed from prose. Each table is written twice: a CSV carrying the full-precision
values, and a markdown rendering at the precision the manuscript quotes, for pasting.

    python MuriVerse/build_publication_tables.py

Two tables are derived rather than selected, and both follow the audit's definitions:

  * The route-restriction comparison (Table 2) scores the SAME 213 intravenous rows under
    the route-pooled fit and the intravenous fit, which is the only like-for-like reading.
  * The per-source polymyxin references (Table 5) take the fold error as the ratio of the
    two MEDIANS, not the median of the per-variant fold errors. The two differ wherever a
    drug's variants straddle the observed value.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TABS = ROOT / "figures" / "tabs"
STRICT = HERE / "rodent" / "strict"

SEED_COLS = [f"y_pred_cluster_PMID_seed{i}" for i in range(3)]
CUT_H = 6.0
PREFIX = "pub-tab__"

ANALYSIS = ROOT / "results" / "revision_analysis"
CAPTION_TABLE1 = "Table 1. Post-audit model: cross-validated performance on log₁₀ reported intravenous half-life for 178 source-verified records. Point estimates and 95% group-bootstrap confidence intervals summarize the mean of three per-seed scores over 50 connected homology-cluster and source groups, using 500 resamples of whole groups. The fixed configuration is the 71-feature configuration used in the final distributed model. Nested feature-configuration selection chooses among eight configurations of at most 100 columns using inner three-fold grouped cross-validation within each outer training fold. Only outer held-out records are scored. Other development decisions were not nested, so residual optimism may remain. These are cross-validation estimates, not scores on the final all-record refit or an independent external test set. AUC for the pre-stated six-hour screen is exploratory. The 213-record pre-audit benchmark is retained separately for reproducibility; it is not the final model."
CAPTION_TABLE2 = "Table 2. Ranking performance of compact and embedding-based feature configurations on the 213 intravenous-labelled records. Spearman ρ under five-fold grouped cross-validation, each value the mean of three seeds, for eleven of the twenty-nine configurations scored in the ablation (Figure 3g). Three grouping schemes are shown: held-out source identifiers, held-out homology clusters, and the headline scheme withholding connected groups of both. \"Retained\" is the ratio of the clusters-withheld to the sources-withheld value, the fraction of apparent ranking performance that survives withholding the peptide family. The shipped configuration is listed first; the others follow in descending order of the headline value. The SMILES-embedding configuration in the sixth row (ρ = 0.247) is the one whose post-audit refit (ρ = 0.292 on the 178 source-verified records) is applied to the polymyxin congeners"

# Markdown column formats, by table. A column absent here is written as-is.
FMT: dict[str, str] = {}


def rho(y: np.ndarray, p: np.ndarray) -> float:
    return float(spearmanr(y, p).statistic)


def auc_at_cut(y_log: np.ndarray, p: np.ndarray, cut_h: float = CUT_H) -> float:
    return float(roc_auc_score((10 ** y_log >= cut_h).astype(int), p))


def write(name: str, frame: pd.DataFrame, caption: str, round_map: dict) -> str:
    """Write one publication table as CSV plus a markdown rendering, and return the markdown."""
    csv_path = TABS / f"{PREFIX}{name}.csv"
    frame.to_csv(csv_path, index=False)

    shown = frame.copy()
    for col, places in round_map.items():
        if col in shown:
            signed = isinstance(places, str) and places.startswith("+")
            n = int(str(places).lstrip("+"))
            shown[col] = shown[col].map(
                lambda v, n=n, signed=signed:
                "" if pd.isna(v) else (f"{v:+.{n}f}" if signed else f"{v:.{n}f}"))
    header = "| " + " | ".join(shown.columns) + " |"
    rule = "|" + "|".join("---" for _ in shown.columns) + "|"
    body = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
            for row in shown.itertuples(index=False)]
    md = "\n".join([f"**{caption}**", "", header, rule, *body])
    print(f"[write] {csv_path.name}  ({len(frame)} rows)")
    return md


# ---------------------------------------------------------------------------------------
# Table 1. Species x route transfer
# ---------------------------------------------------------------------------------------
def species_route_transfer() -> tuple[pd.DataFrame, str]:
    t = pd.read_csv(TABS / "table36_muriverse_halflife_species_route_cells.csv")
    label = {("rat", "IV"): "Rat intravenous", ("rat", "EV"): "Rat extravascular",
             ("mouse", "IV"): "Mouse intravenous", ("mouse", "EV"): "Mouse extravascular"}
    rows = []
    for (sp, adm), cell in t.groupby(["species", "admin"], sort=False):
        by = cell.set_index("train")
        rows.append({
            "Test cell": label[(sp, adm)],
            "n": int(cell.n.iloc[0]),
            "Own-cell rho": by.loc["same cell", "spearman"],
            "Pooled rho": by.loc["pooled", "spearman"],
            "Cross-species rho": by.loc["cross-species", "spearman"],
            "Cross-route rho": by.loc["cross-route", "spearman"],
            "Cross-species bias": by.loc["cross-species", "bias"],
            "Cross-route bias": by.loc["cross-route", "bias"],
        })
    order = ["Rat intravenous", "Rat extravascular", "Mouse intravenous", "Mouse extravascular"]
    frame = pd.DataFrame(rows).set_index("Test cell").loc[order].reset_index()
    cap = ("Source table, not numbered in the manuscript. Cross-route transfer is offset, cross-species transfer is not. The same "
           "ensemble refit on four training pools per species x route cell; Spearman rho and "
           "signed bias in log10 hours, three-seed means.")
    return frame, cap


# ---------------------------------------------------------------------------------------
# Table 2. What the route restriction is worth, like for like
# ---------------------------------------------------------------------------------------
def route_restriction() -> tuple[pd.DataFrame, str]:
    oof = pd.read_csv(STRICT / "oof_predictions_iv.csv")
    pooled = pd.read_csv(STRICT / "oof_predictions.csv")
    # The pooled model restricted to, and aligned on, the same 213 intravenous rows.
    pi = pooled[pooled.id.isin(set(oof.id))].set_index("id").loc[oof.id].reset_index()
    y = oof.y_true.to_numpy(float)
    mouse = (oof.species == "mouse").to_numpy()

    def rmse(d: pd.DataFrame) -> float:
        return float(np.mean([np.sqrt(np.mean((y - d[c].to_numpy(float)) ** 2))
                              for c in SEED_COLS]))

    rows = [
        ("Spearman rho",
         float(np.mean([rho(y, pi[c]) for c in SEED_COLS])),
         float(np.mean([rho(y, oof[c]) for c in SEED_COLS]))),
        ("AUC (>= 6 h)",
         float(np.mean([auc_at_cut(y, pi[c].to_numpy(float)) for c in SEED_COLS])),
         float(np.mean([auc_at_cut(y, oof[c].to_numpy(float)) for c in SEED_COLS]))),
        ("RMSE (log10 h)", rmse(pi), rmse(oof)),
        ("Mouse rho",
         float(np.mean([rho(y[mouse], pi[c].to_numpy(float)[mouse]) for c in SEED_COLS])),
         float(np.mean([rho(y[mouse], oof[c].to_numpy(float)[mouse]) for c in SEED_COLS]))),
    ]
    frame = pd.DataFrame(rows, columns=["Metric", "Pooled-route fit", "Intravenous fit"])
    cap = ("Source table, not numbered in the manuscript. What restricting to intravenous records is worth, measured like for "
           f"like: the same {len(oof)} intravenous rows, the same shipped configuration and "
           "the same grouped protocol, scored under the route-pooled fit and the "
           "intravenous fit.")
    return frame, cap


# ---------------------------------------------------------------------------------------
# Table 3. Held-out performance, over both species and by species
# ---------------------------------------------------------------------------------------
def performance_by_stratum() -> tuple[pd.DataFrame, str]:
    t = pd.read_csv(TABS / "table10_muriverse_halflife_branch_metrics.csv")
    frame = t.loc[[0, 2, 1], ["stratum", "n_val", "spearman", "pearson", "auroc",
                              "r2", "rmse", "mae"]].reset_index(drop=True)
    frame.columns = ["Stratum", "n", "Spearman rho", "Pearson r", "AUC (>= 4 h)",
                     "R2", "RMSE (log10 h)", "MAE (log10 h)"]
    cap = ("Source table, not numbered in the manuscript. Held-out performance on the intravenous cohort, over both species and by species. "
           "Out-of-fold predictions under held-out homology clusters and sources, three-seed "
           "means, on log10 half-life. The AUC is at the 4 h screening threshold; the manuscript reports AUC at 6 h (Table 1, Figure 3e).")
    return frame, cap


# ------------------------------------------------------------------------------
# Manuscript Table 1. Original/pre-audit benchmark and the nested estimate
# ------------------------------------------------------------------------------
def benchmark_table() -> tuple[pd.DataFrame, str]:
    c = pd.read_csv(ANALYSIS / "clean_cohort_metrics.csv")
    c = c[(c.scenario == "clean_cohort") & (c.estimator == "ensemble")].set_index("metric")
    n = pd.read_csv(TABS / "table31b_muriverse_halflife_nested_selection_postaudit.csv").iloc[0]

    def cell(v, lo, hi):
        return f"{v:.3f} ({lo:.3f}–{hi:.3f})"

    rows = []
    for label, metric in (("Spearman ρ", "spearman"), ("Pearson r", "pearson"),
                          ("RMSE (log₁₀ h)", "rmse"), ("MAE (log₁₀ h)", "mae"),
                          ("R²", "r2"), ("AUC (t½ ≥ 6 h)", "auc")):
        r = c.loc[metric]
        nested = "not computed"
        if metric in ("spearman", "pearson", "auc"):
            nested = cell(n[metric], n[f"{metric}_lo"], n[f"{metric}_hi"])
        rows.append({"Metric": label,
                     "Post-audit fixed configuration, point estimate (95% CI)": cell(r.value, r.ci_low, r.ci_high),
                     "Nested feature-configuration selection, point estimate (95% CI)": nested})
    return pd.DataFrame(rows), CAPTION_TABLE1


# ------------------------------------------------------------------------------
# Manuscript Table 2. Ranking performance of compact and embedding-based configurations
# ------------------------------------------------------------------------------
TABLE2_ROWS = [
    ("annotation + context + termini + physchem", "Molecule annotation + study context + terminal classes + physicochemical (shipped)", 71),
    ("physchem + study context + termini", "Physicochemical + study context + terminal classes", 62),
    ("physchem + study context", "Physicochemical + study context", 49),
    ("study context only", "Study context alone", 33),
    ("sequence + physchem + study context", "Sequence embedding + physicochemical + study context", 1329),
    ("SMILES + physchem + study context", "SMILES embedding + physicochemical + study context", 817),
    ("all molecular blocks + study context", "Shipped blocks + both embeddings", 2119),
    ("sequence + SMILES", "Both embeddings", 2048),
    ("SMILES (PeptideCLM) only", "SMILES embedding (PeptideCLM) alone", 768),
    ("sequence (ESM-2) only", "Sequence embedding (ESM-2) alone", 1280),
    ("molecule annotation only", "Molecule annotation alone", 9),
]


def ablation_table() -> tuple[pd.DataFrame, str]:
    t = pd.read_csv(TABS / "table31_muriverse_halflife_feature_ablation_strict_iv.csv")
    g = {k: t[t.grouping == k].set_index("features") for k in ("PMID", "cluster", "cluster_PMID")}
    out = []
    for key, label, cols in TABLE2_ROWS:
        srcs, clus, both = (g[k].loc[key] for k in ("PMID", "cluster", "cluster_PMID"))
        if int(srcs.dims) != cols:
            raise ValueError(f"{key}: {int(srcs.dims)} columns, manuscript says {cols}")
        out.append({"Feature configuration": label, "Columns": f"{cols:,}",
                    "Sources withheld": srcs.spearman, "Clusters withheld": clus.spearman,
                    "Clusters + sources withheld": both.spearman,
                    "Retained": f"{100 * clus.spearman / srcs.spearman:.0f}%"})
    return pd.DataFrame(out), CAPTION_TABLE2


# ---------------------------------------------------------------------------------------
# Table 5. Evaluation stringency
# ---------------------------------------------------------------------------------------
def stringency() -> tuple[pd.DataFrame, str]:
    t = pd.read_csv(TABS / "table31c_muriverse_halflife_stringency.csv")
    shown = {
        "sources held out only": "Sources held out only",
        "clusters held out only": "Clusters held out only",
        "strict: clusters + sources held out": "Clusters + sources held out (headline)",
        "strict, configuration chosen inside folds": "Configuration chosen inside folds",
        "within source (source mean removed)": "Within source (source mean removed)",
    }
    t = t.copy()
    t["key"] = t.evaluation.str.lower()
    out = []
    for key, label in shown.items():
        r = t[t.key == key]
        if r.empty:
            raise KeyError(f"stringency row missing from the result table: {key}")
        r = r.iloc[0]
        out.append({"Evaluation": label, "n": int(r.n),
                    "Pearson r": r.pearson, "Spearman rho": r.spearman})
    cap = ("Source table, not numbered in the manuscript. Performance degrades as the holdout gets stricter. The last row removes "
           "each source's mean from both prediction and label, isolating within-source "
           "ordering; the gap between it and the headline is the source-level offset.")
    return pd.DataFrame(out), cap


# ---------------------------------------------------------------------------------------
# Table 6. Polymyxin observed against predicted
# ---------------------------------------------------------------------------------------
def polymyxin_observed_vs_predicted() -> tuple[pd.DataFrame, str]:
    obs = pd.read_csv(TABS / "table49_polymyxin_observed_by_drug_species_iv.csv")
    pred = pd.read_csv(TABS / "table20_polymyxin_predictions.csv")
    base = pd.read_csv(TABS / "table27_polymyxin_fold_error_vs_baseline_iv.csv")
    sh = pred[pred.role == "shipped"]
    # Per drug AND species. Species is a model input, so a per-drug median would average two
    # different predictions and report neither.
    med = sh.groupby(["drug", "species"]).pred_half_life_h.median()
    # Fold error here is the MEDIAN OF THE PER-VARIANT fold errors, which is the quantity
    # the text quotes for this table, and is not the ratio of the two medians used for the
    # per-source table below. A fold error so defined is a magnitude with no direction, so
    # the signed ratio of the medians is carried in its own column rather than collapsed
    # into an "over"/"under" word that the four variants need not agree on.
    fe = base[base.role == "shipped"].groupby(["drug", "species"]).model_fold_error.median()
    out = []
    for r in obs.sort_values(["drug", "species"]).itertuples(index=False):
        key = next(k for k, _ in med.index if k.startswith(r.drug.split()[0]))
        p = float(med[(key, r.species)])
        out.append({
            "Drug": r.drug.split(" (")[0], "Species": r.species.title(),
            "Observed half-life (h)": r.observed_h,
            "Dose groups": int(r.n_groups), "Animals": int(r.n_animals),
            "Studies": int(r.n_studies),
            "Median prediction (h)": p,
            "Fold error": float(fe[(key, r.species)]),
            "Median prediction / observed": p / r.observed_h,
        })
    cap = ("Table 6. Model predictions versus observed intravenous polymyxin half-lives by species "
           "and drug. Prediction, median is the median of the four prediction rows for that drug and "
           "species (two congeners under two equivalent N-terminal encodings, which give identical "
           "model inputs). Fold error is the median of the per-row fold errors, each taken as the larger "
           "of prediction/observed and observed/prediction. Prediction / observed is the median prediction "
           "divided by the observed value. The last two columns are not reciprocals of one another: fold "
           "error is aggregated per row and discards the direction of each error, whereas prediction / "
           "observed is a ratio of medians and retains it. In every cell the per-congener predictions lie "
           "on one side of the intravenous reference, so the ratio of medians and the median fold error "
           "agree in direction.")
    return pd.DataFrame(out), cap


# ---------------------------------------------------------------------------------------
# Table 5. Polymyxin fold error by source of observation
# ---------------------------------------------------------------------------------------
def polymyxin_by_source() -> tuple[pd.DataFrame, str]:
    t = pd.read_csv(TABS / "table53_polymyxin_fold_error_by_source.csv")
    t = t[(t.role == "shipped") & (~t.source_group.str.startswith("pooled"))]
    # Per drug AND species: species is a model input, and each reference is for one animal.
    med_pred = t.groupby(["drug", "species"]).pred_half_life_h.median()
    # Mean residence times are reported beside the half-lives but are not compared with the
    # model: a mean residence time is not a half-life. They are the rows with no value
    # basis, and they are the reason this table has ten rows rather than eleven.
    t = t[t.value_basis.notna() & (t.source_group != "mean residence time")]
    keys = t[["drug", "species", "source_group", "value_basis"]].drop_duplicates()
    out = []
    for r in keys.itertuples(index=False):
        rows = t[(t.drug == r.drug) & (t.species == r.species)
                 & (t.source_group == r.source_group) & (t.value_basis == r.value_basis)]
        obs = float(rows.observed_h.iloc[0])
        p = float(med_pred[(r.drug, r.species)])
        # Ratio of the two medians, which is what the manuscript quotes.
        fold = p / obs if p >= obs else obs / p
        out.append({
            "Drug": r.drug.split(" (")[0], "Species": r.species.title(),
            "Reference": f"{r.source_group}, {r.value_basis}",
            "Observed half-life (h)": obs,
            "Fold error": fold,
            "Direction": "over" if p >= obs else "under",
            "Animals": int(rows.n_animals.iloc[0]),
            "Studies": int(rows.n_studies.iloc[0]),
        })
    frame = pd.DataFrame(out).sort_values(
        ["Drug", "Species", "Observed half-life (h)"], ascending=[True, True, False]
    ).reset_index(drop=True)
    cap = ("Table 5. Polymyxin half-life estimates by study group. Observed half-lives are "
           "animal-count-weighted geometric means of the contributing dose groups, across all "
           "administration routes. Fold error is the ratio of the drug-and-species median "
           "prediction to that group's observed value, taken in whichever direction exceeds "
           "one. Table 6 instead summarizes the median of the four per-congener fold errors, so "
           "the two conventions can differ for the same cell.")
    return frame, cap


# ---------------------------------------------------------------------------------------
# Table 7. Sensitivity of the pooled observed reference
# ---------------------------------------------------------------------------------------
def pooled_sensitivity() -> tuple[pd.DataFrame, str]:
    t = pd.read_csv(TABS / "table55_polymyxin_pooled_sensitivity.csv")
    pred = pd.read_csv(TABS / "table20_polymyxin_predictions.csv")
    sh = pred[pred.role == "shipped"]

    def fold(variants: np.ndarray, obs: float) -> float:
        """Median of the per-variant fold errors, as for Table 6."""
        if pd.isna(obs):
            return float("nan")
        return float(np.median(np.where(variants >= obs, variants / obs, obs / variants)))

    out = []
    for r in t.sort_values(["drug", "species"]).itertuples(index=False):
        stem = r.drug.split()[0]
        v = sh[(sh.drug.str.startswith(stem)) & (sh.species == r.species)]
        p = v.pred_half_life_h.to_numpy(float)
        out.append({
            "Drug": r.drug.split(" (")[0], "Species": r.species.title(),
            "Pooled, animal-weighted (h)": r.pooled_animal_weighted_h,
            "Pooled fold error": fold(p, r.pooled_animal_weighted_h),
            "Each study once (h)": r.each_study_once_h,
            "Each study once fold error": fold(p, r.each_study_once_h),
            "Reported values only (h)": r.reported_only_h,
            "Reported only fold error": fold(p, r.reported_only_h),
            "Intravenous only (h)": r.iv_only_h,
            "Intravenous only fold error": fold(p, r.iv_only_h),
            "Leave-one-out min (h)": r.leave_one_out_min_h,
            "Leave-one-out max (h)": r.leave_one_out_max_h,
        })
    cap = ("Table 7. Robustness of pooled polymyxin half-life estimates to aggregation method. "
           "Each fold error is the median of the four per-congener values for that drug and species, each "
           "taken as the larger of prediction/observed and observed/prediction. Leave-one-out columns give "
           "the smallest and largest pooled estimate when each contributing study is omitted in turn. "
           "IV = intravenous.")
    return pd.DataFrame(out), cap


# ---------------------------------------------------------------------------------------
def main() -> None:
    TABS.mkdir(parents=True, exist_ok=True)
    r3 = {c: 3 for c in ("Spearman rho", "Pearson r",
                         "AUC", "AUC (>= 4 h)", "R2", "RMSE (log10 h)", "MAE (log10 h)",
                         "Pooled-route fit", "Intravenous fit")}
    r2 = {c: 2 for c in ("Observed half-life (h)", "Median prediction (h)", "Fold error",
                         "Pooled, animal-weighted (h)", "Pooled fold error",
                         "Each study once (h)", "Each study once fold error",
                         "Reported values only (h)", "Reported only fold error",
                         "Intravenous only (h)", "Intravenous only fold error",
                         "Median prediction / observed",
                         "Leave-one-out min (h)", "Leave-one-out max (h)")}
    specs = [
        ("table1_benchmark_metrics", benchmark_table, {}),
        ("table2_feature_ablation", ablation_table,
         {c: 3 for c in ("Sources withheld", "Clusters withheld", "Clusters + sources withheld")}),
        ("1_species_route_transfer", species_route_transfer,
         {**{c: 2 for c in ("Own-cell rho", "Pooled rho", "Cross-species rho",
                            "Cross-route rho")},
          **{c: "+2" for c in ("Cross-species bias", "Cross-route bias")}}),
        ("2_route_restriction", route_restriction, r3),
        ("3_performance_by_stratum", performance_by_stratum, r3),
        ("5_evaluation_stringency", stringency, r3),
        ("table5_polymyxin_by_source", polymyxin_by_source, r2),
        ("table6_polymyxin_observed_vs_predicted", polymyxin_observed_vs_predicted, r2),
        ("table7_polymyxin_pooled_sensitivity", pooled_sensitivity, r2),
    ]
    blocks = []
    for name, fn, fmt in specs:
        frame, caption = fn()
        blocks.append(write(name, frame, caption, fmt))

    index = TABS / f"{PREFIX}all.md"
    index.write_text(
        "# Publication tables\n\n"
        "Generated by MuriVerse/build_publication_tables.py from the pipeline's own result "
        "tables. Tables 1, 2 and 5 to 7 (prefixed table) are the manuscript's tables of the same number. The other tables are source "
        "tables for Supplementary Figure 2 and for reference; they are not numbered in the "
        "manuscript. Manuscript Tables 3 and 4 are not generated here. Values are shown at the precision the manuscript quotes; the companion CSV "
        "beside each carries full precision.\n\n" + "\n\n".join(blocks) + "\n")
    print(f"[write] {index.name}  ({len(specs)} tables)")


if __name__ == "__main__":
    main()
