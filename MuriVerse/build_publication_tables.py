"""Render the publication tables from the pipeline's own result tables.

Emits standalone publication tables under figures/tabs/ with the `pub-tab__` prefix, so a table
quoted in the text has one generated source and cannot drift from it.

Every value here is read or re-derived from a committed result table, or from the
out-of-fold predictions. Nothing is
transcribed from prose. Each table is written twice: a CSV carrying the full-precision
values, and a markdown rendering at the precision the manuscript quotes, for pasting.

    python MuriVerse/build_publication_tables.py

Two tables are derived rather than selected, and both follow the audit's definitions:

  * The route-restriction comparison (Table 2) scores the SAME 213 intravenous rows under
    the route-pooled fit and the intravenous fit, which is the only like-for-like reading.
  * The per-source polymyxin references (Table 7) take the fold error as the ratio of the
    two MEDIANS, not the median of the per-variant fold errors. The two differ wherever a
    drug's variants straddle the observed value (Monash mouse: 1.04 against 1.31).
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
    cap = ("Table 1. Cross-route transfer is offset, cross-species transfer is not. The same "
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
    cap = ("Table 2. What restricting to intravenous records is worth, measured like for "
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
    cap = ("Table 3. Held-out performance on the intravenous cohort, over both species and by species. "
           "Out-of-fold predictions under held-out homology clusters and sources, three-seed "
           "means, on log10 half-life.")
    return frame, cap


# ---------------------------------------------------------------------------------------
# Table 4. Feature-block ablation
# ---------------------------------------------------------------------------------------
def feature_ablation() -> tuple[pd.DataFrame, str]:
    t = pd.read_csv(TABS / "table31_muriverse_halflife_feature_ablation_strict_iv.csv")
    g = t[t.grouping == "cluster_PMID"].set_index("features")
    # The configurations the manuscript tabulates, in the order it tabulates them, with the
    # plain-language block names used in the text.
    rows = [
        ("physchem + study context + termini", "physicochemical + study context + terminal classes (top-scoring, not shipped)"),
        ("annotation + context + termini + physchem", "annotation + study context + terminal classes + physicochemical (shipped)"),
        ("physchem + study context", "physicochemical + study context"),
        ("study context + termini", "study context + terminal classes"),
        ("molecule annotation + study context", "molecule annotation + study context"),
        ("study context only", "study context only"),
        ("SMILES + physchem + study context", "SMILES embedding + physicochemical + study context"),
        ("sequence + SMILES", "sequence + SMILES embeddings"),
        ("SMILES (PeptideCLM) only", "SMILES embedding only"),
        ("sequence (ESM-2) only", "sequence embedding only"),
        ("molecule annotation only", "molecule annotation only"),
        ("baseline (train mean)", "baseline (train mean)"),
    ]
    out = []
    for key, shown in rows:
        if key not in g.index:
            raise KeyError(f"ablation configuration missing from the result table: {key}")
        r = g.loc[key]
        out.append({"Feature block": shown, "Dims": int(r.dims),
                    "Spearman rho": r.spearman, "Pearson r": r.pearson, "AUC": r.auc})
    cap = (f"Table 4. Feature-block ablation under the headline grouping, held-out homology "
           f"clusters and sources. {int((t.grouping == 'cluster_PMID').sum())} configurations "
           "were scored; the twelve the text discusses are shown. Most usable signal is study "
           "context rather than molecular structure, and both pretrained embeddings score "
           "below the compact set.")
    return pd.DataFrame(out), cap


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
    cap = ("Table 5. Performance degrades as the holdout gets stricter. The last row removes "
           "each source's mean from both prediction and label, isolating within-source "
           "ordering; the gap between it and the headline is the source-level offset.")
    return pd.DataFrame(out), cap


# ---------------------------------------------------------------------------------------
# Table 6. Polymyxin observed against predicted
# ---------------------------------------------------------------------------------------
def polymyxin_observed_vs_predicted() -> tuple[pd.DataFrame, str]:
    obs = pd.read_csv(TABS / "table49_polymyxin_observed_by_drug_species.csv")
    pred = pd.read_csv(TABS / "table20_polymyxin_predictions.csv")
    base = pd.read_csv(TABS / "table27_polymyxin_fold_error_vs_baseline.csv")
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
    cap = ("Table 6. Held-out application to colistin and polymyxin B. Observed values are "
           "animal-weighted geometric means across dose groups; predictions are the median "
           "shipped-model value per drug and species, taken over the two congeners of each "
           "drug. Fold error is the median of the per-congener fold errors, so it is a "
           "magnitude; the final column gives the signed ratio. Species is a model input, "
           "so the mouse and rat predictions differ; the difference the model applies is "
           "small, and on these molecules it runs opposite to the observed one. The "
           "C-terminus is annotated only as unmodified, the macrolactam being carried by "
           "the topology block, so each cell is one prediction and not a range over "
           "terminal-class annotations.")
    return pd.DataFrame(out), cap


# ---------------------------------------------------------------------------------------
# Table 7. Polymyxin fold error by source of observation
# ---------------------------------------------------------------------------------------
def polymyxin_by_source() -> tuple[pd.DataFrame, str]:
    t = pd.read_csv(TABS / "table53_polymyxin_fold_error_by_source.csv")
    t = t[(t.role == "shipped") & (~t.source_group.str.startswith("pooled"))]
    med_pred = t.groupby("drug").pred_half_life_h.median()
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
        p = float(med_pred[r.drug])
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
    cap = ("Table 7. The verdict depends on which published half-life is the reference. Fold "
           "error is the ratio of the median shipped prediction for that drug to the "
           "source's observed value, so the same prediction is short against "
           "the longest published values and long against the shortest.")
    return frame, cap


# ---------------------------------------------------------------------------------------
# Table 8. Sensitivity of the pooled observed reference
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
    cap = ("Table 8. The pooled observed reference is a property of the weighting. Leaving "
           "any single study out moves the pooled polymyxin B mouse value only between 0.99 "
           "and 1.33 h, but restricting to values the studies themselves report more than "
           "doubles its fold error.")
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
        ("1_species_route_transfer", species_route_transfer,
         {**{c: 2 for c in ("Own-cell rho", "Pooled rho", "Cross-species rho",
                            "Cross-route rho")},
          **{c: "+2" for c in ("Cross-species bias", "Cross-route bias")}}),
        ("2_route_restriction", route_restriction, r3),
        ("3_performance_by_stratum", performance_by_stratum, r3),
        ("4_feature_ablation", feature_ablation, r3),
        ("5_evaluation_stringency", stringency, r3),
        ("6_polymyxin_observed_vs_predicted", polymyxin_observed_vs_predicted, r2),
        ("7_polymyxin_by_source", polymyxin_by_source, r2),
        ("8_polymyxin_pooled_sensitivity", pooled_sensitivity, r2),
    ]
    blocks = []
    for name, fn, fmt in specs:
        frame, caption = fn()
        blocks.append(write(name, frame, caption, fmt))

    index = TABS / f"{PREFIX}all.md"
    index.write_text(
        "# Publication tables\n\n"
        "Generated by MuriVerse/build_publication_tables.py from the pipeline's own result "
        "tables. Values are shown at the precision the manuscript quotes; the companion CSV "
        "beside each carries full precision.\n\n" + "\n\n".join(blocks) + "\n")
    print(f"[write] {index.name}  ({len(specs)} tables)")


if __name__ == "__main__":
    main()
