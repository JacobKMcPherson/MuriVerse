"""Spearman and paired MAE differences (ensemble minus comparator) from the saved identical-fold predictions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ANALYSIS = Path(__file__).resolve().parents[2] / "results" / "revision_analysis"
PREDICTIONS = ANALYSIS / "model_member_baselines_oof.csv"
OUT = ANALYSIS / "model_member_paired.csv"
SEEDS = (0, 1, 2)
N_BOOT = 500
COMPARATORS = ("ridge", "svr", "xgboost", "species_median")


def seed_matrix(frame: pd.DataFrame, name: str) -> np.ndarray:
    return np.vstack([frame[f"{name}_seed{seed}"].to_numpy(float) for seed in SEEDS])


def mean_mae(y: np.ndarray, pred: np.ndarray, idx: np.ndarray) -> float:
    return float(np.mean([np.mean(np.abs(y[idx] - row[idx])) for row in pred]))


def mean_rho(y: np.ndarray, pred: np.ndarray) -> float:
    return float(np.mean([spearmanr(y, row).statistic for row in pred]))


def main() -> None:
    frame = pd.read_csv(PREDICTIONS)
    y = frame["y_true"].to_numpy(float)
    groups = frame["group"].astype(str).to_numpy()
    unique = np.unique(groups)
    members = {g: np.flatnonzero(groups == g) for g in unique}
    preds = {n: seed_matrix(frame, n) for n in ("ensemble", *COMPARATORS)}

    rng = np.random.default_rng(0)
    draws = []
    for _ in range(N_BOOT):
        chosen = rng.choice(unique, size=len(unique), replace=True)
        draws.append(np.concatenate([members[g] for g in chosen]))

    everything = np.arange(len(y))
    rows = []
    for name in ("ensemble", *COMPARATORS):
        row = {"estimator": name, "spearman": mean_rho(y, preds[name]),
               "mae_log10": mean_mae(y, preds[name], everything)}
        if name != "ensemble":
            delta = [mean_mae(y, preds["ensemble"], idx) - mean_mae(y, preds[name], idx)
                     for idx in draws]
            row.update({
                "mae_delta_vs_ensemble": mean_mae(y, preds["ensemble"], everything)
                - mean_mae(y, preds[name], everything),
                "delta_ci_low": float(np.percentile(delta, 2.5)),
                "delta_ci_high": float(np.percentile(delta, 97.5)),
                "ensemble_better_share": float(np.mean(np.asarray(delta) < 0)),
            })
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
