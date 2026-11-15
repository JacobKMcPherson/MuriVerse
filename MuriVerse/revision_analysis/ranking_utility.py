"""Tertile recovery of the shipped model's out-of-fold ranking (reviewer 2, major 1)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "MuriVerse"))
from evaluate_halflife_model import union_groups  # noqa: E402

A = Path(__file__).resolve().parents[2] / "results" / "revision_analysis"
SNAPSHOT = REPO / "results/revision_analysis/snapshot/oof_predictions_iv.csv"
OUT = A / "ranking_utility.csv"
N_BOOT = 500


def tertile_scores(y: np.ndarray, seed_predictions: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Mean over seeds of (top recovered, bottom recovered, predicted-top that is truly bottom)."""
    rows = []
    for pred in seed_predictions:
        yy, pp = y[idx], pred[idx]
        k = len(yy) // 3
        true_top, pred_top = set(np.argsort(-yy)[:k]), set(np.argsort(-pp)[:k])
        true_bottom, pred_bottom = set(np.argsort(yy)[:k]), set(np.argsort(pp)[:k])
        rows.append((len(true_top & pred_top) / k, len(true_bottom & pred_bottom) / k,
                     len(pred_top & true_bottom) / k))
    return np.mean(rows, axis=0)


def main() -> None:
    d = pd.read_csv(SNAPSHOT)
    y = d["y_true"].to_numpy(float)
    seeds = np.vstack([d[f"y_pred_cluster_PMID_seed{i}"].to_numpy(float) for i in range(3)])
    groups = union_groups(d["cluster"].to_numpy(), d["PMID"].astype(str).to_numpy())
    members = {g: np.flatnonzero(groups == g) for g in np.unique(groups)}
    keys = list(members)
    rng = np.random.default_rng(0)
    draws = np.array([
        tertile_scores(y, seeds, np.concatenate([members[keys[i]] for i in
                                                  rng.choice(len(keys), len(keys))]))
        for _ in range(N_BOOT)])
    low, high = np.percentile(draws, [2.5, 97.5], axis=0)
    point = tertile_scores(y, seeds, np.arange(len(y)))
    names = ["true_top_tertile_in_predicted_top", "true_bottom_tertile_in_predicted_bottom",
             "predicted_top_that_is_truly_bottom"]
    out = pd.DataFrame({"measure": names, "estimate": point, "ci_low": low, "ci_high": high,
                        "chance": 1 / 3, "n": len(y), "groups": len(keys)})
    out.to_csv(OUT, index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
