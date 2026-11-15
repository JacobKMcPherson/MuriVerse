"""Build the six-metric IV table with group-bootstrap intervals from saved OOF predictions.

Both point estimates and bootstrap draws use the mean of the three seed-specific scores.
Predictions are not averaged before scoring: that would estimate a different statistic and
produce intervals inconsistent with the reported point estimates. This post-hoc script refits
nothing and uses the same group-resampling stream as the evaluator.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "MuriVerse"))

from evaluate_halflife_model import (  # noqa: E402
    N_BOOT,
    _auc,
    _pearson,
    _rho,
    _score,
    union_groups,
)

SNAPSHOT = REPO / "results/revision_analysis/snapshot/oof_predictions_iv.csv"
OUT = (Path(__file__).resolve().parents[2] / "results" / "revision_analysis" / "metrics_iv.csv")

# The published values this script must reproduce before its new columns can be trusted.
EXPECTED = {
    "spearman": (0.574, 0.322, 0.655),
    "pearson": (0.568, 0.381, 0.651),
    "auc": (0.818, 0.620, 0.877),
}
SEED_COLS = [f"y_pred_cluster_PMID_seed{i}" for i in range(3)]


def bootstrap_all(y, predictions, groups):
    """Group-bootstrap the mean of seed-specific scores for all six metrics."""
    rng = np.random.default_rng(0)
    # Groups are drawn in np.unique (sorted-label) order; the reported intervals depend on it.
    members = {g: np.flatnonzero(groups == g) for g in np.unique(groups)}
    keys = list(members)
    stats = {k: [] for k in ("spearman", "pearson", "auc", "rmse", "mae", "r2")}
    for _ in range(N_BOOT):
        pick = rng.choice(len(keys), size=len(keys), replace=True)
        idx = np.concatenate([members[keys[i]] for i in pick])
        if np.ptp(y[idx]) > 0:
            scorers = (
                ("spearman", _rho),
                ("pearson", _pearson),
                ("rmse", lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))),
                ("mae", lambda a, b: float(np.mean(np.abs(a - b)))),
                ("r2", lambda a, b: float(1 - ((a - b) ** 2).sum() / ((a - a.mean()) ** 2).sum())),
                ("auc", _auc),
            )
            for name, scorer in scorers:
                per_seed = [scorer(y[idx], pred[idx]) for pred in predictions]
                valid = [value for value in per_seed if not np.isnan(value)]
                if valid:
                    stats[name].append(float(np.mean(valid)))
    return {
        k: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
        for k, v in stats.items()
    }


def main():
    d = pd.read_csv(SNAPSHOT)
    print(f"input: {SNAPSHOT.name}, {len(d)} rows, {d.shape[1]} columns")

    y = d["y_true"].to_numpy(float)
    predictions = d[SEED_COLS].to_numpy(float).T
    groups = union_groups(d["cluster"].to_numpy(), d["PMID"].to_numpy())
    print(f"grouping: {len(np.unique(groups))} union groups over cluster and PMID")
    print(f"point estimate: mean of {len(SEED_COLS)} per-seed scores")
    print(f"interval: group bootstrap of the mean per-seed score; N_BOOT = {N_BOOT}")

    per_seed = [_score(y, d[c].to_numpy(float)) for c in SEED_COLS]
    point = {k: float(np.mean([s[k] for s in per_seed])) for k in per_seed[0]}
    point["n"] = per_seed[0]["n"]
    ci = bootstrap_all(y, predictions, groups)

    rows = []
    for key, label in [
        ("spearman", "Spearman rho"),
        ("pearson", "Pearson r"),
        ("rmse", "RMSE (log10 h)"),
        ("mae", "MAE (log10 h)"),
        ("r2", "R-squared"),
        ("auc", "AUC (t-half >= 6 h)"),
    ]:
        rows.append(
            {
                "metric": label,
                "value": round(point[key], 3),
                "ci_low": round(ci[key][0], 3),
                "ci_high": round(ci[key][1], 3),
                "n": point["n"],
                "n_groups": len(np.unique(groups)),
            }
        )
    out = pd.DataFrame(rows)

    print()
    print(out.to_string(index=False))

    print()
    ok = True
    for key, (v, lo, hi) in EXPECTED.items():
        got = (round(point[key], 3), round(ci[key][0], 3), round(ci[key][1], 3))
        match = got == (v, lo, hi)
        ok &= match
        print(f"  check {key:9s} expected {v} [{lo}, {hi}]  got {got[0]} "
              f"[{got[1]}, {got[2]}]  {'OK' if match else 'MISMATCH'}")
    if not ok:
        print("\nSTOP: the resampling stream differs from the published one. "
              "Do not use the new columns.")
        return 1

    out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
