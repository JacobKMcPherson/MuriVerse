"""Route-stratified performance for the shipped half-life model.

The headline numbers for this model are pooled over route, and that pooling hides two
things that matter to anyone applying it.

First, route is not a nuisance variable here, it is part of the measurement. After an
intravenous dose the terminal half-life is an elimination half-life; after an extravascular
dose it is one only if absorption is faster than elimination, and when it is not the
terminal slope reports absorption instead (flip-flop kinetics). `evaluate_route_flipflop.py`
shows the consequence directly: IV and SC labels do not share an absolute scale, and a model
trained on one and applied to the other is offset by roughly half a log unit in the
direction flip-flop predicts.

Second, the pooled score is carried by the routes with the most records. IV supplies 213 of
the 320 rows; IM supplies 12. A single rho for the cohort says nothing about whether the
model is usable for a route the reader cares about.

This script reports per-route performance from the SAVED out-of-fold predictions written by
`evaluate_halflife_model.py`, so it costs nothing to run and the numbers are by construction
the same predictions the headline scores come from. Nothing is refitted here.

Bias - mean(predicted - observed) in log10 hours - is reported alongside the ranking
metrics because the ranking metrics cannot see a scale error: rho and AUC are invariant to
a monotone shift and RMSE/MAE are symmetric, so a route that is consistently mispredicted in
one direction looks acceptable on every other column.

Written to figures/tabs/table34_muriverse_halflife_performance_by_route.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
TABS = ROOT.parent / "figures/tabs"
OUT = TABS / "table34_muriverse_halflife_performance_by_route.csv"

# The headline grouping: a held-out peptide shares neither a homology cluster nor a source
# identifier with the training rows.
GROUPING = "cluster_PMID"
SEEDS = (0, 1, 2)
HALF_LIFE_CUTOFF_H = 6.0


def _auc(y: np.ndarray, pred: np.ndarray) -> float:
    """Rank AUC at the 6 h screening cut; NaN when a stratum is all-positive or all-negative."""
    positive = y >= np.log10(HALF_LIFE_CUTOFF_H) - 1e-9
    n1, n0 = int(positive.sum()), int((~positive).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = pd.Series(pred).rank().to_numpy()
    return float((r[positive].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def score(y: np.ndarray, pred: np.ndarray) -> dict:
    fold_error = 10.0 ** np.abs(y - pred)
    return {
        "n": int(len(y)),
        "n_positive_6h": int((y >= np.log10(HALF_LIFE_CUTOFF_H) - 1e-9).sum()),
        "spearman": float(spearmanr(y, pred).correlation),
        "auc_6h": _auc(y, pred),
        "bias": float(np.mean(pred - y)),
        "mdafe": float(np.median(fold_error)),
        "within_2fold": float(np.mean(fold_error <= 2.0)),
        "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
    }


def main() -> None:
    meta_dir = ROOT / "rodent" / "strict"
    oof = pd.read_csv(meta_dir / "oof_predictions.csv")

    # `route_class` is written into the out-of-fold table by the current evaluator, but an
    # older file will not carry it; join it from the metadata rather than failing, so this
    # report still runs against predictions saved before that change.
    if "route_class" not in oof.columns:
        meta = pd.read_csv(meta_dir / "halflife_meta_with_split.csv", low_memory=False)
        oof = oof.merge(meta[["id", "route_class"]], on="id", how="left", validate="1:1")
        if oof["route_class"].isna().any():
            raise SystemExit("route_class missing for some ids; rebuild oof_predictions.csv")

    seed_cols = [f"y_pred_{GROUPING}_seed{k}" for k in SEEDS]
    missing = [c for c in seed_cols if c not in oof.columns]
    if missing:
        raise SystemExit(f"out-of-fold table lacks {missing}; rerun evaluate_halflife_model.py")

    y = oof["y_true"].to_numpy(float)
    route = oof["route_class"].astype(str).to_numpy()

    # One row per (stratum, seed), then averaged, so the reported spread is across complete
    # out-of-fold passes rather than across records.
    strata = [("all routes", np.ones(len(y), dtype=bool))]
    strata += [(r, route == r) for r in sorted(set(route))]

    rows = []
    for name, mask in strata:
        for k, col in zip(SEEDS, seed_cols):
            rows.append({"stratum": name, "seed": k,
                         **score(y[mask], oof[col].to_numpy(float)[mask])})
    per_seed = pd.DataFrame(rows)

    metrics = ["spearman", "auc_6h", "bias", "mdafe", "within_2fold", "rmse"]
    summary = (per_seed.groupby("stratum", as_index=False)
               .agg(n=("n", "first"), n_positive_6h=("n_positive_6h", "first"),
                    **{c: (c, "mean") for c in metrics},
                    **{f"{c}_sd": (c, "std") for c in metrics}))
    # Report in descending size: the pooled row first, then the routes that actually carry it.
    summary = pd.concat([summary[summary.stratum == "all routes"],
                         summary[summary.stratum != "all routes"].sort_values("n",
                                                                              ascending=False)])

    TABS.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT, index=False)

    print(f"Shipped model, held-out {GROUPING}, mean of {len(SEEDS)} out-of-fold passes\n")
    show = summary.copy()
    show["auc_6h"] = show["auc_6h"].map(lambda v: "n/a" if pd.isna(v) else f"{v:.3f}")
    print(show[["stratum", "n", "n_positive_6h", "spearman", "auc_6h", "bias",
                "mdafe", "within_2fold"]].to_string(index=False,
                                                    float_format=lambda v: f"{v:+.3f}"))
    print(f"\n[Wrote] {OUT}")

    # The warning the pooled row cannot give on its own.
    weak = summary[(summary.stratum != "all routes") & (summary.spearman < 0.2)]
    if len(weak):
        names = ", ".join(f"{r.stratum} (n={int(r.n)}, rho={r.spearman:+.2f})"
                          for r in weak.itertuples())
        print(f"\n[Warning] at or below chance ranking: {names}\n"
              f"          The pooled score does not transfer to these routes.")


if __name__ == "__main__":
    main()
