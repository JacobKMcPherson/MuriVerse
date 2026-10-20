"""Cross-species and cross-route transfer, crossed: mouse/rat against IV/extravascular.

Two separate transfer experiments already exist. `evaluate_route_flipflop.py` compares IV
against SC pooled over both animals, and `evaluate_species_transfer.py` compares mouse against
rat pooled over every route. Neither crosses the two, so neither can answer the question an
applier actually has - "is this model usable for MY animal by MY route" - and neither can tell
a species effect from a route effect when the two are confounded, which here they are: rat is
71% intravenous while mouse is only 46%.

This crosses them into four cells and asks, for each, which kind of borrowed data helps.

Routes are grouped as IV against EV (extravascular: subcutaneous, intraperitoneal,
intramuscular). The grouping is mechanistic, not convenience. Only an intravenous dose has no
absorption phase, so only IV guarantees the terminal slope is elimination; every extravascular
route can be absorption-rate-limited (flip-flop kinetics). It also makes the design possible,
since IP (22 rows) and IM (12) are far too small to stand alone. The marginal medians support
the pooling: EV sits above IV by a near-identical factor in both animals independently
(mouse 0.49 -> 1.80 h, rat 0.21 -> 0.83 h).

For each target cell, four models are scored on the SAME held-out rows under ONE set of
grouped folds. The training pools form a 2x2 of what is shared with the target:

    same cell      same species, same route   - the specialist
    cross-species  other species, SAME route  - does the animal transfer?
    cross-route    same species, OTHER route  - does the route transfer?
    pooled         everything                 - the shipped model's training set

Reading `cross-species` against `cross-route` for one target separates the two axes directly:
whichever scores higher is the dimension along which this cohort's label is more portable.

Groups are connected components of (homology cluster, source PMID) over the full cohort, so no
cluster or source leaks from any training pool into the evaluated rows. Cells hold 34-179 rows,
so five seeds are used and the spread across them is reported alongside every mean - several of
these differences are smaller than their own seed noise and must not be read as real.

Writes figures/tabs/table36_muriverse_halflife_species_route_cells.csv (+ _by_seed).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))

from evaluate_halflife_model import (  # noqa: E402
    N_SPLITS, block_slices, fit_predict, load_matrix, union_groups,
)

TABS = Path(__file__).resolve().parent.parent / "figures/tabs"
OUT = TABS / "table36_muriverse_halflife_species_route_cells.csv"
OUT_SEEDS = TABS / "table36b_muriverse_halflife_species_route_cells_by_seed.csv"

SHIPPED_BLOCKS = ["molecule_annotation", "study_context", "termini", "physchem"]
SEEDS = (0, 1, 2, 3, 4)
CUT_H = 6.0
TRAIN_SOURCES = ("same cell", "cross-species", "cross-route", "pooled")


def _auc(y: np.ndarray, pred: np.ndarray) -> float:
    pos = y >= np.log10(CUT_H) - 1e-9
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = pd.Series(pred).rank().to_numpy()
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def score(y: np.ndarray, pred: np.ndarray) -> dict:
    fold_error = 10.0 ** np.abs(y - pred)
    return {
        "n": int(len(y)),
        "spearman": float(pd.Series(y).corr(pd.Series(pred), method="spearman")),
        "auc_6h": _auc(y, pred),
        "mdafe": float(np.median(fold_error)),
        "within_2fold": float(np.mean(fold_error <= 2.0)),
        "bias": float(np.mean(pred - y)),
    }


def transfer(X, y, groups, species, admin, tgt_sp, tgt_ad, seed):
    """Out-of-fold predictions for one cell, from each of the four training pools."""
    other_sp = "rat" if tgt_sp == "mouse" else "mouse"
    other_ad = "EV" if tgt_ad == "IV" else "IV"
    cell = (species == tgt_sp) & (admin == tgt_ad)
    pools = {
        "same cell": cell,
        "cross-species": (species == other_sp) & (admin == tgt_ad),
        "cross-route": (species == tgt_sp) & (admin == other_ad),
        "pooled": np.ones(len(y), dtype=bool),
    }

    eval_idx = np.flatnonzero(cell)
    eval_groups = groups[eval_idx]
    out = {k: np.full(len(eval_idx), np.nan) for k in pools}

    splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
    for _, test_pos in splitter.split(eval_idx, y[eval_idx], eval_groups):
        held = set(eval_groups[test_pos])
        test_rows = eval_idx[test_pos]
        outside = ~np.isin(groups, list(held))
        for name, pool in pools.items():
            train_rows = np.flatnonzero(outside & pool)
            # A held-out group can take every row a small pool owns. Leave the fold NaN
            # rather than fitting on nothing; it is dropped from that pool's score below.
            if len(train_rows) < 5:
                continue
            out[name][test_pos] = fit_predict(X[train_rows], y[train_rows], X[test_rows], seed)
    return eval_idx, out


def main() -> None:
    root = Path(__file__).resolve().parent
    X_all, meta, layout = load_matrix("rodent_strict", root / "rodent" / "strict")
    sl = block_slices(layout)
    X = np.hstack([X_all[:, sl[b]] for b in SHIPPED_BLOCKS])
    y = meta["log_label"].to_numpy(np.float64)
    species = meta["species"].astype(str).str.lower().to_numpy()
    admin = np.where(meta["route_class"].astype(str).str.lower() == "iv", "IV", "EV")
    groups = union_groups(meta["cluster"].to_numpy(), meta["PMID"].to_numpy())

    print("Cell sizes (rows / groups):")
    for sp in ("mouse", "rat"):
        for ad in ("IV", "EV"):
            k = (species == sp) & (admin == ad)
            print(f"  {sp:5} {ad}: {k.sum():3d} rows, {len(set(groups[k])):2d} groups, "
                  f"median {10 ** np.median(y[k]):.2f} h")

    rows = []
    for sp in ("mouse", "rat"):
        for ad in ("IV", "EV"):
            print(f"\n[{sp.upper()} {ad}]")
            for seed in SEEDS:
                idx, preds = transfer(X, y, groups, species, admin, sp, ad, seed)
                for src, pred in preds.items():
                    keep = ~np.isnan(pred)
                    if keep.sum() < 10:
                        continue
                    rows.append({"species": sp, "admin": ad, "train": src, "seed": seed,
                                 "coverage": float(keep.mean()),
                                 **score(y[idx][keep], pred[keep])})
            cur = pd.DataFrame(rows).query("species == @sp and admin == @ad")
            for src in TRAIN_SOURCES:
                s = cur[cur.train == src]
                if len(s):
                    print(f"    {src:<14s} rho={s.spearman.mean():+.3f}+/-{s.spearman.std():.3f}  "
                          f"bias={s.bias.mean():+.3f}  MdAFE={s.mdafe.mean():.2f}x")

    per_seed = pd.DataFrame(rows)
    metrics = ["spearman", "auc_6h", "mdafe", "within_2fold", "bias"]
    summary = (per_seed.groupby(["species", "admin", "train"], as_index=False)
               .agg(n=("n", "max"), n_seeds=("seed", "nunique"),
                    **{c: (c, "mean") for c in metrics},
                    **{f"{c}_sd": (c, "std") for c in metrics}))

    TABS.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT, index=False)
    per_seed.to_csv(OUT_SEEDS, index=False)

    print("\n=== Which axis transfers better? (cross-species vs cross-route, by cell) ===")
    for sp in ("mouse", "rat"):
        for ad in ("IV", "EV"):
            s = summary.query("species == @sp and admin == @ad").set_index("train")
            if {"cross-species", "cross-route"} <= set(s.index):
                cs, cr = s.loc["cross-species"], s.loc["cross-route"]
                gap = cs.spearman - cr.spearman
                noise = max(cs.spearman_sd, cr.spearman_sd)
                verdict = ("species transfers better" if gap > noise else
                           "route transfers better" if -gap > noise else
                           "indistinguishable at seed noise")
                print(f"  {sp:5} {ad}: cross-species {cs.spearman:+.3f} vs "
                      f"cross-route {cr.spearman:+.3f} (seed SD {noise:.3f}) -> {verdict}")
    for p in (OUT, OUT_SEEDS):
        print(f"[Wrote] {p}")


if __name__ == "__main__":
    main()
