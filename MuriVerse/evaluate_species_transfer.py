"""Does the pooled model's species imbalance justify separate mouse and rat models?

Figure 3 reports the pooled model scoring rho = 0.64 in rat but 0.38 in mouse, and an AUC of
0.86 in rat against 0.52 - chance - in mouse. Rat supplies 246 of the 320 rows. Two readings
fit those numbers equally well and they imply opposite actions:

  1. The label means something different in the two animals, and the model is being dragged
     by the majority species. A mouse-only model would then do better on mouse.
  2. Mouse simply has 74 rows, and no model fitted on 74 rows of this cohort would rank
     better. Specialising would then throw away the rat data for nothing.

This script separates them the same way `evaluate_route_flipflop.py` separates route effects.
For each target species, three models are scored on the SAME held-out rows under ONE set of
grouped folds:

    same species  - trained on that species' rows outside the held-out groups
    other species - trained on the other species' rows outside the held-out groups
    pooled        - trained on every cohort row outside the held-out groups

Cross-species transfer is the discriminating test. If the two species' labels are measuring
the same underlying property, a rat-trained model should rank mouse peptides about as well as
a mouse-trained one. Signed bias says whether the two sit on a common absolute scale.

A sample-size control is included for the same reason it was needed for route: mouse's deficit
could be sample size rather than biology, and the two are not distinguishable from the
within-species score alone. Rat is subsampled to mouse's row count and rescored.

Groups are connected components of (homology cluster, source PMID) over the full cohort, so no
cluster or source leaks from a training pool into the evaluated rows.

Writes figures/tabs/table35_muriverse_halflife_species_transfer.csv (+ _by_seed, _size_control).
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
OUT = TABS / "table35_muriverse_halflife_species_transfer.csv"
OUT_SEEDS = TABS / "table35b_muriverse_halflife_species_transfer_by_seed.csv"
OUT_CTL = TABS / "table35c_muriverse_halflife_species_size_control.csv"

SHIPPED_BLOCKS = ["molecule_annotation", "study_context", "termini", "physchem"]
# The shipped block set carries `species=mouse` and `species=rat` in study_context. Within a
# single-species training set those two columns are constant, so they carry no variance for
# the model to fit: ridge leaves a zero-variance column at zero weight and no tree can split
# on it. They are kept rather than dropped so the transfer model is literally the shipped
# feature set and not a near-copy of it, exactly as the route columns are kept under the
# route restriction in evaluate_halflife_model.py.
SPECIES = ("mouse", "rat")
SEEDS = (0, 1, 2, 3, 4)
N_SIZE_DRAWS = 10
HALF_LIFE_CUTOFF_H = 6.0


def _auc(y: np.ndarray, pred: np.ndarray) -> float:
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
        "spearman": float(pd.Series(y).corr(pd.Series(pred), method="spearman")),
        "auc_6h": _auc(y, pred),
        "mdafe": float(np.median(fold_error)),
        "within_2fold": float(np.mean(fold_error <= 2.0)),
        "bias": float(np.mean(pred - y)),
    }


def transfer(X, y, groups, species, target, seed):
    other = [s for s in SPECIES if s != target][0]
    pools = {"same species": species == target,
             "other species": species == other,
             "pooled": np.ones(len(y), dtype=bool)}
    eval_idx = np.flatnonzero(species == target)
    eval_groups = groups[eval_idx]
    out = {k: np.full(len(eval_idx), np.nan) for k in pools}

    splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
    for _, test_pos in splitter.split(eval_idx, y[eval_idx], eval_groups):
        held = set(eval_groups[test_pos])
        test_rows = eval_idx[test_pos]
        outside = ~np.isin(groups, list(held))
        for name, pool in pools.items():
            train_rows = np.flatnonzero(outside & pool)
            out[name][test_pos] = fit_predict(X[train_rows], y[train_rows], X[test_rows], seed)

    if any(np.isnan(v).any() for v in out.values()):
        raise SystemExit(f"some {target} rows never held out; check grouping")
    return eval_idx, out


def size_control(X, y, groups, species, n_target):
    """Rat cut down to mouse's row count, so within-species scores compare at equal n."""
    rat_groups = np.array(sorted(set(groups[species == "rat"])))
    rows = []
    for draw in range(N_SIZE_DRAWS):
        rng = np.random.default_rng(2000 + draw)
        picked, n_rows = [], 0
        for g in rng.permutation(rat_groups):
            if n_rows >= n_target:
                break
            picked.append(g)
            n_rows += int(((groups == g) & (species == "rat")).sum())
        mask = (species == "rat") & np.isin(groups, picked)
        idx = np.flatnonzero(mask)
        if len(set(groups[idx])) < N_SPLITS:
            continue
        pred = np.full(len(idx), np.nan)
        for _, te in GroupKFold(n_splits=N_SPLITS, shuffle=True,
                                random_state=draw).split(idx, y[idx], groups[idx]):
            held = set(groups[idx][te])
            tr = np.flatnonzero(mask & ~np.isin(groups, list(held)))
            pred[te] = fit_predict(X[tr], y[tr], X[idx[te]], draw)
        rows.append({"draw": draw, "n_groups": len(picked), **score(y[idx], pred)})
    return rows


def main() -> None:
    root = Path(__file__).resolve().parent
    X_all, meta, layout = load_matrix("rodent_strict", root / "rodent" / "strict")
    sl = block_slices(layout)
    X = np.hstack([X_all[:, sl[b]] for b in SHIPPED_BLOCKS])
    y = meta["log_label"].to_numpy(np.float64)
    species = meta["species"].astype(str).str.lower().to_numpy()
    groups = union_groups(meta["cluster"].to_numpy(), meta["PMID"].to_numpy())

    rows = []
    for target in SPECIES:
        n_rows = int((species == target).sum())
        print(f"\n[Transfer] target {target.upper()}: {n_rows} rows, "
              f"{len(set(groups[species == target]))} groups")
        for seed in SEEDS:
            idx, preds = transfer(X, y, groups, species, target, seed)
            for src, pred in preds.items():
                rows.append({"target": target, "train": src, "seed": seed,
                             **score(y[idx], pred)})
        avg = pd.DataFrame(rows).query("target == @target").groupby("train")[
            ["spearman", "auc_6h", "mdafe", "bias"]].mean()
        for src, r in avg.iterrows():
            print(f"    train={src:<14s} rho={r.spearman:+.3f}  AUC={r.auc_6h:.3f}  "
                  f"MdAFE={r.mdafe:.2f}x  bias={r.bias:+.3f}")

    per_seed = pd.DataFrame(rows)
    metrics = ["spearman", "auc_6h", "mdafe", "within_2fold", "bias"]
    summary = (per_seed.groupby(["target", "train"], as_index=False)
               .agg(n=("n", "first"), n_seeds=("seed", "nunique"),
                    **{c: (c, "mean") for c in metrics},
                    **{f"{c}_sd": (c, "std") for c in metrics}))

    n_mouse = int((species == "mouse").sum())
    ctl = pd.DataFrame(size_control(X, y, groups, species, n_mouse))
    mouse_own = summary.query("target == 'mouse' and train == 'same species'").iloc[0]
    print(f"\n[Size control] rat subsampled to mouse's {n_mouse} rows, {len(ctl)} draws")
    print(f"  rat at mouse size : rho {ctl.spearman.mean():+.3f} +/- {ctl.spearman.std():.3f} "
          f"(median n={int(ctl.n.median())})")
    print(f"  mouse own species : rho {mouse_own.spearman:+.3f} +/- {mouse_own.spearman_sd:.3f} "
          f"(n={int(mouse_own.n)})")
    print("  => " + ("rat still ranks better at equal n: mouse deficit is NOT just sample size"
                     if ctl.spearman.mean() - ctl.spearman.std() > mouse_own.spearman
                     else "mouse's deficit is NOT separable from its smaller sample size"))

    TABS.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT, index=False)
    per_seed.to_csv(OUT_SEEDS, index=False)
    ctl.to_csv(OUT_CTL, index=False)
    for p in (OUT, OUT_SEEDS, OUT_CTL):
        print(f"[Wrote] {p}")


if __name__ == "__main__":
    main()
