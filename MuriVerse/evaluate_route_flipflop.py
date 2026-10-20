"""Are subcutaneous half-life labels contaminated by flip-flop kinetics?

After intravenous dosing the terminal half-life is an elimination half-life. After
subcutaneous dosing it is only an elimination half-life when absorption is fast relative to
elimination; when absorption is the slower step the terminal phase reports the ABSORPTION
rate instead. That is flip-flop kinetics, and it means an SC label and an IV label can be
measuring two different physical quantities while sharing a column name.

This matters for the model because the cohort pools both routes into one `half_life_h`
target. If the two routes measure different things, the pooled label is heterogeneous and
molecular features cannot explain the SC portion, because absorption is set largely by
formulation and injection site rather than by the molecule.

Two analyses, written to two tables.

1. Paired flip-flop evidence (table32). Peptides measured by BOTH routes in the SAME
   species, so the molecule and the animal are held fixed and only the route varies. A
   systematic SC > IV shift is the flip-flop signature. This is the only comparison here
   that can separate flip-flop from composition: the marginal SC/IV half-life gap is roughly
   4x, but most of that is simply that the two routes were used on different peptides.

2. Route transfer (table33). For each target route, three models scored on the SAME held-out
   rows under one set of grouped folds:

     same route   - trained on the target route's own rows outside the held-out groups
     other route  - trained on the other route's rows outside the held-out groups
     pooled       - trained on every cohort row outside the held-out groups

   Cross-route transfer is the sharp test. If IV and SC labels measure the same property, a
   model trained on one should rank the other about as well as its own-route model does. If
   flip-flop contaminates SC, transfer should degrade, and the SIGNED bias should go
   negative when training on IV and predicting SC - predicting an elimination half-life for
   a label that is really absorption-limited, and therefore longer.

Groups are connected components of (homology cluster, source PMID) over the FULL cohort, so
no cluster or source leaks from a training set into the evaluated rows, in either direction.
SC carries only 73 records in 20 groups, so this uses five seeds rather than the pipeline's
usual three and reports the spread across them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))

from evaluate_halflife_model import (  # noqa: E402
    N_SPLITS, block_slices, fit_predict, load_matrix, union_groups,
)

TABS = Path(__file__).resolve().parent.parent / "figures/tabs"
OUT_PAIRS = TABS / "table32_muriverse_halflife_route_flipflop_pairs.csv"
OUT_TRANSFER = TABS / "table33_muriverse_halflife_route_transfer.csv"
OUT_TRANSFER_SEEDS = TABS / "table33b_muriverse_halflife_route_transfer_by_seed.csv"
# Sample-size control: SC's own-route score is both low and very unstable, which could mean
# its labels are contaminated OR simply that 73 rows in 19 groups is too few to fit on. This
# table scores IV cut down to SC's size so the two are compared at equal n.
OUT_SIZECTL = TABS / "table33c_muriverse_halflife_route_size_control.csv"
N_SIZE_DRAWS = 10

SHIPPED_BLOCKS = ["molecule_annotation", "study_context", "termini", "physchem"]
ROUTES = ("iv", "sc")
# SC is small (73 rows, 20 groups), so average over more seeds than the main pipeline's three.
SEEDS = (0, 1, 2, 3, 4)


def score(y: np.ndarray, pred: np.ndarray) -> dict:
    """Ranking, absolute error, and SIGNED bias.

    Bias is the mean of (predicted - observed) in log10 hours. It is the metric that
    diagnoses flip-flop directly: a model trained on IV elimination half-lives should
    UNDER-predict SC labels that are absorption-limited and therefore longer.
    """
    fold_error = 10.0 ** np.abs(y - pred)
    return {
        "n": int(len(y)),
        "spearman": float(pd.Series(y).corr(pd.Series(pred), method="spearman")),
        "mdafe": float(np.median(fold_error)),
        "within_2fold": float(np.mean(fold_error <= 2.0)),
        "bias": float(np.mean(pred - y)),
    }


def paired_flipflop(meta: pd.DataFrame) -> pd.DataFrame:
    """Same sequence, same species, both routes: the route is the only thing that varies."""
    d = meta[meta["route_class"].isin(ROUTES)]
    per_cell = (d.groupby(["sequence", "species", "route_class"], as_index=False)
                  .agg(log_label=("log_label", "mean"), n=("id", "size")))
    wide = per_cell.pivot_table(index=["sequence", "species"], columns="route_class",
                                values="log_label").dropna(subset=list(ROUTES)).reset_index()
    wide["delta"] = wide["sc"] - wide["iv"]
    return wide


def transfer(X, y, groups, route, target, seed) -> tuple[np.ndarray, dict]:
    """Out-of-fold predictions on the target route's rows, from three training pools."""
    target_mask = route == target
    other = [r for r in ROUTES if r != target][0]
    pools = {
        "same route": route == target,
        "other route": route == other,
        "pooled": np.ones(len(y), dtype=bool),
    }

    eval_idx = np.flatnonzero(target_mask)
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
        raise SystemExit(f"some {target} records were never held out; check the grouping key")
    return eval_idx, out


def size_control(X, y, groups, route, n_target_rows, n_target_groups) -> list[dict]:
    """Score IV cut down to SC's size, so own-route performance is compared at equal n.

    Whole groups are drawn (never partial groups), accumulating until the row budget is
    met, so the subsample keeps the same leak-free structure as the full evaluation.
    """
    iv_groups = np.array(sorted(set(groups[route == "iv"])))
    rows = []
    for draw in range(N_SIZE_DRAWS):
        rng = np.random.default_rng(1000 + draw)
        # Match the ROW budget, not the group count. Model capacity scales with rows, and IV
        # groups are smaller on average than SC's, so capping on groups would hand the
        # control far fewer rows than SC has and rig the comparison against IV.
        picked, n_rows = [], 0
        for g in rng.permutation(iv_groups):
            if n_rows >= n_target_rows:
                break
            picked.append(g)
            n_rows += int(((groups == g) & (route == "iv")).sum())
        mask = (route == "iv") & np.isin(groups, picked)
        idx = np.flatnonzero(mask)
        if len(set(groups[idx])) < N_SPLITS:
            continue
        splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=draw)
        pred = np.full(len(idx), np.nan)
        for _, test_pos in splitter.split(idx, y[idx], groups[idx]):
            held = set(groups[idx][test_pos])
            train_rows = np.flatnonzero(mask & ~np.isin(groups, list(held)))
            pred[test_pos] = fit_predict(X[train_rows], y[train_rows], X[idx[test_pos]], draw)
        rows.append({"draw": draw, "n_groups": len(picked), **score(y[idx], pred)})
    return rows


def main() -> None:
    root = Path(__file__).resolve().parent
    X_all, meta, layout = load_matrix("rodent_strict", root / "rodent" / "strict")
    sl = block_slices(layout)
    missing = [b for b in SHIPPED_BLOCKS if b not in sl]
    if missing:
        raise SystemExit(f"saved dataset lacks blocks {missing}; rebuild the ablation dataset")
    X = np.hstack([X_all[:, sl[b]] for b in SHIPPED_BLOCKS])

    y = meta["log_label"].to_numpy(np.float64)
    route = meta["route_class"].astype(str).str.lower().to_numpy()
    groups = union_groups(meta["cluster"].to_numpy(), meta["PMID"].to_numpy())

    # --- 1. paired flip-flop evidence ---------------------------------------------------
    pairs = paired_flipflop(meta)
    delta = pairs["delta"].to_numpy()
    wilcox = stats.wilcoxon(delta)
    ttest = stats.ttest_1samp(delta, 0.0)
    print(f"[Flip-flop] {len(pairs)} sequence x species pairs measured by both routes "
          f"({pairs['species'].value_counts().to_dict()})")
    print(f"  SC longer in {(delta > 0).sum()}/{len(delta)} pairs ({100 * (delta > 0).mean():.0f}%)")
    print(f"  median delta {np.median(delta):+.3f} log10 (SC = {10 ** np.median(delta):.2f}x IV); "
          f"mean {delta.mean():+.3f} (={10 ** delta.mean():.2f}x)")
    print(f"  paired Wilcoxon p = {wilcox.pvalue:.2g}; paired t-test p = {ttest.pvalue:.2g}")

    pairs = pairs.assign(
        n_pairs=len(pairs), median_delta=np.median(delta), mean_delta=delta.mean(),
        frac_sc_longer=(delta > 0).mean(), wilcoxon_p=wilcox.pvalue, ttest_p=ttest.pvalue)

    # Marginal (unpaired) route medians, to show how much of the raw gap is composition.
    for r in ROUTES:
        m = meta[meta["route_class"] == r]
        print(f"  [marginal] {r.upper()}: n={len(m)} median t1/2 "
              f"{10 ** m['log_label'].median():.2f} h")

    # --- 2. route transfer ---------------------------------------------------------------
    rows = []
    for target in ROUTES:
        n_groups = len(set(groups[route == target]))
        print(f"\n[Transfer] target {target.upper()}: {(route == target).sum()} rows, "
              f"{n_groups} groups")
        for seed in SEEDS:
            eval_idx, preds = transfer(X, y, groups, route, target, seed)
            for src, pred in preds.items():
                rows.append({"target": target, "train": src, "seed": seed,
                             **score(y[eval_idx], pred)})
        last = pd.DataFrame(rows).query("target == @target").groupby("train")[
            ["spearman", "mdafe", "bias"]].mean()
        for src, r in last.iterrows():
            print(f"    train={src:<12s} rho={r.spearman:+.3f}  "
                  f"MdAFE={r.mdafe:.2f}x  bias={r.bias:+.3f} log10")

    per_seed = pd.DataFrame(rows)
    metrics = ["spearman", "mdafe", "within_2fold", "bias"]
    summary = (per_seed.groupby(["target", "train"], as_index=False)
               .agg(n=("n", "first"), n_seeds=("seed", "nunique"),
                    **{c: (c, "mean") for c in metrics},
                    **{f"{c}_sd": (c, "std") for c in metrics}))

    # --- 3. sample-size control ----------------------------------------------------------
    n_sc = int((route == "sc").sum())
    n_sc_groups = len(set(groups[route == "sc"]))
    ctl = pd.DataFrame(size_control(X, y, groups, route, n_sc, n_sc_groups))
    sc_own = summary.query("target == 'sc' and train == 'same route'").iloc[0]
    print(f"\n[Size control] IV subsampled to SC's size ({n_sc} rows / {n_sc_groups} groups), "
          f"{len(ctl)} draws")
    print(f"  IV at SC size : rho {ctl.spearman.mean():+.3f} +/- {ctl.spearman.std():.3f}  "
          f"(median n={int(ctl.n.median())})")
    print(f"  SC own route  : rho {sc_own.spearman:+.3f} +/- {sc_own.spearman_sd:.3f}  "
          f"(n={int(sc_own.n)})")
    print("  => " + ("SC still underperforms IV at equal n: consistent with label contamination"
                     if ctl.spearman.mean() - ctl.spearman.std() > sc_own.spearman
                     else "SC's deficit is NOT separable from its smaller sample size"))

    TABS.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(OUT_PAIRS, index=False)
    summary.to_csv(OUT_TRANSFER, index=False)
    per_seed.to_csv(OUT_TRANSFER_SEEDS, index=False)
    ctl.to_csv(OUT_SIZECTL, index=False)
    for p in (OUT_PAIRS, OUT_TRANSFER, OUT_TRANSFER_SEEDS, OUT_SIZECTL):
        print(f"[Wrote] {p}")


if __name__ == "__main__":
    main()
