"""Evaluate the rodent half-life model with grouped cross-validation and block ablation.

Every record receives an out-of-fold prediction under three groupings: HMMER homology
cluster, source identifier, and the union of the two (a held-out group shares neither a
cluster nor a source with training), which is the headline. The saved 80/20 split is an artifact for
training data consumers, not the reported score.

The model is `halflife_model.HalfLifeEnsemble` (ridge + RBF SVR + XGBoost, equal weight) with
every hyperparameter fixed in advance, so nothing is tuned on the records it is scored on.
Fold membership is reshuffled per seed (GroupKFold with shuffle), so the score does not hang on
one lopsided group, and each seed is a complete out-of-fold pass. Reported are Pearson r,
Spearman rho and the AUC of the continuous prediction against a half-life of at least 6 h, each with
a 95% group bootstrap over the union groups. After evaluation the selected and
structure-reading configurations are refitted on all rows for application. Training-set
predictions are never scores.

`--nested` additionally estimates the cost of CHOOSING the configuration: the configuration is
picked by inner cross-validation inside each outer training fold, and only the outer fold's
held-out records are scored.

  python MuriVerse/evaluate_halflife_model.py --cohort rodent --tier strict
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from datasets import load_from_disk
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

from halflife_model import MEMBERS, HalfLifeEnsemble

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
DATA = REPO / "results/half_life_muriverse"
TABS = REPO / "figures/tabs"

N_SPLITS = 5
INNER_SPLITS = 3
SEEDS = (0, 1, 2)
N_BOOT = 500
HALF_LIFE_CUTOFF_H = 6.0   # the screening label for the AUC: a half-life at or above this is positive

# Configurations that only price a block; never offered as the shipped model.
DIAGNOSTIC = {"molecule annotation only", "study context only"}

# The SHIPPED configuration is the one that gets applied to molecules the cohort does not
# contain (predict_polymyxins.py), so it has two required blocks: `molecule_annotation`
# and `study_context`. Each is required for a reason that stands on its own.
#
# `study_context` is required because it is the only block that carries the ANIMAL. Species
# is a one-hot pair in it (`species=mouse`, `species=rat`), and a shipped configuration
# without it would return one number for a mouse and a rat dosed with the same molecule -
# which is exactly the defect the species fix removed, and it is not a defect that should
# be reachable by a 0.01 rho difference in the ablation. The block is also at least as
# definable for an outside molecule as `molecule_annotation` is: every column in it is a
# study-design choice known whenever a study is being planned.
#
# `molecule_annotation` is required so the shipped model always carries the engineering
# class, topology and chirality the curator recorded; without it the model would commit to
# predicting a molecule it has no bucketed description of. On the IV cohort the
# top-scoring configuration that drops this block, `physchem + study context + termini`,
# scores rho = 0.587 against the shipped set's 0.574 - a difference of 0.013, well inside
# the seed spread and the bootstrap interval. Both configurations give physically plausible
# predictions on the Figure 4 polymyxins (shipped 1.0-2.5 h, unconstrained 0.5-1.7 h), so
# the choice is between carrying the annotation block and shedding it for 0.013 in rho.
#
# A previous build reported this constraint as the difference between applicable and
# meaningless predictions: it quoted the unconstrained configuration returning polymyxin
# half-lives of about 10^-120 h, and the SMILES-embedding configuration returning about
# 10^-131 h. Both of those extreme values were a block-order bug in save_production_models
# and predict_polymyxins.assemble - the saved block list was in layout order while training
# had concatenated the slices in configuration-list order, so inference assembled polymyxin
# features in a different column order than training. The scaler stats for a sparse one-hot
# column (mean 0.005, scale 0.068) were then applied to polymyxin's molecular weight
# (1203), producing z-scores in the thousands and ridge predictions of -400 log10 h.
# The fix is in `_blocks_of` here (preserve concatenation order) and
# `predict_polymyxins.assemble` (do not reorder the saved block list). Post-fix the
# applicability-pathology argument is gone and the argument above - species via
# study_context and carrying the molecular annotations - is what the constraint rests on.
#
# This constrains only which configuration is SHIPPED. The ablation still scores and reports
# every configuration, including the unconstrained winner, so the figure shows what the
# constraint cost.
REQUIRED_SHIPPED_BLOCKS = {"molecule_annotation", "study_context"}


def _shippable(name: str, blocks: dict) -> bool:
    """Is this configuration one we could apply to a molecule outside the cohort?"""
    return (name not in DIAGNOSTIC
            and REQUIRED_SHIPPED_BLOCKS.issubset(set(blocks.get(name, ()))))
THREADS = 16
MAX_NESTED_DIMS = 100


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    return float(spearmanr(a, b).correlation)


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a, b)[0, 1])


def _auc(y: np.ndarray, pred: np.ndarray) -> float:
    positive = y >= np.log10(HALF_LIFE_CUTOFF_H) - 1e-9   # half-life >= the cutoff
    if positive.all() or not positive.any():
        return float("nan")
    return float(roc_auc_score(positive, pred))


def fit_predict(X_train, y_train, X_test, seed: int) -> np.ndarray:
    model = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(X_train, y_train)
    return model.predict(X_test)


def union_groups(*keys: np.ndarray) -> np.ndarray:
    """Connected components over several grouping keys: rows sharing ANY key share a group."""
    parent: dict = {}

    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for row in zip(*keys):
        roots = [find((i, v)) for i, v in enumerate(row)]
        for r in roots[1:]:
            parent[r] = roots[0]
    return np.array([str(find((0, v))) for v in keys[0]])


def bootstrap_ci(y: np.ndarray, pred: np.ndarray, groups: np.ndarray
                 ) -> dict[str, tuple[float, float]]:
    """95% intervals for Spearman, Pearson and AUC, resampling whole groups with replacement."""
    rng = np.random.default_rng(0)
    members = {g: np.flatnonzero(groups == g) for g in np.unique(groups)}
    keys = list(members)
    stats = {"spearman": [], "pearson": [], "auc": []}
    for _ in range(N_BOOT):
        pick = rng.choice(len(keys), size=len(keys), replace=True)
        idx = np.concatenate([members[keys[i]] for i in pick])
        if np.ptp(y[idx]) > 0:
            stats["spearman"].append(_rho(y[idx], pred[idx]))
            stats["pearson"].append(_pearson(y[idx], pred[idx]))
            auc = _auc(y[idx], pred[idx])
            if not np.isnan(auc):
                stats["auc"].append(auc)
    return {k: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
            for k, v in stats.items()}


def load_matrix(cohort: str, meta_dir: Path | None = None
                ) -> tuple[np.ndarray, pd.DataFrame, dict]:
    """The full-block feature matrix in metadata row order, plus its block layout.

    Reads the `_ablation_` dataset, which always carries all three blocks, rather than the
    shipped one, which by default carries only the two the ablation selected - otherwise
    this script could not measure the block it caused to be dropped.

    `save_dataset` writes the rows in metadata order and then splits them, and both
    `Dataset.from_dict` and `.filter` preserve order, but the two splits are saved
    separately - so the rows are re-joined on `id` rather than by concatenating, which would
    silently interleave train and val wrongly if that ever stopped holding.
    """
    ds = load_from_disk(str(DATA / f"halflife_{cohort}_ablation_with_embeddings"))
    vectors = {i: np.asarray(e, dtype=np.float32)
               for split in ds for i, e in zip(ds[split]["id"], ds[split]["embedding"])}
    meta = pd.read_csv((meta_dir or ROOT / cohort) / "halflife_meta_with_split.csv",
                       low_memory=False)
    missing = set(meta["id"]) - set(vectors)
    extra = set(vectors) - set(meta["id"])
    if missing or extra:
        raise SystemExit(f"Feature dataset and metadata IDs differ: {len(missing)} missing "
                         f"vectors and {len(extra)} extra vectors; rebuild both with "
                         f"build_halflife_training_set.py")
    layout = json.loads((DATA / f"halflife_{cohort}_ablation_with_embeddings"
                         / "feature_names.json").read_text())
    return np.stack([vectors[i] for i in meta["id"]]), meta, layout


def block_slices(layout: dict) -> dict[str, slice]:
    return {b["name"]: slice(b["start"], b["start"] + b["dim"]) for b in layout["blocks"]}


def _folds(X, y, groups, seeds=SEEDS):
    """One complete out-of-fold pass per seed; each seed reshuffles which groups share a fold."""
    preds = np.full((len(seeds), len(y)), np.nan)
    for row, seed in enumerate(seeds):
        splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        for train_idx, test_idx in splitter.split(X, y, groups):
            preds[row, test_idx] = fit_predict(X[train_idx], y[train_idx], X[test_idx], seed)
    if np.isnan(preds).any():
        raise SystemExit("some records were never held out; check the grouping key")
    return preds


STRATUM_METRICS = ("spearman", "pearson", "auc", "bias")


def cross_validate(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                   strata: dict[str, np.ndarray]) -> tuple[list[dict], np.ndarray]:
    """Per-seed scores and the seed-averaged out-of-fold prediction.

    `strata` maps a stratification name to a per-row label array; every distinct value gets
    its own set of columns. Species was the only stratification for a long time, which hid
    a route effect: IV and SC terminal half-lives do not share an absolute scale, because an
    SC terminal slope can report absorption rather than elimination (flip-flop kinetics), so
    a route-blind summary averages over two regimes. Stratum values are assumed distinct
    across stratifications - they are here, {mouse, rat} against {im, ip, iv, sc} - so the
    column suffix is the value alone and existing per-species column names are unchanged.
    """
    preds = _folds(X, y, groups)
    rows = []
    for seed, oof in zip(SEEDS, preds):
        row = {"seed": seed, **_score(y, oof)}
        for values in strata.values():
            for v in sorted(set(values)):
                m = values == v
                scored = _score(y[m], oof[m])
                for metric in STRATUM_METRICS:
                    row[f"{metric}_{v}"] = scored[metric]
        rows.append(row)
    return rows, preds.mean(axis=0)


def _score(y: np.ndarray, pred: np.ndarray) -> dict:
    """Ranking, error and SIGNED bias.

    `bias` is mean(predicted - observed) in log10 hours, and it is reported because the
    error metrics above it cannot see a systematic offset: RMSE and MAE are symmetric, and
    rho and AUC are invariant to any monotone shift. A stratum whose labels sit on a
    different absolute scale from the training mixture therefore looks fine on every other
    column here while being consistently wrong in one direction. Route is the case that
    forced this - see `evaluate_route_flipflop.py` - so bias travels with every stratum.
    """
    return {
        "n": int(len(y)),
        "spearman": _rho(y, pred),
        "pearson": _pearson(y, pred),
        "auc": _auc(y, pred),
        "bias": float(np.mean(pred - y)),
        "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
        "mae": float(np.mean(np.abs(y - pred))),
        "r2": float(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()),
    }


def _init_worker(nthread: int) -> None:
    """Give each worker its share of the thread budget instead of all of it."""
    global THREADS
    THREADS = nthread


def _score_configuration(job: tuple) -> dict:
    """Score one (grouping, configuration) pair. Pure in its arguments, so it parallelises.

    The folds are fixed by the grouping key and the seeds are fixed constants, so this
    returns the same summary whether it is called in this process or another one.
    """
    group_name, groups, boot_groups, name, X_sub, y, strata = job
    scored, mean_pred = cross_validate(X_sub, y, groups, strata)
    frame = pd.DataFrame(scored)
    ci = bootstrap_ci(y, mean_pred, boot_groups)
    summary = {
        "grouping": group_name,
        "features": name,
        "dims": int(X_sub.shape[1]),
        "n": int(len(y)),
        "model": "ridge + SVR + XGBoost",
        "r2": frame["r2"].mean(),
        "rmse": frame["rmse"].mean(),
        "mae": frame["mae"].mean(),
        "bias": frame["bias"].mean(),
    }
    for metric in ("spearman", "pearson", "auc"):
        summary[metric] = frame[metric].mean()
        summary[f"{metric}_sd"] = frame[metric].std()
        summary[f"{metric}_lo"], summary[f"{metric}_hi"] = ci[metric]
    for values in strata.values():
        for v in sorted(set(values)):
            for metric in STRATUM_METRICS:
                summary[f"{metric}_{v}"] = frame[f"{metric}_{v}"].mean()
    return summary


# Output namespace for a route-restricted run. Empty for the default pooled-route cohort, so
# an unrestricted run writes exactly the paths it always wrote; set in main() from --route.
TAG = ""


def main() -> None:
    global TAG
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=1,
                    help="Score this many (grouping x configuration) pairs at once. The "
                         "pairs are independent over fixed folds, so this changes the wall "
                         "clock only; the thread budget is divided among the workers.")
    ap.add_argument("--nested", action="store_true",
                    help="Also estimate the optimism of choosing the configuration on the "
                         "same folds that score it (writes table31b).")
    ap.add_argument("--models-only", action="store_true",
                    help="Skip scoring: refit and save the production models using the "
                         "existing ablation table to pick the configurations.")
    ap.add_argument("--cohort", default="rodent")
    ap.add_argument("--tier", default="strict", choices=["strict", "pooled"],
                    help="Which tier's saved dataset to score. 'strict' is the default, "
                         "matching build_halflife_training_set.py.")
    # Restricting the cohort to one route is a measurement, not a new shipped model, so it
    # writes to its own namespace (see `suffix` below) and leaves the pooled tier's
    # oof_predictions.csv, table31 and production models exactly where they were. The
    # motivation is physical: only an intravenous dose has no absorption phase, so IV
    # terminal half-lives are elimination half-lives throughout, while extravascular ones
    # are a mixture of elimination and absorption rates (evaluate_route_flipflop.py). A
    # route-restricted cohort therefore has ONE physical quantity as its target.
    ap.add_argument("--route", default=None,
                    help="Score only the records with this route_class (e.g. 'iv'), or a "
                         "comma-separated set. Route becomes constant under the restriction, "
                         "so the route columns of the molecule_annotation and study_context "
                         "blocks carry no information and the model is fitted without their "
                         "variance - StandardScaler leaves a zero-variance column at zero.")
    ap.add_argument("--models-suffix", default=None,
                    help="Override the output namespace. Defaults to the route restriction.")
    args = ap.parse_args()

    strict = args.tier == "strict"
    dataset = f"{args.cohort}_strict" if strict else args.cohort
    meta_dir = ROOT / args.cohort / "strict" if strict else ROOT / args.cohort
    X, meta, layout = load_matrix(dataset, meta_dir)

    # Row restriction happens here, before anything derives folds, groups or strata from the
    # metadata, so every downstream object describes the restricted cohort and none of them
    # can be left describing the full one.
    route_tag = ""
    if args.route:
        keep_routes = [r.strip() for r in args.route.split(",") if r.strip()]
        mask = meta["route_class"].astype(str).isin(keep_routes).to_numpy()
        if not mask.any():
            raise SystemExit(f"No records with route_class in {keep_routes}; the cohort has "
                             f"{sorted(set(meta['route_class'].astype(str)))}")
        before = len(meta)
        X, meta = X[mask], meta.loc[mask].reset_index(drop=True)
        route_tag = "_" + "".join(sorted(keep_routes))
        print(f"[Route] restricted to {'/'.join(keep_routes)}: {len(meta)} of {before} records")
    TAG = args.models_suffix if args.models_suffix is not None else route_tag
    y = meta["log_label"].to_numpy(np.float64)
    species = meta["species"].to_numpy()
    # Species and route are reporting strata here AND model inputs: species is the
    # `species=mouse`/`species=rat` pair in the study_context block and route is
    # `route_model` in the same block. Stratifying by them changes what is reported, not
    # what is fitted - but a route-blind summary averages IV and SC, whose terminal
    # half-lives do not share an absolute scale, and a species-blind one hides the offset
    # the species columns exist to carry.
    route = meta["route_class"].astype(str).to_numpy()
    strata = {"species": species, "route": route}
    sl = block_slices(layout)
    print(f"[Data] {X.shape[0]} records x {X.shape[1]} features  "
          f"({', '.join(f'{k} {v.stop - v.start}' for k, v in sl.items())})")

    # Molecular identity and experimental context are priced separately. This directly tests
    # whether the model transfers molecular signal or mainly recognizes animal/study setup.
    configurations = {
        "molecule annotation only": ["molecule_annotation"],
        "study context only": ["study_context"],
        "termini only": ["termini"],
        "sequence (ESM-2) only": ["sequence"],
        "SMILES (PeptideCLM) only": ["smiles"],
        "sequence + SMILES": ["sequence", "smiles"],
        "molecule annotation + study context": ["molecule_annotation", "study_context"],
        "study context + termini": ["study_context", "termini"],
        "sequence + termini": ["sequence", "termini"],
        "SMILES + termini": ["smiles", "termini"],
        "sequence + SMILES + termini": ["sequence", "smiles", "termini"],
        "sequence + molecule annotation": ["sequence", "molecule_annotation"],
        "SMILES + molecule annotation": ["smiles", "molecule_annotation"],
        "sequence + SMILES + molecule annotation": ["sequence", "smiles", "molecule_annotation"],
        "sequence + study context": ["sequence", "study_context"],
        "SMILES + study context": ["smiles", "study_context"],
        "sequence + SMILES + study context": ["sequence", "smiles", "study_context"],
        "sequence + SMILES + both annotations": ["sequence", "smiles", "molecule_annotation",
                                                 "study_context"],
        "sequence + SMILES + context + termini": ["sequence", "smiles", "study_context",
                                                   "termini"],
        "sequence + SMILES + molecule + termini": ["sequence", "smiles",
                                                    "molecule_annotation", "termini"],
        "all molecular blocks + study context": ["sequence", "smiles", "molecule_annotation",
                                                 "study_context", "termini", "physchem"],
        # Interpretable descriptors read from the verified molecular structure.
        "physchem only": ["physchem"],
        "physchem + study context": ["physchem", "study_context"],
        "physchem + molecule annotation": ["physchem", "molecule_annotation"],
        "physchem + study context + termini": ["physchem", "study_context", "termini"],
        "SMILES + physchem + study context": ["smiles", "physchem", "study_context"],
        "sequence + physchem + study context": ["sequence", "physchem", "study_context"],
        # The compact set: every low-dimensional block and no embedding.
        "annotation + context + termini + physchem": ["molecule_annotation", "study_context",
                                                      "termini", "physchem"],
    }
    # A block can be absent in older saved datasets; configurations that need one are
    # skipped rather than silently scored on a different feature set.
    skipped = {name for name, blocks in configurations.items()
               if any(b not in sl for b in blocks)}
    if skipped:
        print(f"[Blocks] {', '.join(sorted(sl))} present; "
              f"skipping {len(skipped)} configuration(s) needing a block this tier has "
              f"not got: {', '.join(sorted(skipped))}")
    # The block NAMES are kept beside the slice lists, because the shipped-configuration
    # constraint below is stated in terms of blocks and the slice lists cannot be read back.
    config_blocks = {name: tuple(blocks) for name, blocks in configurations.items()
                     if name not in skipped}
    configurations = {name: [sl[b] for b in blocks]
                      for name, blocks in configurations.items()
                      if name not in skipped}

    key_cluster = meta["cluster"].to_numpy()
    key_pmid = meta["PMID"].astype(str).to_numpy()
    grouping = {"cluster": key_cluster, "PMID": key_pmid,
                "cluster_PMID": union_groups(key_cluster, key_pmid)}
    boot_groups = grouping["cluster_PMID"]
    print(f"[Groups] cluster {len(set(key_cluster))}, PMID {len(set(key_pmid))}, "
          f"cluster+PMID {len(set(boot_groups))} (largest "
          f"{pd.Series(boot_groups).value_counts().iloc[0]} records)")

    if args.models_only:
        table = pd.read_csv(TABS / f"table31_muriverse_halflife_feature_ablation"
                            f"{('_strict' if strict else '') + TAG}.csv")
        cluster_scores = (table[table["grouping"] == "cluster_PMID"]
                          .set_index("features")["spearman"])
        shipped = max((c for c in configurations
                       if c in cluster_scores and _shippable(c, config_blocks)),
                      key=lambda c: cluster_scores[c])
        save_production_models(X, y, sl, configurations, cluster_scores, shipped, meta_dir,
                               strata)
        return

    # The label's own mean is the floor every configuration has to clear; a tree given a
    # single constant column reproduces exactly that, so it is run through the same folds
    # rather than computed analytically.
    constant = np.zeros((len(y), 1), dtype=np.float32)
    jobs = [(group_name, groups, boot_groups, name, X_sub, y, strata)
            for group_name, groups in grouping.items()
            for name, X_sub in [("baseline (train mean)", constant)] + [
                (n, np.concatenate([X[:, s] for s in parts], axis=1))
                for n, parts in configurations.items()]]

    # Each (grouping, configuration) pair is an independent computation over the same fixed
    # folds, so running several at once changes the wall clock and nothing else. It is opt-in
    # (--jobs) and off by default so the script's default behaviour is the sequential one.
    # The worker count is spent on `nthread` when it is 1, so the total thread budget is the
    # same either way; `nthread` is already excluded from the saved params as an environment
    # detail rather than part of the configuration being reported.
    if args.jobs > 1:
        print(f"\n[Parallel] {len(jobs)} (grouping x configuration) pairs over "
              f"{args.jobs} workers, {max(1, THREADS // args.jobs)} threads each")
        with ProcessPoolExecutor(max_workers=args.jobs, initializer=_init_worker,
                                 initargs=(max(1, THREADS // args.jobs),)) as pool:
            records = list(pool.map(_score_configuration, jobs))
        for summary in records:
            print(f"  [{summary['grouping']:12s}] {summary['features']:42s} "
                  f"dims={summary['dims']:5d}  r={summary['pearson']:.3f} "
                  f"rho={summary['spearman']:.3f} auc={summary['auc']:.3f}")
    else:
        records = []
        for job in jobs:
            summary = _score_configuration(job)
            records.append(summary)
            print(f"  [{summary['grouping']:12s}] {summary['features']:42s} "
                  f"dims={summary['dims']:5d}  r={summary['pearson']:.3f} "
                  f"rho={summary['spearman']:.3f} auc={summary['auc']:.3f}")

    out = pd.DataFrame(records)

    # Out-of-fold predictions for the configuration that scores best on held-out homology
    # clusters, which is the honest grouping on this cohort. Written per record so the
    # figures plot predictions the model made without having seen the record OR any
    # relative of it, rather than a single small validation draw.
    cluster_scores = out[out["grouping"] == "cluster_PMID"].set_index("features")["spearman"]
    shipped = max((c for c in configurations
                   if c in cluster_scores and _shippable(c, config_blocks)),
                  key=lambda c: cluster_scores[c])
    # What the applicability constraint cost, printed so it is visible in the run log rather
    # than only discoverable by reading the ablation table.
    unconstrained = max((c for c in configurations
                         if c in cluster_scores and c not in DIAGNOSTIC),
                        key=lambda c: cluster_scores[c])
    print(f"\n[Shipped] best APPLICABLE on held-out clusters and sources: {shipped} "
          f"(rho = {cluster_scores[shipped]:.3f})")
    if unconstrained != shipped:
        print(f"[Shipped] unconstrained best was {unconstrained} "
              f"(rho = {cluster_scores[unconstrained]:.3f}); it omits "
              f"{'/'.join(sorted(REQUIRED_SHIPPED_BLOCKS))} and so cannot be applied to a "
              f"molecule outside the cohort - see REQUIRED_SHIPPED_BLOCKS")
    blocks = np.concatenate([X[:, s] for s in configurations[shipped]], axis=1)
    # `route_class` travels with the saved out-of-fold predictions so any downstream
    # report can stratify by route without re-joining the metadata.
    predictions = meta[["id", "species", "route_class", "cluster", "PMID",
                        "log_label", "half_life_h"]].copy()
    predictions["features"] = shipped
    # `y_pred_<grouping>` is the average of the per-seed passes, convenient for plotting. A
    # deployed model is ONE fit, so the headline scores are the mean over the per-seed passes
    # kept in `y_pred_<grouping>_seed<k>`; averaging three fold shuffles is a small bagging
    # ensemble and scores a little higher than any single pass.
    for name, groups in grouping.items():
        passes = _folds(blocks, y, groups)
        predictions[f"y_pred_{name}"] = passes.mean(axis=0)
        for seed, row in zip(SEEDS, passes):
            predictions[f"y_pred_{name}_seed{seed}"] = row
    predictions = predictions.rename(columns={"log_label": "y_true"})
    oof_path = meta_dir / f"oof_predictions{TAG}.csv"
    predictions.to_csv(oof_path, index=False)
    print(f"[Saved] {oof_path}")

    TABS.mkdir(parents=True, exist_ok=True)
    # One table per tier. The two score different cohorts over different block sets, so
    # writing both to one path silently replaced the published pooled-tier numbers with
    # strict-tier ones that are not comparable to them.
    suffix = ("_strict" if strict else "") + TAG
    path = TABS / f"table31_muriverse_halflife_feature_ablation{suffix}.csv"
    out.to_csv(path, index=False)
    print(f"\n[Saved] {path}")

    if args.nested:
        nested_confirmation(X, y, species, configurations, boot_groups)

    save_production_models(X, y, sl, configurations, cluster_scores, shipped, meta_dir,
                           strata)


def _nested_fold(job: tuple) -> tuple:
    """Choose a configuration by inner grouped CV on the training rows, then score the fold."""
    seed, train_idx, test_idx, groups, names, features, y = job
    inner = GroupKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=seed)
    scores = []
    for name in names:
        Xf = features[name]
        pred = np.zeros(len(train_idx))
        for a, b in inner.split(Xf[train_idx], y[train_idx], groups[train_idx]):
            pred[b] = fit_predict(Xf[train_idx][a], y[train_idx][a], Xf[train_idx][b], seed)
        scores.append(_pearson(y[train_idx], pred))
    chosen = names[int(np.argmax(scores))]
    Xf = features[chosen]
    return test_idx, fit_predict(Xf[train_idx], y[train_idx], Xf[test_idx], seed), chosen


def nested_confirmation(X, y, species, configurations, groups) -> None:
    """The headline score when the feature configuration is itself chosen inside each fold."""
    # The candidate set is fixed in advance as every non-diagnostic configuration of at most
    # MAX_NESTED_DIMS columns. The embedding configurations are left out of the SELECTION (not
    # out of the ablation): every one of them scores below the compact set in the ablation,
    # and searching 2,000-column configurations inside each fold costs hours.
    features = {n: np.concatenate([X[:, s] for s in configurations[n]], axis=1)
                for n in configurations if n not in DIAGNOSTIC}
    features = {n: f for n, f in features.items() if f.shape[1] <= MAX_NESTED_DIMS}
    names = list(features)
    preds = np.full((len(SEEDS), len(y)), np.nan)
    picks = []
    global THREADS
    THREADS = 1
    jobs = []
    for row, seed in enumerate(SEEDS):
        splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        for train_idx, test_idx in splitter.split(X, y, groups):
            jobs.append((row, (seed, train_idx, test_idx, groups, names, features, y)))
    # `spawn`, not fork: the parent has already run multithreaded XGBoost, and forking after
    # OpenMP threads exist can deadlock the children.
    with ProcessPoolExecutor(max_workers=16, mp_context=multiprocessing.get_context("spawn"),
                             initializer=_init_worker, initargs=(1,)) as pool:
        results = list(pool.map(_nested_fold, [j for _, j in jobs]))
    for (row, _), (test_idx, pred, chosen) in zip(jobs, results):
        preds[row, test_idx] = pred
        picks.append(chosen)
    mean_pred = preds.mean(axis=0)
    ci = bootstrap_ci(y, mean_pred, groups)
    rows = []
    for row, seed in enumerate(SEEDS):
        rows.append(_score(y, preds[row]))
    frame = pd.DataFrame(rows)
    summary = {"selection": "configuration chosen inside each training fold",
               "n": int(len(y)), "n_configurations": len(names), "max_dims": MAX_NESTED_DIMS}
    for metric in ("pearson", "spearman", "auc"):
        summary[metric] = frame[metric].mean()
        summary[f"{metric}_lo"], summary[f"{metric}_hi"] = ci[metric]
    for sp in sorted(set(species)):
        m = species == sp
        per = pd.DataFrame([_score(y[m], p[m]) for p in preds])
        for metric in ("pearson", "spearman", "auc"):
            summary[f"{metric}_{sp}"] = per[metric].mean()
    summary["configurations_chosen"] = json.dumps(
        pd.Series(picks).value_counts().to_dict())
    path = TABS / f"table31b_muriverse_halflife_nested_selection{TAG}.csv"
    pd.DataFrame([summary]).to_csv(path, index=False)
    print(f"\n[Nested] pearson {summary['pearson']:.3f} spearman {summary['spearman']:.3f} "
          f"auc {summary['auc']:.3f}  chosen: {summary['configurations_chosen']}")
    print(f"[Saved] {path}")


def _blocks_of(parts: list, sl: dict[str, slice]) -> list[str]:
    """Which named blocks a configuration is built from, in the ORDER THEY ARE CONCATENATED.

    The order matters: training builds `features = np.concatenate([X[:, s] for s in parts])`,
    so the saved scaler, ridge weights and SVR support vectors were all fit against this
    concrete column order. Returning the blocks in `sl.items()` order instead - which this
    function used to do - silently re-ordered the saved blocks list for any configuration
    whose own ordering disagreed with the layout's, so an applying script that trusted the
    saved list concatenated its features in a different order than training did. For
    `SMILES + physchem + study context` that put the physchem columns at the positions
    the scaler had learned were sparse one-hot context columns, producing z-scores in the
    thousands and polymyxin predictions of about 10^-131 h.
    """
    slice_to_name = {id(s): name for name, s in sl.items()}
    return [slice_to_name[id(p)] for p in parts]


def save_production_models(X: np.ndarray, y: np.ndarray,
                           sl: dict[str, slice], configurations: dict,
                           cluster_scores: pd.Series, shipped: str,
                           meta_dir: Path, strata: dict[str, np.ndarray]) -> None:
    """Fit on EVERY row and save the result, for applying the model to new molecules.

    Everything above this point is measurement: each number comes from a model that never saw
    the record it is scoring. A model you intend to *apply* should instead use all the data,
    so these are refit on all rows with the same fixed hyperparameters. They are deliberately
    not used for any reported metric - `predict_polymyxins.py` is their only consumer.

    TWO configurations are saved, and the second is not an afterthought:

    * `shipped` - whatever won on held-out homology clusters among the configurations that
      carry both required blocks (see REQUIRED_SHIPPED_BLOCKS). On the strict tier that is
      the compact set of every low-dimensional block and no embedding.
    * `structural` - the best-scoring configuration that reads the molecule's STRUCTURE and
      not its sequence, and `structural_embedding` - the best that does so through the
      LEARNED representation (PeptideCLM) rather than the physicochemical descriptors.
      Both are saved because they fail differently and the contrast is the point: the
      descriptors rank better and land far closer in absolute terms, while only the
      embedding returns a different number for each congener. The exclusion of the sequence block is deliberate rather than
      incidental: this model exists to be applied to molecules the cohort does not contain,
      and some of those have no faithful single-letter sequence at all - six of the ten
      residues of a polymyxin are 2,4-diaminobutyric acid, which has no code. Any ESM-2
      input for such a molecule is a fiction, so the configuration offered for application
      is one that never asks for it.

    The ablation selects `shipped` on cross-validated ranking, but the low-dimensional
    blocks encode only bucketed study-design and modification class: `hl_engineering`,
    `topology`, `chirality` in `molecule_annotation`, and `route_model`, `species`,
    `sample_matrix`, `assay_method`, `strain`, `sex`, formulation and dose in
    `study_context`. None of those distinguishes one molecule from another within a class.
    Two different drugs with the same modification class, topology, chirality, species and
    route therefore receive the IDENTICAL feature vector and so, necessarily, the identical
    prediction - which is exactly what happens to the four polymyxin congeners. A
    configuration that can see the molecule is required for the application to say anything
    congener-specific at all, so it is saved alongside and the gap between the two is
    reported rather than papered over. What the shipped set DOES separate is the animal:
    since the species fix put `species` in `study_context`, a mouse and a rat row are no
    longer the same vector, so the degeneracy is across congeners only and not across
    species as well.
    """
    # A route-restricted run is a measurement, so its production models go to their own
    # directory and predict_polymyxins.py keeps loading the pooled-route ones it was built on.
    model_dir = meta_dir / f"model{TAG}"
    model_dir.mkdir(parents=True, exist_ok=True)

    # "Can see the molecule" now means the PeptideCLM embedding of the structure OR the
    # physicochemical descriptors computed from it - both are read off the SMILES, so both
    # are defined for a molecule with no faithful single-letter sequence. The sequence block
    # stays excluded for the reason given above.
    def best_reading(block: str) -> str | None:
        """Best configuration that reads `block` and never the ESM-2 sequence block."""
        return max((c for c in configurations
                    if c in cluster_scores
                    and block in _blocks_of(configurations[c], sl)
                    and "sequence" not in _blocks_of(configurations[c], sl)),
                   key=lambda c: cluster_scores[c], default=None)

    # TWO structure-reading roles are saved, not one, because on this cohort they disagree
    # about something worth showing. `structural` is the best of them on held-out ranking
    # and is the physicochemical block; `structural_embedding` is the best that uses the
    # LEARNED structural representation. The second is kept even though it scores lower
    # because it is the only one that returns a different number for each polymyxin
    # congener, and the comparison between "separates the molecules but ranks worse" and
    # "ranks better but returns one number for four molecules" is the result of Figure 4.
    structural = max((c for c in configurations
                      if c in cluster_scores
                      and any(b in ("smiles", "physchem")
                              for b in _blocks_of(configurations[c], sl))
                      and not any(b == "sequence"
                                  for b in _blocks_of(configurations[c], sl))),
                     key=lambda c: cluster_scores[c], default=None)
    structural_embedding = best_reading("smiles")
    if structural_embedding == structural:
        structural_embedding = None
    # When the selected configuration already reads the structure (it carries the physchem
    # block), a second role with the same features would just duplicate it.
    if structural == shipped:
        structural = None

    saved = {}
    for role, name in (("shipped", shipped), ("structural", structural),
                       ("structural_embedding", structural_embedding)):
        if name is None:
            continue
        blocks = _blocks_of(configurations[name], sl)
        features = np.concatenate([X[:, s] for s in configurations[name]], axis=1)
        model = HalfLifeEnsemble(seed=SEEDS[0], nthread=THREADS).fit(features, y)
        joblib.dump(model, model_dir / f"{role}.joblib")
        saved[role] = {
            "features": name,
            "blocks": blocks,
            "dims": int(features.shape[1]),
            "cluster_pmid_spearman": float(cluster_scores[name]),
            "n_train": int(len(y)),
        }
        print(f"[Model] {role:10s} {name:30s} dims={features.shape[1]:5d} -> "
              f"{model_dir / f'{role}.joblib'}")

    # The block layout and the member list travel with the saved models, so an
    # applying script cannot silently assemble its feature vector in a different order.
    #
    # The training composition travels with them too. Route and species are MODEL INPUTS,
    # not just reporting strata, so a prediction is only defined for a combination the
    # cohort actually contains - and the per-route row counts below are wildly uneven. An
    # applying script that leaves route unset, or sets it to a route with a handful of
    # training rows, gets a confidently wrong number rather than an error, so the counts
    # and the caveat are recorded here where that script can read them.
    composition = {name: {str(v): int((values == v).sum())
                          for v in sorted(set(values))}
                   for name, values in strata.items()}
    (model_dir / "model.json").write_text(json.dumps({
        "tier_rows": int(len(y)),
        "metadata_sha256": hashlib.sha256(
            (meta_dir / "halflife_meta_with_split.csv").read_bytes()).hexdigest(),
        "model": "equal-weight ensemble",
        "members": list(MEMBERS),
        "seed": SEEDS[0],
        "block_order": list(sl),
        "block_slices": {b: [s.start, s.stop] for b, s in sl.items()},
        "training_composition": composition,
        "application_caveats": [
            "Route of administration is a model input. A prediction is route-specific and "
            "must not be read as a prediction for another route: IV and SC terminal "
            "half-lives do not share an absolute scale, because an SC terminal slope can "
            "report absorption rather than elimination (flip-flop kinetics). Transferring "
            "across the two is offset by roughly half a log unit - see "
            "evaluate_route_flipflop.py and figures/tabs/table33.",
            "Per-route performance is very uneven; see "
            "figures/tabs/table34_muriverse_halflife_performance_by_route.csv. The pooled "
            "score is carried by IV. Routes with few training rows rank at or below chance "
            "and the pooled number does not transfer to them.",
        ],
        "models": saved,
    }, indent=2))
    print(f"[Model] layout -> {model_dir / 'model.json'}")


if __name__ == "__main__":
    main()
