"""Compare ensemble members and a fold-local species-median baseline on identical IV folds."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "MuriVerse"))

from evaluate_halflife_model import (  # noqa: E402
    N_SPLITS,
    SEEDS,
    THREADS,
    block_slices,
    load_matrix,
    union_groups,
)
from halflife_model import HalfLifeEnsemble  # noqa: E402

ANALYSIS = Path(__file__).resolve().parents[2] / "results" / "revision_analysis"
SNAPSHOT = REPO / "results/revision_analysis/snapshot/oof_predictions_iv.csv"
OUT_PREDICTIONS = ANALYSIS / "model_member_baselines_oof.csv"
OUT_SUMMARY = ANALYSIS / "model_member_baselines.csv"
FEATURE_BLOCKS = ("molecule_annotation", "study_context", "termini", "physchem")
N_BOOT = 500


def metrics(y: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    error = y - pred
    return float(np.mean(np.abs(error))), float(np.sqrt(np.mean(error**2)))


def main() -> int:
    full_x, meta, layout = load_matrix(
        "rodent_strict", meta_dir=REPO / "MuriVerse/rodent/strict"
    )
    reference = pd.read_csv(SNAPSHOT)
    meta = meta.reset_index(drop=True)
    iv_mask = meta["route_class"].eq("iv").to_numpy()
    iv_meta = meta.loc[iv_mask].reset_index(drop=True)
    reference = reference.set_index("id").loc[iv_meta["id"]].reset_index()

    slices = block_slices(layout)
    missing_blocks = set(FEATURE_BLOCKS) - set(slices)
    if missing_blocks:
        raise SystemExit(f"Feature layout is missing blocks: {sorted(missing_blocks)}")
    X = np.concatenate([full_x[iv_mask, slices[name]] for name in FEATURE_BLOCKS], axis=1)
    y = reference["y_true"].to_numpy(float)
    groups = union_groups(iv_meta["cluster"].to_numpy(), iv_meta["PMID"].to_numpy())
    species = iv_meta["species"].to_numpy()

    names = ("ridge", "svr", "xgboost", "ensemble", "species_median")
    predictions = {name: np.full((len(SEEDS), len(y)), np.nan) for name in names}
    for seed_index, seed in enumerate(SEEDS):
        splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        for train_idx, test_idx in splitter.split(X, y, groups):
            fitted = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(X[train_idx], y[train_idx])
            members = fitted.member_predictions(X[test_idx])
            for name in ("ridge", "svr", "xgboost"):
                predictions[name][seed_index, test_idx] = members[name]
            predictions["ensemble"][seed_index, test_idx] = np.mean(
                [members[name] for name in ("ridge", "svr", "xgboost")], axis=0
            )
            train_species = species[train_idx]
            for animal in np.unique(species[test_idx]):
                test_species = test_idx[species[test_idx] == animal]
                same_species = y[train_idx][train_species == animal]
                fallback = float(np.median(y[train_idx]))
                predictions["species_median"][seed_index, test_species] = (
                    float(np.median(same_species)) if len(same_species) else fallback
                )

    if any(np.isnan(values).any() for values in predictions.values()):
        raise SystemExit("At least one record has no out-of-fold prediction")

    seed_columns = [f"y_pred_cluster_PMID_seed{seed}" for seed in SEEDS]
    saved_ensemble = reference[seed_columns].to_numpy(float).T
    if not np.allclose(predictions["ensemble"], saved_ensemble, rtol=0, atol=1e-6):
        max_diff = float(np.max(np.abs(predictions["ensemble"] - saved_ensemble)))
        raise SystemExit(f"Replayed ensemble does not match saved OOF predictions (max diff {max_diff:.3g})")

    prediction_frame = pd.DataFrame(
        {"id": iv_meta["id"], "species": species, "group": groups, "y_true": y}
    )
    for name, values in predictions.items():
        for seed_index, seed in enumerate(SEEDS):
            prediction_frame[f"{name}_seed{seed}"] = values[seed_index]
    prediction_frame.to_csv(OUT_PREDICTIONS, index=False)

    unique_groups = np.unique(groups)
    members_by_group = {group: np.flatnonzero(groups == group) for group in unique_groups}
    rng = np.random.default_rng(0)
    bootstrap = {name: [] for name in names}
    for _ in range(N_BOOT):
        chosen = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        indices = np.concatenate([members_by_group[group] for group in chosen])
        for name, values in predictions.items():
            per_seed = [metrics(y[indices], values[seed_index, indices])
                        for seed_index in range(len(SEEDS))]
            bootstrap[name].append(np.mean(per_seed, axis=0))

    rows = []
    for name, values in predictions.items():
        per_seed = [metrics(y, values[seed_index]) for seed_index in range(len(SEEDS))]
        estimate = np.mean(per_seed, axis=0)
        draws = np.asarray(bootstrap[name])
        limits = np.percentile(draws, [2.5, 97.5], axis=0)
        rows.append({
            "estimator": name,
            "n": len(y),
            "groups": len(unique_groups),
            "mae_log10": estimate[0],
            "mae_ci_low": limits[0, 0],
            "mae_ci_high": limits[1, 0],
            "rmse_log10": estimate[1],
            "rmse_ci_low": limits[0, 1],
            "rmse_ci_high": limits[1, 1],
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUMMARY, index=False)
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print(f"\nvalidated ensemble OOF against {SNAPSHOT.name}; outputs: {OUT_SUMMARY.name}, {OUT_PREDICTIONS.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())