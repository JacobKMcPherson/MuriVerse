"""Cross-fitted field-level permutation importance for the strict IV ensemble."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

import sys

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
MANIFEST = REPO / "results/revision_analysis/snapshot/feature_manifest.csv"
OUT = ANALYSIS / "feature_permutation_iv.csv"
SUPPLEMENTARY_OUT = REPO / "figures/tabs/supp_table6_feature_permutation_importance.csv"
FEATURE_BLOCKS = ("molecule_annotation", "study_context", "termini", "physchem")
N_PERMUTATIONS = 5
N_BOOT = 500


def main() -> int:
    full_x, meta, layout = load_matrix(
        "rodent_strict", meta_dir=REPO / "MuriVerse/rodent/strict"
    )
    reference = pd.read_csv(SNAPSHOT).set_index("id")
    iv_meta = meta.loc[meta["route_class"].eq("iv")].reset_index(drop=True)
    reference = reference.loc[iv_meta["id"]].reset_index()

    slices = block_slices(layout)
    missing_blocks = set(FEATURE_BLOCKS) - set(slices)
    if missing_blocks:
        raise SystemExit(f"Feature layout is missing blocks: {sorted(missing_blocks)}")
    X = np.concatenate([full_x[meta["route_class"].eq("iv").to_numpy(), slices[name]]
                        for name in FEATURE_BLOCKS], axis=1)
    y = reference["y_true"].to_numpy(float)
    groups = union_groups(iv_meta["cluster"].to_numpy(), iv_meta["PMID"].astype(str).to_numpy())

    manifest = pd.read_csv(MANIFEST).sort_values("feature_index")
    if len(manifest) != X.shape[1] or len(manifest) != 71:
        raise SystemExit(f"Expected 71 compact features, got {X.shape[1]} and {len(manifest)} manifest rows")
    feature_groups = []
    for (block, field), rows in manifest.groupby(["block", "field"], sort=False):
        indices = rows["feature_index"].to_numpy(int)
        feature_groups.append((block, field, indices))
    deltas = np.full((len(feature_groups), len(SEEDS), len(y)), np.nan)
    saved = np.stack([reference[f"y_pred_cluster_PMID_seed{s}"].to_numpy(float)
                      for s in SEEDS])
    rng = np.random.default_rng(0)

    for seed_index, seed in enumerate(SEEDS):
        splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        for train_idx, test_idx in splitter.split(X, y, groups):
            model = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(X[train_idx], y[train_idx])
            test_x = X[test_idx]
            base_pred = model.predict(test_x)
            if not np.allclose(base_pred, saved[seed_index, test_idx], rtol=0, atol=1e-6):
                raise SystemExit(f"Refit OOF predictions differ from saved seed {seed} predictions")
            base_error = np.abs(y[test_idx] - base_pred)

            for feature_index, (_, _, columns) in enumerate(feature_groups):
                repeat_deltas = []
                for _ in range(N_PERMUTATIONS):
                    order = rng.permutation(len(test_idx))
                    permuted = test_x.copy()
                    permuted[:, columns] = test_x[order][:, columns]
                    pred = model.predict(permuted)
                    repeat_deltas.append(np.abs(y[test_idx] - pred) - base_error)
                deltas[feature_index, seed_index, test_idx] = np.mean(repeat_deltas, axis=0)

    if np.isnan(deltas).any():
        raise SystemExit("At least one record has no cross-fitted permutation estimate")

    unique_groups = np.unique(groups)
    members = {group: np.flatnonzero(groups == group) for group in unique_groups}
    bootstrap_rng = np.random.default_rng(0)
    draws = np.empty((N_BOOT, len(feature_groups)))
    for draw_index in range(N_BOOT):
        chosen = bootstrap_rng.choice(unique_groups, size=len(unique_groups), replace=True)
        indices = np.concatenate([members[group] for group in chosen])
        draws[draw_index] = deltas[:, :, indices].mean(axis=(1, 2))

    rows = []
    for index, (block, field, columns) in enumerate(feature_groups):
        low, high = np.percentile(draws[:, index], [2.5, 97.5])
        rows.append({
            "block": block,
            "field": field,
            "n_columns": len(columns),
            "delta_mae_log10": float(deltas[index].mean()),
            "ci_low": float(low),
            "ci_high": float(high),
        })
    result = pd.DataFrame(rows).sort_values("delta_mae_log10", ascending=False)
    result.insert(0, "rank", np.arange(1, len(result) + 1))
    result.to_csv(OUT, index=False)
    result.to_csv(SUPPLEMENTARY_OUT, index=False)
    baseline_mae = np.mean([np.mean(np.abs(y - pred)) for pred in saved])
    print(f"{len(y)} IV records, {len(unique_groups)} connected groups, "
          f"{len(feature_groups)} permuted fields, {N_PERMUTATIONS} repeats")
    print(f"baseline mean per-seed MAE: {baseline_mae:.3f} log10 h")
    print(result.head(20).to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print(f"\nwrote {OUT}")
    print(f"wrote {SUPPLEMENTARY_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())