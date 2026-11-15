"""Refit the IV model after label-certainty and route-evidence restrictions."""

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
OUT = ANALYSIS / "cohort_sensitivity_metrics.csv"
LEDGER = REPO / "MuriVerse/data/curation/halflife_label_verification.csv"
N_BOOT = 500
FEATURE_BLOCKS = ("molecule_annotation", "study_context", "termini", "physchem")


def score(y: np.ndarray, pred: np.ndarray) -> tuple[float, float, float]:
    error = y - pred
    rho = float(pd.Series(y).corr(pd.Series(pred), method="spearman"))
    return rho, float(np.mean(np.abs(error))), float(np.sqrt(np.mean(error**2)))


def main() -> int:
    full_x, meta, layout = load_matrix(
        "rodent_strict", meta_dir=REPO / "MuriVerse/rodent/strict"
    )
    meta = meta.reset_index(drop=True)
    iv_mask = meta["route_class"].eq("iv").to_numpy()
    iv_meta = meta.loc[iv_mask].reset_index(drop=True)
    slices = block_slices(layout)
    X_all = np.concatenate([full_x[iv_mask, slices[name]] for name in FEATURE_BLOCKS], axis=1)
    y_all = iv_meta["log_label"].to_numpy(float)
    groups_all = union_groups(iv_meta["cluster"].to_numpy(), iv_meta["PMID"].astype(str).to_numpy())

    exact_label = (iv_meta["half_life_qualifier"].eq("exact")
                   & ~iv_meta["half_life_provenance"].str.contains("ambiguous", na=False))
    ledger = pd.read_csv(LEDGER).set_index("id").loc[iv_meta["id"]]
    tier = ledger["evidence_tier"].to_numpy()
    corrected = ledger["label_verdict"].eq("corrected").to_numpy()
    scenarios = {
        "all_iv": np.ones(len(iv_meta), dtype=bool),
        "exact_unambiguous_labels": exact_label.to_numpy(),
        "source_contradicted_removed": tier != "discrepant",
        "source_corrected": (tier != "discrepant") | corrected,
        "verified_body_or_abstract": np.isin(tier, ("verified_body", "verified_abstract")),
        "verified_body": tier == "verified_body",
    }
    scenarios["verified_body_excluding_largest_patent"] = (
        scenarios["verified_body"] & (iv_meta["PMID"].astype(str) != "EP2277890A1").to_numpy())
    y_corrected = y_all.copy()
    y_corrected[corrected] = np.log10(ledger.loc[iv_meta["id"][corrected], "corrected_hours"].to_numpy(float))

    saved = pd.read_csv(SNAPSHOT).set_index("id")
    all_predictions = {}
    summary_rows = []
    for scenario, keep in scenarios.items():
        X, y = X_all[keep], (y_corrected if scenario == "source_corrected" else y_all)[keep]
        meta_subset = iv_meta.loc[keep].reset_index(drop=True)
        groups = groups_all[keep]
        species = meta_subset["species"].to_numpy()
        group_count = len(np.unique(groups))
        if group_count < N_SPLITS:
            raise SystemExit(f"{scenario} has only {group_count} groups for {N_SPLITS}-fold CV")

        estimators = {name: np.full((len(SEEDS), len(y)), np.nan)
                      for name in ("ensemble", "species_median")}
        for seed_index, seed in enumerate(SEEDS):
            splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
            for train_idx, test_idx in splitter.split(X, y, groups):
                fitted = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(X[train_idx], y[train_idx])
                estimators["ensemble"][seed_index, test_idx] = fitted.predict(X[test_idx])
                for animal in np.unique(species[test_idx]):
                    test_species = test_idx[species[test_idx] == animal]
                    train_species = y[train_idx][species[train_idx] == animal]
                    estimators["species_median"][seed_index, test_species] = (
                        float(np.median(train_species)) if len(train_species)
                        else float(np.median(y[train_idx]))
                    )

        if any(np.isnan(values).any() for values in estimators.values()):
            raise SystemExit(f"{scenario} has rows without out-of-fold predictions")
        if scenario == "all_iv":
            reference = saved.loc[meta_subset["id"]]
            saved_predictions = reference[
                [f"y_pred_cluster_PMID_seed{seed}" for seed in SEEDS]
            ].to_numpy(float).T
            if not np.allclose(estimators["ensemble"], saved_predictions, rtol=0, atol=1e-6):
                raise SystemExit("Full-cohort refit does not reproduce the saved ensemble OOF predictions")

        all_predictions[scenario] = (meta_subset, y, estimators)
        members_by_group = {group: np.flatnonzero(groups == group) for group in np.unique(groups)}
        group_keys = np.array(list(members_by_group), dtype=object)
        rng = np.random.default_rng(0)
        draws = {name: [] for name in estimators}
        for _ in range(N_BOOT):
            chosen = rng.choice(group_keys, size=len(group_keys), replace=True)
            indices = np.concatenate([members_by_group[group] for group in chosen])
            for name, predictions in estimators.items():
                per_seed = np.array([score(y[indices], predictions[i, indices])
                                     for i in range(len(SEEDS))])
                draws[name].append(per_seed.mean(axis=0))

        for name, predictions in estimators.items():
            per_seed = np.array([score(y, predictions[i]) for i in range(len(SEEDS))])
            estimate = per_seed.mean(axis=0)
            limits = np.percentile(np.asarray(draws[name]), [2.5, 97.5], axis=0)
            summary_rows.append({
                "scenario": scenario,
                "estimator": name,
                "n": len(y),
                "groups": group_count,
                "spearman": estimate[0],
                "spearman_ci_low": limits[0, 0],
                "spearman_ci_high": limits[1, 0],
                "mae_log10": estimate[1],
                "mae_ci_low": limits[0, 1],
                "mae_ci_high": limits[1, 1],
                "rmse_log10": estimate[2],
                "rmse_ci_low": limits[0, 2],
                "rmse_ci_high": limits[1, 2],
            })

    pd.DataFrame(summary_rows).to_csv(OUT, index=False)
    print(pd.DataFrame(summary_rows).to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print("\nLabel restriction: qualifier is exact and provenance is not an ambiguous full-text match.")
    print(f"\nvalidated all_iv ensemble predictions against {SNAPSHOT.name}; wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())