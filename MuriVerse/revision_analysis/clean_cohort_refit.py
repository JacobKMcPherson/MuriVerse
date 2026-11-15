"""Refit the IV model on the cohort whose labels and provenance are confirmed in the primary source.

Inclusion rule, fixed from the verification ledger before any model was scored:
  1. verified_body or verified_abstract: value, species and an intravenous route agree with the source.
  2. corrected: the source prints a different value for the same arm; the source value replaces the label.
  3. species_corrected: the source confirms the value and an intravenous route but names the other
     species; the species label and its feature columns are replaced by the source's.
Everything else is excluded: value or route not confirmed, non-intravenous route, non-elimination
phase, treated or non-plasma arm, value cited from earlier work, ambiguous arm, implausible value.
"""

from __future__ import annotations

import re
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
    _score,
    block_slices,
    load_matrix,
    union_groups,
)
from halflife_model import HalfLifeEnsemble  # noqa: E402
from metrics_iv import bootstrap_all  # noqa: E402

ANALYSIS = Path(__file__).resolve().parents[2] / "results" / "revision_analysis"
LEDGER = REPO / "MuriVerse/data/curation/halflife_label_verification.csv"
SNAPSHOT = REPO / "results/revision_analysis/snapshot/oof_predictions_iv.csv"
OUT_MEMBERSHIP = ANALYSIS / "clean_cohort_membership.csv"
OUT_METRICS = ANALYSIS / "clean_cohort_metrics.csv"
OUT_CONTROL = ANALYSIS / "clean_cohort_size_control.csv"
N_DRAWS = 30
FEATURE_BLOCKS = ("molecule_annotation", "study_context", "termini", "physchem")
LARGEST_PATENT = "EP2277890A1"
FLAGGED = re.compile(r"\s*(CITED_NOT_MEASURED|NON_BASELINE_ARM|NON_PLASMA)")
METRICS = ("spearman", "pearson", "rmse", "mae", "r2", "auc")


def classify(row: pd.Series) -> tuple[str, str]:
    """Return (rule, species_in_source) for one ledger row; rule 'excluded:*' drops the record."""
    flagged = bool(FLAGGED.match(str(row.notes or "")))
    tier = row.evidence_tier
    if tier in ("verified_body", "verified_abstract"):
        return tier, row.species
    # A correction replaces the label with the source's measured value, so flags on the old label no longer apply.
    if row.label_verdict == "corrected" and row.route_verdict in ("iv_confirmed_arm", "iv_confirmed_study_level"):
        return "corrected", row.species
    if flagged or row.route_verdict == "not_iv":
        return "excluded:arm_or_route", row.species
    source_species = re.match(r"\s*(mouse|rat)\b", str(row.species_in_source or ""))
    if (not row.species_matches and row.label_verdict == "confirmed" and source_species
            and row.route_verdict in ("iv_confirmed_arm", "iv_confirmed_study_level")):
        return "species_corrected", source_species.group(1)
    return f"excluded:{tier}", row.species


def main() -> int:
    full_x, meta, layout = load_matrix("rodent_strict", meta_dir=REPO / "MuriVerse/rodent/strict")
    meta = meta.reset_index(drop=True)
    iv_mask = meta["route_class"].eq("iv").to_numpy()
    iv_meta = meta.loc[iv_mask].reset_index(drop=True)
    slices = block_slices(layout)
    X_all = np.concatenate([full_x[iv_mask, slices[name]] for name in FEATURE_BLOCKS], axis=1)

    context = next(b for b in layout["blocks"] if b["name"] == "study_context")
    offset = sum(b["dim"] for b in layout["blocks"]
                 if b["name"] in FEATURE_BLOCKS[:1])  # molecule_annotation precedes study_context
    species_col = {name.split("=")[1]: offset + i for i, name in enumerate(context["columns"])
                   if name.startswith("species=")}

    ledger = pd.read_csv(LEDGER, dtype={"PMID": str}).set_index("id").loc[iv_meta["id"]]
    ledger["species_matches"] = ledger["species_matches"].astype(str).eq("True")
    verdicts = [classify(row) for row in ledger.itertuples()]
    ledger["rule"] = [v[0] for v in verdicts]
    ledger["species_used"] = [v[1] for v in verdicts]
    include = ~ledger["rule"].str.startswith("excluded").to_numpy()

    y_all = iv_meta["log_label"].to_numpy(float).copy()
    corrected = (ledger["rule"] == "corrected").to_numpy()
    y_all_clean = y_all.copy()
    y_all_clean[corrected] = np.log10(ledger.loc[corrected, "corrected_hours"].to_numpy(float))

    X_clean = X_all.copy()
    species_clean = iv_meta["species"].to_numpy(object).copy()
    for position in np.flatnonzero((ledger["rule"] == "species_corrected").to_numpy()):
        new = ledger["species_used"].iloc[position]
        for animal, column in species_col.items():
            X_clean[position, column] = float(animal == new)
        species_clean[position] = new

    membership = pd.DataFrame({
        "id": iv_meta["id"], "PMID": iv_meta["PMID"].astype(str),
        "included": include, "rule": ledger["rule"].to_numpy(),
        "label_hours_original": 10 ** iv_meta["log_label"].to_numpy(float),
        "label_hours_used": np.where(include, 10 ** y_all_clean, np.nan),
        "species_original": iv_meta["species"].to_numpy(), "species_used": species_clean,
    })
    membership.to_csv(OUT_MEMBERSHIP, index=False)
    print(membership["rule"].value_counts().to_string(), "\n")
    print(f"clean cohort: {include.sum()} of {len(include)} IV records, "
          f"{iv_meta.loc[include, 'PMID'].nunique()} sources\n")

    groups_all = union_groups(iv_meta["cluster"].to_numpy(), iv_meta["PMID"].astype(str).to_numpy())
    patent = (iv_meta["PMID"].astype(str) == LARGEST_PATENT).to_numpy()
    scenarios = {
        "original_all_iv": (np.ones(len(include), bool), X_all, y_all, iv_meta["species"].to_numpy()),
        "clean_cohort": (include, X_clean, y_all_clean, species_clean),
        "clean_cohort_excluding_largest_patent": (include & ~patent, X_clean, y_all_clean, species_clean),
    }
    saved = pd.read_csv(SNAPSHOT).set_index("id")

    rows = []
    for scenario, (keep, X_full, y_full, species_full) in scenarios.items():
        X, y, species, groups = X_full[keep], y_full[keep], species_full[keep], groups_all[keep]
        ids = iv_meta["id"].to_numpy()[keep]
        names = ("ridge", "svr", "xgboost", "ensemble", "species_median")
        predictions = {name: np.full((len(SEEDS), len(y)), np.nan) for name in names}
        fold_assignments = np.full((len(SEEDS), len(y)), -1, dtype=int)
        for seed_index, seed in enumerate(SEEDS):
            splitter = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
            for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y, groups)):
                fold_assignments[seed_index, test_idx] = fold
                fitted = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(X[train_idx], y[train_idx])
                members = fitted.member_predictions(X[test_idx])
                for name in ("ridge", "svr", "xgboost"):
                    predictions[name][seed_index, test_idx] = members[name]
                predictions["ensemble"][seed_index, test_idx] = fitted.predict(X[test_idx])
                for animal in np.unique(species[test_idx]):
                    in_test = test_idx[species[test_idx] == animal]
                    same = y[train_idx][species[train_idx] == animal]
                    predictions["species_median"][seed_index, in_test] = (
                        float(np.median(same)) if len(same) else float(np.median(y[train_idx])))
        if any(np.isnan(v).any() for v in predictions.values()):
            raise SystemExit(f"{scenario}: rows without out-of-fold predictions")
        if scenario == "clean_cohort":
            output = iv_meta.loc[keep, ["id", "PMID", "cluster", "route_class"]].reset_index(drop=True)
            output["species"] = species
            output["y_true"] = y
            output["group"] = groups.astype(str)
            for seed_index, seed in enumerate(SEEDS):
                output[f"fold_seed{seed}"] = fold_assignments[seed_index]
                for name, values in predictions.items():
                    output[f"{name}_seed{seed}"] = values[seed_index]
            output.to_csv(ANALYSIS / "clean_cohort_oof_predictions.csv", index=False)
        if scenario == "original_all_iv":
            reference = saved.loc[ids][[f"y_pred_cluster_PMID_seed{s}" for s in SEEDS]].to_numpy(float).T
            if not np.allclose(predictions["ensemble"], reference, rtol=0, atol=1e-6):
                raise SystemExit("Original-cohort refit does not reproduce the saved ensemble predictions")

        n_groups = len(np.unique(groups))
        for name, values in predictions.items():
            per_seed = [_score(y, values[i]) for i in range(len(SEEDS))]
            interval = bootstrap_all(y, values, groups)
            for metric in METRICS:
                rows.append({
                    "scenario": scenario, "estimator": name, "metric": metric,
                    "value": float(np.mean([s[metric] for s in per_seed])),
                    "ci_low": interval[metric][0], "ci_high": interval[metric][1],
                    "n": len(y), "groups": n_groups,
                })
        print(f"{scenario}: n = {len(y)}, groups = {n_groups}")

    # Predictions from the model trained on all 213 records, scored on the clean records only.
    saved_ids = saved.loc[iv_meta["id"]]
    saved_predictions = np.array([saved_ids[f"y_pred_cluster_PMID_seed{s}"].to_numpy(float) for s in SEEDS])
    kept = np.flatnonzero(include)
    per_seed = [_score(y_all_clean[kept], p[kept]) for p in saved_predictions]
    interval = bootstrap_all(y_all_clean[kept], saved_predictions[:, kept], groups_all[kept])
    for metric in METRICS:
        rows.append({"scenario": "clean_cohort_scored_with_original_model", "estimator": "ensemble",
                     "metric": metric, "value": float(np.mean([s[metric] for s in per_seed])),
                     "ci_low": interval[metric][0], "ci_high": interval[metric][1],
                     "n": len(kept), "groups": len(np.unique(groups_all[kept]))})

    # Size control: refit on random record subsets as large as the clean cohort, original labels.
    rng = np.random.default_rng(11)
    controls = []
    for draw in range(N_DRAWS):
        pick = np.sort(rng.choice(len(y_all), size=len(kept), replace=False))
        Xs, ys, gs = X_all[pick], y_all[pick], groups_all[pick]
        ens = np.full((len(SEEDS), len(pick)), np.nan)
        for seed_index, seed in enumerate(SEEDS):
            for train_idx, test_idx in GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed).split(Xs, ys, gs):
                ens[seed_index, test_idx] = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(
                    Xs[train_idx], ys[train_idx]).predict(Xs[test_idx])
        scores = [_score(ys, ens[i]) for i in range(len(SEEDS))]
        controls.append({"draw": draw, "n": len(pick), "groups": len(np.unique(gs)),
                         "spearman": np.mean([s["spearman"] for s in scores]),
                         "mae": np.mean([s["mae"] for s in scores])})
    control = pd.DataFrame(controls)
    control.to_csv(OUT_CONTROL, index=False)
    print(f"size control ({N_DRAWS} random draws of {len(kept)}): Spearman mean {control.spearman.mean():.3f}, "
          f"5th-95th percentile {control.spearman.quantile(.05):.3f}-{control.spearman.quantile(.95):.3f}; "
          f"MAE mean {control.mae.mean():.3f}")

    result = pd.DataFrame(rows)
    result.to_csv(OUT_METRICS, index=False)
    wide = result[result.estimator == "ensemble"].pivot(index="scenario", columns="metric", values="value")
    print(wide[list(METRICS)].to_string(float_format=lambda v: f"{v:.3f}"))
    print(f"\nwrote {OUT_MEMBERSHIP.name} and {OUT_METRICS.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
