"""Two views of the source-verification sensitivity, and the size control that separates them.

A. Score on the verified subset: the saved out-of-fold predictions (trained on all 213 records),
   re-scored on the records whose labels a source supports. Asks whether the headline depends on
   labels the sources contradict or cannot support.
B. Refit on the verified subset: the model is trained only on those records. Asks whether a
   model learned from verified labels alone works. This also shrinks the training set, so
C. controls for size: the same refit on random group-wise subsets of the same number of records.

Reads MuriVerse/data/curation/halflife_label_verification.csv. Writes
verification_scored_on_subset.csv and verification_size_control.csv beside this script.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "MuriVerse"))

from evaluate_halflife_model import (  # noqa: E402
    N_SPLITS, SEEDS, THREADS, _rho, block_slices, load_matrix, union_groups,
)
from halflife_model import HalfLifeEnsemble  # noqa: E402

ANALYSIS = Path(__file__).resolve().parents[2] / "results" / "revision_analysis"
SNAPSHOT = REPO / "results/revision_analysis/snapshot/oof_predictions_iv.csv"
LEDGER = REPO / "MuriVerse/data/curation/halflife_label_verification.csv"
BLOCKS = ("molecule_annotation", "study_context", "termini", "physchem")
N_BOOT = 2000
N_DRAWS = 30


def per_seed(y, preds, idx):
    rho = np.mean([_rho(y[idx], p[idx]) for p in preds])
    mae = np.mean([np.mean(np.abs(y[idx] - p[idx])) for p in preds])
    return rho, mae


def boot(y, preds, groups, idx, seed=0):
    rng = np.random.default_rng(seed)
    members = {g: idx[groups[idx] == g] for g in np.unique(groups[idx])}
    keys = list(members)
    out = []
    for _ in range(N_BOOT):
        pick = rng.choice(len(keys), size=len(keys), replace=True)
        sub = np.concatenate([members[keys[i]] for i in pick])
        if np.ptp(y[sub]) > 0:
            out.append(per_seed(y, preds, sub))
    lo, hi = np.percentile(np.asarray(out), [2.5, 97.5], axis=0)
    return lo, hi


def refit_scores(X, y, groups, species):
    """Mean per-seed Spearman and MAE of the ensemble and a fold-local species median."""
    ens = np.full((len(SEEDS), len(y)), np.nan)
    med = np.full_like(ens, np.nan)
    for s_i, seed in enumerate(SEEDS):
        for tr, te in GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed).split(X, y, groups):
            ens[s_i, te] = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(X[tr], y[tr]).predict(X[te])
            for sp in np.unique(species[te]):
                t = te[species[te] == sp]
                pool = y[tr][species[tr] == sp]
                med[s_i, t] = np.median(pool) if len(pool) else np.median(y[tr])
    allidx = np.arange(len(y))
    return per_seed(y, ens, allidx), per_seed(y, med, allidx)


def main() -> int:
    full_x, meta, layout = load_matrix("rodent_strict", meta_dir=REPO / "MuriVerse/rodent/strict")
    iv = meta["route_class"].eq("iv").to_numpy()
    meta = meta.loc[iv].reset_index(drop=True)
    sl = block_slices(layout)
    X = np.concatenate([full_x[iv, sl[b]] for b in BLOCKS], axis=1)
    y = meta["log_label"].to_numpy(float)
    groups = union_groups(meta["cluster"].to_numpy(), meta["PMID"].astype(str).to_numpy())
    species = meta["species"].to_numpy()
    led = pd.read_csv(LEDGER).set_index("id").loc[meta["id"]]
    tier = led["evidence_tier"].to_numpy()
    corrected = led["label_verdict"].eq("corrected").to_numpy()
    saved = pd.read_csv(SNAPSHOT).set_index("id").loc[meta["id"]]
    preds = np.array([saved[f"y_pred_cluster_PMID_seed{s}"].to_numpy(float) for s in SEEDS])

    y_corr = y.copy()
    y_corr[corrected] = np.log10(led.loc[corrected, "corrected_hours"].to_numpy(float))
    subsets = {
        "all_213": (np.ones(len(y), bool), y),
        "contradicted_removed": (tier != "discrepant", y),
        "contradicted_removed_corrected": ((tier != "discrepant") | corrected, y_corr),
        "verified_body_or_abstract": (np.isin(tier, ("verified_body", "verified_abstract")), y),
        "verified_body": (tier == "verified_body", y),
        "verified_body_excluding_largest_patent": ((tier == "verified_body") & (meta["PMID"].astype(str) != "EP2277890A1").to_numpy(), y),
    }
    rows = []
    for name, (keep, target) in subsets.items():
        idx = np.flatnonzero(keep)
        rho, mae = per_seed(target, preds, idx)
        (rlo, mlo), (rhi, mhi) = boot(target, preds, groups, idx)
        rows.append({"subset": name, "n": len(idx), "groups": len(np.unique(groups[idx])),
                     "spearman": rho, "spearman_lo": rlo, "spearman_hi": rhi,
                     "mae_log10": mae, "mae_lo": mlo, "mae_hi": mhi})
    scored = pd.DataFrame(rows)
    scored.to_csv(ANALYSIS / "verification_scored_on_subset.csv", index=False)
    print("A. saved predictions scored on each subset")
    print(scored.round(3).to_string(index=False))

    rng = np.random.default_rng(11)
    controls = []
    not_patent = (meta["PMID"].astype(str) != "EP2277890A1").to_numpy()
    for target_name, pool in (("verified_body", np.ones(len(y), bool)),
                              ("verified_body_excluding_largest_patent", not_patent)):
        size = int(subsets[target_name][0].sum())
        candidates = np.flatnonzero(pool)
        for draw in range(N_DRAWS):
            keep = np.zeros(len(y), bool)
            keep[rng.choice(candidates, size=size, replace=False)] = True
            n_groups = len(np.unique(groups[keep]))
            (e_rho, e_mae), (m_rho, m_mae) = refit_scores(X[keep], y[keep], groups[keep], species[keep])
            controls.append({"matched_to": target_name, "draw": draw, "n": int(keep.sum()),
                             "groups": n_groups, "ensemble_spearman": e_rho,
                             "ensemble_mae": e_mae, "median_mae": m_mae})
    ctrl = pd.DataFrame(controls)
    ctrl.to_csv(ANALYSIS / "verification_size_control.csv", index=False)
    print("\nC. refit on random group-wise subsets of matching size")
    print(ctrl.groupby("matched_to").agg(n=("n", "mean"), draws=("draw", "size"),
          rho_mean=("ensemble_spearman", "mean"), rho_sd=("ensemble_spearman", "std"),
          rho_p05=("ensemble_spearman", lambda v: np.percentile(v, 5)),
          rho_p95=("ensemble_spearman", lambda v: np.percentile(v, 95)),
          mae_mean=("ensemble_mae", "mean"), median_mae=("median_mae", "mean")).round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
