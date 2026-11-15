"""Sensitivity: the fixed-configuration ensemble on the intravenous cohort without the records of six residues or fewer.

Run: python MuriVerse/revision_analysis/short_peptide_sensitivity.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from sklearn.model_selection import GroupKFold

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "MuriVerse"))
sys.path.insert(0, str(Path(__file__).parent))

from clean_cohort_refit import FEATURE_BLOCKS  # noqa: E402
from evaluate_halflife_model import N_SPLITS, SEEDS, THREADS, _score, block_slices, load_matrix, union_groups  # noqa: E402
from halflife_model import HalfLifeEnsemble  # noqa: E402
from metrics_iv import bootstrap_all  # noqa: E402


def main() -> int:
    full_x, meta, layout = load_matrix("rodent_strict", meta_dir=REPO / "MuriVerse/rodent/strict")
    iv = meta["route_class"].eq("iv").to_numpy()
    sl = block_slices(layout)
    X = np.concatenate([full_x[iv, sl[b]] for b in FEATURE_BLOCKS], axis=1)
    m = meta.loc[iv].reset_index(drop=True)
    y = m["log_label"].to_numpy(float)
    groups = union_groups(m["cluster"].to_numpy(), m["PMID"].astype(str).to_numpy())
    keep = (m["n_residues"].astype(float) > 6).to_numpy()
    for name, mask in (("all 213", np.ones(len(y), bool)), (f"excluding {int((~keep).sum())} short peptides", keep)):
        Xs, ys, gs = X[mask], y[mask], groups[mask]
        pred = np.full((len(SEEDS), len(ys)), np.nan)
        for i, seed in enumerate(SEEDS):
            for tr, te in GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed).split(Xs, ys, gs):
                pred[i, te] = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(Xs[tr], ys[tr]).predict(Xs[te])
        rho = np.mean([_score(ys, p)["spearman"] for p in pred])
        lo, hi = bootstrap_all(ys, pred, gs)["spearman"]
        print(f"{name}: n = {len(ys)}, groups = {len(set(gs))}, Spearman {rho:.3f} ({lo:.3f}-{hi:.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
