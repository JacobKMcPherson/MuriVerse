"""Fit the post-audit models on the source-verified intravenous cohort and save them beside the pre-audit models.

The pre-audit models (MuriVerse/rodent/strict/model_iv/*.joblib, trained on all 213 intravenous-labelled
records) are left untouched as the reproducibility artifact. The inclusion rule, label corrections and species
corrections are those of clean_cohort_refit.py. All three configurations the polymyxin application uses are refit
with the same fixed hyperparameters: shipped (71 columns), structural (physicochemical descriptors + study context
+ terminal classes) and structural_embedding (PeptideCLM SMILES embedding + descriptors + study context). Each is
also scored by the grouped five-fold protocol (three seeds, cluster and source groups) on the 178 records.

Run: python MuriVerse/revision_analysis/refit_postaudit_model.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "MuriVerse"))
sys.path.insert(0, str(Path(__file__).parent))

from clean_cohort_refit import LEDGER, classify  # noqa: E402
from evaluate_halflife_model import (  # noqa: E402
    N_SPLITS,
    SEEDS,
    THREADS,
    block_slices,
    load_matrix,
    union_groups,
)
from halflife_model import HalfLifeEnsemble  # noqa: E402

STRICT = REPO / "MuriVerse/rodent/strict"
OUT = STRICT / "model_iv_postaudit"
PRE = STRICT / "model_iv"


def main() -> int:
    full_x, meta, layout = load_matrix("rodent_strict", meta_dir=STRICT)
    iv = meta["route_class"].eq("iv").to_numpy()
    iv_meta = meta.loc[iv].reset_index(drop=True)
    slices = block_slices(layout)
    Xfull = full_x[iv].copy()

    context = next(b for b in layout["blocks"] if b["name"] == "study_context")
    species_col = {n.split("=")[1]: slices["study_context"].start + i
                   for i, n in enumerate(context["columns"]) if n.startswith("species=")}

    ledger = pd.read_csv(LEDGER, dtype={"PMID": str}).set_index("id").loc[iv_meta["id"]]
    ledger["species_matches"] = ledger["species_matches"].astype(str).eq("True")
    verdicts = [classify(r) for r in ledger.itertuples()]
    rule = np.array([v[0] for v in verdicts])
    include = ~pd.Series(rule).str.startswith("excluded").to_numpy()

    y = iv_meta["log_label"].to_numpy(float).copy()
    fixed = rule == "corrected"
    y[fixed] = np.log10(ledger.loc[fixed, "corrected_hours"].to_numpy(float))
    for pos in np.flatnonzero(rule == "species_corrected"):
        for animal, col in species_col.items():
            Xfull[pos, col] = float(animal == verdicts[pos][1])

    pre_layout = json.loads((PRE / "model.json").read_text())
    groups = union_groups(iv_meta["cluster"].to_numpy(), iv_meta["PMID"].astype(str).to_numpy())
    OUT.mkdir(exist_ok=True)
    models = {}
    for role, spec in pre_layout["models"].items():
        X = np.concatenate([Xfull[:, slices[b]] for b in spec["blocks"]], axis=1)
        if X.shape[1] != spec["dims"]:
            raise SystemExit(f"{role}: {X.shape[1]} columns, expected {spec['dims']}")
        Xk, yk, gk = X[include], y[include], groups[include]
        model = HalfLifeEnsemble(seed=SEEDS[0], nthread=THREADS).fit(Xk, yk)
        joblib.dump(model, OUT / f"{role}.joblib")
        rhos = []
        for seed in SEEDS:
            oof = np.full(len(yk), np.nan)
            for tr, te in GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed).split(Xk, yk, gk):
                oof[te] = HalfLifeEnsemble(seed=seed, nthread=THREADS).fit(Xk[tr], yk[tr]).predict(Xk[te])
            rhos.append(float(spearmanr(yk, oof)[0]))
        models[role] = {**{k: spec[k] for k in ("features", "blocks", "dims")},
                        "cluster_pmid_spearman": float(np.mean(rhos)), "n_train": int(include.sum())}
        print(f"{role}: {include.sum()} records, {spec['dims']} columns, grouped-CV Spearman {np.mean(rhos):.3f}")

    pd.DataFrame({"id": iv_meta["id"], "included": include, "rule": rule,
                  "label_hours_used": np.where(include, 10 ** y, np.nan),
                  "species_used": [v[1] for v in verdicts]}).to_csv(
        OUT / "training_membership.csv", index=False)
    layout_out = {k: v for k, v in pre_layout.items() if k not in ("models", "tier_rows", "training_composition")}
    species_used = np.array([v[1] for v in verdicts])[include]
    layout_out.update({
        "status": "post-audit models, trained on source-verified intravenous records",
        "pre_audit_models": "../model_iv/ (213 intravenous-labelled records, retained for reproducibility)",
        "tier_rows": int(include.sum()),
        "n_sources": int(iv_meta.loc[include, "PMID"].nunique()),
        "training_composition": {"species": {s: int((species_used == s).sum()) for s in sorted(set(species_used))},
                                 "route": {"iv": int(include.sum())}},
        "rules": {k: int((rule == k).sum()) for k in sorted(set(rule))},
        "seed": SEEDS[0],
        "models": models,
        "cross_validated_performance": "grouped five-fold, three seeds, fixed configuration (see cluster_pmid_spearman); "
                                       "clean_cohort_metrics.csv for the shipped configuration",
    })
    (OUT / "model.json").write_text(json.dumps(layout_out, indent=2) + "\n")
    print(pd.Series(rule).value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
