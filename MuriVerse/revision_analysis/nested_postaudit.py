"""Nested feature-configuration selection on the post-audit (source-verified) intravenous cohort.

Same candidate set, inner folds and outer protocol as evaluate_halflife_model.nested_confirmation; only the
records, the three corrected labels and the two corrected species differ (rules of clean_cohort_refit.py).
Writes figures/tabs/table31b_muriverse_halflife_nested_selection_postaudit.csv.

Run: python MuriVerse/revision_analysis/nested_postaudit.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "MuriVerse"))
sys.path.insert(0, str(Path(__file__).parent))

import evaluate_halflife_model as ehm  # noqa: E402
from clean_cohort_refit import LEDGER, classify  # noqa: E402


def main() -> int:
    full_x, meta, layout = ehm.load_matrix("rodent_strict", meta_dir=REPO / "MuriVerse/rodent/strict")
    iv = meta["route_class"].eq("iv").to_numpy()
    X = full_x[iv].copy()
    iv_meta = meta.loc[iv].reset_index(drop=True)
    sl = ehm.block_slices(layout)

    context = next(b for b in layout["blocks"] if b["name"] == "study_context")
    species_col = {n.split("=")[1]: context["start"] + i for i, n in enumerate(context["columns"])
                   if n.startswith("species=")}
    ledger = pd.read_csv(LEDGER, dtype={"PMID": str}).set_index("id").loc[iv_meta["id"]]
    ledger["species_matches"] = ledger["species_matches"].astype(str).eq("True")
    verdicts = [classify(r) for r in ledger.itertuples()]
    rule = np.array([v[0] for v in verdicts])
    include = ~pd.Series(rule).str.startswith("excluded").to_numpy()

    y = iv_meta["log_label"].to_numpy(float).copy()
    fixed = rule == "corrected"
    y[fixed] = np.log10(ledger.loc[fixed, "corrected_hours"].to_numpy(float))
    species = iv_meta["species"].to_numpy(object).copy()
    for pos in np.flatnonzero(rule == "species_corrected"):
        for animal, col in species_col.items():
            X[pos, col] = float(animal == verdicts[pos][1])
        species[pos] = verdicts[pos][1]

    groups = ehm.union_groups(iv_meta["cluster"].to_numpy(), iv_meta["PMID"].astype(str).to_numpy())
    # The eight non-embedding, non-diagnostic candidates of at most 100 columns in the original analysis.
    candidates = {
        "termini only": ["termini"],
        "molecule annotation + study context": ["molecule_annotation", "study_context"],
        "study context + termini": ["study_context", "termini"],
        "physchem only": ["physchem"],
        "physchem + study context": ["physchem", "study_context"],
        "physchem + molecule annotation": ["physchem", "molecule_annotation"],
        "physchem + study context + termini": ["physchem", "study_context", "termini"],
        "annotation + context + termini + physchem": ["molecule_annotation", "study_context", "termini", "physchem"],
    }
    configurations = {n: [sl[b] for b in parts] for n, parts in candidates.items()}

    ehm.TAG = "_postaudit"
    ehm.nested_confirmation(X[include], y[include], species[include], configurations, groups[include])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
