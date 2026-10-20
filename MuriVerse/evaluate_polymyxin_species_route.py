"""Does the cross-route finding hold on polymyxin B and colistin?

Figure 3 established on the training cohort that the two administration routes do not
share an absolute half-life scale (cross-route transfer is offset 0.45-0.65 log10 with the
sign flip-flop kinetics predicts). That is an internal, held-out claim. The polymyxin
literature is an external test of it, because it reports polymyxin B and colistin in BOTH
animals by BOTH routes - one of the few places in this project where the four cells are
observed for the same molecules.

The test is not "is the prediction close", because an absolute prediction for these
molecules is expected to be poor: they sit far outside the training chemistry. The scatter
of observed-against-predicted across the four species x route cells is what gets drawn in
Figure 4c.

Routes are grouped IV against EV (extravascular) as in Figure 3. The model's route levels
are iv and sc, so EV is predicted as sc.

Writes figures/tabs/table56_polymyxin_species_route_cells.csv  (observed vs predicted per cell)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import predict_polymyxins as pp  # noqa: E402

# THIS SCRIPT DELIBERATELY USES THE ROUTE-POOLED MODEL, not the shipped intravenous one.
#
# Its job is to test the route claim externally: does a model that has seen both routes
# mis-transfer across them by the 0.45-0.65 log10 that Figure 3's panels a and b measure
# internally? That question cannot be asked of a model fitted on intravenous records alone.
# Such a model has no subcutaneous level to predict under - route is constant in its training
# matrix - so an EV prediction from it is an extrapolation along an axis with zero variance,
# and the route effect it appeared to report would be an artefact of the scaler.
#
# So the instrument for this diagnostic is the model whose behaviour is being diagnosed. The
# shipped model is the route-restricted one and it is what `predict_polymyxins.py` applies;
# this script is the measurement that justifies restricting it, in the same relation as
# Figure 3's panel a to Figure 3's panel c onwards.
pp.MODEL_TAG = ""
pp.MODEL_DIR = pp.STRICT / "model"
pp.ROUTE_VARIANTS = ("iv", "sc")

TABS = HERE.parent / "figures/tabs"
OUT_CELLS = TABS / "table56_polymyxin_species_route_cells.csv"

DRUGS = {"Polymyxin B": "Polymyxin B", "Colistin (polymyxin E)": "Colistin"}


def observed_cells() -> pd.DataFrame:
    """Observed half-lives collapsed to drug x species x route, one value per study first.

    Averaging within a study before averaging across studies stops a study that reported many
    dose arms from dominating the cell; the cells are small enough for that to matter.

    Only the studies flagged `in_figure4` in data/polymyxin_studies.csv are used, the same
    eighteen as the intravenous Figure 4 (panels a, b, d), Supplementary Figure 3 and Tables 2-5. Without that filter the cells also held
    studies the comparison deliberately excludes (prodrug arms, a renal-function paper,
    the intratracheal arm of a colistin study), so panel (c) was drawn from a different set of observations
    than the rest of the figure.
    """
    studies = pd.read_csv(pp.STUDIES_CSV, dtype=str, keep_default_na=False)
    chosen = set(studies.loc[studies["in_figure4"] == "True", "study_id"])
    obs = pd.read_csv(HERE / "data" / "polymyxin_observed_half_lives.csv")
    obs = obs[obs["study_id"].isin(chosen)]
    obs = obs[obs["use_in_band"] == True]  # noqa: E712
    obs = obs[obs["drug"].isin(DRUGS)].copy()
    obs["drug_short"] = obs["drug"].map(DRUGS)
    obs["admin"] = np.where(obs["route"].astype(str).str.upper().str.startswith("IV"),
                            "IV", "EV")
    obs["log10_h"] = np.log10(obs["half_life_h"])
    per_study = (obs.groupby(["drug_short", "species", "admin", "study_id"], as_index=False)
                 .agg(log10_h=("log10_h", "mean")))
    return (per_study.groupby(["drug_short", "species", "admin"], as_index=False)
            .agg(obs_log10=("log10_h", "mean"),
                 obs_sd=("log10_h", "std"),
                 n_studies=("study_id", "nunique")))


def predict_grid(models, layout, mol_cols, ctx_cols, term_cols) -> pd.DataFrame:
    """Shipped-model prediction for every congener x species x route combination.

    The canonical terminal variant is used - the literal N-terminal acyl and a free
    C-terminus - rather than averaging the interpretations, because the question is how
    the prediction MOVES with study design, and every variant moves identically. Strain
    is left unknown, which is how the shipped model is actually applied.
    """
    rows = []
    blocks = layout["models"]["shipped"]["blocks"]
    physchem_cols = pp.block_columns("physchem")
    for congener in pp.CONGENERS:
        drug = congener["drug"]
        if drug not in DRUGS:
            continue
        phys, _ = pp.physchem_vector(congener["smiles"], physchem_cols)
        term = pp.termini_vector(congener["acyl"], "Free", term_cols)
        for admin, route in (("IV", "iv"), ("EV", "sc")):
            for species in ("mouse", "rat"):
                fields = {"hl_engineering": "lipidated", "route_class": route,
                          "species": species, "topology": "cyclic", "chirality": "mix",
                          "smiles_coverage": "complete"}
                mol, ctx, _ = pp.annotation_vector(fields, mol_cols, ctx_cols)
                X = pp.assemble(blocks, mol, ctx, None, term, layout, phys)
                pred = float(models["shipped"].predict(X)[0])
                rows.append({"drug_short": DRUGS[drug], "congener": congener["name"],
                             "species": species, "admin": admin, "pred_log10": pred})
    return pd.DataFrame(rows)


def main() -> None:
    models, layout = pp.load_models()
    if "shipped" not in models:
        raise SystemExit("no shipped model saved; run evaluate_halflife_model.py first")
    mol_cols, ctx_cols = pp.annotation_columns()
    term_cols = pp.block_columns("termini")

    grid = predict_grid(models, layout, mol_cols, ctx_cols, term_cols)
    obs = observed_cells()

    # Observed vs predicted per cell. The shipped configuration as actually applied leaves
    # strain unknown, so that is the prediction compared against observation.
    pred_cells = (grid.groupby(["drug_short", "species", "admin"], as_index=False)
                  .agg(pred_log10=("pred_log10", "mean"),
                       pred_lo=("pred_log10", "min"), pred_hi=("pred_log10", "max")))
    cells = obs.merge(pred_cells, on=["drug_short", "species", "admin"], how="left")
    cells["error_log10"] = cells["pred_log10"] - cells["obs_log10"]
    cells["fold_error"] = 10 ** cells["error_log10"].abs()
    cells.to_csv(OUT_CELLS, index=False)
    print("Observed vs predicted per cell (log10 h):")
    for _, r in cells.iterrows():
        print(f"  {r.drug_short:12} {r.species:5} {r.admin}: obs {r.obs_log10:+.2f} "
              f"pred {r.pred_log10:+.2f}  error {r.error_log10:+.2f} "
              f"({r.fold_error:.1f}x, {r.n_studies} studies)")

    print(f"[Wrote] {OUT_CELLS}")


if __name__ == "__main__":
    main()
