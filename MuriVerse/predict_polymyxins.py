"""Apply the trained MuriVerse half-life model to colistin (polymyxin E) and polymyxin B.

This is a **held-out application**, not a validation run. Nothing is retrained here. The
models in `rodent/strict/model/` are loaded exactly as `evaluate_halflife_model.py` fitted
them - on all 470 strict rows, with fixed hyperparameters - and asked for a prediction on
molecules that were never in the cohort.

Neither drug is in the strict cohort. `colistin` appears once in the eligible rodent set, in the
`curated_other_agents` free-text column of the ZY4 paper (PMID 31843919), as a comparator
whose dosing was transcribed - never as a modelled row with a half-life label. `polymyxin`
appears nowhere. Both are asserted at runtime below rather than trusted.

TWO saved configurations are applied:

* **shipped** - the ensemble on the 71 compact columns (molecule annotation, study context,
  termini, and the 16 RDKit physicochemical descriptors read off the SMILES). Because the
  descriptors are computed from the structure alone they are fully defined for a molecule
  with no faithful single-letter sequence, and the four congeners get four different vectors.
  This is the model.
* **structural_embedding** - the best configuration that reads the PeptideCLM embedding of the
  SMILES. It is applied for completeness; on molecules this far outside the training chemistry
  it extrapolates to ~1e-8 h, so it is reported in a table and not drawn in Figure 4.

Species is a model input, carried by the `species=mouse` / `species=rat` columns of the
study_context block, and it is varied here rather than assumed: every congener is predicted
once as a mouse and once as a rat, under otherwise identical conditions, so the two numbers
differ by whatever the model learned about the animal. Strain is still left unknown, because
it is a vendor detail the application does not know; it is no longer the only way the animal
reaches the model. Dose is left missing.

Why these molecules are hard, and why that is the point:

* Six of their ten residues are L-2,4-diaminobutyric acid (Dab), which has no single-letter
  code. This no longer affects any prediction — neither saved configuration uses the ESM-2
  sequence block — but it still affects the homology diagnostic, which is why the four
  surrogate spellings are retained there and reported as a spread.
* Residue 6 is D-Phe (polymyxin B) or D-Leu (colistin), so `Chirality` is `Mix`.
* The acyl tail and the Dab4-Thr10 macrolactam are real chemistry that only the structural
  configuration can see at all — and under the physicochemical block it sees them as
  quantities (mass, lipophilicity, ring count, hydrogen-bond donors and the cationic charge
  the five free Dab amines carry) rather than as a position in a learned embedding space. Both are present in the PubChem isomeric SMILES used here,
  which is a different provenance from the RDKit reconstructions the strict set is built from
  and is labelled as such.

Congener-level, not drug-level: both drugs are marketed as mixtures. Polymyxin B is dominated
by B1 and B2; colistin by colistin A (polymyxin E1) and colistin B (polymyxin E2). Each is
predicted separately, from its own structure.

  python MuriVerse/predict_polymyxins.py
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import argparse
import joblib
import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(REPO))

import peptide_physchem  # noqa: E402
from build_halflife_training_set import (  # noqa: E402
    bucket_mod,
    describe_termini,
    embed_peptideclm,
)

STRICT = ROOT / "rodent" / "strict"

# The shipped model is fitted on the INTRAVENOUS cohort only (both species). Figure 3's
# panels a and b are the measurement behind that: every cross-route transfer is offset by
# 0.45-0.65 log10 (a factor of 3-4) with the sign flip-flop kinetics predicts, while every
# cross-species transfer stays inside +/-0.21. A route-pooled target is therefore a mixture
# of two different physical quantities, and the route-restricted model is better on every
# axis over the same held-out IV rows, compared like for like at the shipped configuration
# (rho 0.507 -> 0.574, AUC 0.785 -> 0.818), though that rho gain alone is inside its own
# paired bootstrap interval (+0.063 [-0.005, +0.157]); the bias result is the load-bearing one.
#
# MODEL_TAG selects which saved production models to apply. "_iv" is the shipped set, written
# by `evaluate_halflife_model.py --route iv`; "" is the older route-pooled set, which is left
# on disk for comparison. Override with --model-suffix.
MODEL_TAG = "_iv"
MODEL_DIR = STRICT / f"model{MODEL_TAG}"
DATA = REPO / "results" / "half_life_muriverse"
ABLATION = DATA / "halflife_rodent_strict_ablation_with_embeddings"
TABS = REPO / "figures" / "tabs"

ROLES = ("shipped", "structural", "structural_embedding")

# ---------------------------------------------------------------------------------------
# The molecules. Residue sequences and acyl tails are from the PubChem records cited; the
# curation-style fields are written the way a MuriVerse curator would fill them in, so that
# the pipeline's own bucket_mod decides the features rather than this script.
# ---------------------------------------------------------------------------------------
DAB = "X"  # placeholder for L-2,4-diaminobutyric acid in the residue template

CONGENERS = [
    {
        "name": "Polymyxin B1",
        "drug": "Polymyxin B",
        "pubchem_cid": 49800004,
        "formula": "C56H98N16O13",
        "mw": 1203.5,
        "inchikey": "WQVJHHACXVLGBL-BPJDFBQWSA-N",
        "residues": f"{DAB}T{DAB}{DAB}{DAB}FL{DAB}{DAB}T",
        "d_residue": "D-Phe6",
        "acyl": "(S)-6-methyloctanoyl",
        "smiles": "CCC(C)CCCCC(=O)N[C@@H](CCN)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCN)C(=O)"
                  "N[C@H]1CCNC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H]"
                  "(NC(=O)[C@H](NC1=O)CCN)CC2=CC=CC=C2)CC(C)C)CCN)CCN)[C@@H](C)O",
        "share": "~70-75% of polymyxin B drug substance",
    },
    {
        "name": "Polymyxin B2",
        "drug": "Polymyxin B",
        "pubchem_cid": 11982455,
        "formula": "C55H96N16O13",
        "mw": 1189.5,
        "inchikey": "SGPYLFWAQBAXCZ-RUDZPDEXSA-N",
        "residues": f"{DAB}T{DAB}{DAB}{DAB}FL{DAB}{DAB}T",
        "d_residue": "D-Phe6",
        "acyl": "6-methylheptanoyl",
        "smiles": "C[C@H]([C@H]1C(=O)NCC[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N"
                  "[C@H](C(=O)N[C@H](C(=O)N1)CCN)CCN)CC(C)C)CC2=CC=CC=C2)CCN)NC(=O)[C@H]"
                  "(CCN)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCN)NC(=O)CCCCC(C)C)O",
        "share": "~10-15% of polymyxin B drug substance",
    },
    {
        "name": "Colistin A (polymyxin E1)",
        "drug": "Colistin (polymyxin E)",
        "pubchem_cid": 24825758,
        "formula": "C53H100N16O13",
        "mw": 1169.5,
        "inchikey": "XDJYMJULXQKGMM-HHAJOKTESA-N",
        "residues": f"{DAB}T{DAB}{DAB}{DAB}LL{DAB}{DAB}T",
        "d_residue": "D-Leu6",
        "acyl": "(S)-6-methyloctanoyl",
        "smiles": "CC[C@H](C)CCCCC(=O)N[C@@H](CCN)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCN)"
                  "C(=O)N[C@H]1CCNC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)"
                  "[C@H](NC(=O)[C@@H](NC1=O)CCN)CC(C)C)CC(C)C)CCN)CCN)[C@@H](C)O",
        "share": "~60-70% of colistin drug substance",
    },
    {
        "name": "Colistin B (polymyxin E2)",
        "drug": "Colistin (polymyxin E)",
        "pubchem_cid": 25138298,
        "formula": "C52H98N16O13",
        "mw": 1155.4,
        "inchikey": "KNIWPHSUTGNZST-SSWRVQTPSA-N",
        "residues": f"{DAB}T{DAB}{DAB}{DAB}LL{DAB}{DAB}T",
        "d_residue": "D-Leu6",
        "acyl": "6-methylheptanoyl",
        "smiles": "C[C@H]([C@H]1C(=O)NCC[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N"
                  "[C@H](C(=O)N[C@H](C(=O)N1)CCN)CCN)CC(C)C)CC(C)C)CCN)NC(=O)[C@H](CCN)NC(=O)"
                  "[C@H]([C@@H](C)O)NC(=O)[C@H](CCN)NC(=O)CCCCC(C)C)O",
        "share": "~25-35% of colistin drug substance",
    },
]

# Dab has no single-letter code. These are used ONLY by the homology diagnostic - no saved
# configuration reads the ESM-2 sequence block - and running all four keeps the dependence on
# an arbitrary spelling visible rather than hidden.
DAB_SURROGATES = {
    "K": "Lys - cationic, side chain two carbons longer",
    "R": "Arg - cationic, guanidinium rather than amine",
    "Q": "Gln - same side-chain length, neutral",
    "A": "Ala - side chain removed entirely",
}

# How the N-terminal acyl tail is written decides which bucket it lands in. `literal` is the
# chemical name as PubChem gives it; `normalised` is what a curator following MuriVerse's own
# bucket vocabulary would write. bucket_mod resolves both.
NTER_VARIANTS = {
    "literal": lambda acyl: acyl,
    "normalised": lambda acyl: f"Lipidation ({acyl})",
}

# The C-terminal Thr10 carboxyl is consumed by the macrolactam to Dab4's side-chain amine,
# so the molecule has neither a free terminal carboxyl nor a C-terminal amide. PEPlife2 has
# no vocabulary for "in a ring", and `topology` already carries the ring, so `Free` is used
# in its sense of "no terminal modification recorded" and the cyclisation is carried by
# `topology` alone.
#
# `Amidation` was previously run alongside it, on the reasoning that the macrolactam bond is
# itself an amide. That reasoning was wrong and the variant has been removed. A C-terminal
# amide is -C(=O)NH2; the macrolactam is Thr10-C(=O)-NH-Dab4, an isopeptide bond to a
# side-chain amine. Both are amides, but only one is a C-TERMINAL amide, and colistin and
# polymyxin B have the ring, not the terminal amide. Setting `Amidation` therefore switched
# on an annotation for a modification these molecules do not carry, while the SMILES - which
# already encodes the macrocycle - stayed identical. It did not describe an alternative
# chemistry, it described a false one, and because it moved the prediction further than the
# drug, the congener or the species did, it dominated every real contrast in the analysis.
CTER_VARIANTS = ("Free",)

# Route IS one of the shipped model's columns, and the retrieved observations span several
# routes, so it is varied rather than fixed at one value and hidden. It is a property of the
# STUDY, not of the molecule, which is the comparison the figure is built on: a study-design
# choice moves this model's prediction, and swapping one drug for another does not.
# Only the routes the applied model was actually FITTED on. Under the intravenous cohort the
# route columns are constant in training, so a subcutaneous input is not a different point in
# a dimension the model has learned - it is an extrapolation along an axis with zero training
# variance, and the number that comes back is an artefact of the scaler, not a prediction.
# Panel b is the measurement that says what such a prediction would be worth: a cross-route
# read is off by 3-4 fold. So the SC variant is dropped rather than drawn with a caveat.
ROUTE_VARIANTS = ("iv",) if MODEL_TAG == "_iv" else ("iv", "sc")
SPECIES_LIST = ("mouse", "rat")   # the two species the model was trained on

# ---------------------------------------------------------------------------------------
# Experimental mouse half-lives, for the comparison axis. Retrieved from PubMed / PMC; each
# is a total-drug plasma terminal half-life in mice, with its route and dose recorded
# because for these drugs both matter and the literature spread is wide.
# ---------------------------------------------------------------------------------------
# Observed half-lives live in data/polymyxin_observed_half_lives.csv, one row per value, with
# the study they come from (data/polymyxin_studies.csv) and whether the value is a plasma or
# serum elimination half-life that can be compared with the model (`use_in_band`). Citations
# were checked against PubMed: an earlier table gave the wrong PMID for two studies.
OBSERVED_CSV = ROOT / "data" / "polymyxin_observed_half_lives.csv"
STUDIES_CSV = ROOT / "data" / "polymyxin_studies.csv"
GROUPS_CSV = ROOT / "data" / "polymyxin_study_groups.csv"   # lab of origin and how t1/2 was computed


def observed_groups(obs: pd.DataFrame) -> pd.DataFrame:
    """One log10 half-life per dose group (study, drug, species, route, dose) with its animals.

    Components of one dose group (B1, B2, ...) come from the same animals, so they are averaged
    and the group's animals counted once."""
    band = obs[obs["use_in_band"]].copy()
    band["n_animals"] = pd.to_numeric(band["n_animals"], errors="coerce")
    band["sd_h"] = pd.to_numeric(band["sd_h"], errors="coerce")
    if band["n_animals"].isna().any():
        missing = band.loc[band["n_animals"].isna(), ["study_id", "dose"]].drop_duplicates()
        raise SystemExit(f"n_animals is missing for in-band values:\n{missing}")
    keys = ["drug", "species", "study_id", "route", "dose"]
    band["derived"] = band["detail"].astype(str).str.startswith("derived")
    return (band.groupby(keys)
            .agg(group_log10_h=("log10_half_life_h", "mean"), n_animals=("n_animals", "first"),
                 half_life_h=("half_life_h", "mean"), sd_h=("sd_h", "first"),
                 derived=("derived", "all"))
            .reset_index())


def centre_and_sd(g: pd.DataFrame) -> tuple[float, float, str]:
    """The animal-weighted mean log10 half-life of dose groups, its SD and where the SD comes from.

    Several groups: the animal-weighted SD across them. One group: the SD the study reports,
    converted to log10 units by the delta method (sd_h / (half_life * ln 10)). Otherwise NaN."""
    w = g["n_animals"].to_numpy(float)
    x = g["group_log10_h"].to_numpy(float)
    centre = float(np.average(x, weights=w))
    if len(g) > 1:
        return centre, float(np.sqrt(np.average((x - centre) ** 2, weights=w))), \
            "animal-weighted SD across dose groups"
    if g["sd_h"].notna().all():
        return centre, float(g["sd_h"].iloc[0] / (g["half_life_h"].iloc[0] * np.log(10))), \
            "SD reported by the study"
    return centre, np.nan, "single group, no SD"


def observed_summary(obs: pd.DataFrame) -> pd.DataFrame:
    """The observed half-life per drug and species, weighted by the animals dosed.

    The centre is the animal-weighted mean of the dose-group log10 half-lives. The band is that
    mean +/- one animal-weighted SD across dose groups. A drug and species with a single dose
    group has no between-group spread, so the band uses the SD the paper reports for it,
    converted to log10 units by the delta method (sd_h / (half_life * ln 10))."""
    rows = []
    for (drug, species), g in observed_groups(obs).groupby(["drug", "species"]):
        centre, sd, basis = centre_and_sd(g)
        rows.append({"drug": drug, "species": species, "observed_log10_h": centre,
                     "band_low_log10_h": centre - sd, "band_high_log10_h": centre + sd,
                     "band_basis": basis, "n_studies": g["study_id"].nunique(),
                     "n_groups": len(g), "n_animals": int(g["n_animals"].sum())})
    out = pd.DataFrame(rows)
    out["observed_h"] = 10 ** out["observed_log10_h"]
    return out


def observed_by_source(obs: pd.DataFrame, studies: pd.DataFrame) -> pd.DataFrame:
    """The observed half-life per drug, species and source of the value.

    The sources are the lab groups in data/polymyxin_study_groups.csv (the Tam lab studies, Jiao
    2023, all others), each split into the values the studies report and the values we derived
    from their published parameters (`value_basis`). Each is the animal-weighted geometric mean
    of its dose groups, as for the pooled value, with its SD (`sd_log10_h`).

    Three sensitivity versions of the pooled value are added (`pooled: ...`): every study counted
    once, whatever its number of dose groups or animals; reported values only; intravenous
    doses only. A mean residence time is not a half-life, so it is never pooled; it gets a row
    of its own."""
    lab = studies.set_index("study_id")["source_group"]
    g = observed_groups(obs)
    g["source_group"] = g["study_id"].map(lab)
    g["value_basis"] = np.where(g["derived"], "derived", "reported")
    rows = []

    def add(drug, species, source, basis, s, centre, sd, n_studies):
        rows.append({"drug": drug, "species": species, "source_group": source, "study_id": "",
                     "observed_h": 10 ** centre, "n_studies": n_studies,
                     "n_animals": int(s["n_animals"].sum()), "quantity": "elimination half-life",
                     "value_basis": basis, "sd_log10_h": sd})

    for (drug, species, source, basis), s_ in g.groupby(["drug", "species", "source_group", "value_basis"]):
        centre, sd, _ = centre_and_sd(s_)
        add(drug, species, source, basis, s_, centre, sd, s_["study_id"].nunique())
    for (drug, species), s_ in g.groupby(["drug", "species"]):
        # each study once: the animal-weighted mean of each study's own dose groups, then an
        # equal-weight mean of those study means
        study_means = s_.groupby("study_id").apply(
            lambda d: np.average(d["group_log10_h"], weights=d["n_animals"]), include_groups=False)
        sd = float(study_means.std(ddof=0)) if len(study_means) > 1 else np.nan
        add(drug, species, "pooled: each study once", "", s_, float(study_means.mean()), sd,
            len(study_means))
        for label, keep in (("pooled: reported only", ~s_["derived"]),
                            ("pooled: IV only", s_["route"].astype(str).str.startswith("IV"))):
            t = s_[keep]
            if len(t):
                centre, sd, _ = centre_and_sd(t)
                add(drug, species, label, "", t, centre, sd, t["study_id"].nunique())
    for _, r in obs[obs["kind"] == "mrt"].iterrows():
        rows.append({"drug": r["drug"], "species": r["species"],
                     "source_group": "mean residence time", "study_id": r["study_id"],
                     "observed_h": float(r["half_life_h"]), "n_studies": 1,
                     "n_animals": int(float(r["n_animals"])), "quantity": "mean residence time",
                     "value_basis": "", "sd_log10_h": np.nan})
    return pd.DataFrame(rows)


def pooled_sensitivity(obs: pd.DataFrame) -> pd.DataFrame:
    """How far the pooled observed half-life moves under other reasonable choices.

    Per drug and species: the animal-weighted pooled value used in Figure 4, the same with every
    study counted once, with reported values only, with intravenous doses only, and the lowest
    and highest value when each study is left out in turn. `top_study` is the study with the
    largest share of the animal weight."""
    g = observed_groups(obs)
    rows = []
    for (drug, species), s_ in g.groupby(["drug", "species"]):
        def pooled(t):
            return 10 ** float(np.average(t["group_log10_h"], weights=t["n_animals"]))
        study_means = s_.groupby("study_id").apply(
            lambda d: np.average(d["group_log10_h"], weights=d["n_animals"]), include_groups=False)
        loso = {sid: pooled(s_[s_["study_id"] != sid]) for sid in s_["study_id"].unique()
                if s_["study_id"].nunique() > 1}
        share = s_.groupby("study_id")["n_animals"].sum() / s_["n_animals"].sum()
        iv = s_[s_["route"].astype(str).str.startswith("IV")]
        rep_ = s_[~s_["derived"]]
        rows.append({
            "drug": drug, "species": species, "n_studies": s_["study_id"].nunique(),
            "pooled_animal_weighted_h": pooled(s_),
            "each_study_once_h": 10 ** float(study_means.mean()),
            "reported_only_h": pooled(rep_) if len(rep_) else np.nan,
            "iv_only_h": pooled(iv) if len(iv) else np.nan,
            "leave_one_out_min_h": min(loso.values()) if loso else np.nan,
            "leave_one_out_max_h": max(loso.values()) if loso else np.nan,
            "leave_out_for_min": min(loso, key=loso.get) if loso else "",
            "leave_out_for_max": max(loso, key=loso.get) if loso else "",
            "top_study": share.idxmax(), "top_study_weight_share": float(share.max()),
            "derived_weight_share": float(s_.loc[s_["derived"], "n_animals"].sum() / s_["n_animals"].sum()),
        })
    return pd.DataFrame(rows)


def fold_error_by_source(err: pd.DataFrame, by_source: pd.DataFrame) -> pd.DataFrame:
    """Every prediction against every source's observed value. `pred_over_observed` is signed
    (above 1 = the model is too long); `fold_error` is the same ratio folded to at least 1."""
    out = err[["role", "drug", "species", "name", "pred_half_life_h"]].merge(
        by_source, on=["drug", "species"], how="inner")
    out["pred_over_observed"] = out["pred_half_life_h"] / out["observed_h"]
    out["fold_error"] = np.maximum(out["pred_over_observed"], 1 / out["pred_over_observed"])
    return out


def parse_dose(text: str) -> tuple[float, float]:
    """The dose in mg/kg (a range gives its low and high), or NaN when none is stated."""
    m = re.search(r"(\d+(?:\.\d+)?)(?:\s*-\s*(\d+(?:\.\d+)?))?\s*mg(?: base)?/kg", str(text))
    if not m:
        return np.nan, np.nan
    lo = float(m.group(1))
    return lo, float(m.group(2)) if m.group(2) else lo


def compare_with_observed(pred: pd.DataFrame, baseline: dict[str, float], iv_only: bool = False
                          ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The observed half-lives, and the fold error of every prediction against them.

    `iv_only` keeps only intravenous dose groups. The shipped model is fitted on intravenous
    records and predicts intravenous inputs, so the intravenous comparison is the primary one;
    the all-route comparison (`iv_only=False`) is the secondary analysis.

    Only the studies flagged `in_figure4` in data/polymyxin_studies.csv are used. A value is
    compared with the model when it is a plasma or serum elimination half-life of the drug
    (`use_in_band`), and is weighted by the animals dosed in its group (`n_animals`).
    `baseline` is each species' training-set median half-life, which every configuration has
    to beat."""
    studies = pd.read_csv(STUDIES_CSV, dtype=str, keep_default_na=False)
    groups = pd.read_csv(GROUPS_CSV, dtype=str, keep_default_na=False)
    studies = studies.merge(groups, on="study_id", how="left")
    studies[["source_group", "half_life_method"]] = studies[["source_group", "half_life_method"]].fillna("")
    chosen = set(studies.loc[studies["in_figure4"] == "True", "study_id"])
    obs = pd.read_csv(OBSERVED_CSV, keep_default_na=False)
    obs = obs[obs["study_id"].isin(chosen)].copy()
    obs["use_in_band"] = obs["use_in_band"].astype(str).eq("True")
    obs[["dose_lo_mg_kg", "dose_hi_mg_kg"]] = [parse_dose(d) for d in obs["dose"]]
    obs["log10_half_life_h"] = np.log10(obs["half_life_h"].astype(float))
    if iv_only:
        obs = obs[obs["route"].astype(str).str.startswith("IV")].copy()

    obs_sum = observed_summary(obs)
    err = pred.merge(obs_sum[["drug", "species", "observed_h", "n_studies"]],
                     on=["drug", "species"], how="left")
    err["baseline_half_life_h"] = err["species"].map(baseline)

    def fold(a, b):
        return np.maximum(a / b, b / a)

    err["model_fold_error"] = fold(err["pred_half_life_h"], err["observed_h"])
    err["baseline_fold_error"] = fold(err["baseline_half_life_h"], err["observed_h"])
    err["model_beats_baseline"] = err["model_fold_error"] < err["baseline_fold_error"]
    err = err.dropna(subset=["observed_h"]).reset_index(drop=True)
    return studies, obs, obs_sum, err


def write_observed_tables(studies: pd.DataFrame, obs: pd.DataFrame, obs_sum: pd.DataFrame,
                          err: pd.DataFrame, suffix: str = "") -> None:
    """Tables 21, 27, 49, 50 and 52: everything that depends on which studies are drawn.

    With `suffix="_iv"` (observations from `compare_with_observed(..., iv_only=True)`) only the
    four tables that carry a comparison value are written, under an `_iv` suffix; the study
    lists (21, 50) and the pooling sensitivity (55) are properties of the all-route set."""
    TABS.mkdir(parents=True, exist_ok=True)
    by_source = observed_by_source(obs, studies)
    if suffix:
        by_source.to_csv(TABS / f"table52_polymyxin_observed_by_source{suffix}.csv", index=False)
        fold_error_by_source(err, by_source).to_csv(
            TABS / f"table53_polymyxin_fold_error_by_source{suffix}.csv", index=False)
        obs_sum.to_csv(TABS / f"table49_polymyxin_observed_by_drug_species{suffix}.csv", index=False)
        err.to_csv(TABS / f"table27_polymyxin_fold_error_vs_baseline{suffix}.csv", index=False)
        return
    by_source.to_csv(TABS / "table52_polymyxin_observed_by_source.csv", index=False)
    pooled_sensitivity(obs).to_csv(TABS / "table55_polymyxin_pooled_sensitivity.csv", index=False)
    fold_error_by_source(err, by_source).to_csv(
        TABS / "table53_polymyxin_fold_error_by_source.csv", index=False)
    obs_out = obs.merge(studies[["study_id", "citation", "pmid", "pmcid", "doi", "source_group"]],
                        on="study_id", how="left")
    obs_out.to_csv(TABS / "table21_polymyxin_observed_half_lives.csv", index=False)
    studies[studies["in_figure4"] == "True"].to_csv(
        TABS / "table50_polymyxin_studies_used.csv", index=False)
    obs_sum.to_csv(TABS / "table49_polymyxin_observed_by_drug_species.csv", index=False)
    err.to_csv(TABS / "table27_polymyxin_fold_error_vs_baseline.csv", index=False)


# ---------------------------------------------------------------------------------------
# The saved model
# ---------------------------------------------------------------------------------------
def load_models() -> tuple[dict, dict]:
    """The saved models plus the layout that says how to assemble their inputs."""
    if not (MODEL_DIR / "model.json").exists():
        raise SystemExit(
            f"{MODEL_DIR / 'model.json'} is missing. Fit and save the production models "
            f"first:\n  python MuriVerse/evaluate_halflife_model.py --cohort rodent "
            f"--tier strict")
    layout = json.loads((MODEL_DIR / "model.json").read_text())
    metadata_path = STRICT / "halflife_meta_with_split.csv"
    expected_hash = layout.get("metadata_sha256")
    if expected_hash:
        actual_hash = hashlib.sha256(metadata_path.read_bytes()).hexdigest()
        if actual_hash != expected_hash:
            raise SystemExit("Saved model was fitted against a different strict metadata "
                             "snapshot; rerun evaluate_halflife_model.py before prediction")
    models = {}
    for role in ROLES:
        if role not in layout["models"]:
            continue
        models[role] = joblib.load(MODEL_DIR / f"{role}.joblib")
    return models, layout


def block_columns(name: str) -> list[str]:
    """One block's column names, in the order the training matrix used them.

    Read from the saved dataset rather than re-derived: the one-hot levels are taken in
    sorted order from the training data, so re-deriving them here from a different frame
    would be a second implementation that could silently disagree about column order.
    """
    blocks = json.loads((ABLATION / "feature_names.json").read_text())["blocks"]
    return next(b for b in blocks if b["name"] == name)["columns"]


def termini_vector(n_text: str, c_text: str, columns: list[str]) -> np.ndarray:
    """The partner block for a congener.

    Neither drug carries a macromolecular partner: the acyl tail is a C8/C9 fatty acyl on
    the N-terminus, which is a terminal modification rather than something bolted on, so
    `partner_class` is `none` and both partner masses are zero. The block still has to be
    built, because the selected model reads it and a block of zeros would assert "no
    N-terminal modification", which is false.
    """
    described = describe_termini(n_text, c_text)
    values = {
        "n_ter_class": described.n_class,
        "c_ter_class": described.c_class,
        "partner_class": "none",
    }
    vector = np.zeros(len(columns), dtype=np.float32)
    for i, column in enumerate(columns):
        if column == "log1p_peg_kda":
            vector[i] = float(np.log1p(described.peg_kda))
        elif column == "log1p_partner_kda":
            vector[i] = float(np.log1p(described.partner_kda))
        elif column == "has_n_mod":
            vector[i] = float(described.n_class != "free")
        elif column == "has_c_mod":
            vector[i] = float(described.c_class != "free")
        else:
            covariate, _, level = column.partition("=")
            if values.get(covariate) == level:
                vector[i] = 1.0
    return vector


# The model also reads study conditions (see study_conditions.py). The polymyxin pharmacokinetic
# literature measures drug in plasma by LC-MS/MS; strain and sex are left at the unknown level
# rather than guessed, so a prediction is for "plasma, LC-MS, route as varied below".
# Species is NOT in here, because it is not an assumption: it is varied per prediction, in
# the loop over SPECIES_LIST, and passed through `annotation_vector`. The model carries
# `species=mouse` / `species=rat` as columns of the study_context block, so a mouse and a rat
# prediction for the same molecule differ even with strain unknown. Dose is left missing (as
# it is for most training rows).
ASSUMED_STUDY_CONDITIONS = {"sample_matrix": "plasma", "assay_method": "LC-MS",
                            "strain": "unknown", "sex": "unknown",
                            "formulation_class": "not reported",
                            "pk_dose_unit": "not reported"}


def annotation_columns() -> tuple[list[str], list[str]]:
    """The molecule-annotation and study-context column names, in training-matrix order."""
    return block_columns("molecule_annotation"), block_columns("study_context")


def annotation_vector(fields: dict, molecule_columns: list[str], context_columns: list[str]
                      ) -> tuple[np.ndarray, np.ndarray, dict]:
    """One-hot the congener's curation fields into the two training blocks.

    A level the training data never saw simply leaves every column of its covariate at zero,
    which is the honest encoding of "outside the training vocabulary" - it is not folded into
    a neighbouring bucket. Continuous study-context columns (dose) stay missing (NaN).
    """
    molecule_values = {
        "hl_engineering": fields["hl_engineering"],
        "topology": fields["topology"],
        "chirality": fields["chirality"],
        "smiles_coverage": fields["smiles_coverage"],
    }
    context_values = {"route_model": fields["route_class"], "species": fields["species"],
                      **ASSUMED_STUDY_CONDITIONS}

    def one_hot(values: dict, columns: list[str]) -> np.ndarray:
        vector = np.zeros(len(columns), dtype=np.float32)
        for i, column in enumerate(columns):
            if "=" not in column:
                vector[i] = np.nan
                continue
            covariate, _, level = column.partition("=")
            if values.get(covariate) == level:
                vector[i] = 1.0
        return vector

    unseen = {k: v for k, v in {**molecule_values, **context_values}.items()
              if not any(c == f"{k}={v}" for c in molecule_columns + context_columns)}
    return (one_hot(molecule_values, molecule_columns),
            one_hot(context_values, context_columns), unseen)


def physchem_vector(smiles: str, columns: list[str]) -> tuple[np.ndarray, int]:
    """The physicochemical block for one congener, plus how many columns are NOT defined.

    The RDKit half of the block is read off the PubChem SMILES and is exactly as defined
    here as it is for a cohort record. The sequence-composition half is NOT: these molecules
    have no faithful single-letter sequence - six of the ten residues of a polymyxin are
    2,4-diaminobutyric acid, which has no code - so GRAVY and the residue fractions would be
    computed over the four residues that happen to have one. That is a fabricated number,
    not a measured one, so those columns are emitted as NaN, which XGBoost treats as missing
    and routes down the default branch it learned.

    The count is returned rather than swallowed because a structural configuration built on
    columns this script cannot fill is not the configuration that was cross-validated, and
    the caller refuses to ship such a model silently.
    """
    values = peptide_physchem.structure_descriptors(smiles)
    vector = np.array([values.get(c, np.nan) for c in columns], dtype=np.float32)
    return vector, int(np.isnan(vector).sum())


def assemble(blocks_needed: list[str], molecule: np.ndarray, context: np.ndarray,
             smiles_embedding: np.ndarray, termini: np.ndarray,
             layout: dict, physchem: np.ndarray | None = None) -> np.ndarray:
    """Concatenate the blocks a given role needs, in the EXACT order training used.

    `blocks_needed` is the saved per-role block list from `model.json`, written by
    `evaluate_halflife_model.py` in the same order the training matrix was concatenated.
    Reordering by `layout["block_order"]` here would silently diverge from that order
    whenever a configuration list disagrees with the layout, and the saved scaler + ridge
    weights would then be applied to the wrong columns - which is how a previous build
    produced polymyxin predictions of about 10^-131 h. The saved order is the truth.
    """
    available = {"molecule_annotation": molecule, "study_context": context,
                 "smiles": smiles_embedding, "termini": termini}
    if physchem is not None:
        available["physchem"] = physchem
    missing = [b for b in blocks_needed if b not in available]
    if missing:
        raise SystemExit(
            f"this script can build {sorted(available)} but the saved model needs {missing}. "
            f"The ESM-2 sequence block is deliberately not built here: these molecules have "
            f"no faithful single-letter sequence (6 of 10 residues are Dab), and no saved "
            f"configuration used it. Re-run evaluate_halflife_model.py if that has changed.")
    return np.concatenate([available[b] for b in blocks_needed]).reshape(1, -1)


# ---------------------------------------------------------------------------------------
# Leakage check
# ---------------------------------------------------------------------------------------
def assert_not_in_cohort() -> pd.DataFrame:
    """The two drugs must not appear as modelled rows. Asserted, not assumed."""
    raw = pd.read_csv(STRICT / "in_vivo_half_life.csv", low_memory=False)
    meta = pd.read_csv(STRICT / "halflife_meta_with_split.csv", low_memory=False)

    pat = r"polymyxin|colistin|colisti"
    name_hits = raw[raw["Name"].astype(str).str.contains(pat, case=False, na=False)]
    meta_hits = meta[meta["Name"].astype(str).str.contains(pat, case=False, na=False)]
    if len(name_hits) or len(meta_hits):
        raise RuntimeError(
            f"Leakage: {len(name_hits)} cohort rows and {len(meta_hits)} modelled rows are "
            "named for a polymyxin. This is no longer a held-out prediction.")

    # Widest surface this search has: the eligible set, which is the strict set's own parent
    # and carries the curated free-text columns (dose notes, extraction evidence) that a
    # comparator mention hides in. Reading it rather than the earlier per-species pooled
    # table keeps the whole check inside one nested chain - eligible -> strict -> modelled -
    # so a row counted here is guaranteed to be an ancestor of the rows counted above.
    # Verified equivalent: both surfaces return the single ZY4 comparator mention.
    eligible = pd.read_csv(ROOT / "rodent" / "eligible_in_vivo_half_life.csv", low_memory=False)
    anywhere = eligible.apply(
        lambda r: r.astype(str).str.contains(pat, case=False, na=False).any(), axis=1)
    return pd.DataFrame([
        {"check": "rows whose Name is a polymyxin, in the strict set", "count": len(name_hits)},
        {"check": "rows whose Name is a polymyxin, in the modelled cohort",
         "count": len(meta_hits)},
        {"check": "rows mentioning a polymyxin in any column of the ELIGIBLE set "
                  "(comparator free text)", "count": int(anywhere.sum())},
    ])


# ---------------------------------------------------------------------------------------
# Applicability domain: sequence homology to the training set
# ---------------------------------------------------------------------------------------
def _stockholm_rows(path: Path) -> tuple[str, dict[str, str]]:
    query, hits = "", {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        name, _, seq = line.partition(" ")
        seq = seq.strip()
        if not seq:
            continue
        if not query:
            query = seq
        hits[name] = hits.get(name, "") + seq
    return query, hits


def phmmer_identity(queries: dict[str, str], targets: list[str]) -> pd.DataFrame:
    """Best percent identity of each query against the strict training sequences.

    The same tool the training split's clusters come from, so "related to the training set"
    means the same thing here as it does there.
    """
    if shutil.which("phmmer") is None:
        return pd.DataFrame()
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        target_fa = tmp / "targets.fasta"
        target_fa.write_text("".join(
            f">train{i}\n{s}\n" for i, s in enumerate(targets) if s))
        for label, seq in queries.items():
            query_fa = tmp / "q.fasta"
            query_fa.write_text(f">{label}\n{seq}\n")
            out = tmp / "q.sto"
            result = subprocess.run(
                ["phmmer", "--max", "-A", str(out), "--noali", "-E", "1000",
                 str(query_fa), str(target_fa)],
                capture_output=True, text=True)
            if result.returncode != 0 or not out.exists():
                rows.append({"query": label, "best_percent_identity": np.nan, "n_hits": 0})
                continue
            qseq, hits = _stockholm_rows(out)
            best, n = 0.0, 0
            for name, hseq in hits.items():
                if name == label:
                    continue
                n += 1
                pairs = [(a, b) for a, b in zip(qseq, hseq) if a != "-" and b != "-"]
                if pairs:
                    best = max(best, 100.0 * sum(a.upper() == b.upper()
                                                 for a, b in pairs) / len(pairs))
            rows.append({"query": label, "best_percent_identity": best, "n_hits": n})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# Dose sensitivity
# ---------------------------------------------------------------------------------------
DOSE_GRID_MG_KG = (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 64)
DOSE_COLUMN = "log10_dose_mg_kg"


def dose_sensitivity(models, layout, molecule_columns, context_columns, termini_columns,
                     smiles_by_name, physchem_by_name) -> pd.DataFrame:
    """The shipped model's prediction for every congener and encoding variant when a dose is
    supplied, and when it is left missing (as in the main predictions).

    Dose is the one study condition the observed half-lives differ in (0.5 to 64 mg/kg), so this
    asks whether the model could be moved by it at all. Only 3 of the 320 training rows carry a
    dose, so the answer describes the model, not the pharmacology.

    Both animals are swept, because species is a model input: the fold change a dose produces
    is not guaranteed to be the same in a mouse and a rat, and reporting it from one animal
    would hide that. The consumer medians over the rows it selects, so the extra column
    changes the table's shape and not its summary."""
    role = "shipped"
    spec = layout["models"][role]
    dose_at = context_columns.index(DOSE_COLUMN)
    rows = []
    for species, congener in [(sp, c) for sp in SPECIES_LIST for c in CONGENERS]:
        for nter_name, nter_fn in NTER_VARIANTS.items():
            for cter in CTER_VARIANTS:
                for route in ROUTE_VARIANTS:
                    fields = {"hl_engineering": "lipidated", "route_class": route,
                              "species": species, "topology": "cyclic", "chirality": "mix",
                              "smiles_coverage": "complete"}
                    molecule_vec, context_vec, _ = annotation_vector(
                        fields, molecule_columns, context_columns)
                    termini = termini_vector(nter_fn(congener["acyl"]), cter, termini_columns)
                    for dose in (None, *DOSE_GRID_MG_KG):
                        context = context_vec.copy()
                        context[dose_at] = np.nan if dose is None else np.log10(dose)
                        X = assemble(spec["blocks"], molecule_vec, context,
                                     smiles_by_name[congener["name"]], termini, layout,
                                     physchem_by_name[congener["name"]])
                        pred = float(models[role].predict(X)[0])
                        rows.append({"species": species,
                                     "name": congener["name"], "drug": congener["drug"],
                                     "variant": f"{nter_name} N-term / {cter} C-term / {route}",
                                     "route_class": route, "dose_mg_kg": dose,
                                     "pred_half_life_h": 10 ** pred})
    out = pd.DataFrame(rows)
    missing = (out[out["dose_mg_kg"].isna()]
               .set_index(["species", "name", "variant"])["pred_half_life_h"])
    key = pd.MultiIndex.from_frame(out[["species", "name", "variant"]])
    out["fold_vs_dose_missing"] = out["pred_half_life_h"].to_numpy() / missing.reindex(key).to_numpy()
    return out


# ---------------------------------------------------------------------------------------
def main() -> None:
    global MODEL_TAG, MODEL_DIR, ROUTE_VARIANTS
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-suffix", default=MODEL_TAG,
                    help="Which saved production models to apply: '_iv' (the shipped, "
                         "route-restricted set) or '' (the older route-pooled set).")
    args = ap.parse_args()
    if args.model_suffix != MODEL_TAG:
        MODEL_TAG = args.model_suffix
        MODEL_DIR = STRICT / f"model{MODEL_TAG}"
        ROUTE_VARIANTS = ("iv",) if MODEL_TAG == "_iv" else ("iv", "sc")
    print(f"[Model] {MODEL_DIR}  routes applied: {'/'.join(ROUTE_VARIANTS)}")
    models, layout = load_models()
    molecule_columns, context_columns = annotation_columns()
    termini_columns = block_columns("termini")
    leak = assert_not_in_cohort()
    print(f"[Leakage] clear:\n{leak.to_string(index=False)}\n")

    meta = pd.read_csv(STRICT / "halflife_meta_with_split.csv", low_memory=False)
    # The baseline this application has to beat is "a constant, taken from the training data,
    # that needs no model". It must therefore be taken from the SAME training data the applied
    # model saw - which, under the intravenous shipped model, is the IV rows of the training
    # split and not all 257 of them. Taking it from the route-pooled split let the constant
    # borrow the extravascular labels, which are longer (the flip-flop offset Figure 3 panels a
    # and b measure, 0.45-0.65 log10). That pulled the mouse constant up from 0.49 h to 0.89 h,
    # onto the observed polymyxin values of 1.1-1.2 h, and produced the conclusion that a
    # constant beats the model in mouse. It beat it by holding labels the model was not allowed
    # to use. Restricted to the routes the model was fitted on, the model is the better
    # predictor in both species.
    train = meta[meta["split"] == "train"]
    if ROUTE_VARIANTS:
        train = train[train["route_class"].isin(ROUTE_VARIANTS)]
    print(f"[Baseline] species medians from the {'/'.join(ROUTE_VARIANTS)} training split: "
          f"{len(train)} of {int((meta['split'] == 'train').sum())} training rows")
    train_labels_by_species = {sp: train.loc[train["species"] == sp, "log_label"].to_numpy(float)
                               for sp in SPECIES_LIST}

    # One PeptideCLM embedding per congener; the SMILES is PubChem's, not a reconstruction.
    device = torch.device("cpu")
    smiles_list = [c["smiles"] for c in CONGENERS]
    print(f"[Encode] PeptideCLM over {len(smiles_list)} congener SMILES")
    smiles_emb = embed_peptideclm(smiles_list, device, batch_size=4,
                                  cache_path=ROOT / "peptideclm_embedding_cache.npz")
    smiles_by_name = {c["name"]: smiles_emb[i] for i, c in enumerate(CONGENERS)}

    # The physicochemical block, read off the same PubChem SMILES. Any column this script
    # cannot define for a molecule with no single-letter sequence comes back NaN, and a role
    # whose feature set depends on those columns is reported rather than quietly shipped.
    physchem_columns = (block_columns("physchem")
                        if any("physchem" in layout["models"][r]["blocks"] for r in models)
                        else [])
    physchem_by_name, undefined_physchem = {}, {}
    for congener in CONGENERS:
        vector, n_missing = physchem_vector(congener["smiles"], physchem_columns)
        physchem_by_name[congener["name"]] = vector
        undefined_physchem[congener["name"]] = n_missing
    if physchem_columns:
        worst = max(undefined_physchem.values())
        print(f"[Physchem] {len(physchem_columns)} columns; {worst} undefined per congener "
              f"(sequence-composition columns: these molecules have no faithful sequence)")
        if worst:
            print("  WARNING: the structural configuration reads columns this script cannot "
                  "fill, so the applied model is not exactly the cross-validated one.")

    rows, unseen_levels = [], []
    for species, congener in [(sp, c) for sp in SPECIES_LIST for c in CONGENERS]:
        for nter_name, nter_fn in NTER_VARIANTS.items():
          for cter in CTER_VARIANTS:
            for route in ROUTE_VARIANTS:
                nter_text = nter_fn(congener["acyl"])
                # The acyl tail is a C8/C9 fatty acyl on the N-terminus, so the engineering
                # class is `lipidated` however the tail is spelled. bucket_mod is still run,
                # because WHICH n-terminal bucket the spelling lands in is the thing the two
                # NTER_VARIANTS exist to vary, and it is recorded per row below.
                nter_bucket = bucket_mod(nter_text)
                fields = {
                    "hl_engineering": "lipidated",
                    "route_class": route,
                    "species": species,
                    "topology": "cyclic",
                    "chirality": "mix",
                    "smiles_coverage": "complete",
                }
                molecule_vec, context_vec, unseen = annotation_vector(
                    fields, molecule_columns, context_columns)
                termini = termini_vector(nter_text, cter, termini_columns)
                if unseen:
                    unseen_levels.append({"species": species, "congener": congener["name"], **unseen})
                variant = f"{nter_name} N-term / {cter} C-term / {route}"
                for role, model in models.items():
                    spec = layout["models"][role]
                    X = assemble(spec["blocks"], molecule_vec, context_vec,
                                 smiles_by_name[congener["name"]], termini, layout,
                                 physchem_by_name[congener["name"]])
                    pred = float(model.predict(X)[0])
                    rows.append({
                        "species": species,
                        "role": role,
                        "features": spec["features"],
                        "name": congener["name"],
                        "drug": congener["drug"],
                        "variant": variant,
                        "n_ter_text": nter_text,
                        "n_ter_bucket": nter_bucket,
                        "c_ter_text": cter,
                        "route_class": route,
                        "pred_log10_half_life_h": pred,
                        "pred_half_life_h": float(10 ** pred),
                    })

    pred = pd.DataFrame(rows)

    # --- the degeneracy, measured ------------------------------------------------------
    # How many DISTINCT predictions does each configuration make across the four congeners,
    # holding the encoding variant fixed? One means the configuration cannot tell the
    # molecules apart at all.
    degeneracy = (pred.groupby(["species", "role", "features", "variant"])
                  .agg(n_congeners=("name", "nunique"),
                       n_distinct_predictions=("pred_half_life_h",
                                               lambda s: s.round(9).nunique()),
                       min_h=("pred_half_life_h", "min"),
                       max_h=("pred_half_life_h", "max"))
                  .reset_index())
    degeneracy["congener_fold_spread"] = degeneracy["max_h"] / degeneracy["min_h"]

    # --- observed, and error against it -------------------------------------------------
    # The baseline every configuration has to beat: that species' training-set median half-life.
    baseline = {sp: float(10 ** np.median(v)) for sp, v in train_labels_by_species.items()}
    studies, obs, obs_sum, err = compare_with_observed(pred, baseline)
    _, obs_iv, obs_sum_iv, err_iv = compare_with_observed(pred, baseline, iv_only=True)

    # --- applicability domain -----------------------------------------------------------
    queries = {}
    for congener in CONGENERS:
        for letter in DAB_SURROGATES:
            queries[f"{congener['name']} [Dab->{letter}]"] = \
                congener["residues"].replace(DAB, letter)
    hom = phmmer_identity(queries, meta.loc[meta["split"] == "train", "sequence"]
                          .dropna().astype(str).tolist())

    ctx = []
    for role, spec in layout["models"].items():
        if role not in models:
            continue
        for sp, labels in train_labels_by_species.items():
            ctx.append({
                "role": role, "species": sp, "features": spec["features"], "dims": spec["dims"],
                "cluster_pmid_spearman": spec["cluster_pmid_spearman"],
                "baseline_split_n": int(len(labels)),
                "train_median_half_life_h": baseline[sp],
                "train_min_half_life_h": float(10 ** labels.min()),
                "train_max_half_life_h": float(10 ** labels.max()),
            })
    ctx = pd.DataFrame(ctx)

    struct = pd.DataFrame([{k: c[k] for k in
                            ("name", "drug", "pubchem_cid", "formula", "mw", "inchikey",
                             "residues", "d_residue", "acyl", "share")} for c in CONGENERS])

    TABS.mkdir(parents=True, exist_ok=True)
    pred.to_csv(TABS / "table20_polymyxin_predictions.csv", index=False)
    write_observed_tables(studies, obs, obs_sum, err)
    write_observed_tables(studies, obs_iv, obs_sum_iv, err_iv, suffix="_iv")
    struct.to_csv(TABS / "table22_polymyxin_structures.csv", index=False)
    degeneracy.to_csv(TABS / "table23_polymyxin_congener_degeneracy.csv", index=False)
    leak.to_csv(TABS / "table24_polymyxin_leakage_check.csv", index=False)
    if not hom.empty:
        hom.to_csv(TABS / "table25_polymyxin_surrogate_homology.csv", index=False)
    ctx.to_csv(TABS / "table26_polymyxin_training_label_context.csv", index=False)
    if "shipped" in models and DOSE_COLUMN in context_columns:
        dose_sensitivity(models, layout, molecule_columns, context_columns, termini_columns,
                         smiles_by_name, physchem_by_name).to_csv(
            TABS / "table54_polymyxin_prediction_vs_dose.csv", index=False)

    print("[Context]")
    print(ctx.to_string(index=False))
    print("\n[Degeneracy] distinct predictions across the four congeners, per variant")
    print(degeneracy.to_string(index=False))
    print("\n[Predictions] median over encoding variants, per congener")
    print(pred.groupby(["species", "role", "drug", "name"])["pred_half_life_h"]
          .median().round(3).to_string())
    print(f"\n[Observed] geometric mean of study values, per drug and species (h):\n"
          f"{obs_sum.round(3).to_string(index=False)}")
    if unseen_levels:
        print(f"\n[Vocabulary] {len(unseen_levels)} field value(s) outside the training "
              f"levels, encoded as all-zero: {unseen_levels[:3]}")
    print(f"\n[Saved] {TABS}/table20..table27")


if __name__ == "__main__":
    main()
