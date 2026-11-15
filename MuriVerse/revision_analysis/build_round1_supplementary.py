"""Write the round-1 supplementary tables (ensemble members, cohort sensitivity) from the analysis outputs."""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
A = Path(__file__).resolve().parents[2] / "results" / "revision_analysis"
TABS = REPO / "figures" / "tabs"

SCENARIOS = {
    "all_iv": "All intravenous records (shipped cohort)",
    "exact_unambiguous_labels": "Exact, unambiguous half-life labels only",
    "source_contradicted_removed": "Records whose source contradicts the label removed",
    "source_corrected": "Contradicted records removed, three corrected values restored",
    "verified_body_or_abstract": "Verified in full primary-source text or abstract",
    "verified_body": "Verified in full primary-source text only",
    "verified_body_excluding_largest_patent": "Verified in full primary-source text, largest patent excluded",
}
ESTIMATORS = {"ensemble": "ensemble", "species_median": "species median"}
BLOCKS = {"molecule_annotation": "Molecule annotation", "study_context": "Study context",
          "termini": "Terminal classes", "physchem": "Physicochemical"}


def importance() -> pd.DataFrame:
    out = pd.read_csv(A / "feature_permutation_iv.csv")
    out["block"] = out["block"].map(BLOCKS)
    return out.rename(columns={"delta_mae_log10": "increase_in_mae_log10", "ci_low": "ci_low",
                               "ci_high": "ci_high"})


def members() -> pd.DataFrame:
    base = pd.read_csv(A / "model_member_baselines.csv")
    paired = pd.read_csv(A / "model_member_paired.csv")
    out = base.merge(paired[["estimator", "spearman", "mae_delta_vs_ensemble", "delta_ci_low",
                             "delta_ci_high"]], on="estimator")
    out = out.rename(columns={
        "mae_delta_vs_ensemble": "ensemble_minus_comparator_mae",
        "delta_ci_low": "difference_ci_low", "delta_ci_high": "difference_ci_high"})
    order = ["ensemble", "ridge", "svr", "xgboost", "species_median"]
    return out.set_index("estimator").loc[order].reset_index()


def sensitivity() -> pd.DataFrame:
    s = pd.read_csv(A / "cohort_sensitivity_metrics.csv")
    s = s[s["scenario"].isin(SCENARIOS)].copy()
    s["scenario"] = s["scenario"].map(SCENARIOS).fillna(s["scenario"])
    s["estimator"] = s["estimator"].map(ESTIMATORS).fillna(s["estimator"])
    return s


if __name__ == "__main__":
    importance().to_csv(TABS / "supp_table6_feature_permutation_importance.csv", index=False)
    members().to_csv(TABS / "supp_table7_ensemble_member_comparison.csv", index=False)
    sensitivity().to_csv(TABS / "supp_table8_cohort_sensitivity.csv", index=False)
    print(members().round(3).to_string(index=False))
