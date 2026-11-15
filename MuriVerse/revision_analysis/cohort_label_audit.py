"""Editor item 4: audit of the PK labels and route assignments of the 213 intravenous records.

"Verify the PK labels and route assignments, prioritizing ambiguous records, extreme values, and
heavily represented sources."

What this script can and cannot do. Verifying a half-life against its source paper needs the
paper; the retrieval ledger marks most source identifiers as not retrieved, and 207 of the 213
labels are PEPlife2 transcriptions with no curated full text. So this script does three things that
do not need the papers:

1. HARD CHECKS of internal consistency (route fields agree, the label is the log10 of the hours,
   the curated full-text rows agree with the label, route evidence quotes name the assigned route);
2. a RISK-RANKED WORKLIST - every record scored for being ambiguous, extreme, from a heavily
   represented source, or having its route assigned from the database column alone - to say which
   records to check against their papers first;
3. the IMPACT of each flagged group on the headline result, from the saved out-of-fold predictions,
   so the worklist can be read against what is at stake.

It writes analysis/cohort_label_route_audit_iv.csv and analysis/cohort_label_route_summary.json and
exits non-zero if a hard check fails.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ANALYSIS = Path(__file__).resolve().parents[2] / "results" / "revision_analysis"
sys.path.insert(0, str(REPO / "MuriVerse"))
import evaluate_halflife_model as E  # noqa: E402

META = REPO / "results" / "revision_analysis" / "snapshot" / "halflife_meta_with_split.csv"
OOF = REPO / "results" / "revision_analysis" / "snapshot" / "oof_predictions_iv.csv"
SEEDS = [f"y_pred_cluster_PMID_seed{i}" for i in range(3)]

IV_TERMS = re.compile(r"intravenous|\bi\.?v\.?\b|tail vein|jugular|bolus|infus|injected into the (tail )?vein|vein", re.I)
Z_EXTREME = 2.5
SHORT_H, LONG_H = 0.01, 24.0
HEAVY_SOURCE_N = 5     # a source contributing this many of the 213 records or more

FAILS: list[str] = []
REVIEW: list[str] = []   # findings that need a human look at the source, not script errors


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"  {'ok  ' if ok else 'FAIL'}  {label}{('  ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(label)


def rho_by_seed(y, preds, mask):
    return float(np.mean([E._rho(y[mask], p[mask]) for p in preds]))


def main() -> int:
    m = pd.read_csv(META, low_memory=False)
    iv = m[m["route_class"] == "iv"].reset_index(drop=True)
    oof = pd.read_csv(OOF)
    assert (oof["id"].to_numpy() == iv["id"].to_numpy()).all()
    y = oof["y_true"].to_numpy(float)
    preds = [oof[c].to_numpy(float) for c in SEEDS]
    iv["PMID"] = iv["PMID"].astype(str)

    print("1. hard consistency checks")
    check("route_class, route_model and the database route agree on every record",
          bool((iv.route_class == iv.route_model).all() and (iv.route_class == iv.database_route_model).all()))
    check("every record is intravenous", bool((iv.route_class == "iv").all()))
    check("log_label equals log10 of the half-life in hours",
          bool(np.allclose(iv.log_label, np.log10(iv.half_life_h), atol=1e-9)))
    check("the label (hours) equals half_life_h", bool(np.allclose(iv.label.astype(float), iv.half_life_h, atol=1e-9)))
    check("half-lives are positive and finite", bool(np.isfinite(iv.half_life_h).all() and (iv.half_life_h > 0).all()))
    check("no record carries the pipeline's own extreme-value QC flag", not iv.half_life_extreme_qc.astype(str).eq("True").any())
    curated = iv[iv.evidence_half_life.notna()]
    for _, r in curated.iterrows():
        num, unit = re.match(r"([\d.]+)\s*(\w+)", str(r.evidence_half_life)).groups()
        hours = float(num) * {"h": 1, "min": 1 / 60, "d": 24}.get(unit, np.nan)
        check(f"curated full-text half-life ({r.evidence_half_life}) equals the label", abs(hours - r.half_life_h) < 0.005,
              f"{r.id}: label {r.half_life_h:.3f} h")
    quoted = iv[iv.route_evidence_quote.notna()]
    not_iv = quoted[~quoted.route_evidence_quote.astype(str).apply(lambda q: bool(IV_TERMS.search(q)))]
    iv["flag_route_quote_not_iv"] = iv.id.isin(not_iv.id)
    print(f"  {len(quoted)} records quote the passage that fixes their route; {len(not_iv)} quote does not name an intravenous route")
    for _, r in not_iv.iterrows():
        REVIEW.append(f"{r.id} (PMID {r.PMID}, {r.species}, {r.Name}, {r.half_life_h:.3g} h): assigned iv, but the evidence quote reads "
                      f"\"{r.route_evidence_quote}\"")
    check("the 2 records whose curated full text is an ambiguous match are not auto-applied",
          int(iv.half_life_provenance.str.contains("ambiguous").sum()) == 2)

    print("2. flags")
    z = (iv.log_label - iv.log_label.mean()) / iv.log_label.std()
    iv["flag_ambiguous"] = (iv.half_life_provenance.str.contains("ambiguous")
                            | (iv.half_life_qualifier != "exact")
                            | iv.exact_duplicate_candidate.astype(str).eq("True")
                            | iv.repeat_observation_candidate.astype(str).eq("True"))
    iv["flag_extreme"] = (z.abs() > Z_EXTREME) | (iv.half_life_h < SHORT_H) | (iv.half_life_h > LONG_H)
    n_src = iv.PMID.map(iv.PMID.value_counts())
    iv["source_n_records"] = n_src
    iv["source_share"] = (n_src / len(iv)).round(3)
    iv["flag_heavy_source"] = n_src >= HEAVY_SOURCE_N
    iv["flag_route_from_database_only"] = iv.route_evidence_level == "database_route_column"
    iv["flag_phase_not_stated"] = iv.pk_phase == "unlabelled"
    iv["z_log_label"] = z.round(2)
    iv["priority_score"] = (3 * iv.flag_ambiguous + 3 * iv.flag_extreme + 2 * iv.flag_heavy_source
                            + 2 * iv.flag_route_from_database_only + 3 * iv.flag_route_quote_not_iv)
    for f in ("flag_ambiguous", "flag_extreme", "flag_heavy_source", "flag_route_from_database_only",
              "flag_route_quote_not_iv", "flag_phase_not_stated"):
        print(f"  {f:34s} {int(iv[f].sum()):4d} of 213")
    print(f"  label provenance: {iv.half_life_provenance.value_counts().to_dict()}")
    print(f"  route evidence:   {iv.route_evidence_level.value_counts().to_dict()}")

    top = iv.groupby("PMID").size().sort_values(ascending=False)
    print("  heaviest sources:", {k: f"{v} ({v / 213:.0%})" for k, v in top.head(4).items()})

    print("3. impact on the headline (Spearman rho of the saved out-of-fold predictions, mean over 3 seeds)")
    base = rho_by_seed(y, preds, np.ones(len(y), bool))
    impact = [{"group": "all 213 records", "n_removed": 0, "n_left": 213, "rho": base}]
    groups = {"the heaviest source (EP2277890A1, a patent)": iv.PMID == top.index[0],
              "the four heavy sources (>= 5 records)": iv.flag_heavy_source,
              "route from the database column only": iv.flag_route_from_database_only,
              "ambiguous / non-exact / repeat candidates": iv.flag_ambiguous,
              "extreme values": iv.flag_extreme,
              "any of the above": iv.flag_ambiguous | iv.flag_extreme | iv.flag_heavy_source | iv.flag_route_from_database_only}
    for name, mask in groups.items():
        keep = ~mask.to_numpy()
        impact.append({"group": f"without {name}", "n_removed": int(mask.sum()), "n_left": int(keep.sum()),
                       "rho": rho_by_seed(y, preds, keep) if keep.sum() > 10 else None})
        print(f"  without {name:52s} n={int(keep.sum()):3d}  rho = {impact[-1]['rho']:.3f}" if impact[-1]["rho"] is not None
              else f"  without {name}: too few records left")

    out = iv[["id", "PMID", "species", "Name", "sequence", "half_life_h", "half_life_qualifier", "half_life_provenance",
              "pk_phase", "assay_method", "route_class", "route_evidence_level", "route_evidence_quote",
              "source_n_records", "source_share", "z_log_label", "flag_ambiguous", "flag_extreme", "flag_heavy_source",
              "flag_route_from_database_only", "flag_route_quote_not_iv", "flag_phase_not_stated", "priority_score", "split"]].copy()
    out["oof_pred_log10_h"] = np.mean(preds, axis=0)
    out["oof_residual_log10"] = out["oof_pred_log10_h"] - np.log10(out.half_life_h)
    out = out.sort_values(["priority_score", "source_n_records"], ascending=False)
    out.to_csv(ANALYSIS / "cohort_label_route_audit_iv.csv", index=False)

    summary = {
        "n": 213, "flags": {f: int(iv[f].sum()) for f in iv.columns if f.startswith("flag_")},
        "label_provenance": iv.half_life_provenance.value_counts().to_dict(),
        "route_evidence": iv.route_evidence_level.value_counts().to_dict(),
        "heaviest_sources": {k: int(v) for k, v in top.head(6).items()},
        "heaviest_source_share": float(top.iloc[0] / 213),
        "four_heavy_sources_share": float(top[top >= HEAVY_SOURCE_N].sum() / 213),
        "extreme_records": iv[iv.flag_extreme][["id", "PMID", "half_life_h", "z_log_label"]].to_dict("records"),
        "impact_on_rho": impact,
        "needs_human_review": REVIEW,
        "thresholds": {"z": Z_EXTREME, "short_h": SHORT_H, "long_h": LONG_H, "heavy_source_n": HEAVY_SOURCE_N},
        "verified_against_source_papers": 0,
        "note": "No record was re-verified against its paper here; 207 of 213 labels are PEPlife2 "
                "transcriptions with no curated full text. The worklist ranks which to check first.",
    }
    (ANALYSIS / "cohort_label_route_summary.json").write_text(json.dumps(summary, indent=2, default=float))

    print("\ntop of the worklist:")
    print(out.head(12)[["id", "PMID", "species", "half_life_h", "source_n_records", "flag_ambiguous", "flag_extreme",
                        "flag_route_from_database_only", "priority_score"]].to_string(index=False))
    if REVIEW:
        print("\nREVIEW (needs a look at the source paper):")
        for r in REVIEW:
            print("  -", r)
    if FAILS:
        print(f"\n{len(FAILS)} HARD CHECK(S) FAILED:")
        for f in FAILS:
            print("  -", f)
        return 1
    print("\nHard checks passed. The worklist is a ranking, not a verification against the papers.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
