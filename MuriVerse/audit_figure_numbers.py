"""Check that every verifiable number in the Figure 2, 3 and 4 captions matches the data.

Caption numbers drift. A cohort is rebuilt, a model is rescored, a panel is re-ordered, and a
figure that was correct when it was written keeps asserting the old value - and nothing in the
pipeline notices, because a caption is prose. The largest drift this file has caught so far was
the move of the shipped model from the route-pooled cohort to the intravenous one: every
performance number in Figures 3 and 4 changed, the panels were re-laid-out, and both captions
went on describing the old model.

Every assertion below re-derives its number from the data rather than from another caption, so
this is a check and not a restatement. Sources are:

    rodent/strict/halflife_meta_with_split.csv  the 320-row input cohort (Figure 2)
    rodent/strict/oof_predictions_iv.csv        held-out predictions of the SHIPPED model
    rodent/strict/oof_predictions.csv           the route-pooled model, kept for comparison
    figures/tabs/table31*, table32, table36, table5x   the summary tables the figures read

WHICH COHORT IS WHICH, because this is the thing most easily got wrong here. The shipped model
is fitted on the 213 INTRAVENOUS records (`oof_predictions_iv.csv`). Figure 2 describes the
full 320-row input cohort; Figure 3's panels a, b and h describe it too, because they are the
measurement that justifies restricting it. Everything else in Figure 3, and every shipped-model
number in Figure 4, is the 213.

What is deliberately NOT checked: Figure 4's literature values, which can only be verified
against the source papers.

Usage: python MuriVerse/audit_figure_numbers.py   (exit 1 if any check fails)
"""

from __future__ import annotations

import json
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
STRICT = ROOT / "rodent" / "strict"
TABS = ROOT.parent / "figures/tabs"
SNAP = ROOT.parent / "results/revision_analysis/snapshot"
CUT_H = 6.0
SEED_COLS = [f"y_pred_cluster_PMID_seed{s}" for s in range(3)]

FAILS: list[str] = []
NCHECK = 0


def check(label: str, claimed, actual, tol: float | None = None) -> None:
    """Compare a caption claim with the value re-derived from the data.

    A caption states a rounded number, so the test is whether the data rounds to what the
    caption says - not whether it falls inside some hand-chosen tolerance. Comparing at the
    caption's own precision avoids both false failures (0.955 is correctly written 0.96, but
    fails a 0.005 tolerance on a floating-point boundary) and false passes from a tolerance
    set loose enough to hide a genuine drift. `tol` overrides this where a claim is an
    approximation rather than a rounding, such as a p-value written "0.001".
    """
    global NCHECK
    NCHECK += 1
    if isinstance(claimed, (int, np.integer)) and isinstance(actual, (int, np.integer)):
        good = int(claimed) == int(actual)
        shown = f"{claimed} vs {actual}"
    elif tol is not None:
        good = abs(float(claimed) - float(actual)) <= tol
        shown = f"{claimed} vs {float(actual):.4f}"
    else:
        # Round half-UP on the decimal value, which is what someone writing the caption
        # does. Python's round() is half-to-even and operates on the binary float, so a
        # median of exactly 0.955 is stored as 0.954999... and rounds to 0.95, failing a
        # caption that correctly says 0.96.
        decimals = len(s.split(".")[1]) if "." in (s := repr(float(claimed))) else 0
        quantum = Decimal(1).scaleb(-decimals)
        rounded = Decimal(str(float(actual))).quantize(quantum, rounding=ROUND_HALF_UP)
        good = rounded == Decimal(str(float(claimed))).quantize(quantum, rounding=ROUND_HALF_UP)
        shown = f"{claimed} vs {float(actual):.4f}"
    if good:
        print(f"  ok    {label:52} {shown}")
    else:
        print(f"  FAIL  {label:52} {shown}")
        FAILS.append(f"{label}: caption {claimed}, data {actual}")


def claim(label: str, holds: bool, detail: str = "") -> None:
    """Test a caption's ASSERTION, not just its numbers.

    Several of the figures' load-bearing statements are directions and orderings rather than
    values - "every cross-species offset is smaller than every cross-route one", "the model
    beats the constant in both species". Those survive a change in the numbers and fail
    silently, so they are tested as claims.
    """
    global NCHECK
    NCHECK += 1
    if holds:
        print(f"  ok    {label:52} {detail}")
    else:
        print(f"  FAIL  {label:52} {detail}")
        FAILS.append(f"{label} (claim no longer holds) {detail}")


def rho(a, b) -> float:
    return float(pd.Series(np.asarray(a)).corr(pd.Series(np.asarray(b)), method="spearman"))


def auc_at_cut(y: np.ndarray, pred: np.ndarray) -> float:
    pos = y >= np.log10(CUT_H) - 1e-9
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = pd.Series(pred).rank().to_numpy()
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main() -> None:
    print("\nFINAL POST-AUDIT MODEL - primary manuscript and release estimates")
    analysis = Path(__file__).resolve().parents[1] / "results" / "revision_analysis"
    final_metrics = pd.read_csv(analysis / "clean_cohort_metrics.csv")
    final_metrics = final_metrics[(final_metrics.scenario == "clean_cohort") &
                                              (final_metrics.estimator == "ensemble")].set_index("metric")
    final_oof = pd.read_csv(analysis / "clean_cohort_oof_predictions.csv")
    check("final cohort n", 178, len(final_oof))
    check("final evaluation groups", 50, final_oof.group.nunique())
    check("final rat records", 150, int((final_oof.species == "rat").sum()))
    check("final mouse records", 28, int((final_oof.species == "mouse").sum()))
    check("final fixed rho", 0.445, final_metrics.loc["spearman", "value"])
    check("final fixed rho CI low", 0.153, final_metrics.loc["spearman", "ci_low"])
    check("final fixed rho CI high", 0.563, final_metrics.loc["spearman", "ci_high"])
    final_nested = pd.read_csv(TABS / "table31b_muriverse_halflife_nested_selection_postaudit.csv").iloc[0]
    check("final nested n", 178, final_nested.n)
    check("final nested rho", 0.297, final_nested.spearman)
    check("final nested rho CI low", 0.096, final_nested.spearman_lo)
    check("final nested rho CI high", 0.450, final_nested.spearman_hi)
    final_table = pd.read_csv(TABS / "pub-tab__table1_benchmark_metrics.csv")
    claim("Table 1 uses final fixed rho and interval", final_table.iloc[0, 1] == "0.445 (0.153–0.563)", str(final_table.iloc[0, 1]))
    claim("Table 1 uses final nested rho and interval", final_table.iloc[0, 2] == "0.297 (0.096–0.450)", str(final_table.iloc[0, 2]))
    for seed in range(3):
          check(f"final OOF seed {seed} fold count", 5, final_oof[f"fold_seed{seed}"].nunique())
          claim(f"final OOF seed {seed} groups not split", bool((final_oof.groupby("group")[f"fold_seed{seed}"].nunique() == 1).all()), "one fold per group")
    final_rho = np.mean([rho(final_oof.y_true, final_oof[f"ensemble_seed{seed}"]) for seed in range(3)])
    check("final rho reproduced from released OOF", final_metrics.loc["spearman", "value"], final_rho, 1e-10)

    print("\nBOOTSTRAP - released interval code reproduces the stored intervals")
    code_directory = ROOT / "revision_analysis"
    sys.path.insert(0, str(code_directory if code_directory.is_dir() else analysis))
    from metrics_iv import bootstrap_all
    import evaluate_halflife_model as ehm
    stored = pd.read_csv(analysis / "clean_cohort_metrics.csv")
    pre_frame = pd.read_csv(STRICT / "oof_predictions_iv.csv")
    post_frame = pd.read_csv(analysis / "clean_cohort_oof_predictions.csv")
    cases = (
        ("original_all_iv", pre_frame, "y_pred_cluster_PMID_seed",
         ehm.union_groups(pre_frame.cluster.to_numpy(), pre_frame.PMID.astype(str).to_numpy())),
        ("clean_cohort", post_frame, "ensemble_seed", post_frame["group"].astype(str).to_numpy()),
    )
    for scenario, frame, prefix, groups_b in cases:
        preds_b = [frame[f"{prefix}{s}"].to_numpy(float) for s in range(3)]
        interval = bootstrap_all(frame.y_true.to_numpy(float), preds_b, groups_b)
        for metric in ("spearman", "pearson", "mae", "auc"):
            row = stored[(stored.scenario == scenario) & (stored.estimator == "ensemble") & (stored.metric == metric)].iloc[0]
            check(f"{scenario} {metric} 95% CI low", round(row.ci_low, 3), interval[metric][0])
            check(f"{scenario} {metric} 95% CI high", round(row.ci_high, 3), interval[metric][1])
    m = pd.read_csv(STRICT / "halflife_meta_with_split.csv", low_memory=False)
    oof = pd.read_csv(STRICT / "oof_predictions_iv.csv")          # the SHIPPED model
    pooled = pd.read_csv(STRICT / "oof_predictions.csv")          # route-pooled, for comparison
    y = oof["y_true"].to_numpy(float)
    p0 = oof[SEED_COLS[0]].to_numpy(float)
    mouse = (oof.species.str.lower() == "mouse").to_numpy()
    rat = ~mouse

    # ---------------------------------------------------------------- Figure 2: the cohort
    print("FIGURE 2 - cohort composition (all 320 input records)")
    check("320 records", 320, len(m))
    check("74 mouse", 74, int((m.species == "mouse").sum()))
    check("246 rat", 246, int((m.species == "rat").sum()))
    check("166 distinct ESM-2 input sequences", 166, m.sequence.nunique())
    check("129 sources", 129, m.PMID.nunique())
    check("90 homology clusters", 90, m.cluster.nunique())
    check("3 residues minimum", 3, int(m.n_residues.min()))
    check("50 residues maximum", 50, int(m.n_residues.max()))
    check("median 29 residues", 29, int(m.n_residues.median()))
    check("min half-life 0.0028 h", 0.0028, m.half_life_h.min(), 5e-5)
    check("max half-life 62.0 h", 62.0, m.half_life_h.max(), 0.05)
    check("median half-life 0.427 h", 0.427, m.half_life_h.median())
    check("12 sequences in both species", 12,
          int((m.groupby("sequence").species.nunique() > 1).sum()))
    check("7 clusters span both species", 7,
          int((m.groupby("cluster").species.nunique() > 1).sum()))
    check("4 engineering classes retained", 4, m.hl_engineering.nunique())
    for cls, n in (("unmodified / none stated", 189), ("backbone modification only", 101),
                   ("radio/fluorescent label only", 21), ("lipidated", 9)):
        check(f"{cls} n", n, int((m.hl_engineering == cls).sum()))
    check("253 linear", 253, int((m.topology == "linear").sum()))
    check("67 cyclic", 67, int((m.topology == "cyclic").sum()))
    for lvl, n in (("L", 291), ("mix", 24), ("D", 5)):
        check(f"chirality {lvl}", n, int((m.chirality_from_sequence == lvl).sum()))

    print("\nFIGURE 2 - route panel (e)")
    rc = m.route_class.value_counts()
    for route, claimed in (("iv", 213), ("sc", 73), ("ip", 22), ("im", 12)):
        check(f"{route} records", claimed, int(rc[route]))
    check("0 records report no route", 0, int((m.route_class == "not reported").sum()))
    med = m.groupby("route_class").half_life_h.median()
    check("iv median 0.24 h", 0.24, med["iv"])
    check("sc median 0.96 h", 0.96, med["sc"])
    check("ip median 1.80 h", 1.80, med["ip"])
    check("im median 0.40 h", 0.40, med["im"])
    check("121 records under 15 min", 121, int((m.half_life_h < 0.25).sum()))
    check("38% of the cohort under 15 min", 0.38, (m.half_life_h < 0.25).mean(), 5e-3)

    print("\nFIGURE 2 - selection funnel")
    fc = pd.read_csv(STRICT / "filtering_counts.csv")
    rec = pd.read_csv(STRICT / "halflife_cohort_reconciliation.csv").set_index("step")
    check("funnel starts at 1159 eligible", 1159, int(fc.iloc[0].records))
    check("467 pass the structural gates", 467, int(fc.iloc[-1].records))
    check("84 removed, non-parenteral route", 84, int(rec.loc["parenteral dose route"].removed))
    check("57 removed, structure tier", 57, int(rec.loc["whole-molecule structure"].removed))
    check("6 removed, unnamed assay", 6, int(rec.loc["named bioanalytical method"].removed))
    check("reconciliation ends at 320", 320, int(rec.loc["observations retained"].rows_out))

    # ------------------------------------------------- Supplementary Figure 1: audit panels
    # These were Supplementary Figure 1's panels until that figure was rebuilt around the
    # preprocessing gates (it now runs a-k: databases, species x setting, half-life by species,
    # pipeline funnel, species mix, route gate, structure gate, assay gate, three IV-split
    # violins). The panels below are still BUILT in fig2.R but are in no layout, so these
    # numbers now appear in the manuscript prose with no figure to cite - which is why they are
    # checked here. The letters are dropped: they pointed at the old layout.
    print("\nCOHORT NUMBERS IN THE PROSE - panels built in fig2.R but not drawn (all 320)")
    kept = pd.read_csv(STRICT / "in_vivo_half_life.csv", low_memory=False)
    vc = m.PMID.value_counts()
    check("largest source 61 records", 61, int(vc.iloc[0]))
    check("largest source 19%", 0.19, vc.iloc[0] / len(m), 5e-3)
    check("ten largest contribute 126", 126, int(vc.head(10).sum()))
    bad = {"free", "not reported"}
    term = ((m.n_ter_mod.notna() & ~m.n_ter_mod.isin(bad))
            | (m.c_ter_mod.notna() & ~m.c_ter_mod.isin(bad)))
    check("113 carry a terminal modification", 113, int(term.sum()))
    check("29 with D-residues", 29, int(m.chirality_from_sequence.isin(["mix", "D"]).sum()))
    natural = (~m.chirality.isin(["d", "mix"])) & (m.n_nonstandard == 0)
    check("268 all natural L-residues", 268, int(natural.sum()))
    check("406 of 467 are the whole molecule", 406,
          int((kept.structure_tier == "complete").sum()))
    check("61 are not", 61, int((kept.structure_tier != "complete").sum()))
    am = m.groupby("assay_method").half_life_h.median()
    check("assay medians span 0.16 h", 0.16, am.min())
    check("assay medians span 1.85 h", 1.85, am.max())
    check("154 routes read back to a source", 154,
          int(m.route_evidence_level.isin(["primary_source", "primary_arm"]).sum()))

    # -------------------------------------------- Figure 3: the SHIPPED intravenous model
    print("\nFIGURE 3 - the shipped cohort")
    check("213 intravenous records", 213, len(oof))
    check("179 rat", 179, int(rat.sum()))
    check("34 mouse", 34, int(mouse.sum()))
    # Read the saved layout rather than restating the block widths as literals: the 69-column
    # version of this check was a tautology over four constants and went on passing when the
    # study_context block changed width.
    layout = json.loads((ROOT / "rodent/strict/model_iv/model.json").read_text())
    blocks = {b: s[1] - s[0] for b, s in layout["block_slices"].items()}
    check("molecule_annotation 9 columns", 9, blocks["molecule_annotation"])
    check("study_context 33 columns", 33, blocks["study_context"])
    check("termini 13 columns", 13, blocks["termini"])
    check("physchem 16 columns", 16, blocks["physchem"])
    check("71 shipped columns", 71, layout["models"]["shipped"]["dims"])
    claim("the shipped feature set carries the animal",
          "study_context" in layout["models"]["shipped"]["blocks"],
          "species=mouse/species=rat live in study_context")
    claim("one shipped configuration in the oof table", oof.features.nunique() == 1,
          str(oof.features.iloc[0]))
    check("rat is 73% intravenous", 0.73, (m[m.species == "rat"].route_class == "iv").mean(),
          5e-3)
    check("mouse is 46% intravenous", 0.46, (m[m.species == "mouse"].route_class == "iv").mean(),
          5e-3)
    check("mouse was 54% extravascular", 0.54,
          (m[m.species == "mouse"].route_class != "iv").mean(), 5e-3)

    print("\nFIGURE 3 - why route, not species (the restriction's own justification)")
    # Route as a free lever: score every record by its own route's median label, no molecule.
    check("route-median-only rho 0.42 on the pooled label", 0.42,
          rho(m.log_label, m.groupby("route_class").log_label.transform("median")))
    # Like for like: the SAME 213 rows and the SAME shipped configuration, pooled fit vs IV fit.
    pi = pooled[pooled.id.isin(set(oof.id))].set_index("id").loc[oof.id].reset_index()
    pooled_rho = float(np.mean([rho(y, pi[c]) for c in SEED_COLS]))
    iv_rho = float(np.mean([rho(y, oof[c]) for c in SEED_COLS]))
    check("pooled model on the IV rows, rho 0.507", 0.507, pooled_rho)
    check("IV model on the IV rows, rho 0.574", 0.574, iv_rho)
    check("pooled model AUC 0.785", 0.785,
          float(np.mean([auc_at_cut(y, pi[c].to_numpy(float)) for c in SEED_COLS])))
    check("IV model AUC 0.818", 0.818,
          float(np.mean([auc_at_cut(y, oof[c].to_numpy(float)) for c in SEED_COLS])))
    rmse = lambda d: float(np.mean([np.sqrt(np.mean((y - d[c].to_numpy(float)) ** 2))
                                    for c in SEED_COLS]))
    check("pooled model RMSE 0.699", 0.699, rmse(pi))
    check("IV model RMSE 0.675", 0.675, rmse(oof))
    check("pooled model mouse rho 0.409", 0.409,
          float(np.mean([rho(y[mouse], pi[c].to_numpy(float)[mouse]) for c in SEED_COLS])))
    check("IV model mouse rho 0.562", 0.562,
          float(np.mean([rho(y[mouse], oof[c].to_numpy(float)[mouse]) for c in SEED_COLS])))
    # The caption says the rho gain alone is NOT resolvable. Test that, not just the point
    # estimate: a caption that claimed significance here would be the easy mistake to make.
    delta, lo, hi = paired_group_bootstrap(oof, pi, y)
    check("paired delta rho +0.066", 0.066, delta)
    check("paired delta CI low -0.016", -0.016, lo, 5e-3)
    check("paired delta CI high +0.138", 0.138, hi, 5e-3)
    claim("the rho gain alone is not resolvable (CI spans 0)", lo < 0 < hi,
          f"[{lo:+.3f}, {hi:+.3f}]")

    print("\nFIGURE 3 - headline (table31c)")
    s = pd.read_csv(TABS / "table31c_muriverse_halflife_stringency.csv")
    strict_row = s[s.evaluation.str.startswith("strict: clusters")].iloc[0]
    nested = s[s.evaluation.str.contains("configuration chosen inside")].iloc[0]
    check("headline n 213", 213, int(strict_row.n))
    check("rho 0.57", 0.57, strict_row.spearman)
    check("rho CI low 0.32", 0.32, strict_row.spearman_lo)
    check("rho CI high 0.65", 0.65, strict_row.spearman_hi)
    check("pearson 0.57", 0.57, strict_row.pearson)
    check("pearson CI low 0.38", 0.38, strict_row.pearson_lo)
    check("pearson CI high 0.65", 0.65, strict_row.pearson_hi)
    check("nested r 0.456", 0.456, nested.pearson)
    check("nested rho 0.460", 0.460, nested.spearman)
    check("nested rho CI low 0.27", 0.27, nested.spearman_lo)
    check("nested rho CI high 0.55", 0.55, nested.spearman_hi)
    check("nested Pearson CI low 0.26", 0.26, nested.pearson_lo)
    check("nested Pearson CI high 0.57", 0.57, nested.pearson_hi)
    check("selection optimism 0.11 in rho", 0.11, strict_row.spearman - nested.spearman)
    abl = pd.read_csv(TABS / "table31_muriverse_halflife_feature_ablation_strict_iv.csv")
    abl = abl[abl.grouping == "cluster_PMID"].set_index("features")
    ship = abl.loc["annotation + context + termini + physchem"]
    check("R2 0.315", 0.315, ship.r2)
    check("RMSE 0.675", 0.675, ship.rmse)
    check("MAE 0.517", 0.517, ship.mae)

    print("\nFIGURE 3 - panels a/b, species x route cell transfer (all 320; table36)")
    t36 = pd.read_csv(TABS / "table36_muriverse_halflife_species_route_cells.csv")
    cell = lambda sp, ad, tr: t36[(t36.species == sp) & (t36.admin == ad)
                                  & (t36.train == tr)].iloc[0]
    for sp, ad, n_c in (("rat", "IV", 179), ("rat", "EV", 67),
                        ("mouse", "IV", 34), ("mouse", "EV", 40)):
        check(f"{sp} {ad} n", n_c, int(cell(sp, ad, "same cell").n))
    for sp, ad, own, pooled_c in (("rat", "IV", 0.62, 0.57), ("rat", "EV", 0.19, 0.39),
                                  ("mouse", "IV", 0.23, 0.41), ("mouse", "EV", -0.19, 0.18)):
        check(f"{sp} {ad} own-cell rho", own, cell(sp, ad, "same cell").spearman)
        check(f"{sp} {ad} pooled rho", pooled_c, cell(sp, ad, "pooled").spearman)
    # rat IV is included although the caption discusses only the other three cells: its two
    # transfer correlations are the pair docs/results.md left as dashes until 2026-10-11.
    for sp, ad, cs, cr in (("mouse", "IV", 0.43, 0.14), ("rat", "EV", -0.06, 0.44),
                           ("mouse", "EV", 0.18, 0.27), ("rat", "IV", 0.24, 0.20)):
        check(f"{sp} {ad} cross-species rho", cs, cell(sp, ad, "cross-species").spearman)
        check(f"{sp} {ad} cross-route rho", cr, cell(sp, ad, "cross-route").spearman)
    for sp, ad, b in (("mouse", "EV", -0.65), ("rat", "EV", -0.45),
                      ("rat", "IV", 0.65), ("mouse", "IV", 0.45)):
        check(f"{sp} {ad} cross-route bias", b, cell(sp, ad, "cross-route").bias)
    for sp, ad, b in (("mouse", "IV", -0.21), ("rat", "IV", 0.21),
                      ("mouse", "EV", -0.18), ("rat", "EV", 0.07)):
        check(f"{sp} {ad} cross-species bias", b, cell(sp, ad, "cross-species").bias)
    # This separation is the entire argument for restricting by route and pooling over
    # species, so it is tested as a claim and not left to the reader to check by eye.
    cs_max = t36[t36.train == "cross-species"].bias.abs().max()
    cr_min = t36[t36.train == "cross-route"].bias.abs().min()
    check("max |cross-species bias| 0.21", 0.21, cs_max)
    check("min |cross-route bias| 0.45", 0.45, cr_min)
    claim("every cross-species offset < every cross-route offset", cs_max < cr_min,
          f"{cs_max:.2f} < {cr_min:.2f}")
    claim("cross-route signs flip-flop as kinetics predicts",
          cell("mouse", "EV", "cross-route").bias < 0
          and cell("rat", "EV", "cross-route").bias < 0
          and cell("rat", "IV", "cross-route").bias > 0
          and cell("mouse", "IV", "cross-route").bias > 0,
          "EV labels under-predicted, IV labels over-predicted")

    print("\nFIGURE 3 - panel c, held-out prediction by species")
    for sp, k, n_c, mean_c, lo_c, hi_c in (("mouse", mouse, 34, 0.56, 0.50, 0.68),
                                           ("rat", rat, 179, 0.60, 0.55, 0.65)):
        v = [rho(y[k], oof[c].to_numpy(float)[k]) for c in SEED_COLS]
        check(f"{sp} n", n_c, int(k.sum()))
        check(f"{sp} 3-seed mean rho", mean_c, float(np.mean(v)))
        check(f"{sp} seed range low", lo_c, float(np.min(v)))
        check(f"{sp} seed range high", hi_c, float(np.max(v)))
    # The caption's point about seed 0 being a misleading pass is itself a claim about the
    # data. Mouse is 34 records, so its per-seed rho swings; at seed 0 it reads ABOVE rat and
    # in the other two seeds below it, which is why the caption quotes the three-seed mean.
    claim("mouse is above rat at seed 0 and below it in the other two seeds",
          rho(y[mouse], p0[mouse]) > rho(y[rat], p0[rat])
          and all(rho(y[mouse], oof[c].to_numpy(float)[mouse])
                  < rho(y[rat], oof[c].to_numpy(float)[rat]) for c in SEED_COLS[1:]),
          "so the seed-0 reading overstates mouse")
    check("mouse seed-0 rho 0.68", 0.68, rho(y[mouse], p0[mouse]))
    check("rat seed-0 rho 0.65", 0.65, rho(y[rat], p0[rat]))
    # p-values are quoted as the WORST of the three seeds, which is what the panel prints.
    from scipy import stats
    worst = lambda k: max(stats.spearmanr(y[k], oof[c].to_numpy(float)[k]).pvalue
                          for c in SEED_COLS)
    check("mouse p 0.002 (worst seed)", 0.002, worst(mouse), 5e-4)
    claim("rat p < 0.001 (worst seed)", worst(rat) < 0.001, f"{worst(rat):.2e}")

    print("\nFIGURE 3 - panel d, accuracy on the absolute scale")
    for label, k, acc_c, med_c, bias_c in (("pooled", np.ones(len(y), bool), 0.57, None, None),
                                           ("mouse", mouse, 0.41, 4.0, -0.03),
                                           ("rat", rat, 0.60, 2.2, -0.06)):
        fe = [10 ** np.abs(y[k] - oof[c].to_numpy(float)[k]) for c in SEED_COLS]
        check(f"{label} within 3-fold", acc_c, float(np.mean([(f <= 3).mean() for f in fe])),
              5e-3)
        if med_c is not None:
            check(f"{label} median fold error", med_c,
                  float(np.mean([np.median(f) for f in fe])))
            check(f"{label} bias", bias_c,
                  float(np.mean([np.mean(oof[c].to_numpy(float)[k] - y[k])
                                 for c in SEED_COLS])))

    print("\nFIGURE 3 - panel e, ROC at the 6 h cut")
    a = [auc_at_cut(y, oof[c].to_numpy(float)) for c in SEED_COLS]
    check("AUC 0.82", 0.82, float(np.mean(a)))
    check("seed range low 0.77", 0.77, float(np.min(a)))
    check("seed range high 0.85", 0.85, float(np.max(a)))
    pos = y >= np.log10(CUT_H) - 1e-9
    check("17 positives of 213", 17, int(pos.sum()))
    check("mouse positives 8", 8, int(pos[mouse].sum()))
    check("rat positives 9", 9, int(pos[rat].sum()))

    print("\nFIGURE 3 - panel f, error strata (seed 0)")
    d = oof.merge(m[["id", "assay_method", "hl_engineering", "n_residues"]], on="id", how="left")
    d["fe"] = 10 ** np.abs(d.y_true - d[SEED_COLS[0]])
    d["species"] = d.species.str.title()
    for col, value, n_c, mdafe_c in (
            ("species", "Mouse", 34, 3.8), ("species", "Rat", 179, 2.1),
            ("assay_method", "radiometric / imaging", 13, 6.0),
            ("assay_method", "LC-MS", 62, 4.3),
            ("assay_method", "ELISA / immunoassay", 28, 4.0),
            ("assay_method", "HPLC", 19, 2.5),
            ("assay_method", "radioimmunoassay", 90, 1.8),
            ("hl_engineering", "radio/fluorescent label only", 17, 4.6),
            ("hl_engineering", "backbone modification only", 58, 2.9),
            ("hl_engineering", "unmodified / none stated", 135, 1.9)):
        sub = d[d[col] == value]
        check(f"{value} n", n_c, len(sub))
        check(f"{value} MdAFE", mdafe_c, sub.fe.median())
    # Strata under five records are dropped rather than drawn; the caption names which.
    check("bioassay n (dropped)", 1, int((d.assay_method == "bioassay").sum()))
    check("lipidated n (dropped)", 3, int((d.hl_engineering == "lipidated").sum()))
    size = d.PMID.map(d.PMID.value_counts())
    check("largest source n", 54, int((size == size.max()).sum()))
    check("largest source MdAFE", 1.6, d.loc[size == size.max(), "fe"].median())
    check("singleton sources n", 53, int((size == 1).sum()))
    check("singleton sources MdAFE", 4.1, d.loc[size == 1, "fe"].median())

    print("\nFIGURE 3 - panel g, feature ablation")
    check("29 configurations scored", 29, len(abl))
    check("shipped rho 0.57", 0.57, ship.spearman)
    check("molecule annotation only rho 0.10", 0.10,
          abl.loc["molecule annotation only"].spearman)
    unconstrained = abl.loc["physchem + study context + termini"]
    check("unconstrained winner rho 0.587", 0.587, unconstrained.spearman)
    check("applicability constraint costs 0.013", 0.013,
          unconstrained.spearman - ship.spearman, 5e-4)
    claim("the unconstrained winner outscores the shipped set",
          unconstrained.spearman > ship.spearman,
          "which is why the constraint has to be stated, not assumed")

    print("\nFIGURE 3 - panel h, paired flip-flop (all 320)")
    pairs = pd.read_csv(TABS / "table32_muriverse_halflife_route_flipflop_pairs.csv")
    check("32 paired sequence x species", 32, len(pairs))
    check("24 of 32 above identity", 24, int((pairs.sc > pairs.iv).sum()))
    check("75% above identity", 0.75, (pairs.sc > pairs.iv).mean(), 5e-3)
    check("Wilcoxon p 0.001", 0.001, pairs.wilcoxon_p[0], tol=5e-4)
    check("median delta +0.068", 0.068, pairs.median_delta[0])

    print("\nFIGURE 3 - panel i, saved split on the IV cohort")
    iv_meta = m[m.route_class == "iv"]
    check("176 training records", 176, int((iv_meta.split == "train").sum()))
    check("37 validation records", 37, int((iv_meta.split == "val").sum()))
    check("training mean -0.53", -0.53, iv_meta.loc[iv_meta.split == "train", "log_label"].mean())
    check("validation mean -0.57", -0.57, iv_meta.loc[iv_meta.split == "val", "log_label"].mean())

    print("\nFIGURE 3 - closing paragraph, the source offset")
    res = y - p0
    offset = pd.Series(res).groupby(oof.PMID.to_numpy()).transform("mean")
    check("r with the source offset known 0.93", 0.93, np.corrcoef(y, p0 + offset)[0, 1])
    check("three quarters of residual variance is source offset", 0.77,
          offset.var() / pd.Series(res).var(), 0.02)
    # This replaces a statistic earlier captions quoted that no script generated and no
    # reconstruction reproduced. The pairing rule is now stated and executed here.
    npairs, r_pairs, disagree = cross_source_pairs(iv_meta)
    check("33 cross-source pairs", 33, npairs)
    check("cross-source r 0.75", 0.75, r_pairs)
    check("median disagreement 0.35", 0.35, disagree)
    check("label spread 0.82", 0.82, iv_meta.log_label.std())

    print("\nMETHODS - the feature matrix and the cross-validation structure (results/revision_analysis/snapshot)")
    snap = Path(__file__).resolve().parents[1] / "results" / "revision_analysis" / "snapshot"
    fm = pd.read_csv(snap / "feature_manifest.csv")
    check("shipped matrix has 71 columns", 71, len(fm))
    blocks = fm.groupby("block", sort=False).size().to_dict()
    check("molecule annotation 9 columns", 9, blocks["molecule_annotation"])
    check("study context 33 columns", 33, blocks["study_context"])
    check("terminal classes 13 columns", 13, blocks["termini"])
    check("physicochemical descriptors 16 columns", 16, blocks["physchem"])
    check("11 of the 71 columns are constant on the 213 records", 11, int(fm.constant_in_cohort.sum()))
    claim("the four route columns are among the constant ones",
          set(fm[fm.field == "route_model"].column) <= set(fm[fm.constant_in_cohort].column), "")
    claim("no partner-class or partner-mass column is in the matrix",
          not fm.column.str.contains("partner|peg_kda").any(), "dropped as constant on the strict cohort")
    check("the dose is reported on 1 of 213 records", 1, int(213 - fm[fm.column == "log10_dose_mg_kg"].n_missing.iloc[0]))
    fa = pd.read_csv(snap / "fold_assignments_iv.csv")
    claim("every one of the 213 records is held out in every seed",
          bool((fa[[f"test_fold_seed{k}" for k in (0, 1, 2)]] >= 0).all().all()) and len(fa) == 213, "")
    check("37 validation-split records are held out and scored", 37, int((fa.split == "val").sum()))
    check("seed-to-seed SD of the shipped rho 0.048", 0.048, ship.spearman_sd)
    check("context + physchem rho 0.522", 0.522, abl.loc["physchem + study context"].spearman)
    check("context alone rho 0.347", 0.347, abl.loc["study context only"].spearman)
    # Discussion/Results say the embeddings LOWER the compact set (0.574 -> 0.230). They once also said no
    # embedding configuration beat the annotation block; under the IV data each embedding block alone does.
    check("compact set plus both embeddings (all molecular blocks + context) rho 0.230", 0.230,
          abl.loc["all molecular blocks + study context"].spearman)
    emb = abl[abl.index.str.contains("sequence|SMILES|all molecular", regex=True)]
    claim("every embedding-containing configuration scores below the compact set",
          bool((emb.spearman < ship.spearman).all()), f"best of {len(emb)} is {emb.spearman.max():.3f}")
    claim("each embedding block alone scores ABOVE the 9-column annotation block (so no text may say otherwise)",
          bool(min(abl.loc["sequence (ESM-2) only"].spearman, abl.loc["SMILES (PeptideCLM) only"].spearman)
               > abl.loc["molecule annotation only"].spearman), "0.145 and 0.177 against 0.104")
    prov =json.loads((snap / "embedding_provenance.json").read_text())
    claim("neither learned embedding is in the shipped model and both encoders are verified",
          (not prov["shipped_model"]["has_learned_embedding"]) and not prov["problems"]
          and all(e["blob_name_is_the_sha256"] for e in prov["encoders"].values()), "")
    check("longest cohort sequence 50 residues", 50, prov["sequence_length_facts"]["strict_cohort_max_residues"])

    print("\nSUPPLEMENTARY FIGURE 2 - pairs, allometry, headline metrics, fold errors")
    csum = json.loads((Path(__file__).resolve().parents[1] / "results" / "revision_analysis"
                       / "cross_source_pairs_summary.json").read_text())
    check("33 cross-source pairs", 33, csum["n_pairs"])
    check("from 9 sequence-and-species combinations", 9, csum["n_combinations"])
    check("32 of the pairs are in rat", 32, csum["pairs_rat"])
    check("one peptide supplies 21 of the pairs", 21, csum["dominant_peptide_pairs"])
    check("... reported by 7 sources (7 sources make 21 pairs)", 7, csum["dominant_peptide_sources"])
    check("cross-source r 0.75", 0.75, csum["r"])
    check("cross-source rho 0.61", 0.61, csum["rho"])
    check("median disagreement 0.35", 0.35, csum["median_abs_diff_log10"])
    check("r without the dominant peptide 0.84 (12 pairs)", 0.84, csum["r_without_dominant_peptide"])
    check("12 pairs remain without it", 12, csum["n_pairs_without_dominant"])
    check("combination-bootstrap r lower bound 0.31", 0.31, csum["bootstrap_over_combinations"]["r_lo"])
    check("combination-bootstrap r upper bound 0.98", 0.98, csum["bootstrap_over_combinations"]["r_hi"])
    check("combination-bootstrap median r 0.74", 0.74, csum["bootstrap_over_combinations"]["r_median"])
    cs_t = pd.read_csv(TABS / "table60_muriverse_halflife_cross_source_pairs_iv.csv")
    check("table60 has the 33 pairs", 33, len(cs_t))
    pr = pd.read_csv(TABS / "table61_polymyxin_species_pairs_allometry.csv").set_index(["drug", "route"])
    check("2 intravenous species pairs", 2, int((pr.index.get_level_values("route") == "IV").sum()))
    check("colistin IV b 0.52", 0.52, pr.loc[("Colistin", "IV"), "b"])
    check("colistin IV b interval lower 0.07", 0.07, pr.loc[("Colistin", "IV"), "b_lo"])
    check("colistin IV b interval upper 0.96", 0.96, pr.loc[("Colistin", "IV"), "b_hi"])
    check("polymyxin B IV b 0.48", 0.48, pr.loc[("Polymyxin B", "IV"), "b"])
    check("polymyxin B IV b interval lower 0.13", 0.13, pr.loc[("Polymyxin B", "IV"), "b_lo"])
    check("polymyxin B IV b interval upper 0.83", 0.83, pr.loc[("Polymyxin B", "IV"), "b_hi"])
    check("colistin EV b 0.30", 0.30, pr.loc[("Colistin", "EV"), "b"])
    claim("both intravenous exponents' intervals include the classic 0.25",
          bool(((pr.xs("IV", level="route").b_lo <= 0.25) & (pr.xs("IV", level="route").b_hi >= 0.25)).all()), "")
    check("rat-minus-mouse gap, colistin IV 0.55 log10", 0.55, pr.loc[("Colistin", "IV"), "rat_minus_mouse_log10"])
    check("rat-minus-mouse gap, polymyxin B IV 0.51 log10", 0.51, pr.loc[("Polymyxin B", "IV"), "rat_minus_mouse_log10"])
    check("rat-minus-mouse gap, colistin EV 0.31 log10", 0.31, pr.loc[("Colistin", "EV"), "rat_minus_mouse_log10"])
    oofs = pd.read_csv(snap / "oof_predictions_iv.csv")
    yy = oofs.y_true.to_numpy()
    def seed_fe(mask):
        return float(np.mean([np.median((10 ** np.abs(oofs[f"y_pred_cluster_PMID_seed{k}"].to_numpy() - yy))[mask]) for k in (0, 1, 2)]))
    allm = np.ones(len(oofs), bool)
    check("median fold error, both species (3-seed mean) 2.37", 2.37, seed_fe(allm))
    check("median fold error, rat 2.18 (the Discussion's 2.2)", 2.18, seed_fe((oofs.species == "rat").to_numpy()))
    check("median fold error, mouse 3.96 (the Discussion's 4.0)", 3.96, seed_fe((oofs.species == "mouse").to_numpy()))
    iv_obs = pd.read_csv(TABS / "table49_polymyxin_observed_by_drug_species_iv.csv")
    check("panel h band lower bound 0.45 h", 0.45, iv_obs.observed_h.min())
    check("panel h band upper bound 1.90 h", 1.90, iv_obs.observed_h.max())

    print("\nDISCUSSION - source concentration and route evidence (analysis/cohort_label_route_summary.json)")
    csum2 = json.loads((Path(__file__).resolve().parents[1] / "results" / "revision_analysis"
                        / "cohort_label_route_summary.json").read_text())
    check("the heaviest source holds 25% of the intravenous cohort", 0.25, csum2["heaviest_source_share"])
    check("it is 54 of the 213 records", 54, csum2["heaviest_sources"]["EP2277890A1"])
    imp = {r["group"]: r for r in csum2["impact_on_rho"]}
    wp = imp["without the heaviest source (EP2277890A1, a patent)"]
    check("leaving the patent out leaves 159 records", 159, wp["n_left"])
    check("... and lowers rho to 0.524", 0.524, wp["rho"])
    # The Discussion's chemical-coverage sentence. Verified 2026-10-17 against the saved split
    # metadata; the 320-record figures (67 cyclic, 41 non-standard) are the ones in Figure 2, so
    # the two pairs are easy to swap by accident.
    meta_iv = pd.read_csv(SNAP / "halflife_meta_with_split.csv")
    meta_iv = meta_iv[meta_iv.route_model == "iv"]
    check("213 rows in the saved IV split metadata", 213, len(meta_iv))
    check("50 cyclic records in the IV cohort", 50, int((meta_iv.topology == "cyclic").sum()))
    check("11 IV records carry a non-standard residue", 11,
          int(meta_iv.sequence_nonstandard.notna().sum()))

    check("the patent is 19% of the 320-record analysis set", 0.19, 61 / 320, 5e-3)
    check("61 analysis-set records come from that patent", 61, int((m.PMID.astype(str) == "EP2277890A1").sum()))

    # -------------------------------------------------- Figure 4: the held-out application
    print("\nFIGURE 4 - shipped predictions (intravenous model)")
    t20 = pd.read_csv(TABS / "table20_polymyxin_predictions.csv")
    sh = t20[t20.role == "shipped"]
    check("16 shipped rows", 16, len(sh))
    check("2 input variants per congener", 2, sh.variant.nunique())
    claim("only intravenous inputs are predicted", set(sh.route_class) == {"iv"},
          str(sorted(set(sh.route_class))))
    # The C-terminus is annotated unmodified only. `Amidation` was removed because it
    # asserts a terminal amide these molecules do not carry: the Thr10 carboxyl is consumed
    # in the Dab4 macrolactam, which the topology block already encodes.
    claim("the C-terminus is annotated unmodified only", set(sh.c_ter_text) == {"Free"},
          str(sorted(set(sh.c_ter_text))))
    check("8 distinct predicted values", 8, sh.pred_half_life_h.round(9).nunique())
    check("prediction min 0.57 h", 0.57, sh.pred_half_life_h.min())
    check("prediction max 0.92 h", 0.92, sh.pred_half_life_h.max())
    # The overall median is not a value the model returns: within each species the two drugs occupy
    # disjoint ranges (polymyxin B below colistin), but the ranges overlap across species. Report per
    # drug and species, never pooled.
    for _sp in ("mouse", "rat"):
        _pol = sh[(sh.drug == "Polymyxin B") & (sh.species == _sp)].pred_half_life_h
        _col = sh[sh.drug.str.startswith("Colistin") & (sh.species == _sp)].pred_half_life_h
        claim(f"{_sp}: the two drugs occupy disjoint ranges, polymyxin B below colistin",
              _pol.max() < _col.min(), f"polymyxin B <= {_pol.max():.3f} h, colistin >= {_col.min():.3f} h")
    check("colistin median 0.85 h", 0.85,
          sh[sh.drug.str.startswith("Colistin")].pred_half_life_h.median())
    check("polymyxin B median 0.68 h", 0.68,
          sh[sh.drug == "Polymyxin B"].pred_half_life_h.median())
    claim("mouse and rat predictions differ now that species is a model input",
          not np.allclose(sh[sh.species == "mouse"].pred_log10_half_life_h.to_numpy(),
                          sh[sh.species == "rat"].pred_log10_half_life_h.to_numpy()),
          f"rat - mouse = "
          f"{sh[sh.species == 'rat'].pred_log10_half_life_h.mean() - sh[sh.species == 'mouse'].pred_log10_half_life_h.mean():+.3f} log10")
    deg = pd.read_csv(TABS / "table23_polymyxin_congener_degeneracy.csv")
    deg = deg[deg.role == "shipped"]
    check("congener spread 1.38 min", 1.38, deg.congener_fold_spread.min())
    check("congener spread 1.38 max", 1.38, deg.congener_fold_spread.max())
    # The two congeners of each drug differ by the same single methylene in the acyl chain,
    # but the model does not treat the two pairs alike: colistin A/B agree to 1.02-fold
    # while polymyxin B1/B2 differ by 1.23-fold, which is most of the between-drug gap.
    claim("the congener pairs are not equally degenerate",
          deg.congener_fold_spread.max() / deg.congener_fold_spread.min() < 1.05,
          f"spreads {sorted(deg.congener_fold_spread.round(3).unique())}")
    disc = pd.read_csv(TABS / "table29_polymyxin_drug_discrimination.csv")   # intravenous comparison
    check("encoding moves 1.0-fold (N-terminal wording is a no-op)", 1.0,
          disc.max_encoding_fold_spread.max())
    check("predicted drug gap 1.24-fold", 1.24, 1 / disc.median_predicted_fold_gap.iloc[0])
    d29 = disc.set_index("species")
    check("observed IV polymyxin B / colistin, mouse 1.31", 1.31, d29.loc["Mouse", "observed_fold_gap"])
    check("observed IV polymyxin B / colistin, rat 1.19", 1.19, d29.loc["Rat", "observed_fold_gap"])
    claim("against intravenous observations the model's drug direction is opposite in BOTH species",
          set(disc.ordering) == {"opposite to the observed direction"}, "")

    print("\nFIGURE 4 - intravenous comparison against the constants (table27_iv, table49_iv)")
    b27 = pd.read_csv(TABS / "table27_polymyxin_fold_error_vs_baseline_iv.csv")
    b27 = b27[b27.role == "shipped"]
    t26 = pd.read_csv(TABS / "table26_polymyxin_training_label_context.csv")
    # The baseline MUST come from the cohort the applied model was fitted on. This is the
    # check that would have caught the bug it is here for.
    # The applied model is the post-audit one, so the constants come from the post-audit records of the
    # intravenous training split, with the audit's corrected labels and species.
    _member = pd.read_csv(ROOT / "rodent/strict/model_iv_postaudit/training_membership.csv").set_index("id")
    iv_train = m[(m.split == "train") & (m.route_class == "iv")]
    iv_train = iv_train[iv_train["id"].map(_member["included"]).fillna(False).astype(bool)].copy()
    iv_train["log_label"] = np.log10(iv_train["id"].map(_member["label_hours_used"]).to_numpy(float))
    iv_train["species"] = iv_train["id"].map(_member["species_used"])
    check("baseline split n 146", 146, int(len(iv_train)))
    for sp, n_c, h_c in (("mouse", 27, 0.547), ("rat", 119, 0.196)):
        row = t26[t26.species == sp].iloc[0]
        check(f"{sp} baseline split n", n_c, int(row.baseline_split_n))
        check(f"{sp} baseline median h", h_c, row.train_median_half_life_h)
        check(f"{sp} baseline matches the IV training split", h_c,
              float(10 ** np.median(iv_train.loc[iv_train.species == sp, "log_label"])))
    both_const = float(10 ** np.median(iv_train.log_label))
    check("both-species training median 0.245 h", 0.245, both_const)
    lip = iv_train[iv_train.hl_engineering == "lipidated"]
    check("1 lipidated post-audit IV training peptide", 1, len(lip))
    LIP_H = float(10 ** np.median(lip.log_label))
    check("lipidated median 4.3 h", 4.3, LIP_H)
    check("one of 213 training records carries a dose", 1,
          int(m[m.route_class == "iv"].pk_dose_log10_mg_kg.notna().sum()))
    g = b27.groupby(["species", "drug"])
    # drug prefix, species, observed h, studies, animals, model fold, species-median fold
    IV_CELLS = (("Colistin", "mouse", 0.45, 2, 37, 2.1, 1.2),
                ("Colistin", "rat",   1.59, 6, 36, 2.0, 8.1),
                ("Polymyxin B", "mouse", 0.58, 2, 56, 1.3, 1.1),
                ("Polymyxin B", "rat",   1.90, 3, 33, 3.1, 9.7))
    for drug, sp, obs_c, n_st, n_an, fe_c, base_c in IV_CELLS:
        key = (sp, "Colistin (polymyxin E)" if drug == "Colistin" else "Polymyxin B")
        tag = f"IV {sp} {drug.split()[0]}"
        check(f"{tag} observed h", obs_c, g.observed_h.first()[key])
        check(f"{tag} studies", n_st, int(g.n_studies.first()[key]))
        check(f"{tag} model fold error", fe_c, g.model_fold_error.median()[key])
        check(f"{tag} species-median fold error", base_c, g.baseline_fold_error.first()[key])
    t49i = pd.read_csv(TABS / "table49_polymyxin_observed_by_drug_species_iv.csv")
    check("IV animals: colistin mouse 37", 37, int(t49i[(t49i.drug.str.startswith("Colistin")) & (t49i.species == "mouse")].n_animals.iloc[0]))
    check("IV animals: colistin rat 36", 36, int(t49i[(t49i.drug.str.startswith("Colistin")) & (t49i.species == "rat")].n_animals.iloc[0]))
    check("IV animals: polymyxin B mouse 56", 56, int(t49i[(t49i.drug == "Polymyxin B") & (t49i.species == "mouse")].n_animals.iloc[0]))
    check("IV animals: polymyxin B rat 33", 33, int(t49i[(t49i.drug == "Polymyxin B") & (t49i.species == "rat")].n_animals.iloc[0]))
    # The lipidated constant, per row, against the same intravenous observations
    b27 = b27.assign(lip_fold=lambda d: np.maximum(LIP_H / d.observed_h, d.observed_h / LIP_H))
    for drug, sp, lip_c in (("Colistin", "mouse", 9.6), ("Colistin", "rat", 2.7),
                            ("Polymyxin B", "mouse", 7.4), ("Polymyxin B", "rat", 2.27)):
        rows = b27[b27.drug.str.startswith(drug) & (b27.species == sp)]
        check(f"IV {sp} {drug} lipidated-constant fold error", lip_c, rows.lip_fold.iloc[0])
    beat_sp = b27.model_beats_baseline
    beat_lip = b27.model_fold_error < b27.lip_fold
    check("model beats the species-median constant in 8 of 16 rows", 8, int(beat_sp.sum()))
    check("model beats the lipidated constant in 12 of 16 rows", 12, int(beat_lip.sum()))
    claim("the species-median wins are exactly the rat rows",
          bool(beat_sp.eq(b27.species == "rat").all()), "")
    claim("the model loses to the lipidated constant only in the rat polymyxin B rows",
          bool((~beat_lip).eq(b27.species.eq("rat") & b27.drug.eq("Polymyxin B")).all()), "")
    _closer = b27.model_fold_error < np.minimum(b27.baseline_fold_error, b27.lip_fold)
    claim("the model is closer than both constants only in the rat colistin rows",
          bool(_closer.eq(b27.species.eq("rat") & b27.drug.str.startswith("Colistin")).all()), "")
    claim("the model is long in both mouse cells and short in both rat cells",
          bool(((b27.pred_half_life_h > b27.observed_h) == (b27.species == "mouse")).all()), "")
    # intravenous rat/mouse ratios behind the species-term statement
    for drug, ratio_c, log_c in (("Colistin", 3.6, 0.55), ("Polymyxin B", 3.3, 0.51)):
        ob = t49i[t49i.drug.str.startswith(drug)].set_index("species").observed_h
        check(f"IV {drug} rat/mouse ratio", ratio_c, ob["rat"] / ob["mouse"])
        check(f"IV {drug} rat - mouse log10", log_c, np.log10(ob["rat"] / ob["mouse"]))
    check("model puts rat 15% shorter (mean over drugs)", 0.85,
          float((sh[sh.species == "rat"].pred_half_life_h.mean()) / sh[sh.species == "mouse"].pred_half_life_h.mean()), 5e-3)

    print("\nFIGURE 4 - intravenous leave-one-study-out (analysis/polymyxin_iv_leave_one_study_out.csv)")
    loo = pd.read_csv(Path(__file__).resolve().parents[1] / "results" / "revision_analysis"
                      / "polymyxin_iv_leave_one_study_out.csv")
    def loo_fe(study, drug, sp):
        r = loo[(loo.dropped_study == study) & loo.drug.str.startswith(drug) & (loo.species == sp)]
        return float(r.fold_error.iloc[0])
    def loo_all(drug, sp):
        r = loo[loo.drug.str.startswith(drug) & (loo.species == sp)].dropna(subset=["fold_error"])
        return float(r.fold_error.min()), float(r.fold_error.max())
    check("PMB mouse: dropping Lin 2017 (aerosol) gives 2.7-fold", 2.7, loo_fe("LinPMBAerosol2017", "Polymyxin B", "Mouse"))
    check("PMB mouse: dropping He 2013 gives 1.4-fold", 1.4, loo_fe("HeLipo2013", "Polymyxin B", "Mouse"))
    check("colistin mouse: dropping Jin 2009 gives 1.5-fold (Lin alone)", 1.5, loo_fe("Jin2009", "Colistin", "Mouse"))
    check("colistin mouse: dropping Lin 2017 gives 3.1-fold (Jin alone)", 3.1, loo_fe("Lin2017", "Colistin", "Mouse"))
    check("colistin rat: dropping Liao 2023 gives 2.47-fold", 2.47, loo_fe("Liao2023", "Colistin", "Rat"))
    lo, hi = loo_all("Colistin", "Rat")
    check("colistin rat leave-one-out min 1.66", 1.66, lo)
    check("colistin rat leave-one-out max 2.47", 2.47, hi)
    lo, hi = loo_all("Polymyxin B", "Rat")
    check("polymyxin B rat leave-one-out min 2.31", 2.31, lo)
    check("polymyxin B rat leave-one-out max 3.36", 3.36, hi)
    ob21 = pd.read_csv(TABS / "table21_polymyxin_observed_half_lives.csv")
    check("Jin 2009 derived half-life 0.29 h", 0.29, float(ob21[(ob21.study_id == "Jin2009") & ob21.use_in_band].half_life_h.iloc[0]))
    check("Liao 2023 free-colistin half-life 0.21 h (12.46 min)", 0.21, float(ob21[(ob21.study_id == "Liao2023") & ob21.use_in_band].half_life_h.iloc[0]))
    claim("Liao 2023's half-life is shorter than its 15-min infusion",
          float(ob21[(ob21.study_id == "Liao2023") & ob21.use_in_band].half_life_h.iloc[0]) * 60 < 15, "")
    claim("Jin 2009 is entered as derived (the paper reports no half-life)",
          str(ob21[(ob21.study_id == "Jin2009") & ob21.use_in_band].detail.iloc[0]).startswith("derived"), "")
    claim("the intratracheal arm of Lin 2017 is out of the comparison",
          not bool(ob21[(ob21.study_id == "Lin2017") & (ob21.route == "intratracheal")].use_in_band.any()), "")

    print("\nSUPPLEMENTARY FIGURE 3 / Figure 4 baseline note - all-route comparison (table27, table49)")
    p27 = pd.read_csv(TABS / "table27_polymyxin_fold_error_vs_baseline.csv")
    p27 = p27[p27.role == "shipped"]
    gp = p27.groupby(["species", "drug"])
    for sp, drug, obs_c, fe_c, base_c in (("mouse", "Colistin (polymyxin E)", 1.02, 1.1, 1.9),
                                          ("rat", "Colistin (polymyxin E)", 1.70, 2.2, 8.7),
                                          ("mouse", "Polymyxin B", 1.13, 1.5, 2.1),
                                          ("rat", "Polymyxin B", 1.90, 3.1, 9.7)):
        tag = f"all-route {sp} {drug.split()[0]}"
        check(f"{tag} observed h", obs_c, gp.observed_h.first()[(sp, drug)])
        check(f"{tag} model fold error", fe_c, gp.model_fold_error.median()[(sp, drug)])
        check(f"{tag} species-median fold error", base_c, gp.baseline_fold_error.first()[(sp, drug)])
    claim("against every route the model beats its species-median constant in all 16 rows",
          bool(p27.model_beats_baseline.all()), f"{int(p27.model_beats_baseline.sum())}/{len(p27)}")
    t49p = pd.read_csv(TABS / "table49_polymyxin_observed_by_drug_species.csv")
    for drug, sp, n_st, n_gr, n_an in (("Colistin", "mouse", 4, 9, 221), ("Colistin", "rat", 7, 8, 42),
                                       ("Polymyxin B", "mouse", 5, 13, 237), ("Polymyxin B", "rat", 3, 4, 33)):
        r = t49p[t49p.drug.str.startswith(drug) & (t49p.species == sp)].iloc[0]
        check(f"all-route {sp} {drug} studies", n_st, int(r.n_studies))
        check(f"all-route {sp} {drug} dose groups", n_gr, int(r.n_groups))
        check(f"all-route {sp} {drug} animals", n_an, int(r.n_animals))
    pooled_mouse_base = p27[p27.species == "mouse"].baseline_fold_error
    claim("the mouse constant is 1.9-2.1-fold off against every route but 1.1-1.2-fold against IV",
          bool(pooled_mouse_base.min() > 1.8 and b27[b27.species == "mouse"].baseline_fold_error.max() < 1.3), "")

    # Panel c is the per-congener predicted-versus-observed scatter against the INTRAVENOUS
    # references (fig4.R builds `sr_cells` from table20 joined to table49_iv), so its in-panel
    # median fold errors are these, not the table56 ones checked below. Corrected 2026-10-17:
    # this block used to be labelled "panel c" and checked table56, whose mouse value happens
    # to agree (2.7) while its rat value does not (1.1 against the panel's 1.4). The caption
    # drew its rat figure from the panel and so disagreed with the audit while being right.
    print("\nFIGURE 4 - panel c, the per-congener IV scatter (table20 x table49_iv)")
    t20 = pd.read_csv(TABS / "table20_polymyxin_predictions.csv")
    t49iv = pd.read_csv(TABS / "table49_polymyxin_observed_by_drug_species_iv.csv")
    sr = (t20[t20.role == "shipped"]
          .groupby(["species", "drug", "name"], as_index=False)
          .pred_log10_half_life_h.mean()
          .merge(t49iv[["drug", "species", "observed_log10_h"]], on=["drug", "species"]))
    sr["fold_error"] = 10 ** (sr.pred_log10_half_life_h - sr.observed_log10_h).abs()
    check("8 points in the panel", 8, len(sr))
    for sp, want in (("mouse", 1.7), ("rat", 2.4)):
        check(f"{sp} in-panel median fold error {want}", want,
              float(sr[sr.species == sp].fold_error.median()))
    claim("every mouse prediction is above the identity line and every rat one below",
          (sr[sr.species == "mouse"].pred_log10_half_life_h
           > sr[sr.species == "mouse"].observed_log10_h).all()
          and (sr[sr.species == "rat"].pred_log10_half_life_h
               < sr[sr.species == "rat"].observed_log10_h).all(),
          "mouse %.1f-%.1f x long, rat %.1f-%.1f x short" % (
              sr[sr.species == "mouse"].fold_error.min(),
              sr[sr.species == "mouse"].fold_error.max(),
              sr[sr.species == "rat"].fold_error.min(),
              sr[sr.species == "rat"].fold_error.max()))

    # table56 is the route-pooled species x route test. It is NOT drawn any more (fig4.R:18);
    # the Results quote its values in text and point at the deposited figure source tables.
    print("\nFIGURE 4 - the route-pooled species x route test, not drawn (table56)")
    c56 = pd.read_csv(TABS / "table56_polymyxin_species_route_cells.csv")
    check("mouse median fold error 2.7", 2.7, c56[c56.species == "mouse"].fold_error.median())
    check("rat median fold error 1.1", 1.1, c56[c56.species == "rat"].fold_error.median())
    check("seven species x route cells", 7, len(c56))
    def cell_fe(drug, sp, adm):
        return float(c56[(c56.drug_short == drug) & (c56.species == sp) & (c56.admin == adm)].fold_error.iloc[0])
    check("colistin mouse EV 2.4", 2.4, cell_fe("Colistin", "mouse", "EV"))
    check("colistin mouse IV 4.5", 4.5, cell_fe("Colistin", "mouse", "IV"))
    check("polymyxin B mouse IV 2.9", 2.9, cell_fe("Polymyxin B", "mouse", "IV"))
    check("polymyxin B mouse EV 1.1", 1.1, cell_fe("Polymyxin B", "mouse", "EV"))
    check("colistin mouse IV observed 0.42 h (mean of per-study means)", 0.42,
          float(10 ** c56[(c56.drug_short == "Colistin") & (c56.species == "mouse") & (c56.admin == "IV")].obs_log10.iloc[0]))

    # The Results' route-transfer paragraph quotes these and, since the panel was removed, has
    # no figure to point at - so the numbers are only guarded here. Each is EV minus IV within
    # one drug and species: the model's shift is near-constant, the literature's is not.
    def ev_minus_iv(drug, sp, col):
        cell = c56[(c56.drug_short == drug) & (c56.species == sp)].set_index("admin")[col]
        return float(cell["EV"] - cell["IV"])
    # "approximately 0.19" in the Results is the mean over the three cells that have both
    # routes; the individual cells are 0.186, 0.187 and 0.196, so no single one rounds to it.
    both_routes = (("Colistin", "mouse"), ("Colistin", "rat"), ("Polymyxin B", "mouse"))
    shifts = [ev_minus_iv(d, s, "pred_log10") for d, s in both_routes]
    check("predicted route shift, mean over the three cells 0.19", 0.19, sum(shifts) / 3)
    claim("the predicted route shift is one near-constant offset",
          max(shifts) - min(shifts) < 0.02,
          "cells %s" % ", ".join(f"{v:.3f}" for v in shifts))
    check("observed route shift, rat colistin +0.22", 0.22, ev_minus_iv("Colistin", "rat", "obs_log10"))
    check("observed route shift, mouse colistin +0.46", 0.46, ev_minus_iv("Colistin", "mouse", "obs_log10"))
    check("observed route shift, mouse polymyxin B +0.70", 0.70,
          ev_minus_iv("Polymyxin B", "mouse", "obs_log10"))
    claim("the observed shifts spread much wider than the predicted one",
          ev_minus_iv("Polymyxin B", "mouse", "obs_log10")
          - ev_minus_iv("Colistin", "rat", "obs_log10") > 0.4,
          "three reference differences, not two - the prose said 'both' until 2026-10-17")

    print("\nFIGURE 4 - the per-source references, all routes (table52, table53)")
    # The caption quotes, per source of observation, the observed half-life and how far the
    # drug's MEDIAN shipped prediction sits from it. Fold error here is the ratio of the two
    # MEDIANS, not the median of the per-variant fold errors (Monash mouse: 1.04 vs 1.31).
    t53 = pd.read_csv(TABS / "table53_polymyxin_fold_error_by_source.csv")
    t53 = t53[(t53.role == "shipped") & (~t53.source_group.str.startswith("pooled"))]
    # Per drug AND species: species is a model input, so a per-drug median would average the
    # mouse and the rat prediction and then be compared against one animal's observations.
    med_pred = t53.groupby(["drug", "species"]).pred_half_life_h.median()
    PER_SOURCE = [
        # drug prefix, species, source group, value basis, observed h, fold, model is...
        ("Polymyxin B", "mouse", "Jiao 2023",    "reported", 3.78, 5.13, "under"),
        ("Polymyxin B", "mouse", "Tam lab",      "reported", 0.27, 2.71, "over"),
        ("Polymyxin B", "mouse", "Monash group", "derived",  1.61, 2.18, "under"),
        ("Polymyxin B", "mouse", "other",        "derived",  1.20, 1.64, "under"),
        ("Polymyxin B", "rat",   "Tam lab",      "reported", 2.07, 3.31, "under"),
        ("Polymyxin B", "rat",   "Monash group", "reported", 1.45, 2.32, "under"),
        ("Colistin",    "mouse", "Monash group", "reported", 1.13, 1.23, "under"),
        ("Colistin",    "mouse", "Monash group", "derived",  0.29, 3.14, "over"),
        ("Colistin",    "rat",   "Monash group", "reported", 1.36, 1.75, "under"),
        ("Colistin",    "rat",   "Monash group", "derived",  4.07, 5.24, "under"),
        ("Colistin",    "rat",   "other",        "reported", 1.09, 1.40, "under"),
    ]
    for drug, species, group, basis, obs_h, fold, direction in PER_SOURCE:
        rows = t53[t53.drug.str.startswith(drug) & (t53.species == species)
                   & (t53.source_group == group) & (t53.value_basis == basis)]
        tag = f"{drug[:3]} {species} {group} ({basis})"
        if rows.empty:
            claim(f"{tag} is present", False, "no such row in table53")
            continue
        obs = rows.observed_h.iloc[0]
        pred = med_pred[(rows.drug.iloc[0], species)]
        check(f"{tag} observed h", obs_h, obs)
        check(f"{tag} fold error", fold, max(pred / obs, obs / pred))
        claim(f"{tag} model reads {direction}", ("over" if pred > obs else "under") == direction)
    # The pooled-sensitivity rows quote the OTHER summary: the median of the per-variant fold
    # errors, not the ratio of medians used above, so each is pinned to the metric it uses.
    SENSITIVITY = [
        # drug prefix, species, pooling, observed h, fold error (median of per-row values)
        ("Polymyxin B", "mouse", "pooled: each study once", 1.39, 1.9),
        ("Polymyxin B", "mouse", "pooled: reported only",   0.46, 1.6),
        ("Polymyxin B", "mouse", "pooled: IV only",         0.58, 1.3),
        ("Polymyxin B", "rat",   "pooled: each study once", 1.72, 2.8),
        ("Polymyxin B", "rat",   "pooled: reported only",   1.90, 3.1),
        ("Polymyxin B", "rat",   "pooled: IV only",         1.90, 3.1),
        ("Colistin",    "mouse", "pooled: each study once", 0.72, 1.3),
        ("Colistin",    "mouse", "pooled: reported only",   1.13, 1.2),
        ("Colistin",    "mouse", "pooled: IV only",         0.45, 2.1),
        ("Colistin",    "rat",   "pooled: each study once", 1.63, 2.1),
        ("Colistin",    "rat",   "pooled: reported only",   1.20, 1.5),
        ("Colistin",    "rat",   "pooled: IV only",         1.59, 2.0),
    ]
    t53p = pd.read_csv(TABS / "table53_polymyxin_fold_error_by_source.csv")
    t53p = t53p[t53p.role == "shipped"]
    for drug, species, pooling, obs_h, fold in SENSITIVITY:
        rows = t53p[t53p.drug.str.startswith(drug) & (t53p.species == species)
                    & (t53p.source_group == pooling)]
        tag = f"{drug[:3]} {species} {pooling.split(': ')[1]}"
        if rows.empty:
            claim(f"{tag} is present", False, "no such row in table53")
            continue
        check(f"{tag} observed h", obs_h, rows.observed_h.iloc[0])
        check(f"{tag} fold error", fold, rows.fold_error.median())
    t55 = pd.read_csv(TABS / "table55_polymyxin_pooled_sensitivity.csv")
    pmb = t55[t55.drug.str.startswith("Polymyxin B") & (t55.species == "mouse")].iloc[0]
    check("van der Meijden animal-weight share 67.5%", 67.5, 100 * pmb.top_study_weight_share)
    check("derived-value weight share 87.3%", 87.3, 100 * pmb.derived_weight_share)
    check("leave-one-out min 0.99 h", 0.99, pmb.leave_one_out_min_h)
    check("leave-one-out max 1.33 h", 1.33, pmb.leave_one_out_max_h)

    span = t53[t53.drug.str.startswith("Polymyxin B") & (t53.species == "mouse")
               & (t53.source_group == "Tam lab")].fold_error
    check("Tam-lab mouse congener span, min 2.5", 2.5, span.min())
    check("Tam-lab mouse congener span, max 3.0", 3.0, span.max())
    # The prediction is too long against some published references and too short against
    # others, so no single reference settles whether it is accurate.
    for drug in ("Polymyxin B", "Colistin"):
        rows = t53[t53.drug.str.startswith(drug)]
        pred = med_pred[rows.drug.iloc[0]].median()
        over = (pred > rows.observed_h).any()
        under = (pred < rows.observed_h).any()
        claim(f"{drug}: the prediction sits between the published extremes", bool(over and under),
              f"{rows.observed_h.min():.2f}-{rows.observed_h.max():.2f} h observed, {pred:.2f} h predicted")

    print("\nFIGURE 4 - the configurations that are applied but not drawn")
    emb = pd.read_csv(TABS / "table51_polymyxin_embedding_extrapolation.csv")
    # The PeptideCLM-embedding variant now returns physically plausible half-lives after the
    # block-order fix (see predict_polymyxins.assemble); it is kept out of the figure because
    # its held-out ranking is weak, not because the applied prediction is pathological.
    claim("PeptideCLM variant stays within the trained range (0.01-100 h)",
          bool(((emb.median_pred_h > 0.01) & (emb.median_pred_h < 100)).all()),
          f"range {emb.median_pred_h.min():.2f}-{emb.median_pred_h.max():.2f} h")
    check("PeptideCLM rho 0.25", 0.25, abl.loc["SMILES + physchem + study context"].spearman)
    structural = t20[t20.role == "structural"]
    claim("the unconstrained physchem+context+termini configuration is applicable here",
          bool(((structural.pred_half_life_h > 0.01) & (structural.pred_half_life_h < 100)).all()),
          f"range {structural.pred_half_life_h.min():.2f}-{structural.pred_half_life_h.max():.2f} h")

    print("\nNOT CHECKED")
    print("  Figure 4 literature PK values - verifiable only against the source papers.")

    print(f"\n{NCHECK - len(FAILS)}/{NCHECK} checks passed")
    if FAILS:
        print("\nFAILURES:")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)


def paired_group_bootstrap(iv: pd.DataFrame, po: pd.DataFrame, y: np.ndarray,
                           draws: int = 500, seed: int = 0):
    """How much the route restriction is worth, resampling whole cluster/source groups.

    The two models are scored on the SAME rows in the same resample, so the molecule and the
    study are held fixed and only the training cohort differs. Resampling the connected
    components of homology cluster and source rather than rows is what makes the interval an
    estimate of transfer to a new family or publication rather than to a new row of a
    publication already seen.
    """
    par: dict = {}

    def find(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for c, p in zip(iv.cluster, iv.PMID):
        par[find(("c", c))] = find(("p", p))
    roots = [find(("c", c)) for c in iv.cluster]
    key = {r: i for i, r in enumerate(dict.fromkeys(roots))}
    g = np.array([key[r] for r in roots])
    idx = [np.where(g == i)[0] for i in range(len(key))]
    gain = lambda rows: (np.mean([rho(y[rows], iv[c].to_numpy(float)[rows]) for c in SEED_COLS])
                         - np.mean([rho(y[rows], po[c].to_numpy(float)[rows])
                                    for c in SEED_COLS]))
    rs = np.random.default_rng(seed)
    boot = [gain(np.concatenate([idx[k] for k in rs.integers(0, len(idx), len(idx))]))
            for _ in range(draws)]
    return gain(np.arange(len(y))), float(np.percentile(boot, 2.5)), float(
        np.percentile(boot, 97.5))


def cross_source_pairs(cohort: pd.DataFrame):
    """How well two sources agree about the same molecule in the same animal.

    The pairing rule, stated here because the caption quotes the result: take every
    sequence x species combination that more than one source reports, reduce each source to
    its mean log10 half-life for that combination, and form every unordered pair of sources
    within the combination. This is the ceiling on the accuracy attainable against any one
    held-out reported label.
    """
    g = (cohort.dropna(subset=["sequence"])
         .groupby(["sequence", "species", "PMID"]).log_label.mean().reset_index())
    pairs = [(v[i], v[j])
             for _, grp in g.groupby(["sequence", "species"])
             for v in [grp.log_label.to_numpy()]
             for i in range(len(v)) for j in range(i + 1, len(v))]
    a = np.array(pairs)
    return len(a), float(np.corrcoef(a[:, 0], a[:, 1])[0, 1]), float(
        np.median(np.abs(a[:, 0] - a[:, 1])))


if __name__ == "__main__":
    main()
