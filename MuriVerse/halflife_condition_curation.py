"""Primary-source curation of the study conditions for the rodent half-life cohort.

Two of the study-condition covariates read "not reported" or "unknown" far more often than
the sources warrant, because the database column they are derived from is blank. This
module curates both from the sources themselves, each against its own ledger in
`data/curation/`, and each ledger records the verbatim quote and where it came from.

DOSE ROUTE. PEPlife2 leaves `Dose Route` blank on 453 of the 1,159 eligible rodent records
(39%) and on 238 of the 467 strict records, so `route_class` reads "not reported" for those
rows even where the source states the route plainly. The database does often carry it, but only as prose inside `Test Sample`
("Rat plasma (Intravenous)"), which is a transcription rather than evidence.

This module closes the gap from the sources themselves. `data/curation/route_curation.csv`
records, for every source behind a "not reported" record, which routes the PRIMARY source
actually evidences, with the verbatim quote and where it came from (full text, patent text,
abstract, or NLM MeSH indexing). A record's route is then the database-stated route only
when the primary source evidences that route too; a record whose source evidences nothing,
or whose annotation the source contradicts, is left unresolved.

The cohort keeps parenterally dosed records only - iv, sc, ip, im. An oral, intranasal,
rectal or site-unstated record is dropped because its label is an absorption-limited
composite, not the elimination half-life the model is fitted on, and an unresolved record
is dropped because nothing establishes which of the two it is.

Evidence levels, worst last:
  primary_arm     the source names the route for THIS arm - a single-route study, or the
                  record's half-life matched to its row in the source's own PK table
  primary_source  the source names the routes it used and all of them are parenteral, so
                  the record is attributed at source level rather than arm level
  none            no primary source names a route for the PK arm

ASSAY METHOD. `assay_method` is bucketed from the free-text `Assay` column, which for some
records names a separation rather than a detector ("TLC"), names the readout rather than the
method ("Serum decay experiment", "Plasma stability assay"), or is simply "N.A.".
`assay_curation.csv` carries the method the source states, for the records whose `Assay`
text cannot yield one. Unlike the route, an unresolved assay is not a reason to drop a
record: the half-life is still a half-life, and `unknown` is an honest level for it.

Usage: imported by build_halflife_training_set.py; running this module prints both
reconciliation summaries.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROUTE_LEDGER = HERE / "data" / "curation" / "route_curation.csv"
ASSAY_LEDGER = HERE / "data" / "curation" / "assay_curation.csv"

# The routes whose half-life is an elimination half-life of the dosed molecule.
PARENTERAL = ("iv", "sc", "ip", "im")

# Route vocabulary as the sources and the database spell it, including PEPlife2's
# "i/v", "Intravenouly" and "Intarvenous". Ordered: the first match wins, and the
# intravenous patterns come first so "femoral vein" is not read as a bare injection.
ROUTE_PATTERNS = (
    ("iv", r"intra-?venous|intravenou|intarvenous|\bi\.?/?v\.?\b|tail ?vein|jugular|"
           r"retro-?orbital|femoral vein|mesenteric vein"),
    ("sc", r"sub-?cutan|\bs\.?/?c\.?\b"),
    ("ip", r"intra-?periton|intraperitonial|\bi\.?/?p\.?\b"),
    ("im", r"intra-?muscul|\bi\.?/?m\.?\b"),
    ("oral", r"\boral|gavage|per ?os"),
    ("intranasal", r"intra-?nasal|nasal"),
    ("rectal", r"rectal"),
)


def stated_route(text: object) -> str:
    """The route the database text names, or "" when it names none or more than one.

    Two matches is not a route: "(i/v and s/c route of injection)" describes a study,
    not an arm, and the arm is then settled by the ledger instead.
    """
    value = "" if text is None or pd.isna(text) else str(text).lower()
    hits = [label for label, pattern in ROUTE_PATTERNS if re.search(pattern, value)]
    return hits[0] if len(hits) == 1 else ""


def load_ledger(path: Path = ROUTE_LEDGER) -> pd.DataFrame:
    ledger = pd.read_csv(path, dtype=str).fillna("")
    ledger["routes"] = ledger["evidenced_routes"].map(
        lambda v: [r for r in v.split("|") if r])
    return ledger.set_index("source")


def resolve(frame: pd.DataFrame, ledger: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per-record curated route, evidence level and quote, for records needing curation.

    `frame` needs `PMID` and the database text columns `Test Sample` and `Dose Route`.
    Returns a frame aligned to `frame.index`.
    """
    ledger = load_ledger() if ledger is None else ledger
    stated = frame["Test Sample"].map(stated_route)
    blank = stated.eq("")
    stated = stated.mask(blank, frame.get("Dose Route", "").map(stated_route))

    route, level, quote, why = [], [], [], []
    for source, said in zip(frame["PMID"].astype(str), stated):
        if source not in ledger.index:
            route.append(""); level.append("none"); quote.append("")
            why.append("source not in the route curation ledger")
            continue
        entry = ledger.loc[source]
        evidenced = entry["routes"]
        if not evidenced:
            route.append(""); level.append("none"); quote.append("")
            why.append(entry["note"] or "no primary-source route evidence")
        elif len(evidenced) == 1:
            # One route in the source: every arm of it used that route, whatever the
            # database says - a single-route study cannot have a second arm.
            route.append(evidenced[0]); level.append(entry["evidence_level"])
            quote.append(entry["evidence_quote"]); why.append("")
        elif said in evidenced:
            route.append(said); level.append(entry["evidence_level"])
            quote.append(entry["evidence_quote"]); why.append("")
        else:
            # The source used several routes and nothing ties this record to one of them.
            route.append(""); level.append("none"); quote.append("")
            why.append(f"source evidences {'/'.join(evidenced)} but the record's "
                       f"annotation ({said or 'none'}) matches no evidenced arm")

    out = pd.DataFrame({"curated_route": route, "route_evidence_level": level,
                        "route_evidence_quote": quote, "route_unresolved_reason": why},
                       index=frame.index)
    out["route_curation_keep"] = out["curated_route"].isin(PARENTERAL)
    return out


def apply_route(frame: pd.DataFrame) -> pd.DataFrame:
    """Write the curated route and its provenance onto a copy of `frame`.

    The curated route is a fact about the RECORD, so it belongs in the published kept file
    beside the record rather than being derived again downstream - otherwise the kept file
    and the modelling file disagree about the same record's route, and a table built from
    one contradicts a figure built from the other.

    Whether a record is KEPT for its route is a separate, cohort-level decision, and this
    function does not take it: it only fills the column in and says where the value came
    from. `route_evidence_level` reads `database_route_column` for a record that already
    carried a route, which no source was re-read for.
    """
    out = frame.copy()
    for column in ("route_evidence_level", "route_evidence_quote",
                   "route_unresolved_reason"):
        out[column] = ""
    out["route_evidence_level"] = "database_route_column"

    needs = out["route_class"].fillna("not reported").eq("not reported")
    if not needs.any():
        return out
    curated = resolve(out[needs])
    for column in ("route_evidence_level", "route_evidence_quote",
                   "route_unresolved_reason"):
        out.loc[needs, column] = curated[column]
    resolved_rows = curated.index[curated["curated_route"].ne("")]
    out.loc[resolved_rows, "route_class"] = curated.loc[resolved_rows, "curated_route"]

    # A curated route without the sentence it rests on is the one state this module must
    # never ship: it would be indistinguishable from a transcription while claiming to be
    # evidence. Checked rather than trusted, because a ledger row is hand-written.
    claimed = out["route_evidence_level"].isin({"primary_arm", "primary_source"})
    unsupported = claimed & out["route_evidence_quote"].eq("")
    if unsupported.any():
        bad = sorted(out.loc[unsupported, "PMID"].astype(str).unique())
        raise SystemExit(
            f"{int(unsupported.sum())} records claim a primary-source route with no "
            f"evidence quote; add one to data/curation/route_curation.csv for: "
            + ", ".join(bad))
    return out


def summarise(resolved: pd.DataFrame) -> pd.DataFrame:
    return (resolved.groupby(["curated_route", "route_evidence_level"], dropna=False)
            .size().rename("records").reset_index()
            .sort_values("records", ascending=False, kind="stable"))


def resolve_assay(frame: pd.DataFrame, bucketed: pd.Series,
                  path: Path = ASSAY_LEDGER) -> pd.DataFrame:
    """Assay method with the source-curated value filled in where the text yielded none.

    `bucketed` is what `study_conditions` made of the `Assay` column. The ledger is only
    consulted where that came out "unknown": a method the column names is read from the
    column, so the vocabulary stays the single place that decides what "LC/MS" means, and
    the ledger holds only what the column could not say.
    """
    ledger = pd.read_csv(path, dtype=str).fillna("").set_index("source")
    method = bucketed.astype(object).copy()
    quote = pd.Series("", index=frame.index, dtype=object)
    source = pd.Series("", index=frame.index, dtype=object)

    for idx, (said, src) in zip(frame.index, zip(method, frame["PMID"].astype(str))):
        if said != "unknown" or src not in ledger.index:
            continue
        entry = ledger.loc[src]
        if entry["assay_method"]:
            method[idx] = entry["assay_method"]
            quote[idx] = entry["evidence_quote"]
            source[idx] = entry["evidence_source"]

    return pd.DataFrame({"assay_method": method, "assay_evidence_quote": quote,
                         "assay_evidence_source": source}, index=frame.index)


if __name__ == "__main__":
    meta = pd.read_csv(HERE / "rodent" / "strict" / "halflife_meta_with_split.csv",
                       low_memory=False)
    kept = pd.read_csv(HERE / "rodent" / "strict" / "in_vivo_half_life.csv",
                       low_memory=False)
    meta["peplife_id"] = meta["peplife_id"].astype(str)
    kept["ID"] = kept["ID"].astype(str)
    joined = meta.merge(kept.set_index("ID")[["Test Sample", "Dose Route"]],
                        left_on="peplife_id", right_index=True, how="left")
    target = joined[joined["route_class"].eq("not reported")]
    resolved = resolve(target)
    print("ROUTE")
    print(summarise(resolved).to_string(index=False))
    print(f"\nkeep {int(resolved['route_curation_keep'].sum())} "
          f"of {len(resolved)} curated records")

    import sys
    sys.path.insert(0, str(HERE))
    from study_conditions import study_conditions  # noqa: E402

    kept["route_class"] = kept["route_class"].fillna("not reported")
    bucketed = study_conditions(kept, kept["route_class"])["assay_method"]
    assay = resolve_assay(kept, bucketed)
    print("\n\nASSAY")
    print(assay["assay_method"].value_counts().to_string())
    print(f"\ncurated from the sources: "
          f"{int((assay['assay_evidence_quote'] != '').sum())} records")
