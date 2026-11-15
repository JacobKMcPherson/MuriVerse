"""Consolidate the per-source verification verdicts into one ledger and assign evidence tiers.

Reads MuriVerse/data/curation/label_verification/verdicts/*.json, joins the worklist, re-checks
every quote against the source text it claims to come from, and writes
MuriVerse/data/curation/halflife_label_verification.csv (one row per intravenous record).

Tiers (the first that applies):
  discrepant        the source contradicts or does not support the label: a corrected value, a
                    non-intravenous route, a different species, a non-elimination phase, a value
                    only cited from other work, a treated or non-plasma arm, or an arm that
                    cannot be tied to the record
  verified_body     value, species, phase and an intravenous route all agree in the article or
                    patent body
  verified_abstract the same, from the abstract only (no body text could be retrieved)
  value_only        the value is confirmed but the source does not state an intravenous route
  unverifiable      no source text, or the text does not contain the value
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
LV = REPO / "MuriVerse/data/curation/label_verification"
CORPUS = REPO / "MuriVerse/data/curation_cache/label_sources"
OUT = REPO / "MuriVerse/data/curation/halflife_label_verification.csv"
PAGE_FETCHED = {"10377427", "17204551", "19805567", "21244372", "25039358", "22934681"}  # article pages read directly
UNRESOLVED_CITATION: set[str] = set()
# The document actually read where PEPlife2's citation is malformed or wrong (author-supplied).
SOURCE_DOCUMENT = {
    "US 2022/0048569 W": "WO2023081149A1 (PEPlife2 id 4304)",
    "EP 2022074607 W": "WO2023031455A1",
    "30817059": "Xenotransplantation 2019;26:e12505 (PEPlife2 cites a JBC URL)",
    "8542672": "Clin Exp Pharmacol Physiol 1995;22:614-618 (author-supplied page images)",
}

# Reader verdicts that the evidence does not support at face value, with the reason.
OVERRIDES = {
    "hl_00299": ("discrepant", "printed '362 +/- 2.3 min' is implausible next to CNP-A 21.7 min; probable lost decimal (36.2 min)"),
    "hl_00184": ("discrepant", "value is cited from earlier work, not measured in this paper"),
    "hl_00026": ("discrepant", "value is a cited literature figure; free GLP-1 was not dosed in this study"),
}


def norm(text: str) -> str:
    text = text.replace("\u00ad", "").replace("-\n", "")
    text = re.sub(r"(\w)- (\w)", r"\1\2", text)
    return re.sub(r"\s+", " ", text).strip().lower()


def quote_ok(quote: object, haystack: str) -> bool:
    if not isinstance(quote, str) or not quote.strip():
        return True
    return norm(quote).replace("…", "") in haystack or all(
        part.strip() and norm(part) in haystack for part in re.split(r"\.\.\.|…|\|", quote))


def tier(row: pd.Series) -> str:
    note = str(row.get("notes") or "")
    flagged = bool(re.match(r"\s*(CITED_NOT_MEASURED|NON_BASELINE_ARM|NON_PLASMA)", note))
    if row.label_verdict in ("corrected", "phase_not_elimination", "ambiguous_arm"):
        return "discrepant"
    if row.route_verdict == "not_iv" or row.species_matches is False or flagged:
        return "discrepant"
    if row.label_verdict == "not_located":
        return "unverifiable"
    iv = row.route_verdict in ("iv_confirmed_arm", "iv_confirmed_study_level")
    if row.label_verdict == "confirmed" and iv:
        return "verified_abstract" if row.evidence_level == "abstract" else "verified_body"
    if row.label_verdict == "confirmed":
        return "value_only"
    return "unverifiable"


def main() -> int:
    work = pd.read_csv(LV / "worklist_all.csv", dtype={"PMID": str})
    rows, bad_quotes = [], []
    for path in sorted((LV / "verdicts").glob("*.json")):
        source = path.stem
        text_path = CORPUS / f"{source}.txt"
        abstract_path = LV / "abstracts" / f"{source}.txt"
        for rec in json.loads(path.read_text()):
            level = rec.get("evidence_level") or ("body" if text_path.exists() or source in PAGE_FETCHED else "none")
            haystack = ""
            if level == "abstract" and abstract_path.exists():
                haystack = norm(abstract_path.read_text())
            elif text_path.exists():
                haystack = norm(text_path.read_text(errors="replace"))
            checkable = bool(haystack)
            for key in ("source_quote", "route_quote"):
                if checkable and not quote_ok(rec.get(key), haystack):
                    bad_quotes.append((rec["id"], key))
            rows.append({**rec, "evidence_level": level, "source_file": source,
                         "quotes_checked": checkable})
    led = pd.DataFrame(rows)
    led = work[["id", "PMID", "species", "Name", "label_hours", "pk_phase",
                "route_evidence_level", "half_life_provenance"]].merge(
        led.drop(columns=["species_matches"]).assign(species_matches=led["species_matches"]),
        on="id", how="left")
    miss = led.label_verdict.isna()
    led.loc[miss, "label_verdict"] = "not_located"
    led.loc[miss, "route_verdict"] = "iv_not_stated"
    led.loc[miss, "evidence_level"] = "none"
    led.loc[miss, "notes"] = "no retrievable source text"
    led["source_document"] = led.PMID.map(SOURCE_DOCUMENT).fillna("")
    led["species_matches"] = led["species_matches"].map(lambda v: v if isinstance(v, bool) else None)
    led["evidence_tier"] = led.apply(tier, axis=1)
    led["tier_override_reason"] = ""
    for rid, (new_tier, why) in OVERRIDES.items():
        led.loc[led.id == rid, ["evidence_tier", "tier_override_reason"]] = [new_tier, why]
    led.to_csv(OUT, index=False)
    print(f"{len(led)} records; quotes that fail verbatim re-check: {len(bad_quotes)}")
    for item in bad_quotes[:20]:
        print("  ", item)
    print(led.evidence_tier.value_counts().to_string())
    print(pd.crosstab(led.evidence_tier, led.evidence_level).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
