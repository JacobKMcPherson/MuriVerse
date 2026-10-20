"""Build the curation worklist for the ten largest current strict-set sources."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from halflife_data_quality import prepare_training_metadata, priority_source_summary

ROOT = Path(__file__).resolve().parent
STRICT = ROOT / "rodent" / "strict"


def source_files(source: str) -> dict:
    if not source.isdigit():
        return {
            "kind": "non_numeric_source",
            "fulltext": [],
            "extractions": [],
            "note": "Patent or other non-PMID source; requires source-specific primary-text curation.",
        }
    found = {"kind": "PMID", "fulltext": [], "extractions": []}
    for species in ("mouse", "rat"):
        for folder in ("fulltext_txt", "pdf_txt"):
            path = ROOT / species / folder / f"{source}.txt"
            if path.exists():
                found["fulltext"].append(str(path.relative_to(ROOT.parent)))
        extraction = ROOT / species / "extractions" / f"{source}.json"
        if extraction.exists():
            found["extractions"].append(str(extraction.relative_to(ROOT.parent)))
    found["note"] = (
        "Exact-ID pk_observations in the extraction are eligible to promote rows to A/B."
        if found["extractions"]
        else "No structured extraction is currently present; retain rows as C until curated."
    )
    return found


def main() -> None:
    meta = pd.read_csv(STRICT / "halflife_meta_with_split.csv", low_memory=False)
    meta = prepare_training_metadata(meta)
    summary = priority_source_summary(meta, 10)
    for item in summary:
        item["source_files"] = source_files(item["source"])

    output = {
        "selection_rule": "ten largest sources in the current strict model table",
        "records_in_top10": int(
            meta.loc[meta["priority_source_top10"]].shape[0]
        ),
        "sources": summary,
        "curation_rule": (
            "A/B require exact PEPlife2-ID linkage to a primary-source arm and an agreeing "
            "numeric half-life. Metadata-only or narrative matching remains C."
        ),
    }
    path = STRICT / "priority_source_curation_worklist.json"
    path.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    print(f"[Saved] {path}")
    print(
        f"[Top 10] {output['records_in_top10']} rows; "
        f"{sum(x['gold_labels'] for x in summary)} currently A/B"
    )


if __name__ == "__main__":
    main()
