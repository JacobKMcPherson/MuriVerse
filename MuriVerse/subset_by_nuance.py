"""Step 4 - subset <species>/extractions/*.json by study-design nuance.

Collapses each PMID's extraction JSON into one row that flags the paper-level nuances needed
to tell in vivo animal pharmacokinetics apart from ex vivo plasma/serum stability (and from
clinical or pure in vitro work), plus a route/matrix breakdown. PMIDs with no completed
extraction are still emitted (see `extraction_status`) so the subset stays a complete PMID
inventory, not just the classified ones.

  python subset_by_nuance.py --species mouse
"""
import argparse
import glob
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent

ROUTE_PATTERNS = [
    ("iv", r"intravenous|\bi\.?v\.?\b|tail vein"),
    ("sc", r"subcutaneous|\bs\.?c\.?\b"),
    ("ip", r"intraperitoneal|\bi\.?p\.?\b"),
    ("im", r"intramuscular|\bi\.?m\.?\b"),
    ("oral", r"\boral\b|gavage|\bp\.?o\.?\b"),
    ("icv", r"intracerebroventricular|\bicv\b"),
    ("infusion_pump", r"infusion|osmotic pump"),
    ("not_administered", r"not applicable|not administered|in vitro only"),
    ("not_reported", r"not reported"),
]

MATRIX_PATTERNS = [
    ("plasma", r"plasma"),
    ("serum", r"serum"),
    ("blood", r"\bblood\b|circulation|dynamic pet"),
    ("tissue_organ", r"tissue|organ|liver|kidney|spleen|brain|lavage|injection site"),
    ("urine", r"urine"),
    ("csf", r"csf|cerebrospinal"),
    ("not_applicable", r"not applicable"),
    ("not_reported", r"not reported"),
]


def categorize(text: str, patterns: list[tuple[str, str]]) -> set[str]:
    text = (text or "").lower()
    hits = {label for label, pat in patterns if re.search(pat, text)}
    return hits or ({"other"} if text.strip() else set())


def load_status(species_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(species_dir / "fulltext_retrieval_status.csv", dtype={"PMID": str})
    return df.set_index("PMID")


def summarize_extraction(pmid: str, d: dict, n_rows: int) -> dict:
    study_types = set(d.get("study_types") or [])
    agents = d.get("agents") or []

    in_vivo_routes, ex_vivo_matrices = set(), set()
    hl_summaries = []
    for a in agents:
        if a.get("half_life") and a["half_life"] != "not reported":
            hl_summaries.append(f"{a.get('name', '?')}: {a['half_life']} ({a.get('half_life_matrix', '?')})")
        in_vivo_routes |= categorize(a.get("route"), ROUTE_PATTERNS)
        if "ex vivo" in study_types:
            ex_vivo_matrices |= categorize(a.get("half_life_matrix"), MATRIX_PATTERNS)
    ex_vivo_matrices |= categorize(d.get("ex_vivo_detail"), MATRIX_PATTERNS) if "ex vivo" in study_types else set()
    # A route category is only a real in vivo route if the paper actually reports in vivo dosing.
    if "in vivo" not in study_types:
        in_vivo_routes = set()

    return {
        "pmid": pmid,
        "extraction_status": "extracted",
        "n_peplife_rows": n_rows,
        "study_types": "|".join(sorted(study_types)),
        "is_clinical": "clinical" in study_types,
        "is_in_vivo": "in vivo" in study_types,
        "is_ex_vivo": "ex vivo" in study_types,
        "is_in_vitro": "in vitro" in study_types,
        "in_vivo_species": "|".join(s.get("species", "") for s in (d.get("in_vivo_species") or [])),
        "in_vivo_routes": "|".join(sorted(in_vivo_routes)),
        "ex_vivo_origin": d.get("ex_vivo_origin") or "",
        "ex_vivo_matrices": "|".join(sorted(ex_vivo_matrices)),
        "n_agents": len(agents),
        "n_agents_with_half_life": len(hl_summaries),
        "half_life_summaries": " ; ".join(hl_summaries),
        "bioanalytical_method": d.get("bioanalytical_method") or "",
        "concentration_time_available": bool((d.get("concentration_time_figure") or {}).get("available")),
        "concentration_time_locator": (d.get("concentration_time_figure") or {}).get("locator") or "",
        "has_notes": bool((d.get("notes") or "").strip()),
        "notes": d.get("notes") or "",
    }


def unclassified_row(pmid: str, status: str, n_rows: int) -> dict:
    return {
        "pmid": pmid, "extraction_status": status, "n_peplife_rows": n_rows,
        "study_types": "", "is_clinical": False, "is_in_vivo": False, "is_ex_vivo": False,
        "is_in_vitro": False, "in_vivo_species": "", "in_vivo_routes": "", "ex_vivo_origin": "",
        "ex_vivo_matrices": "", "n_agents": 0, "n_agents_with_half_life": 0,
        "half_life_summaries": "", "bioanalytical_method": "",
        "concentration_time_available": False, "concentration_time_locator": "",
        "has_notes": False, "notes": "",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", choices=["mouse", "rat"], default="mouse")
    args = ap.parse_args()
    species_dir = ROOT / args.species

    status = load_status(species_dir)
    extracted_pmids = {Path(f).stem for f in glob.glob(str(species_dir / "extractions" / "*.json"))}
    rendered_pmids = {Path(f).stem for f in glob.glob(str(species_dir / "fulltext_txt" / "*.txt"))}

    rows = []
    for pmid, rec in status.iterrows():
        n_rows = int(rec.get("n_peplife_rows", 0))
        if pmid in extracted_pmids:
            d = json.loads((species_dir / "extractions" / f"{pmid}.json").read_text())
            rows.append(summarize_extraction(pmid, d, n_rows))
        elif pmid in rendered_pmids:
            rows.append(unclassified_row(pmid, "rendered_not_yet_classified", n_rows))
        elif rec.get("fulltext_retrieved"):
            rows.append(unclassified_row(pmid, "retrieved_not_yet_rendered", n_rows))
        else:
            rows.append(unclassified_row(pmid, "fulltext_unavailable", n_rows))

    out = pd.DataFrame(rows)
    out_path = species_dir / "nuance_subset.csv"
    out.to_csv(out_path, index=False)
    print(f"[done] wrote {len(out)} PMID rows to {out_path}")
    print(out["extraction_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
