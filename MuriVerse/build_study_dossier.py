"""Assemble a per-record provenance table and a per-study dossier for a species cohort.

Answers "where did this half-life value come from?" for every record in the cohort by joining,
on one row per PEPlife2 entry:

  * the PEPlife2 fields themselves (sequence, reported half-life, matrix, `Vivo/Vitro` flag,
    dose route, assay),
  * the study design recovered from the `Test Sample` free text, because PEPlife2's flag has
    no ex vivo category and is not internally consistent,
  * whether a live animal was actually dosed (`animal_dosed`), which is not the same question
    as the design bucket - see `dosed_of`,
  * the PubMed record behind the PMID (title, journal, year, DOI, publication types),
  * how far that publication has got through the full-text curation pipeline, and the curated
    per-paper evidence where an extraction exists.

This module owns the design classification; `figures/code/04_data_muriverse.R` reads the
`matrix`/`design` columns from here rather than re-deriving them.

  python build_study_dossier.py --species mouse
"""
import argparse
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
PUBMED = ROOT / "data/pubmed"
PAPERS = ROOT / "data/peplife2/papers"

CELL_LINE = {
    "mouse": r"\b(?:3t3|raw ?264|l929|b16|neuro-?2a|c2c12|hepa|p388|ct26|4t1|mc38)\b",
    "rat": r"\b(?:pc-?12|h9c2|gh3|rin-?m|ins-?1|c6)\b",
}
TISSUE = (r"liver|kidney|lung|brain|intestin|gastric|muscle|skin|homogenate|microsome|"
          r"hepatocyte")
# A concentration of the matrix itself ("25% mouse serum", "50% v/v mouse plasma", "Full mouse
# serum") can only describe peptide added to matrix in a tube, whatever the Vivo/Vitro flag says.
# The percentage must qualify the matrix: "serum (oral delivery in 0.5% Intravail)" is in vivo.
DILUTED_MATRIX = (r"(?:\d+(?:\.\d+)?\s*%|\bfull\b)\s*(?:\(?v\s*/\s*v\)?\s*)?(?:of\s+)?"
                  r"(?:mouse|mice|murine|rat|human|fetal|fbs)?\s*(?:serum|plasma|blood)")
# Evidence of an in-life administration, used only to overrule a wrong `in vitro` flag.
ADMINISTERED = (r"inject|administer|gavage|infus|\bdose[ds]?\b|route of delivery|"
                r"\bi\.\s?[vpm]\.|\bs\.\s?c\.")

DESIGNS = ["in vivo PK", "ex vivo (fluid stability)", "ex vivo (tissue homogenate)",
           "in vitro (other)", "ambiguous / no flag"]
DOSING = ["dosed", "not dosed", "unknown"]


def norm_seq(s: object) -> str:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return ""
    return re.sub(r"[^A-Za-z]", "", str(s)).upper()


def setting_of(flag: object) -> str:
    f = re.sub(r"\s+", " ", str(flag or "")).strip().lower()
    if "vivo/vitro" in f:
        return "other"
    if f.startswith("in vivo"):
        return "in vivo"
    if re.match(r"^in vi(t)?ro", f):
        return "in vitro"
    if f.startswith("ex vivo"):
        return "ex vivo"
    return "other"


def matrix_of(sample: str, cell_line: bool) -> str:
    s = sample.lower()
    if cell_line:
        return "cell line"
    if "plasma" in s:
        return "plasma"
    if "serum" in s:
        return "serum"
    if "blood" in s:
        return "blood"
    if re.search(TISSUE, s):
        return "tissue / homogenate"
    return "whole animal / other"


def design_of(setting: str, matrix: str, sample: str = "") -> str:
    if matrix == "tissue / homogenate":
        return "ex vivo (tissue homogenate)"
    if matrix in {"plasma", "serum", "blood"} and (
            setting == "in vitro" or re.search(DILUTED_MATRIX, sample, re.I)):
        return "ex vivo (fluid stability)"
    if setting == "in vivo":
        return "in vivo PK"
    if setting == "in vitro":
        return "in vitro (other)"
    return "ambiguous / no flag"


def dosed_of(setting: str, design: str, sample: str = "") -> str:
    """Did a live animal receive the peptide, or was matrix spiked in a tube?

    Not recoverable from `design` alone: `design_of` buckets every tissue matrix as ex vivo
    before it reads the flag, so an animal dosed in life whose organ was then harvested lands
    in `ex vivo (tissue homogenate)`. Conversely an `in vivo`-flagged diluted matrix ("25%
    mouse serum") is a spiked tube and must not count as dosed.
    """
    if design == "in vivo PK":
        return "dosed"
    if design == "ex vivo (tissue homogenate)" and setting == "in vivo":
        return "dosed"
    if design == "in vitro (other)" and re.search(ADMINISTERED, sample, re.I):
        return "dosed"
    if design == "ambiguous / no flag":
        return "unknown"
    return "not dosed"


def load_extractions(species_dir: Path) -> pd.DataFrame:
    """One row per curated PMID: the agent-level evidence collapsed to readable summaries."""
    def reported(agents: list[dict], key: str) -> list[dict]:
        return [a for a in agents
                if a.get(key) and "not reported" not in str(a[key]).lower()]

    rows = []
    for f in sorted((species_dir / "extractions").glob("*.json")):
        j = json.loads(f.read_text())
        agents = j.get("agents") or []
        rows.append({
            "PMID": str(j.get("pmid", f.stem)),
            "curated_study_types": "|".join(j.get("study_types") or []),
            "curated_species_dosed": "|".join(
                s.get("species", "") for s in (j.get("in_vivo_species") or []) if isinstance(s, dict)),
            "curated_n_agents": len(agents),
            "curated_agents_with_half_life": len(reported(agents, "half_life")),
            "curated_incomplete": not (j.get("study_types") or agents),
            "curated_half_lives": " ;; ".join(
                f"{a.get('name', '?')}: {a['half_life']}"
                f"{' [' + a['half_life_matrix'] + ']' if a.get('half_life_matrix') else ''}"
                for a in reported(agents, "half_life")),
            "curated_doses": " ;; ".join(
                f"{a.get('name', '?')}: {a['dose']} {a.get('route') or ''}".strip()
                for a in reported(agents, "dose")),
            "curated_bioanalytical_method": j.get("bioanalytical_method") or "",
            "curated_notes": j.get("notes") or "",
        })
    return pd.DataFrame(rows)


def build(species: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sdir = ROOT / species
    rec = pd.read_csv(sdir / "peplife2_rows.csv", low_memory=False)
    rec["PMID"] = rec["PMID"].astype(str).str.strip()
    rec["sample"] = rec["Test Sample"].fillna("").astype(str)
    rec["sequence"] = rec["Sequence"].map(norm_seq)
    rec["setting"] = rec["Vivo/Vitro"].map(setting_of)
    rec["matrix"] = [matrix_of(s, bool(c)) for s, c in zip(rec["sample"], rec["matrix_is_cell_line"])]
    rec["design"] = [design_of(s, m, t)
                     for s, m, t in zip(rec["setting"], rec["matrix"], rec["sample"])]
    rec["animal_dosed"] = [dosed_of(s, d, t)
                           for s, d, t in zip(rec["setting"], rec["design"], rec["sample"])]
    rec["in_cohort"] = ~rec["matrix_is_cell_line"].astype(bool)
    rec["half_life_h"] = pd.to_numeric(rec["half-life(in Sec)"], errors="coerce") / 3600

    pub = pd.read_csv(PUBMED / f"pubmed_{species}.csv", low_memory=False)
    pub["PMID"] = pub["PMID"].astype(str)
    pub = pub.rename(columns={"year": "pub_year", "title": "study_title",
                              "journal": "study_journal", "doi": "study_doi",
                              "first_author": "study_first_author",
                              "publication_types": "study_publication_types",
                              "mesh_terms": "study_mesh_terms"})
    keep = ["PMID", "study_title", "study_journal", "pub_year", "study_doi", "pmcid",
            "study_first_author", "study_publication_types", "study_mesh_terms",
            "abstract_chars"]

    status = pd.read_csv(sdir / "fulltext_retrieval_status.csv", dtype={"PMID": str})
    status = status[["PMID", "fulltext_retrieved", "fulltext_status", "fulltext_chars"]]

    nuance = pd.read_csv(sdir / "nuance_subset.csv", dtype={"pmid": str}).rename(
        columns={"pmid": "PMID"})
    ncols = [c for c in ["PMID", "extraction_status", "is_in_vivo", "is_ex_vivo", "is_in_vitro",
                         "is_clinical", "in_vivo_routes", "ex_vivo_matrices",
                         "half_life_summaries", "bioanalytical_method"] if c in nuance.columns]

    pdfs = PAPERS / "paper_download_status.csv"
    pdf_ok = pd.DataFrame(columns=["PMID", "peplife2_pdf"])
    if pdfs.exists():
        p = pd.read_csv(pdfs, dtype={"PMID": str})
        pdf_ok = p.assign(peplife2_pdf=p["status"].isin(["downloaded", "already_present"])
                          )[["PMID", "peplife2_pdf"]]

    curated = load_extractions(sdir)

    prov = (rec
            .merge(pub[keep], on="PMID", how="left")
            .merge(status, on="PMID", how="left")
            .merge(nuance[ncols], on="PMID", how="left")
            .merge(pdf_ok, on="PMID", how="left")
            .merge(curated, on="PMID", how="left")
            .sort_values("ID"))
    prov["evidence"] = [
        "curated extraction" if e == "extracted" else
        "full text retrieved" if bool(r) else
        "publisher PDF only" if bool(d) else "metadata only"
        for e, r, d in zip(prov["extraction_status"].fillna(""),
                           prov["fulltext_retrieved"].fillna(False),
                           prov["peplife2_pdf"].fillna(False))]
    prov["pub_year"] = pd.to_numeric(prov["pub_year"], errors="coerce").astype("Int64")
    # PEPlife2 repeats identical measurements as separate entries; an unweighted mean or
    # median over the raw rows silently counts those studies many times over.
    dup_key = ["PMID", "sequence", "sample", "Half Life", "Units Half Life"]
    prov["identical_copies"] = prov.groupby(dup_key, dropna=False)["ID"].transform("size")
    prov["is_repeat"] = prov.duplicated(dup_key, keep="first")

    cohort = prov[prov["in_cohort"]]
    studies = (cohort.groupby("PMID", dropna=False).agg(
        records=("ID", "size"),
        distinct_measurements=("is_repeat", lambda s: int((~s).sum())),
        sequences=("sequence", lambda s: s[s != ""].nunique()),
        records_with_half_life=("half_life_h", "count"),
        min_half_life_h=("half_life_h", "min"),
        median_half_life_h=("half_life_h", "median"),
        max_half_life_h=("half_life_h", "max"),
        designs=("design", lambda s: "|".join(sorted(set(s)))),
        records_dosed=("animal_dosed", lambda s: int((s == "dosed").sum())),
        matrices=("matrix", lambda s: "|".join(sorted(set(s)))),
        settings_flagged=("setting", lambda s: "|".join(sorted(set(s)))),
        dose_routes=("Dose Route", lambda s: "|".join(sorted({str(x) for x in s.dropna()}))),
        assays=("Assay", lambda s: "|".join(sorted({str(x) for x in s.dropna()}))),
    ).reset_index())
    for d in DESIGNS:
        studies[f"n_{d}"] = (cohort.assign(hit=cohort["design"] == d)
                             .groupby("PMID")["hit"].sum().reindex(studies["PMID"]).values)
    meta_cols = [c for c in prov.columns if c.startswith(("study_", "curated_")) or c in
                 {"pub_year", "pmcid", "fulltext_retrieved", "fulltext_status", "fulltext_chars",
                  "extraction_status", "is_in_vivo", "is_ex_vivo", "is_in_vitro", "is_clinical",
                  "in_vivo_routes", "ex_vivo_matrices", "peplife2_pdf", "abstract_chars",
                  "evidence"}]
    studies = studies.merge(prov.drop_duplicates("PMID")[["PMID"] + meta_cols], on="PMID", how="left")

    seqs = (cohort[cohort["sequence"] != ""].groupby("sequence").agg(
        records=("ID", "size"),
        studies=("PMID", "nunique"),
        length=("Length", "max"),
        names=("Name", lambda s: "|".join(sorted({str(x) for x in s.dropna()})[:5])),
        designs=("design", lambda s: "|".join(sorted(set(s)))),
        values_with_half_life=("half_life_h", "count"),
        min_half_life_h=("half_life_h", "min"),
        median_half_life_h=("half_life_h", "median"),
        max_half_life_h=("half_life_h", "max"),
        pmids=("PMID", lambda s: "|".join(sorted(set(s), key=lambda p: (not p.isdigit(), p)))),
    ).reset_index())
    # A sequence whose reported half-lives span >1 order of magnitude is worth inspecting before
    # it is collapsed to a single training label.
    seqs["fold_range"] = seqs["max_half_life_h"] / seqs["min_half_life_h"].replace(0, pd.NA)
    seqs["multi_study_conflict"] = (seqs["studies"] > 1) & (seqs["fold_range"] > 10)
    return prov, studies.sort_values("records", ascending=False), seqs.sort_values(
        ["studies", "records"], ascending=False)


def write_dossier(species: str, studies: pd.DataFrame, prov: pd.DataFrame, dest: Path) -> None:
    cohort = prov[prov["in_cohort"]]
    lines = [f"# {species.capitalize()} cohort - study dossier", "",
             f"One section per publication behind the {len(cohort)} {species} records in "
             f"`{species}/peplife2_rows.csv`, ordered by how many records it contributes. "
             f"Generated by `build_study_dossier.py`; do not edit by hand.", ""]
    for _, s in studies.iterrows():
        pmid = s["PMID"]
        head = f"## PMID {pmid}" if str(pmid).isdigit() else f"## {pmid} (non-numeric identifier)"
        lines += [head, ""]
        cite = " ".join(x for x in [
            str(s.get("study_first_author") or ""),
            f"({s['pub_year']})" if pd.notna(s.get("pub_year")) else "",
            str(s.get("study_title") or ""),
            f"*{s['study_journal']}*." if pd.notna(s.get("study_journal")) else "",
        ] if x)
        if cite.strip():
            lines.append(cite)
        if isinstance(s.get("study_doi"), str) and s["study_doi"]:
            lines.append(f"DOI: https://doi.org/{s['study_doi']}")
        lines += ["",
                  f"- Contributes **{s['records']} records**"
                  + (f" ({s['distinct_measurements']} distinct, the rest are exact repeats)"
                     if s["distinct_measurements"] < s["records"] else "")
                  + f", {s['sequences']} unique sequences; "
                  f"{s['records_with_half_life']} carry a numeric half-life.",
                  f"- Recovered design: {s['designs']}; PEPlife2 `Vivo/Vitro`: {s['settings_flagged']}.",
                  f"- Live animal dosed in {s['records_dosed']} of {s['records']} records.",
                  f"- Assay matrices: {s['matrices']}.",
                  f"- Evidence available: **{s['evidence']}**"
                  + (f" (PMC {s['pmcid']})" if isinstance(s.get("pmcid"), str) and s["pmcid"] else "")
                  + "."]
        if pd.notna(s.get("min_half_life_h")):
            lines.append(f"- Reported half-life: {s['min_half_life_h']:.3g}-{s['max_half_life_h']:.3g} h "
                         f"(median {s['median_half_life_h']:.3g} h).")
        if isinstance(s.get("dose_routes"), str) and s["dose_routes"]:
            lines.append(f"- PEPlife2 dose route: {s['dose_routes']}.")
        if s.get("extraction_status") == "extracted":
            if s.get("curated_incomplete") is True:
                lines.append("- Curated extraction: **file present but empty** - re-extract.")
            else:
                lines.append(f"- Curated extraction: study types {s.get('curated_study_types') or 'n/a'}; "
                             f"{int(s.get('curated_n_agents') or 0)} agents, "
                             f"{int(s.get('curated_agents_with_half_life') or 0)} with a body-quoted "
                             f"half-life.")
            if isinstance(s.get("in_vivo_routes"), str) and s["in_vivo_routes"]:
                lines.append(f"- Curated in vivo routes: {s['in_vivo_routes']}.")
            if isinstance(s.get("ex_vivo_matrices"), str) and s["ex_vivo_matrices"]:
                lines.append(f"- Curated ex vivo matrices: {s['ex_vivo_matrices']}.")
            if isinstance(s.get("curated_bioanalytical_method"), str) and s["curated_bioanalytical_method"]:
                lines.append(f"- Bioanalytical method: {s['curated_bioanalytical_method']}.")
            if isinstance(s.get("curated_half_lives"), str) and s["curated_half_lives"]:
                lines.append(f"- Body-quoted half-lives: {s['curated_half_lives']}.")
            if isinstance(s.get("curated_notes"), str) and s["curated_notes"]:
                lines.append(f"- Caveats: {s['curated_notes']}")
        rows = cohort[cohort["PMID"] == pmid]
        lines += ["", "| ID | peptide | half-life | matrix (`Test Sample`) | flag | design |",
                  "|---|---|---|---|---|---|"]
        for _, r in rows.iterrows():
            hl = f"{r['Half Life']} {r['Units Half Life']}" if pd.notna(r["Half Life"]) else "not reported"
            lines.append(f"| {r['ID']} | {str(r['Name'] or '')[:40]} | {hl} | "
                         f"{str(r['sample'])[:60]} | {r['Vivo/Vitro']} | {r['design']} |")
        lines.append("")
    dest.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="mouse")
    args = ap.parse_args()

    sdir = ROOT / args.species
    prov, studies, seqs = build(args.species)
    prov.to_csv(sdir / "records_provenance.csv", index=False)
    studies.to_csv(sdir / "studies.csv", index=False)
    seqs.to_csv(sdir / "sequences.csv", index=False)
    write_dossier(args.species, studies, prov, sdir / "STUDIES.md")

    cohort = prov[prov["in_cohort"]]
    print(f"[records]   {len(cohort)} in cohort ({len(prov) - len(cohort)} cell-line rows excluded)")
    print(cohort["design"].value_counts().to_string())
    print(f"[dosed]     {int((cohort['animal_dosed'] == 'dosed').sum())} records dosed a live animal")
    print(cohort["animal_dosed"].value_counts().to_string())
    print(f"[repeats]   {int(cohort['is_repeat'].sum())} records are exact repeats of another "
          f"record in the same study")
    print(f"[studies]   {len(studies)} publications")
    print(studies["evidence"].value_counts().to_string())
    print(f"[sequences] {len(seqs)} unique sequences, "
          f"{int(seqs['multi_study_conflict'].sum())} with >10x spread across >1 study")
    print(f"[done] records_provenance.csv, studies.csv, sequences.csv, STUDIES.md in {sdir}")


if __name__ == "__main__":
    main()
