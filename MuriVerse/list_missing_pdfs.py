"""List dosed studies whose source PDF is missing, as a manual-retrieval worklist.

`fetch_peplife2_papers.py` only mirrors what PEPlife2 hosts, and the dosed subset skews to
pharma and patent literature that is neither there nor in the PMC OA subset. Those studies
have no figure source at all, so `extract_pk_figures.py` can never reach them.

Emits `<species>/pdf_retrieval_worklist.csv` and a readable `<species>/PDF_WORKLIST.md`
ordered by how much of the cohort each missing paper blocks.

  python list_missing_pdfs.py --species mouse
"""
import argparse
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
PAPERS = ROOT / "data" / "peplife2" / "papers"


def build(species: str) -> pd.DataFrame:
    sdir = ROOT / species
    rec = pd.read_csv(sdir / "records_provenance.csv", low_memory=False)
    dosed = rec[rec["in_cohort"].astype(bool) & rec["animal_dosed"].eq("dosed")].copy()
    dosed["PMID"] = dosed["PMID"].astype(str)

    agg = (dosed.groupby("PMID")
           .agg(dosed_records=("ID", "size"),
                agents=("Name", lambda s: " ;; ".join(sorted({str(v) for v in s if pd.notna(v)}))),
                study_title=("study_title", "first"),
                study_journal=("study_journal", "first"),
                pub_year=("pub_year", "first"),
                study_doi=("study_doi", "first"),
                pmcid=("pmcid", "first"),
                first_author=("study_first_author", "first"),
                fulltext_status=("fulltext_status", "first"),
                extraction_status=("extraction_status", "first"),
                evidence=("evidence", "first"),
                peplife2_year=("YEAR", "first"),
                patent_no=("Patent No.", "first"),
                reference=("Reference", "first"))
           .reset_index())

    agg["has_pdf"] = agg["PMID"].map(lambda p: (PAPERS / f"{p}.pdf").exists())
    miss = agg[~agg["has_pdf"]].drop(columns=["has_pdf"]).copy()

    miss["is_patent"] = ~miss["PMID"].str.fullmatch(r"\d+")
    miss["pubmed_url"] = [f"https://pubmed.ncbi.nlm.nih.gov/{p}/" if not pat else ""
                          for p, pat in zip(miss["PMID"], miss["is_patent"])]
    miss["doi_url"] = ["https://doi.org/" + str(d) if pd.notna(d) and str(d).strip() else ""
                       for d in miss["study_doi"]]
    miss["pmc_url"] = ["https://pmc.ncbi.nlm.nih.gov/articles/" + str(c) + "/"
                       if pd.notna(c) and str(c).strip() else "" for c in miss["pmcid"]]
    miss["route"] = ["already have full text (figures only)"
                     if s == "full_text_retrieved" else
                     "PMC record exists - try the OA package or the publisher site"
                     if str(c).strip() and str(c) != "nan" else
                     "patent - retrieve from the patent office" if pat else
                     "no PMC record - publisher or interlibrary loan"
                     for s, c, pat in zip(miss["fulltext_status"], miss["pmcid"].fillna(""),
                                          miss["is_patent"])]
    return miss.sort_values(["dosed_records", "PMID"], ascending=[False, True])


def write_markdown(miss: pd.DataFrame, dest: Path, species: str, total: int) -> None:
    lines = [f"# {species.capitalize()} dosed studies with no source PDF",
             "",
             f"{len(miss)} of {total} dosed publications are not mirrored by PEPlife2 and have "
             f"no local PDF, so `extract_pk_figures.py` cannot reach them. Between them they "
             f"back **{int(miss['dosed_records'].sum())} dosed records**.",
             "",
             "Regenerate with `python list_missing_pdfs.py --species "
             f"{species}`; drop retrieved files into `data/peplife2/papers/<PMID>.pdf` and "
             "re-run `extract_pk_figures.py`.",
             ""]

    for label, sub in [("Have full text already (need the figures only)",
                        miss[miss["fulltext_status"].eq("full_text_retrieved")]),
                       ("Patents (no PMID)", miss[miss["is_patent"]]),
                       ("No full text and no PDF",
                        miss[~miss["is_patent"] & miss["fulltext_status"].ne("full_text_retrieved")])]:
        if sub.empty:
            continue
        lines += [f"## {label} ({len(sub)})", "",
                  "| Records | PMID | Year | Citation | DOI | PMCID |",
                  "|---:|---|---:|---|---|---|"]
        for _, r in sub.iterrows():
            title = re.sub(r"\s+", " ", str(r["study_title"] or "")).strip()
            if title and title != "nan":
                cite = f"{r['first_author'] or '?'}, *{r['study_journal'] or '?'}* — {title}"
            else:
                # Patents carry no PubMed record; PEPlife2 keeps the number and a Lens link.
                patent = str(r["patent_no"] or "").strip()
                cite = f"Patent {patent}" if patent and patent != "nan" else "—"
                cite += f" (agents: {r['agents']})" if r["agents"] else ""
            doi = f"[{r['study_doi']}]({r['doi_url']})" if r["doi_url"] else "—"
            ref = str(r["reference"] or "").strip()
            if doi == "—" and ref.startswith("http"):
                doi = f"[Lens]({ref})"
            pmc = f"[{r['pmcid']}]({r['pmc_url']})" if r["pmc_url"] else "—"
            pmid = f"[{r['PMID']}]({r['pubmed_url']})" if r["pubmed_url"] else f"`{r['PMID']}`"
            year = r["pub_year"] if pd.notna(r["pub_year"]) else r["peplife2_year"]
            year = "" if pd.isna(year) else int(year)
            lines.append(f"| {r['dosed_records']} | {pmid} | {year} | {cite} | {doi} | {pmc} |")
        lines.append("")
    dest.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", default="mouse", choices=["mouse", "rat"])
    args = ap.parse_args()

    sdir = ROOT / args.species
    rec = pd.read_csv(sdir / "records_provenance.csv", low_memory=False)
    total = rec[rec["in_cohort"].astype(bool) & rec["animal_dosed"].eq("dosed")]["PMID"].nunique()

    miss = build(args.species)
    csv = sdir / "pdf_retrieval_worklist.csv"
    miss.to_csv(csv, index=False)
    write_markdown(miss, sdir / "PDF_WORKLIST.md", args.species, total)

    print(f"[write] {csv}")
    print(f"[write] {sdir / 'PDF_WORKLIST.md'}")
    print(f"[missing] {len(miss)} of {total} dosed publications, "
          f"backing {int(miss['dosed_records'].sum())} dosed records")
    print(miss["route"].value_counts().to_string())


if __name__ == "__main__":
    main()
