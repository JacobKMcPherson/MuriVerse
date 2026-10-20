"""Fetch PubMed metadata for every publication behind a cohort.

The PEPlife2 record table carries only a PMID and a free-text `Reference` string, which is not
enough to say what kind of study a half-life value came from. This pulls the PubMed record for
each PMID (title, journal, year, DOI, PMCID, publication types, MeSH terms, abstract) so the
cohort can be described at the level of the underlying literature.

Raw efetch XML is archived per batch under `data/pubmed/` so the tidy table can be rebuilt
offline.

  python fetch_pubmed_metadata.py --species mouse
  python fetch_pubmed_metadata.py --species mouse --offline
"""
import argparse
import re
import time
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "data/pubmed"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = {"tool": "MuriVerse", "email": "muriverse@example.org"}
BATCH = 100

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "MuriVerse/1.0 (research; muriverse@example.org)"


def text_of(node: ET.Element | None) -> str:
    """Flatten an element's text, including inline markup such as <i>/<sup>."""
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip() if node is not None else ""


def fetch_batch(pmids: list[str], api_key: str = "") -> str:
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml", **TOOL}
    if api_key:
        params["api_key"] = api_key
    r = SESSION.post(EFETCH, data=params, timeout=120)
    r.raise_for_status()
    return r.text


def parse(xml: str) -> list[dict]:
    root = ET.fromstring(xml)
    out = []
    for art in root.findall(".//PubmedArticle"):
        cit = art.find("MedlineCitation")
        a = cit.find("Article")
        # Only the article's own ID list; .//ArticleId would also pick up cited references.
        ids = {i.get("IdType"): (i.text or "").strip()
               for i in art.findall("PubmedData/ArticleIdList/ArticleId")}
        journal = a.find("Journal")
        # PubDate is Year, or MedlineDate free text ("2011 Jan-Feb") when no clean year exists.
        pub = journal.find(".//PubDate") if journal is not None else None
        year = text_of(pub.find("Year")) if pub is not None else ""
        if not year and pub is not None:
            year = (re.search(r"\d{4}", text_of(pub.find("MedlineDate"))) or [""])[0] \
                if re.search(r"\d{4}", text_of(pub.find("MedlineDate"))) else ""
        abstract = " ".join(text_of(t) for t in a.findall(".//Abstract/AbstractText"))
        out.append({
            "PMID": text_of(cit.find("PMID")),
            "title": text_of(a.find("ArticleTitle")),
            "journal": text_of(journal.find("ISOAbbreviation")) if journal is not None else "",
            "journal_full": text_of(journal.find("Title")) if journal is not None else "",
            "year": year,
            "doi": ids.get("doi", ""),
            "pmcid": ids.get("pmc", ""),
            "publication_types": "|".join(
                text_of(p) for p in a.findall(".//PublicationTypeList/PublicationType")),
            "mesh_terms": "|".join(
                text_of(d.find("DescriptorName")) for d in cit.findall(".//MeshHeadingList/MeshHeading")),
            "n_authors": len(a.findall(".//AuthorList/Author")),
            "first_author": next(
                (f"{text_of(au.find('LastName'))} {text_of(au.find('Initials'))}".strip()
                 for au in a.findall(".//AuthorList/Author") if au.find("LastName") is not None), ""),
            "abstract_chars": len(abstract),
            "abstract": abstract,
        })
    return out


def cohort_pmids(species: str) -> list[str]:
    df = pd.read_csv(ROOT / f"{species}/peplife2_rows.csv", low_memory=False)
    return sorted({p for p in df["PMID"].astype(str).str.strip() if p.isdigit()}, key=int)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="mouse")
    ap.add_argument("--offline", action="store_true", help="parse archived XML only")
    ap.add_argument("--api-key", default="")
    args = ap.parse_args()

    ARCHIVE.mkdir(parents=True, exist_ok=True)
    pmids = cohort_pmids(args.species)
    print(f"[pubmed] {len(pmids)} numeric PMIDs in the {args.species} cohort")

    records = []
    for i in range(0, len(pmids), BATCH):
        chunk = pmids[i : i + BATCH]
        cache = ARCHIVE / f"efetch_{args.species}_{i // BATCH:02d}.xml"
        if args.offline or cache.exists():
            if not cache.exists():
                raise SystemExit(f"--offline but {cache} is missing")
            xml = cache.read_text()
        else:
            xml = fetch_batch(chunk, args.api_key)
            cache.write_text(xml)
            time.sleep(0.5)
        got = parse(xml)
        print(f"  [batch {i // BATCH}] {len(got)}/{len(chunk)} records")
        records.extend(got)

    df = pd.DataFrame(records).drop_duplicates("PMID")
    missing = sorted(set(pmids) - set(df["PMID"]), key=int)
    if missing:
        print(f"  [warn] {len(missing)} PMIDs returned no PubMed record: {missing[:10]}")
    dest = ARCHIVE / f"pubmed_{args.species}.csv"
    df.sort_values("PMID", key=lambda s: s.astype(int)).to_csv(dest, index=False)
    print(f"[done] {len(df)} PubMed records -> {dest}")
    print(df["year"].value_counts().sort_index().tail(10).to_string())


if __name__ == "__main__":
    main()
