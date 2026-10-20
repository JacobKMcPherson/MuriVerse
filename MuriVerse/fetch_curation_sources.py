"""Mirror the source documents the curation ledgers were read from.

`data/curation/route_curation.csv` and `data/curation/assay_curation.csv` each quote the
sentence their verdict rests on, so the evidence travels with the repository. This script
re-fetches the documents those quotes came from, into `data/curation_cache/`, so a reader can
check a quote in context or extend a ledger when a paywalled paper becomes reachable.

Four kinds of source, in the order the curation preferred them:

  pubmed/    PubMed efetch XML, one file per PMID: title, abstract and NLM MeSH indexing.
             MeSH route headings ("Injections, Intravenous") are curated by NLM from the
             full text, which is why they count as primary-source evidence here even when
             the article body is paywalled.
  fulltext/  Europe PMC full-text XML, for the open-access subset. The repository's own
             mouse/ and rat/ fulltext_txt directories already hold some of these; this
             fetches the rest.
  probe/     For a non-open-access article whose full text Europe PMC has nonetheless
             INDEXED, a section-restricted term search ("does 'subcutaneous' occur in the
             METHODS of this PMID?"). It returns presence, not the sentence, which is why a
             probe alone never settles a route - but an absence is decisive, and a METHODS
             section carrying "tail vein" and no other dosing term resolves one.
  patent/    Google Patents HTML for the 8 patent sources, which have no PMID.

Nothing here is consulted at build time; the ledgers are. This is for re-reading and for
extending them.

Usage: python MuriVerse/fetch_curation_sources.py [--probe]

`--probe` adds the Europe PMC section probes, which are one request per term per source and
so are off by default.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
LEDGERS = HERE / "data" / "curation"
CACHE = HERE / "data" / "curation_cache"

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
PATENTS = "https://patents.google.com/patent/{}/en"

# The patent sources, with the Google Patents document id each is published under. The
# ledger keys them the way PEPlife2 does, which is not always the publication number.
PATENT_IDS = {
    "EP2277890A1": "EP2277890A1",
    "US6656906": "US6656906B1",
    "US7112567B2": "US7112567B2",
    "WO2009027844": "WO2009027844A2",
    "EP1605897B1": "EP1605897B1",
    "EP2423233A2": "EP2423233A2",
    "WO2012138867A2": "WO2012138867A2",
    "WO2013106273A2": "WO2013106273A2",
}

ROUTE_TERMS = ["intravenous", "intravenously", "tail vein", "jugular", "subcutaneous",
               "subcutaneously", "intraperitoneal", "intraperitoneally", "intramuscular",
               "intramuscularly", "oral", "gavage", "intranasal", "rectal", "infusion",
               "bolus", "femoral"]


def sources() -> list[str]:
    """Every source identifier either ledger names."""
    found: list[str] = []
    for name in ("route_curation.csv", "assay_curation.csv"):
        found += pd.read_csv(LEDGERS / name, dtype=str)["source"].tolist()
    return sorted(dict.fromkeys(found))


def get(url: str, timeout: int = 60) -> str | None:
    for _ in range(3):
        try:
            return urllib.request.urlopen(url, timeout=timeout).read().decode(
                "utf-8", "replace")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            time.sleep(2)
        except Exception:
            time.sleep(2)
    return None


def fetch_pubmed(pmids: list[str]) -> None:
    out = CACHE / "pubmed"
    out.mkdir(parents=True, exist_ok=True)
    todo = [p for p in pmids if not (out / f"{p}.xml").exists()]
    print(f"[pubmed]  {len(pmids)} sources, {len(todo)} to fetch")
    for i in range(0, len(todo), 20):
        batch = todo[i:i + 20]
        xml = get(f"{EFETCH}?db=pubmed&retmode=xml&id=" + ",".join(batch))
        if xml is None:
            continue
        # Split the set back into one file per record, so the cache is addressable by PMID.
        for chunk in xml.split("<PubmedArticle>")[1:]:
            body = "<PubmedArticle>" + chunk.split("</PubmedArticleSet>")[0]
            pmid = body.split('<PMID Version="1">')[1].split("</PMID>")[0]
            (out / f"{pmid}.xml").write_text(body)
        time.sleep(0.5)


def pmc_id(pmid: str) -> str | None:
    """The article's own PMC id, from the record's id list.

    Read from `PubmedData/ArticleIdList` rather than by searching the whole record: a
    PubMed record also carries ArticleId nodes for its comments-and-corrections links, and
    taking the first match anywhere returns another article's PMC id.
    """
    path = CACHE / "pubmed" / f"{pmid}.xml"
    if not path.exists():
        return None
    ids = ET.fromstring(path.read_text()).find("./PubmedData/ArticleIdList")
    if ids is None:
        return None
    return next((i.text for i in ids if i.get("IdType") == "pmc"), None)


def fetch_fulltext(pmids: list[str]) -> None:
    out = CACHE / "fulltext"
    out.mkdir(parents=True, exist_ok=True)
    pairs = [(p, pmc_id(p)) for p in pmids]
    pairs = [(p, c) for p, c in pairs if c]
    got = missing = 0
    for pmid, pmc in pairs:
        if (out / f"{pmid}.xml").exists():
            got += 1
            continue
        body = get(f"{EPMC}/{pmc}/fullTextXML")
        if body and "<body" in body:
            (out / f"{pmid}.xml").write_text(body)
            got += 1
        else:
            # In PMC but outside the open-access subset: the publisher does not allow the
            # body to be downloaded, and only the probes below can see into it.
            missing += 1
        time.sleep(0.3)
    print(f"[fulltext] {len(pairs)} have a PMC id; {got} full texts cached, "
          f"{missing} not open access")


def hits(query: str) -> int | None:
    body = get(f"{EPMC}/search?query={urllib.parse.quote(query)}"
               f"&format=json&resultType=idlist", timeout=45)
    return json.loads(body)["hitCount"] if body else None


def fetch_probes(pmids: list[str]) -> None:
    out = CACHE / "probe"
    out.mkdir(parents=True, exist_ok=True)
    print(f"[probe]   {len(pmids)} sources x {len(ROUTE_TERMS)} terms x 2 sections")
    for pmid in pmids:
        path = out / f"{pmid}.json"
        if path.exists():
            continue
        data: dict = {"has_ft": hits(f"EXT_ID:{pmid} AND HAS_FT:Y")}
        if data["has_ft"]:
            for field in ("METHODS", "BODY"):
                data[f"{field}_route"] = {}
                for term in ROUTE_TERMS:
                    data[f"{field}_route"][term] = hits(
                        f'EXT_ID:{pmid} AND {field}:"{term}"')
                    time.sleep(0.15)
        path.write_text(json.dumps(data, indent=1))


def fetch_patents() -> None:
    out = CACHE / "patent"
    out.mkdir(parents=True, exist_ok=True)
    print(f"[patent]  {len(PATENT_IDS)} documents")
    for key, doc in PATENT_IDS.items():
        html = out / f"{doc}.html"
        if not html.exists():
            body = get(PATENTS.format(doc))
            if body is None:
                print(f"  {key}: not retrievable as {doc}")
                continue
            html.write_text(body)
        # A plain-text rendering beside it, which is what the quotes were read from.
        text = out / f"{doc}.txt"
        if not text.exists():
            raw = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", html.read_text(),
                         flags=re.S | re.I)
            text.write_text(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", raw)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true",
                    help="also run the Europe PMC section probes (many small requests)")
    args = ap.parse_args()

    identifiers = sources()
    pmids = [s for s in identifiers if s.isdigit()]
    CACHE.mkdir(parents=True, exist_ok=True)
    fetch_pubmed(pmids)
    fetch_fulltext(pmids)
    fetch_patents()
    if args.probe:
        fetch_probes(pmids)
    print(f"\ncached under {CACHE.relative_to(HERE.parent)}")


if __name__ == "__main__":
    main()
