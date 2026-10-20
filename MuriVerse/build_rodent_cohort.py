"""Step 1 - build a rodent cohort from PEPlife2 and fetch full text for each PMID.

Cohort membership comes from PEPlife 2.0's own organism API (see fetch_peplife2_api.py): the
mouse cohort is the union of the `Mouse` and `Mice` organism queries, the rat cohort of `Rat`
and `Rats`. If those responses have not been archived yet, the script falls back to matching
the rodent name in the `Test Sample` free text and says so; the two agree exactly for the
organism values verified so far.

Every extraction downstream must come from the article body, so this script also resolves each
PMID to a PMC full-text record and stores it, recording per-PMID whether retrieval succeeded.

  python fetch_peplife2_api.py --species mouse     # defines the cohort
  python build_rodent_cohort.py --species mouse
"""
import argparse
import json
import os
import re
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data/peplife2"
PEPLIFE = DATA / "peplife2_records.csv"
ELINK = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
BIOC = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{}/unicode"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = {"tool": "MuriVerse", "email": "muriverse@example.org"}

# NCBI returns 403 to the default requests user-agent.
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "MuriVerse/1.0 (research; muriverse@example.org)"

MOUSE = re.compile(r"\bmice\b|\bmouse\b", re.I)
RAT = re.compile(r"\brats?\b|wistar|sprague", re.I)
# Matrix strings that name a rodent but are cell lines, not the animal.
MOUSE_CELL_LINE = re.compile(r"\b(?:3t3|raw ?264|l929|b16|neuro-?2a|c2c12|hepa|p388|ct26|4t1|mc38)\b", re.I)
RAT_CELL_LINE = re.compile(r"\b(?:pc-?12|h9c2|gh3|rin-?m|ins-?1|c6)\b", re.I)

SPECIES = {"mouse": (MOUSE, MOUSE_CELL_LINE), "rat": (RAT, RAT_CELL_LINE)}


def rodent_cohort(df: pd.DataFrame, species: str) -> pd.DataFrame:
    """Cohort by PEPlife2 organism API where available, else by `Test Sample` text."""
    animal, cell_line = SPECIES[species]
    api_csv = DATA / f"peplife2_api_{species}.csv"
    text_hit = df["Test Sample"].fillna("").str.contains(animal)

    if api_csv.exists():
        api = pd.read_csv(api_csv, low_memory=False)
        api_ids = set(api["id"])
        selected = df[df["ID"].isin(api_ids)].copy()
        # api_query records which organism spelling (e.g. Mouse vs Mice) returned the entry.
        selected = selected.merge(
            api[["id", "api_query"]].rename(columns={"id": "ID"}), on="ID", how="left")
        selected["cohort_source"] = "peplife2 organism API"

        # A partially archived API must not silently shrink the cohort.
        remainder = df[text_hit & ~df["ID"].isin(api_ids)].copy()
        if len(remainder):
            remainder["api_query"] = pd.NA
            remainder["cohort_source"] = "Test Sample text match (organism API pending)"
            print(f"  [warn] API covers {len(selected)} entries; {len(remainder)} more are "
                  f"text-matched until every organism query is archived")
        out = pd.concat([selected, remainder], ignore_index=True).sort_values("ID")

        unmatched = api_ids - set(df["ID"])
        if unmatched:
            print(f"  [warn] {len(unmatched)} API ids absent from the record table")
    else:
        out = df[text_hit].copy()
        out["api_query"] = pd.NA
        out["cohort_source"] = "Test Sample text match (organism API not archived)"
        print(f"  [warn] {api_csv.name} missing; falling back to text matching. "
              f"Run: python fetch_peplife2_api.py --species {species}")

    out["matrix_is_cell_line"] = out["Test Sample"].fillna("").str.contains(cell_line)
    out["pmid_valid"] = out["PMID"].astype(str).str.fullmatch(r"\d+")
    return out


def resolve_pmcids(pmids: list[str], api_key: str = "") -> dict[str, dict]:
    """PMID -> {pmcid} via elink; an empty linkset means the article has no PMC record."""
    found = {}
    for i in range(0, len(pmids), 100):
        chunk = pmids[i : i + 100]
        # Repeated `id=` params keep one linkset per PMID; a comma list merges them all.
        q = "&".join(f"id={p}" for p in chunk)
        params = {"dbfrom": "pubmed", "db": "pmc", "retmode": "json",
              "linkname": "pubmed_pmc", **TOOL}
        # An empty api_key= is rejected with HTTP 400, so only send it when set.
        if api_key:
            params["api_key"] = api_key
        r = SESSION.get(f"{ELINK}?{q}", params=params, timeout=90)
        r.raise_for_status()
        for ls in r.json().get("linksets", []):
            links = [l for db in ls.get("linksetdbs", []) for l in db.get("links", [])]
            found[ls["ids"][0]] = {"pmcid": f"PMC{links[0]}" if links else None}
        print(f"  [elink] {min(i + 100, len(pmids))}/{len(pmids)}")
        time.sleep(0.5)
    return found


def fetch_bioc(pmcid: str, tries: int = 3) -> tuple[dict | None, str]:
    """Fetch a PMC OA BioC document, returning (doc, reason).

    The service answers both throttles and genuine misses in ways a status check cannot
    separate: a miss is an `[Error] : No result can be found` page served as HTTP 200 whose
    text still starts with `[`. Collapsing every failure to one reason cached transient 429s
    as permanent, so each outcome is reported distinctly and 429/5xx are retried.
    """
    for attempt in range(tries):
        r = SESSION.get(BIOC.format(pmcid), timeout=120)
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(5 * (attempt + 1))
            continue
        if r.status_code != 200:
            return None, f"HTTP {r.status_code} ({pmcid})"
        if r.text.lstrip().startswith("[Error]"):
            return None, f"not in PMC OA subset ({pmcid})"
        try:
            doc = r.json()[0]["documents"][0]
        except (ValueError, KeyError, IndexError):
            return None, f"unparseable BioC ({pmcid})"
        # A body-less record is an abstract-only deposit and is unusable for this task.
        sections = {p.get("infons", {}).get("section_type", "") for p in doc.get("passages", [])}
        if not sections - {"TITLE", "ABSTRACT", "REF", "COMP_INT", "AUTH_CONT", "", "SUPPL"}:
            return None, f"abstract-only deposit ({pmcid})"
        return doc, f"PMC OA BioC ({pmcid})"
    return None, f"throttled, retry needed ({pmcid})"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", choices=sorted(SPECIES), default="mouse")
    ap.add_argument("--outdir", default="")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    ap.add_argument("--skip-existing", action="store_true",
                    help="keep PMIDs already recorded in fulltext_retrieval_status.csv")
    args = ap.parse_args()
    out = Path(args.outdir or ROOT / args.species)
    (out / "fulltext").mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(PEPLIFE, low_memory=False)
    cohort = rodent_cohort(df, args.species)
    cohort.to_csv(out / "peplife2_rows.csv", index=False)
    mix = "; ".join(f"{n} via {src}" for src, n in cohort["cohort_source"].value_counts().items())
    print(f"[cohort] {len(cohort)} of {len(df)} PEPlife2 entries ({mix}), "
          f"{cohort['PMID'].nunique()} distinct PMID strings")

    pmids = sorted({p for p in cohort.loc[cohort["pmid_valid"], "PMID"].astype(str)})
    status_csv = out / "fulltext_retrieval_status.csv"
    done = pd.DataFrame()
    if args.skip_existing and status_csv.exists():
        done = pd.read_csv(status_csv, dtype={"PMID": str})
        pmids = [p for p in pmids if p not in set(done["PMID"])]
    print(f"[cohort] {len(pmids)} numeric PMIDs to retrieve")
    ids = resolve_pmcids(pmids, args.api_key) if pmids else {}

    status = []
    for n, pmid in enumerate(pmids, 1):
        rec = ids.get(pmid, {})
        pmcid, doc, why = rec.get("pmcid"), None, "no PMC record"
        if pmcid:
            doc, why = fetch_bioc(pmcid)
            time.sleep(0.4)
        if doc:
            (out / "fulltext" / f"{pmid}.json").write_text(json.dumps(doc))
            chars = sum(len(p.get("text", "")) for p in doc["passages"])
            state, source = "full_text_retrieved", why
        else:
            chars, state, source = 0, "full_text_unavailable", why
        status.append({"PMID": pmid, "pmcid": pmcid,
                       "fulltext_retrieved": bool(doc), "fulltext_status": state,
                       "fulltext_source": source, "fulltext_chars": chars})
        if n % 25 == 0 or n == len(pmids):
            print(f"  [fetch] {n}/{len(pmids)}")

    st = pd.concat([done, pd.DataFrame(status)], ignore_index=True)
    st = st[st["PMID"].astype(str).isin(set(cohort["PMID"].astype(str)))].copy()
    counts = cohort["PMID"].astype(str).value_counts()
    st["n_peplife_entries"] = st["PMID"].astype(str).map(counts).fillna(0).astype(int)
    st = st.drop(columns=["n_peplife_rows"], errors="ignore")

    # A stored BioC document is the ground truth; a resumed run must not carry a stale "unavailable"
    for i, pmid in st["PMID"].astype(str).items():
        f = out / "fulltext" / f"{pmid}.json"
        if not f.exists():
            continue
        doc = json.loads(f.read_text())
        st.loc[i, ["fulltext_retrieved", "fulltext_status", "fulltext_source", "fulltext_chars"]] = [
            True, "full_text_retrieved", f"PMC OA BioC ({st.loc[i, 'pmcid']})",
            sum(len(p.get("text", "")) for p in doc.get("passages", []))]
    st["fulltext_retrieved"] = st["fulltext_retrieved"].astype(bool)
    st["fulltext_chars"] = st["fulltext_chars"].astype(int)

    st.sort_values("PMID").to_csv(status_csv, index=False)
    print(f"[done] full text retrieved for {st['fulltext_retrieved'].sum()}/{len(st)} PMIDs")
    missed = st.loc[~st["fulltext_retrieved"], "fulltext_source"].str.replace(
        r" \(PMC?\w*\)$", "", regex=True).value_counts()
    if len(missed):
        print("[missed] " + "; ".join(f"{n} {why}" for why, n in missed.items()))
    retry = int(st["fulltext_source"].str.startswith("throttled").sum())
    if retry:
        print(f"[warn] {retry} PMIDs were throttled, not genuinely unavailable - re-run to collect")


if __name__ == "__main__":
    main()
