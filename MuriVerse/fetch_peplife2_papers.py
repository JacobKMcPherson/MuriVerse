"""Fetch the PEPlife 2.0 source PDF for every publication in the record table.

PEPlife 2.0 mirrors one PDF per PubMed ID at
`https://webs.iiitd.edu.in/raghava/peplife2/Papers/<PMID>.pdf`. This downloads them to
`data/peplife2/papers/` and records the outcome per PMID.

Re-runnable: PDFs already on disk are skipped, so an interrupted run resumes where it stopped.

  python fetch_peplife2_papers.py
  python fetch_peplife2_papers.py --species mouse   # cohort subset only
"""
import argparse
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
PEPLIFE = ROOT / "data/peplife2/peplife2_records.csv"
OUTDIR = ROOT / "data/peplife2/papers"
URL = "https://webs.iiitd.edu.in/raghava/peplife2/Papers/{}.pdf"

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "MuriVerse/1.0 (research; muriverse@example.org)"


def fetch(pmid: str, outdir: Path, retries: int = 2, delay: float = 1.0) -> dict:
    dest = outdir / f"{pmid}.pdf"
    if dest.exists() and dest.stat().st_size > 0:
        return {"PMID": pmid, "status": "already_present", "http_code": 200,
                "bytes": dest.stat().st_size}
    for attempt in range(retries + 1):
        try:
            r = SESSION.get(URL.format(pmid), timeout=120)
        except requests.RequestException as e:
            if attempt == retries:
                return {"PMID": pmid, "status": f"error: {type(e).__name__}", "http_code": 0,
                        "bytes": 0}
            time.sleep(5 * (attempt + 1))
            continue
        finally:
            time.sleep(delay)
        if r.status_code == 404:
            return {"PMID": pmid, "status": "not_on_server", "http_code": 404, "bytes": 0}
        # The server answers a missing paper with an HTML error page, not a PDF.
        if r.status_code == 200 and r.content[:5] == b"%PDF-":
            dest.write_bytes(r.content)
            return {"PMID": pmid, "status": "downloaded", "http_code": 200, "bytes": len(r.content)}
        if attempt == retries:
            return {"PMID": pmid, "status": "not_a_pdf", "http_code": r.status_code,
                    "bytes": len(r.content)}
        time.sleep(2 * (attempt + 1))
    return {"PMID": pmid, "status": "unknown", "http_code": 0, "bytes": 0}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", choices=["mouse", "rat"], default="",
                    help="restrict to a rodent cohort instead of the whole database")
    ap.add_argument("--outdir", default=str(OUTDIR))
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--delay", type=float, default=1.0,
                    help="seconds to pause after each request; the host is a small university server")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.species:
        src = ROOT / f"{args.species}/peplife2_rows.csv"
    else:
        src = PEPLIFE
    df = pd.read_csv(src, low_memory=False)
    pmid = df["PMID"].astype(str).str.strip()
    counts = pmid.value_counts()
    pmids = sorted(p for p in pmid.unique() if p.isdigit())
    print(f"[papers] {len(pmids)} numeric PMIDs from {src.name} "
          f"({pmid.nunique() - len(pmids)} non-numeric identifiers skipped)")

    rows, done = [], 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for rec in pool.map(lambda p: fetch(p, outdir, delay=args.delay), pmids):
            rows.append(rec)
            done += 1
            if done % 50 == 0 or done == len(pmids):
                print(f"  [fetch] {done}/{len(pmids)}")

    man = pd.DataFrame(rows)
    man["n_peplife_entries"] = man["PMID"].map(counts).fillna(0).astype(int)
    man.sort_values("PMID").to_csv(outdir / "paper_download_status.csv", index=False)
    print(man["status"].value_counts().to_string())
    print(f"[done] {man['bytes'].sum() / 1e9:.2f} GB in {outdir}")


if __name__ == "__main__":
    main()
