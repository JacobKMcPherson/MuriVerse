"""Select a species cohort from PEPlife 2.0's own organism API.

The database exposes `api/api.php?dataType=organism&dataValue=<organism>`, which is the
authoritative membership test for a species - no free-text matching on our side. PEPlife 2.0
splits the mouse literature over two organism values, `Mouse` and `Mice`, so the mouse cohort
is the union of both queries; rat is `Rat` + `Rats`.

Raw responses are archived verbatim under `data/peplife2/api/` and the merged selection is
written next to the record table. The API also carries `dose_route`, a route-of-administration
field absent from the site's .xlsx exports.

  python fetch_peplife2_api.py --species mouse
  python fetch_peplife2_api.py --species mouse --offline   # merge archived responses only
"""
import argparse
import json
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data/peplife2"
ARCHIVE = DATA / "api"
API = "https://webs.iiitd.edu.in/raghava/peplife2/api/api.php"

ORGANISM_QUERIES = {"mouse": ["Mouse", "Mice"], "rat": ["Rat", "Rats"]}

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "MuriVerse/1.0 (research; muriverse@example.org)"


def fetch_organism(value: str, timeout: int = 120, retries: int = 2) -> dict | None:
    for attempt in range(retries + 1):
        try:
            r = SESSION.get(API, params={"dataType": "organism", "dataValue": value},
                            timeout=timeout)
            r.raise_for_status()
            payload = r.json()
        except (requests.RequestException, ValueError) as e:
            print(f"  [warn] {value}: {type(e).__name__}")
            if attempt == retries:
                return None
            time.sleep(5 * (attempt + 1))
            continue
        if "data" not in payload:
            print(f"  [warn] {value}: response has no 'data' key")
            return None
        return payload
    return None


def load_archive(value: str) -> dict | None:
    f = ARCHIVE / f"organism_{value}.json"
    return json.loads(f.read_text()) if f.exists() else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", choices=sorted(ORGANISM_QUERIES), default="mouse")
    ap.add_argument("--organism", nargs="*", default=None,
                    help="override the organism values queried for this species")
    ap.add_argument("--offline", action="store_true",
                    help="skip the network and merge whatever is already archived")
    ap.add_argument("--allow-partial", action="store_true",
                    help="write the selection even if some organism queries are still missing")
    args = ap.parse_args()

    ARCHIVE.mkdir(parents=True, exist_ok=True)
    queries = args.organism or ORGANISM_QUERIES[args.species]

    frames, missing = [], []
    for value in queries:
        payload = None if args.offline else fetch_organism(value)
        if payload is not None:
            (ARCHIVE / f"organism_{value}.json").write_text(json.dumps(payload, indent=1))
            print(f"[api] {value}: {payload.get('count')} records (fetched)")
        else:
            payload = load_archive(value)
            if payload is None:
                missing.append(value)
                print(f"[api] {value}: unavailable and not archived")
                continue
            print(f"[api] {value}: {payload.get('count')} records (archived copy)")
        df = pd.DataFrame(payload["data"])
        if payload.get("count") is not None and len(df) != int(payload["count"]):
            raise SystemExit(f"{value}: response count {payload['count']} != {len(df)} rows")
        frames.append(df.assign(api_query=value))

    if missing and not args.allow_partial:
        raise SystemExit(f"[abort] no data for {missing}; re-run when the server responds, "
                         f"or pass --allow-partial")
    if missing:
        print(f"[api] PARTIAL selection: {missing} still missing; entries not covered by the "
              f"API stay text-matched downstream")

    out = pd.concat(frames, ignore_index=True)
    out["id"] = out["id"].astype(int)
    dup = out["id"].duplicated()
    if dup.any():
        # Queries are expected to partition the cohort; keep the first and say so.
        print(f"[api] {dup.sum()} ids returned by more than one query, keeping first")
        out = out[~dup]
    out = out.sort_values("id").reset_index(drop=True)

    dest = DATA / f"peplife2_api_{args.species}.csv"
    out.to_csv(dest, index=False)
    (ARCHIVE / f"queries_{args.species}.json").write_text(json.dumps(
        {"queries": queries, "archived": [q for q in queries if q not in missing],
         "missing": missing}, indent=1))
    print(f"[done] {len(out)} entries from {', '.join(q for q in queries if q not in missing)}"
          f" -> {dest}")
    print(out["api_query"].value_counts().to_string())


if __name__ == "__main__":
    main()
