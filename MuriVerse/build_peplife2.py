"""Step 0 - build the PEPlife 2.0 record table from the database's own REST API.

PEPlife 2.0 documents a REST endpoint at `api/api.php` (see `api/rest.html`) with three query
fields: `organism`, `lin_cyc` and `seq`. Issuing `dataType=organism` with an empty `dataValue`
returns the complete record set, so the whole database can be pulled in one request instead of
stitching together the site's paginated .xlsx exports. The raw response is archived verbatim
so every later step is reproducible without the network.

The API response is the source of truth here. It carries `dose_route`, which the .xlsx exports
omit, but not the exports' precomputed `half-life(in Sec)`; that column is re-derived from
`Half Life` + `Units Half Life`.

Keep `ID` as the only join key for this table. Merging on any shared attribute duplicates
entries many times over and silently weights every downstream mean, median and quantile.

  python build_peplife2.py
  python build_peplife2.py --offline   # rebuild the CSV from the archived response
"""
import argparse
import json
import re
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data/peplife2"
ARCHIVE = DATA / "api/records_all.json"
OUT = DATA / "peplife2_records.csv"
API = "https://webs.iiitd.edu.in/raghava/peplife2/api/api.php"

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "MuriVerse/1.0 (research; muriverse@example.org)"

# API field -> the column name used by every downstream consumer of peplife2_records.csv.
COLUMNS = {
    "id": "ID", "pmid": "PMID", "year": "YEAR", "seq": "Sequence", "name": "Name",
    "length": "Length", "nter": "N-ter MOD", "cter": "C-ter MOD", "lin_cyc": "Linear/Cyclic",
    "chiral": "Chirality", "chem_mod": "Chem-MOD", "origin": "Origin", "nature": "Nature",
    "incubation_time": "Incubation Time", "conc": "Concentration", "half_life": "Half Life",
    "units_half": "Units Half Life", "protease": "Protease", "assay": "Assay",
    "test_sample": "Test Sample", "vivo_vitro": "Vivo/Vitro", "ref": "Reference",
    "patent_no": "Patent No.", "activity_sub": "Activity", "dose_route": "Dose Route",
}

UNIT_SECONDS = {"second": 1, "minute": 60, "hour": 3600, "day": 86400,
                "week": 604800, "month": 2592000, "year": 31536000}

# `Half Life` is free text and some entries are mis-encoded in the database itself
# ("5 ·8 hours", "1.5 €“5 minutes", and ranges Excel turned into dates: "03-Apr" is 3-4).
MOJIBAKE = str.maketrans({"·": ".", "–": "-", "—": "-", "\u00a0": " "})
MONTHS = {m: i for i, m in enumerate(
    "jan feb mar apr may jun jul aug sep oct nov dec".split(), start=1)}
# Only a range that opens the string is a range; "11.2 (t1/2 after 10-60 min)" is not.
RANGE = re.compile(r"^\s*[~<>≈]?\s*(\d+(?:\.\d+)?)\s*(?:-|to|€“)\s*(\d+(?:\.\d+)?)")
NUMBER = re.compile(r"\d+(?:\.\d+)?")


def fetch_all(timeout: int = 300) -> dict:
    r = SESSION.get(API, params={"dataType": "organism", "dataValue": ""}, timeout=timeout)
    r.raise_for_status()
    payload = r.json()
    if "data" not in payload:
        raise SystemExit(f"API response has no 'data' key: {payload}")
    return payload


def half_life_seconds(value: object, unit: object) -> float | None:
    """Reported half-life in seconds: midpoint of a leading range, otherwise its first number.

    Stands in for the site's own `half-life(in Sec)` column, which the API does not return.
    A `±` tolerance is dropped rather than read as the far end of a range.
    """
    s = str(value or "").translate(MOJIBAKE)
    s = re.sub(r"(?:±|\+/-)\s*\d+(?:\.\d+)?", " ", s)
    s = re.sub(r"\s+\.\s*(\d)", r".\1", s)  # "14 .9" -> "14.9"
    s = re.sub(r"\b([a-z]{3})[a-z]*-(\d+)\b",
               lambda m: f"{MONTHS[m.group(1).lower()]}-{m.group(2)}"
               if m.group(1).lower() in MONTHS else m.group(0), s, flags=re.I)
    s = re.sub(r"\b(\d+)-([a-z]{3})[a-z]*\b",
               lambda m: f"{m.group(1)}-{MONTHS[m.group(2).lower()]}"
               if m.group(2).lower() in MONTHS else m.group(0), s, flags=re.I)
    u = re.sub(r"[^a-z]", "", str(unit or "").lower()).rstrip("s")
    if u not in UNIT_SECONDS:
        return None
    span = RANGE.match(s)
    if span:
        return (float(span.group(1)) + float(span.group(2))) / 2 * UNIT_SECONDS[u]
    first = NUMBER.search(s)
    return float(first.group(0)) * UNIT_SECONDS[u] if first else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true",
                    help="skip the network and rebuild from the archived API response")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()

    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    if args.offline:
        if not ARCHIVE.exists():
            raise SystemExit(f"--offline but {ARCHIVE} is not archived yet")
        payload = json.loads(ARCHIVE.read_text())
        print(f"[api] {payload.get('count')} records (archived copy)")
    else:
        payload = fetch_all()
        ARCHIVE.write_text(json.dumps(payload, indent=1))
        print(f"[api] {payload.get('count')} records (fetched) -> {ARCHIVE}")

    df = pd.DataFrame(payload["data"])
    if payload.get("count") is not None and len(df) != int(payload["count"]):
        raise SystemExit(f"response count {payload['count']} != {len(df)} rows")
    missing = [k for k in COLUMNS if k not in df.columns]
    if missing:
        raise SystemExit(f"API response is missing fields {missing}")

    out = df[list(COLUMNS)].rename(columns=COLUMNS)
    # The API returns padded text and a few stray C0 control bytes inside sequences.
    out = out.apply(lambda c: c.str.replace(r"[\x00-\x08\x0b-\x1f]", "", regex=True).str.strip()
                    if c.dtype == object else c)
    out["ID"] = out["ID"].astype(int)
    if out["ID"].duplicated().any():
        raise SystemExit("the API returned duplicate IDs")
    out["Length"] = pd.to_numeric(out["Length"], errors="coerce").astype("Int64")
    out["half-life(in Sec)"] = [half_life_seconds(v, u)
                                for v, u in zip(out["Half Life"], out["Units Half Life"])]
    out = out.sort_values("ID").reset_index(drop=True)
    out.to_csv(args.out, index=False)
    print(f"[done] {len(out)} records, IDs {out['ID'].min()}-{out['ID'].max()}, "
          f"{out['PMID'].nunique()} distinct PMIDs -> {args.out}")


if __name__ == "__main__":
    main()
