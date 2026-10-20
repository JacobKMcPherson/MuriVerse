"""Step 2 - render retrieved BioC full text as plain text for body-only extraction.

References are dropped (they add bulk and no evidence) and the abstract is emitted in a clearly
fenced block so the extraction pass can honour the body-only rule. Figure and table captions are
kept because concentration-time profiles and PK tables live there.

  python render_fulltext.py
"""
import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DROP = {"REF", "ACK_FUND", "AUTH_CONT", "COMP_INT", "SUPPL"}


def render(doc: dict) -> str:
    front = doc.get("passages", [{}])[0].get("infons", {})
    head = [f"PMID: {front.get('article-id_pmid', '?')}   PMC: {doc.get('id', '?')}",
            f"JOURNAL: {front.get('journal', '?')}   YEAR: {front.get('year', '?')}"]
    body, abstract, last = [], [], None
    for p in doc.get("passages", []):
        infons, text = p.get("infons", {}), (p.get("text") or "").strip()
        sec, typ = infons.get("section_type", ""), infons.get("type", "")
        if not text or sec in DROP or typ.startswith("ref"):
            continue
        if sec == "ABSTRACT":
            abstract.append(text)
            continue
        if sec == "TITLE":
            head.append(f"TITLE: {text}")
            continue
        if sec != last:
            body.append(f"\n## [{sec}]")
            last = sec
        if typ.startswith("title"):
            body.append(f"\n### {text}")
        elif typ in ("fig_caption", "table_caption"):
            body.append(f"[{typ.upper()} {infons.get('id', '')}] {text}")
        else:
            body.append(text)
    return "\n".join(head) + "\n\n<<<ABSTRACT - DO NOT USE AS AN EXTRACTION SOURCE>>>\n" \
        + "\n".join(abstract) + "\n<<<END ABSTRACT>>>\n\n=== BODY ===" + "\n".join(body) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(ROOT / "mouse"))
    args = ap.parse_args()
    base = Path(args.dir)
    outdir = base / "fulltext_txt"
    outdir.mkdir(exist_ok=True)

    rows = pd.read_csv(base / "peplife2_rows.csv", low_memory=False)
    rows["PMID"] = rows["PMID"].astype(str)
    keep = ["Sequence", "Name", "Half Life", "Units Half Life", "Test Sample", "Vivo/Vitro", "Assay"]
    cohort_pmids = set(rows["PMID"])

    skipped = 0
    for f in sorted((base / "fulltext").glob("*.json")):
        pmid = f.stem
        if pmid not in cohort_pmids:
            skipped += 1
            continue
        text = render(json.loads(f.read_text()))
        sub = rows[rows["PMID"] == pmid][keep].drop_duplicates()
        ctx = sub.to_csv(index=False)
        (outdir / f"{pmid}.txt").write_text(
            f"=== PEPLIFE2 ROWS FOR THIS PMID ({len(sub)} distinct) ===\n{ctx}\n\n{text}")
    if skipped:
        print(f"[skip] {skipped} stored full texts are not in the current cohort")
    print(f"[done] rendered {len(list(outdir.glob('*.txt')))} full texts to {outdir}")


if __name__ == "__main__":
    main()
