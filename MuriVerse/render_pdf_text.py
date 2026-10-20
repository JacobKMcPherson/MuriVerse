"""Render publisher PDFs as body-only extraction text, for studies absent from PMC.

`render_fulltext.py` only covers the PMC OA BioC deposits. 61 dosed mouse publications have a
PEPlife2-mirrored PDF but no BioC record, so they were uncurateable and contributed no Cmax or
tmax. This produces the same layout `render_fulltext.py` does - PEPlife2 context header, fenced
abstract, `=== BODY ===` - from the PDF instead.

Output goes to `<species>/pdf_txt/`, deliberately NOT `fulltext_txt/`, because PDF-derived text
is a weaker source than a BioC deposit and `evidence` / `fulltext_status` must keep meaning what
they say. Two known fidelity limits are stated in each file's header so the curator can see them:
two-column reading order is reconstructed heuristically, and the abstract fence depends on
finding an "Introduction" heading.

  python render_pdf_text.py --species mouse
"""
import argparse
import re
from pathlib import Path

import pandas as pd
import pymupdf

ROOT = Path(__file__).resolve().parent
PAPERS = ROOT / "data" / "peplife2" / "papers"

INTRO = re.compile(r"^\s*(?:1\.?\s*|I\.\s*)?(?:introduction|background)\b\s*$", re.I)
REFS = re.compile(r"^\s*(?:\d+\.?\s*)?(?:references|bibliography|literature cited|"
                  r"references and notes)\s*:?\s*$", re.I)
KEYWORDS = re.compile(r"^\s*(?:key\s?words|abbreviations|graphical abstract)\b", re.I)
COLUMN_SHARE = 0.2  # each side must hold this much of the page's text to count as two columns


def page_text(page: pymupdf.Page) -> str:
    """Page text in reading order, column-aware.

    A plain y-sort interleaves the two columns of a typical journal page, which scrambles every
    sentence and makes the text useless as evidence.
    """
    blocks = [b for b in page.get_text("blocks") if (b[4] or "").strip()]
    if not blocks:
        return ""
    mid = (page.rect.x0 + page.rect.x1) / 2
    left = [b for b in blocks if (b[0] + b[2]) / 2 < mid]
    right = [b for b in blocks if (b[0] + b[2]) / 2 >= mid]
    total = sum(len(b[4]) for b in blocks)
    two_col = (left and right
               and sum(len(b[4]) for b in left) > COLUMN_SHARE * total
               and sum(len(b[4]) for b in right) > COLUMN_SHARE * total
               and max(b[2] for b in left) <= mid + 12)
    order = (sorted(left, key=lambda b: b[1]) + sorted(right, key=lambda b: b[1])
             if two_col else sorted(blocks, key=lambda b: (b[1], b[0])))
    return "\n".join((b[4] or "").strip() for b in order)


def split_abstract(lines: list[str]) -> tuple[list[str], list[str]]:
    """(abstract, body). Everything before the Introduction heading is treated as front matter."""
    for i, ln in enumerate(lines[:400]):
        if INTRO.match(ln):
            return lines[:i], lines[i:]
    for i, ln in enumerate(lines[:400]):
        if KEYWORDS.match(ln):
            return lines[:i], lines[i:]
    return [], lines


def drop_references(lines: list[str]) -> list[str]:
    for i in range(len(lines) - 1, -1, -1):
        if REFS.match(lines[i]) and i > 0.4 * len(lines):
            return lines[:i]
    return lines


def render(pdf: Path) -> tuple[str, str]:
    doc = pymupdf.open(pdf)
    pages = [page_text(doc[i]) for i in range(len(doc))]
    doc.close()

    lines = [ln.rstrip() for ln in "\n".join(pages).splitlines()]
    lines = [ln for ln in lines if ln.strip()]
    abstract, body = split_abstract(lines)
    body = drop_references(body)
    fenced = "\n".join(abstract) if abstract else (
        "[abstract not separable from this PDF - no Introduction heading found. "
        "The body below therefore still contains the abstract; identify and skip it.]")
    note = "abstract fenced" if abstract else "ABSTRACT NOT FENCED"
    return ("<<<ABSTRACT - DO NOT USE AS AN EXTRACTION SOURCE>>>\n" + fenced
            + "\n<<<END ABSTRACT>>>\n\n=== BODY ===\n" + "\n".join(body) + "\n"), note


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", default="mouse", choices=["mouse", "rat"])
    ap.add_argument("--pmid", nargs="*", help="limit to these PMIDs")
    ap.add_argument("--all-dosed", action="store_true",
                    help="include dosed studies that already have BioC text")
    args = ap.parse_args()

    sdir = ROOT / args.species
    rows = pd.read_csv(sdir / "peplife2_rows.csv", low_memory=False)
    rows["PMID"] = rows["PMID"].astype(str)
    keep = ["Sequence", "Name", "Half Life", "Units Half Life", "Test Sample",
            "Vivo/Vitro", "Assay", "Dose Route", "Concentration"]

    prov = pd.read_csv(sdir / "records_provenance.csv", low_memory=False)
    dosed = prov[prov["in_cohort"].astype(bool) & prov["animal_dosed"].eq("dosed")]
    targets = args.pmid or sorted({str(p) for p in dosed["PMID"]})
    if not args.pmid and not args.all_dosed:
        targets = [p for p in targets if not (sdir / "fulltext_txt" / f"{p}.txt").exists()]

    outdir = sdir / "pdf_txt"
    outdir.mkdir(exist_ok=True)

    written, missing, notes = 0, 0, []
    for pmid in targets:
        pdf = PAPERS / f"{pmid}.pdf"
        if not pdf.exists():
            missing += 1
            continue
        try:
            text, note = render(pdf)
        except Exception as exc:
            notes.append((pmid, f"unreadable: {exc}"))
            continue
        sub = rows[rows["PMID"] == pmid][keep].drop_duplicates()
        header = (f"=== SOURCE: publisher PDF (not a PMC deposit) ===\n"
                  f"Text was extracted from {pdf.name}. Reading order for two-column pages is\n"
                  f"reconstructed heuristically and may still be wrong in places; table cells may\n"
                  f"be split across lines. Verify any value that looks garbled before recording it.\n"
                  f"Status: {note}.\n\n"
                  f"=== PEPLIFE2 ROWS FOR THIS PMID ({len(sub)} distinct) ===\n"
                  f"{sub.to_csv(index=False)}\n\n")
        (outdir / f"{pmid}.txt").write_text(header + text)
        written += 1
        if note != "abstract fenced":
            notes.append((pmid, note))

    print(f"[write] {written} PDFs rendered to {outdir.relative_to(ROOT)}/")
    if missing:
        print(f"[skip]  {missing} targets have no publisher PDF")
    unfenced = [p for p, n in notes if "NOT FENCED" in n]
    if unfenced:
        print(f"[warn]  {len(unfenced)} have no separable abstract: {' '.join(unfenced)}")
    for pmid, note in notes:
        if "unreadable" in note:
            print(f"[error] {pmid}: {note}")


if __name__ == "__main__":
    main()
