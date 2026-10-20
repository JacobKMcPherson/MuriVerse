"""Crop each dosed study's concentration-time figure out of its publisher PDF.

Source is `data/peplife2/papers/<PMID>.pdf` (PEPlife2's author-manuscript mirror). The BioC
full text in `<species>/fulltext/` carries captions but no images, so the PDF is the only
figure source we hold locally.

A figure is located by its caption: a text block that both opens like a caption ("Figure 4.")
and mentions a concentration-time readout. The artwork sits above its caption, so the crop is
the union of the raster images and vector drawings that lie above the caption block and
overlap it horizontally, falling back to the whitespace gap between the caption and the
previous body-text block when a page has no extractable graphics.

  python extract_pk_figures.py --species mouse
  python extract_pk_figures.py --species mouse --pmid 36630826 --page-fallback
"""
import argparse
import re
from pathlib import Path

import pandas as pd
import pymupdf

ROOT = Path(__file__).resolve().parent
PAPERS = ROOT / "data" / "peplife2" / "papers"

CAPTION = re.compile(r"^\s*(?:fig(?:ure)?s?\.?\s*\d+|scheme\s*\d+)", re.I)

# A caption qualifies on any one of these. Kept as separate alternatives rather than a
# "concentration near time" proximity rule because captions bury the two far apart
# ("...were s.c. administrated into naive mice, followed by blood collection at specified
# time points"). Drug-concentration wording is required, so PD time courses (blood *glucose*
# over time, tumour volume, body weight) do not match.
CONC_TIME = re.compile(
    r"pharmacokinetics?"
    r"|concentrations?[\s\-]*(?:vs\.?|versus|over|against)?[\s\-]*time"
    r"|(?:plasma|serum|blood|circulating)\s+(?:concentrations?|levels?|profiles?|clearance)"
    r"|concentrations?\s+(?:in|of)\s+(?:the\s+)?(?:plasma|serum|blood|circulation)"
    r"|\bPK\s+(?:profile|stud|analys|paramet)"
    r"|%\s*(?:of\s*)?(?:the\s*)?injected\s+dose|\b%\s?ID\b"
    r"|(?:percent(?:age)?|fraction)\s+(?:of\s+\w+\s+)?remaining"
    r"|time[\s\-]courses?\s+of\s+.{0,40}(?:plasma|serum|blood|concentration)"
    r"|concentration.{0,30}as\s+a\s+function\s+of\s+time"
    r"|blood\s+collection\s+at\s+.{0,20}time",
    re.I)

RENDER_DPI = 200
PAD = 6
MIN_AREA = 8000  # ignore rules, logos and stray glyph boxes when unioning artwork
MIN_FIG_HEIGHT = 80  # pt; below this the "artwork" is a rule or axis fragment, not a plot


def norm(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def artwork_boxes(page: pymupdf.Page) -> list[pymupdf.Rect]:
    boxes = [pymupdf.Rect(i["bbox"]) for i in page.get_image_info()]
    for d in page.get_drawings():
        boxes.append(pymupdf.Rect(d["rect"]))
    return [b for b in boxes if b.is_valid and not b.is_empty and b.get_area() > MIN_AREA]


def union(boxes: list[pymupdf.Rect]) -> pymupdf.Rect | None:
    if not boxes:
        return None
    rect = boxes[0]
    for b in boxes[1:]:
        rect |= b
    return rect


def overlaps(a: pymupdf.Rect, b: pymupdf.Rect) -> bool:
    span = min(a.x1, b.x1) - max(a.x0, b.x0)
    return span > 0.3 * min(a.width, b.width)


def figure_rect(page: pymupdf.Page, cap: pymupdf.Rect) -> tuple[pymupdf.Rect | None, str]:
    """The artwork belonging to `cap`, searched upward within the caption's column."""
    art = union([b for b in artwork_boxes(page)
                 if b.y1 <= cap.y0 + 2 and overlaps(b, cap)])
    if art is not None and art.height >= MIN_FIG_HEIGHT:
        return art, "artwork"

    # Whitespace between the caption and the text above it. Also used to rescue a thin
    # artwork union, which means only a rule or an axis line cleared MIN_AREA.
    tops = [pymupdf.Rect(b[:4]) for b in page.get_text("blocks")
            if norm(b[4]) and pymupdf.Rect(b[:4]).y1 <= cap.y0 - 2]
    top = max((b.y1 for b in tops if overlaps(b, cap)), default=page.rect.y0 + 36)
    gap = pymupdf.Rect(cap.x0, top, cap.x1, cap.y0)

    if art is not None:
        return (art | gap), "artwork+gap"
    if gap.height >= 60:
        return gap, "text gap"
    return None, "no artwork found"


FIG_PAGE_LABEL = re.compile(r"^\s*fig(?:ure)?s?\.?\s*(\d+)\s*[A-Za-z]?\s*$", re.I)


def plate_index(doc: pymupdf.Document) -> dict[str, list[tuple[int, pymupdf.Rect]]]:
    """Figure number -> plate pages.

    Accepted-manuscript PDFs list every legend on one text page and put the artwork on later
    pages that carry nothing but the figure number, so a caption's own page has no graphics.
    """
    index: dict[str, list[tuple[int, pymupdf.Rect]]] = {}
    for pno in range(len(doc)):
        page = doc[pno]
        art = artwork_boxes(page)
        if not art:
            continue
        for blk in page.get_text("blocks"):
            m = FIG_PAGE_LABEL.match(norm(blk[4]))
            if not m:
                continue
            tag = pymupdf.Rect(blk[:4])
            near = [b for b in art if overlaps(b, tag) or b.get_area() > 0.15 * page.rect.get_area()]
            rect = union(near or art)
            if rect is not None and rect.height >= MIN_FIG_HEIGHT:
                index.setdefault(m.group(1), []).append((pno, rect))
    return index



def label_of(caption: str) -> str:
    m = re.match(r"\s*(fig(?:ure)?s?\.?\s*\d+|scheme\s*\d+)", caption, re.I)
    raw = m.group(1) if m else "figure"
    return re.sub(r"[^a-z0-9]+", "", raw.lower())


def harvest(pmid: str, out_dir: Path, page_fallback: bool) -> list[dict]:
    pdf = PAPERS / f"{pmid}.pdf"
    rows: list[dict] = []
    try:
        doc = pymupdf.open(pdf)
    except Exception as exc:
        return [{"PMID": pmid, "status": f"unreadable pdf: {exc}"}]

    plates = plate_index(doc)
    seen: set[str] = set()
    for pno in range(len(doc)):
        page = doc[pno]
        blocks = page.get_text("blocks")
        for i, blk in enumerate(blocks):
            caption = norm(blk[4])
            if not CAPTION.match(caption):
                continue
            # Some layouts split "Figure 4." from its text; look ahead for the match only.
            window = caption
            for nxt in blocks[i + 1:i + 3]:
                if len(window) > 200:
                    break
                window += " " + norm(nxt[4])
            if not CONC_TIME.search(window):
                continue
            label = label_of(caption)
            if label in seen:
                continue
            seen.add(label)

            cap_rect = pymupdf.Rect(blk[:4])
            rect, how = figure_rect(page, cap_rect)
            src = page

            if rect is None:
                num = re.search(r"\d+", label)
                for plate_no, plate_rect in plates.get(num.group() if num else "", []):
                    if plate_no != pno:
                        src, rect, how = doc[plate_no], plate_rect, f"plate p{plate_no + 1}"
                        break
            if rect is None and page_fallback:
                rect, how = page.rect, "whole page"
            if rect is None:
                rows.append({"PMID": pmid, "page": pno + 1, "figure": label,
                             "caption": window[:400], "image": "", "crop": how,
                             "status": "caption found, crop failed"})
                continue

            clip = (rect + (-PAD, -PAD, PAD, PAD)) & src.rect
            out = out_dir / f"{pmid}_{label}.png"
            src.get_pixmap(dpi=RENDER_DPI, clip=clip).save(out)
            rows.append({"PMID": pmid, "page": pno + 1, "figure": label,
                         "caption": window[:400], "image": str(out.relative_to(ROOT)),
                         "crop": how, "status": "ok"})
    doc.close()
    if not rows:
        rows = [{"PMID": pmid, "status": "no concentration-time caption matched"}]
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", default="mouse", choices=["mouse", "rat"])
    ap.add_argument("--pmid", nargs="*", help="limit to these PMIDs")
    ap.add_argument("--page-fallback", action="store_true",
                    help="render the whole page when the artwork cannot be isolated")
    args = ap.parse_args()

    sdir = ROOT / args.species
    prov = pd.read_csv(sdir / "records_provenance.csv", low_memory=False)
    dosed = prov[prov["in_cohort"].astype(bool) & prov["animal_dosed"].eq("dosed")]
    pmids = args.pmid or sorted({str(p) for p in dosed["PMID"]})

    out_dir = sdir / "pk_figures"
    out_dir.mkdir(exist_ok=True)

    rows: list[dict] = []
    for pmid in pmids:
        if not (PAPERS / f"{pmid}.pdf").exists():
            rows.append({"PMID": pmid, "status": "no publisher PDF"})
            continue
        rows.extend(harvest(pmid, out_dir, args.page_fallback))

    man = pd.DataFrame(rows).reindex(
        columns=["PMID", "page", "figure", "caption", "image", "crop", "status"])
    dest = sdir / "pk_figures.csv"
    man.to_csv(dest, index=False)

    ok = man[man["status"].eq("ok")]
    print(f"[write] {dest}")
    print(f"[figures] {len(ok)} images from {ok['PMID'].nunique()} of {len(pmids)} studies "
          f"-> {out_dir.relative_to(ROOT)}/")
    print(man["status"].value_counts().to_string())
    if not ok.empty:
        print(ok["crop"].value_counts().to_string())


if __name__ == "__main__":
    main()
