"""Audit MuriVerse peptide half-life datasets for provenance, quality, and source-origin differences.

The script is intentionally schema-tolerant. It reads CSV/XLSX files under MuriVerse/data,
profiles source tables, and writes reproducible CSV/Markdown summaries without modifying source data.

Usage:
    python MuriVerse/analyze_data_provenance.py
    python MuriVerse/analyze_data_provenance.py --data-root MuriVerse/data --out MuriVerse/reports/data_provenance
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

SOURCE_DIRS = {
    "peplife2": "peplife2",
    "peptherdia2": "peptherdia2",
    "thpdb2": "thpdb2",
    "pubmed": "pubmed",
}

ALIASES = {
    "sequence": ["sequence", "seq", "amino acid sequence", "amino_acid_sequence"],
    "pmid": ["pmid", "pubmed id", "pubmed_id"],
    "half_life": ["half life", "half-life", "half_life", "halflife", "half-life(in sec)", "half_life_seconds"],
    "matrix": ["test sample", "matrix", "sample", "biological matrix"],
    "assay": ["assay", "bioanalytical method", "bioanalytical_method"],
    "route": ["dose route", "route", "route of administration"],
    "species": ["species", "organism"],
    "id": ["id", "record id", "record_id", "peplife id"],
    "reference": ["reference", "citation", "doi"],
}

def norm(s: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(s).strip().lower()).strip()

def pick_col(df: pd.DataFrame, key: str) -> str | None:
    cols = {norm(c): c for c in df.columns}
    for alias in ALIASES[key]:
        if norm(alias) in cols:
            return cols[norm(alias)]
    return None

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def read_table(path: Path) -> pd.DataFrame | None:
    try:
        if path.suffix.lower() == ".csv":
            return pd.read_csv(path, low_memory=False)
        if path.suffix.lower() in {".xlsx", ".xls"}:
            return pd.read_excel(path)
    except Exception as exc:
        print(f"[warn] {path}: {exc}")
    return None

def source_files(root: Path):
    for source, sub in SOURCE_DIRS.items():
        d = root / sub
        if not d.exists():
            continue
        for p in sorted(d.rglob("*")):
            if p.is_file() and p.suffix.lower() in {".csv", ".xlsx", ".xls"}:
                yield source, p

def profile(source: str, path: Path, df: pd.DataFrame) -> dict:
    seq = pick_col(df, "sequence")
    pmid = pick_col(df, "pmid")
    hl = pick_col(df, "half_life")
    matrix = pick_col(df, "matrix")
    rec = {
        "source": source,
        "file": str(path),
        "sha256": sha256(path),
        "rows": len(df),
        "columns": len(df.columns),
        "sequence_column": seq,
        "pmid_column": pmid,
        "half_life_column": hl,
        "matrix_column": matrix,
        "missing_cells_fraction": float(df.isna().sum().sum() / max(df.size, 1)),
    }
    if seq:
        s = df[seq].astype("string").str.strip()
        rec["unique_sequences"] = int(s.dropna().nunique())
        rec["duplicate_sequence_rows"] = int(s.duplicated(keep=False).sum())
        rec["sequence_missing_fraction"] = float(s.isna().mean())
        rec["median_sequence_length"] = float(s.dropna().str.len().median()) if s.notna().any() else None
    if pmid:
        p = df[pmid].astype("string").str.extract(r"(\d+)")[0]
        rec["unique_pmids"] = int(p.dropna().nunique())
        rec["pmid_missing_fraction"] = float(p.isna().mean())
    if hl:
        x = pd.to_numeric(df[hl], errors="coerce")
        rec["half_life_numeric_fraction"] = float(x.notna().mean())
        if x.notna().any():
            rec["half_life_median_numeric"] = float(x.median())
            rec["half_life_iqr_numeric"] = float(x.quantile(.75) - x.quantile(.25))
    return rec

def collect_overlap(frames: list[tuple[str, Path, pd.DataFrame]]) -> pd.DataFrame:
    rows = []
    sets = {}
    for source, path, df in frames:
        seq = pick_col(df, "sequence")
        if seq:
            vals = set(df[seq].dropna().astype(str).str.strip())
            vals.discard("")
            sets[(source, str(path))] = vals
    keys = list(sets)
    for i, a in enumerate(keys):
        for b in keys[i+1:]:
            inter = sets[a] & sets[b]
            union = sets[a] | sets[b]
            rows.append({
                "source_a": a[0], "file_a": a[1], "source_b": b[0], "file_b": b[1],
                "sequence_overlap_n": len(inter),
                "jaccard": len(inter) / len(union) if union else 0.0,
            })
    return pd.DataFrame(rows)

def quality_flags(source: str, path: Path, df: pd.DataFrame) -> list[dict]:
    flags = []
    seq = pick_col(df, "sequence")
    pmid = pick_col(df, "pmid")
    hl = pick_col(df, "half_life")
    required = {"sequence": seq, "pmid": pmid, "half_life": hl}
    for field, col in required.items():
        if col is None:
            flags.append({"source": source, "file": str(path), "severity": "high",
                          "issue": f"missing recognizable {field} column"})
    if seq:
        s = df[seq].astype("string").str.strip()
        if s.duplicated().any():
            flags.append({"source": source, "file": str(path), "severity": "medium",
                          "issue": f"{int(s.duplicated(keep=False).sum())} rows share a sequence; inspect experiment-level duplication/leakage"})
        bad = s.dropna().str.contains(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO\-*()0-9\[\]{}+.,/ ]", regex=True)
        if bad.any():
            flags.append({"source": source, "file": str(path), "severity": "medium",
                          "issue": f"{int(bad.sum())} sequences contain nonstandard characters/modification notation"})
    if hl:
        x = pd.to_numeric(df[hl], errors="coerce")
        if x.notna().mean() < 0.8:
            flags.append({"source": source, "file": str(path), "severity": "medium",
                          "issue": "half-life field is substantially nonnumeric/censored/free-text; preserve reported text and explicit transform rules"})
    return flags

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="MuriVerse/data")
    ap.add_argument("--out", default="MuriVerse/reports/data_provenance")
    args = ap.parse_args()
    root, out = Path(args.data_root), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    frames, profiles, flags = [], [], []
    for source, path in source_files(root):
        df = read_table(path)
        if df is None:
            continue
        frames.append((source, path, df))
        profiles.append(profile(source, path, df))
        flags.extend(quality_flags(source, path, df))

    prof = pd.DataFrame(profiles)
    qf = pd.DataFrame(flags)
    overlap = collect_overlap(frames)
    prof.to_csv(out / "source_file_profile.csv", index=False)
    qf.to_csv(out / "quality_flags.csv", index=False)
    overlap.to_csv(out / "sequence_overlap.csv", index=False)

    manifest = {
        "data_root": str(root),
        "files_profiled": len(frames),
        "sources": sorted({x[0] for x in frames}),
        "outputs": ["source_file_profile.csv", "quality_flags.csv", "sequence_overlap.csv"],
        "note": "Numeric half-life summaries are only directly comparable when units/transforms are harmonized."
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))

    lines = ["# MuriVerse data provenance audit", "",
             f"Files profiled: **{len(frames)}**", ""]
    if not prof.empty:
        by = prof.groupby("source", dropna=False).agg(files=("file","count"), rows=("rows","sum")).reset_index()
        lines += ["## Source inventory", "", by.to_markdown(index=False), ""]
    lines += ["## Interpretation", "",
              "- Source-level differences can reflect database scope, curation, assay context, publication selection, and biology.",
              "- Compare half-life distributions only after harmonizing units and separating in-vivo PK from ex-vivo stability.",
              "- Group model splits by sequence and preferably publication/experiment to reduce leakage.",
              "- Treat missing provenance as a data-quality defect rather than silently imputing source context.", ""]
    (out / "README.md").write_text("\n".join(lines))
    print(f"[done] wrote audit outputs to {out}")

if __name__ == "__main__":
    main()
