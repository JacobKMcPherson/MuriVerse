"""Assemble the per-agent PK parameter table for the dosed mouse cohort.

One row per (publication, agent, matrix, route, dose) arm, carrying the reported dose, route,
Cmax, tmax, absorption half-life and elimination half-life, plus the cropped
concentration-time figure for that publication.

Three sources are merged, and `source` / `needs_curation` record which supplied what:

  * PEPlife2 records - dose (`Concentration`), route (`Dose Route`) and the half-lives. The
    PK phase is only in the free text appended to the value ("2.12 +/- 0.3 (T1/2
    Elimination)"), so it is parsed out here and pivoted into separate ka / ke columns.
    PEPlife2 has no Cmax or tmax field at all.
  * `<species>/extractions/*.json` - curated Cmax/tmax/dose/route, for the papers read so far.
    Matched to PEPlife2 agents by normalised name, then by sequence; unmatched curated agents
    are emitted as their own rows so nothing curated is silently dropped.
  * `<species>/pk_figures.csv` - the cropped concentration-time figure from the publisher PDF.

  python build_pk_parameters.py --species mouse
"""
import argparse
import json
import re
import unicodedata
from pathlib import Path

import pandas as pd

from subset_in_vivo_half_life import CORRECTIONS

ROOT = Path(__file__).resolve().parent

# Checked in order; the first hit wins. Mojibake spellings are PEPlife2's own (alpha and beta
# round-tripped through a bad encoding), so they have to be matched literally.
PHASES = [
    ("absorption", r"ab[os]{2}rption"),
    ("distribution/alpha", r"distribution|alpha|α|Î±|\bt\s?1/2\s?,?\s?a\b|first phase|initial"),
    ("terminal/beta", r"terminal|beta|β|Î²|Î2|Β|\bt\s?1/2\s?,?\s?b\b|second phase|2nd phase|slow"),
    ("elimination", r"elimination"),
    ("other", r"activity|excretion|radioactivity|metabolic|clearance|fast"),
]

# Curators write "not reported as a plasma Cmax", not just "not reported", so the blank test
# has to be a prefix match rather than an exact one.
BLANK = re.compile(r"^\s*(?:n\.?\s?a\.?|na|-+|\?"
                   r"|not\s+(?:mentioned|reported|applicable|available|stated|given)\b.*)\s*$",
                   re.I | re.S)


def tidy(v: object) -> str:
    s = re.sub(r"\s+", " ", unicodedata.normalize("NFKC", str(v or ""))).strip()
    return "" if BLANK.match(s) else s


def norm_name(v: object) -> str:
    return re.sub(r"[^a-z0-9]", "", unicodedata.normalize("NFKC", str(v or "")).lower())


def agent_key(v: object) -> str:
    """PEPlife2 transcribed several tables twice, once with the size in the name."""
    return norm_name(re.sub(r"\([^)]*\)", " ", tidy(v)))


def dose_key(v: object) -> str:
    s = tidy(v).lower()
    s = re.sub(r"\s*(?:kg|ml|l|g)\s*[-\u2212]\s*1\b", r"/\g<0>", s)  # "nmol kg-1" -> "nmol /kg-1"
    s = re.sub(r"[-\u2212]\s*1\b", "", s)
    return re.sub(r"[^a-z0-9./]", "", s)


def phase_of(label: str) -> str:
    text = label.lower()
    for name, pat in PHASES:
        if re.search(pat, text, re.I):
            return name
    return "unspecified"


def load_curated(sdir: Path) -> pd.DataFrame:
    rows = []
    for f in sorted((sdir / "extractions").glob("*.json")):
        j = json.loads(f.read_text())
        for a in j.get("agents") or []:
            rows.append({
                "PMID": str(j.get("pmid", f.stem)),
                "curated_agent": tidy(a.get("name")),
                "curated_key": agent_key(a.get("name")),
                "curated_seq": norm_name(a.get("sequence")),
                "curated_role": tidy(a.get("role")),
                "curated_dose": tidy(a.get("dose")),
                "curated_route": tidy(a.get("route")),
                "cmax": tidy(a.get("cmax")),
                "tmax": tidy(a.get("tmax")),
                "auc": tidy(a.get("auc")),
                "curated_half_life": tidy(a.get("half_life")),
                "curated_half_life_phase": tidy(a.get("half_life_phase")),
                "curated_ka": tidy(a.get("half_life_absorption")),
                "curated_ke": tidy(a.get("half_life_elimination")),
                "curated_half_life_matrix": tidy(a.get("half_life_matrix")),
                "group_n": tidy(a.get("group_n")),
                "pk_quote": tidy(a.get("pk_quote")),
                "bioanalytical_method": tidy(j.get("bioanalytical_method")),
            })
    return pd.DataFrame(rows)


def build(species: str) -> pd.DataFrame:
    sdir = ROOT / species
    rec = pd.read_csv(sdir / "records_provenance.csv", low_memory=False)
    rec = rec[rec["in_cohort"].astype(bool) & rec["animal_dosed"].eq("dosed")].copy()
    rec["Half Life"] = rec["Half Life"].mask(rec["ID"].isin(CORRECTIONS),
                                             rec["ID"].map(CORRECTIONS))

    rec["half_life_text"] = rec["Half Life"].map(tidy)
    rec["pk_label"] = rec["half_life_text"].str.extract(r"\(([^)]*)\)$")[0].fillna("")
    rec["pk_phase"] = rec["pk_label"].map(phase_of)
    rec["value"] = rec["half_life_text"].str.replace(r"\s*\([^)]*\)$", "", regex=True)
    # pivot_table drops any row with a null in the index, which would silently lose the 10
    # sequence-less records (nanoparticle and protein arms) and 4 whole publications.
    rec["sequence"] = rec["sequence"].fillna("")
    for col, out in [("Name", "agent"), ("Test Sample", "matrix"),
                     ("Dose Route", "route"), ("Concentration", "dose")]:
        rec[out] = rec[col].map(tidy)

    # Group on normalised keys so PEPlife2's twice-transcribed tables ("Ex-PEGMw" at
    # "1000 nmol/kg" and "Ex-PEGMw(46.3KDa)" at "1000 nmol kg-1") form one arm, not two
    # half-populated ones; the readable spellings come back via `first` below.
    rec["agent_key"] = rec["agent"].map(agent_key)
    rec["seq_key"] = rec["sequence"].map(norm_name)
    rec["matrix_key"] = rec["matrix"].map(norm_name)
    rec["route_key"] = rec["route"].map(norm_name)
    rec["dose_key"] = rec["dose"].map(dose_key)
    key = ["PMID", "agent_key", "seq_key", "matrix_key", "route_key", "dose_key"]

    wide = (rec.pivot_table(index=key, columns="pk_phase", values="value",
                            aggfunc=lambda s: " ;; ".join(sorted(set(s))))
            .reset_index().rename_axis(columns=None))
    hours = (rec.pivot_table(index=key, columns="pk_phase", values="half_life_h",
                             aggfunc="median").reset_index().rename_axis(columns=None))

    for phase in ["absorption", "distribution/alpha", "terminal/beta", "elimination",
                  "unspecified", "other"]:
        if phase not in wide.columns:
            wide[phase] = pd.NA
        if phase not in hours.columns:
            hours[phase] = pd.NA

    def longest(s: pd.Series) -> str:
        return max((v for v in s if v), key=len, default="")

    display = (rec.groupby(key, dropna=False)
               .agg(agent=("agent", longest), sequence=("sequence", longest),
                    matrix=("matrix", longest), route=("route", longest),
                    dose=("dose", longest), n_records=("ID", "size"))
               .reset_index())

    out = display.copy()
    out["ka_t_half_reported"] = wide["absorption"]
    out["ka_t_half_h"] = hours["absorption"]
    # PEPlife2 labels the phase only when a study fit more than one; a lone unlabelled value is
    # the study's headline (elimination) half-life, so it backs ke rather than being discarded.
    out["ke_t_half_reported"] = (wide["elimination"].fillna(wide["terminal/beta"])
                                 .fillna(wide["unspecified"]))
    out["ke_t_half_h"] = (hours["elimination"].fillna(hours["terminal/beta"])
                          .fillna(hours["unspecified"]))
    out["ke_basis"] = ["elimination" if pd.notna(e) else
                       "terminal/beta" if pd.notna(t) else
                       "unlabelled (assumed elimination)" if pd.notna(u) else ""
                       for e, t, u in zip(wide["elimination"], wide["terminal/beta"],
                                          wide["unspecified"])]
    out["alpha_t_half_reported"] = wide["distribution/alpha"]
    out["alpha_t_half_h"] = hours["distribution/alpha"]
    out["other_t_half_reported"] = wide["other"]

    curated = load_curated(sdir)
    if not curated.empty:
        curated = curated[curated["PMID"].isin(set(rec["PMID"].astype(str)))]

    cur_cols = ["curated_agent", "curated_role", "curated_dose", "curated_route", "cmax",
                "tmax", "auc", "curated_half_life", "curated_half_life_phase",
                "curated_ka", "curated_ke", "curated_half_life_matrix", "group_n",
                "pk_quote", "bioanalytical_method"]
    if curated.empty:
        for c in cur_cols:
            out[c] = pd.NA
        return finalize(out, sdir, rec)

    by_name = curated.drop_duplicates(["PMID", "curated_key"])
    out = out.merge(by_name[["PMID", "curated_key"] + cur_cols],
                    left_on=["PMID", "agent_key"], right_on=["PMID", "curated_key"],
                    how="left").drop(columns=["curated_key"])

    # Second pass for agents PEPlife2 and the curator named differently.
    by_seq = (curated[curated["curated_seq"].ne("")]
              .drop_duplicates(["PMID", "curated_seq"]))
    gap = out["curated_agent"].isna()
    if gap.any():
        fill = (out.loc[gap, ["PMID", "seq_key"]]
                .merge(by_seq[["PMID", "curated_seq"] + cur_cols],
                       left_on=["PMID", "seq_key"], right_on=["PMID", "curated_seq"],
                       how="left"))
        for c in cur_cols:
            out.loc[gap, c] = fill[c].to_numpy()

    # Anything curated that still matched nothing becomes its own row - a curated Cmax/tmax
    # for an arm PEPlife2 never recorded is exactly the evidence this table exists to hold.
    matched = set(zip(out["PMID"], out["curated_agent"].fillna("")))
    extra = curated[[(p, a) not in matched for p, a in
                     zip(curated["PMID"], curated["curated_agent"])]]
    extra = extra[extra[["cmax", "tmax", "curated_half_life"]].ne("").any(axis=1)]
    if not extra.empty:
        add = extra[["PMID"] + cur_cols].copy()
        add["agent"] = extra["curated_agent"]
        add["sequence"] = ""
        add["matrix"] = extra["curated_half_life_matrix"]
        add["route"] = ""
        add["dose"] = ""
        add["ke_basis"] = ""
        out = pd.concat([out, add], ignore_index=True)

    # Curation is the only source of a phase-resolved half-life for papers PEPlife2 recorded
    # as a single unlabelled value, so let it fill ka/ke rather than sitting in its own column.
    def backfill(target: str, source: str, hours: str | None = None,
                 basis: str | None = None) -> None:
        def present(col: str) -> pd.Series:
            return out[col].notna() & out[col].astype("string").fillna("").str.strip().ne("")

        use = ~present(target) & present(source)
        out.loc[use, target] = out.loc[use, source]
        if hours:
            out.loc[use, hours] = pd.NA
        if basis:
            out.loc[use, "ke_basis"] = basis

    backfill("ka_t_half_reported", "curated_ka", hours="ka_t_half_h")
    backfill("ke_t_half_reported", "curated_ke", hours="ke_t_half_h", basis="curated extraction")
    backfill("dose", "curated_dose")
    backfill("route", "curated_route")

    return finalize(out, sdir, rec)


def finalize(out: pd.DataFrame, sdir: Path, rec: pd.DataFrame) -> pd.DataFrame:
    figs = pd.read_csv(sdir / "pk_figures.csv", dtype={"PMID": str})
    figs = figs[figs["status"].eq("ok")]
    figs = (figs.groupby("PMID")
            .agg(figure_image=("image", lambda s: " ;; ".join(s)),
                 figure_caption=("caption", lambda s: " ;; ".join(s)))
            .reset_index())
    out["PMID"] = out["PMID"].astype(str)
    out = out.merge(figs, on="PMID", how="left")

    meta = rec.drop_duplicates("PMID")[["PMID", "study_title", "study_journal", "pub_year",
                                        "study_doi", "pmcid", "evidence", "extraction_status"]]
    meta["PMID"] = meta["PMID"].astype(str)
    out = out.merge(meta, on="PMID", how="left")

    missing = []
    for _, r in out.iterrows():
        gaps = [n for n, v in [("dose", r["dose"]), ("route", r["route"]),
                               ("cmax", r.get("cmax")), ("tmax", r.get("tmax")),
                               ("ka", r["ka_t_half_reported"]), ("ke", r["ke_t_half_reported"]),
                               ("figure", r.get("figure_image"))]
                if pd.isna(v) or not str(v).strip()]
        missing.append("|".join(gaps))
    out["needs_curation"] = missing

    order = ["PMID", "study_title", "study_journal", "pub_year", "study_doi", "pmcid",
             "agent", "curated_agent", "sequence", "matrix",
             "dose", "curated_dose", "route", "curated_route",
             "cmax", "tmax", "auc",
             "ka_t_half_reported", "ka_t_half_h",
             "ke_t_half_reported", "ke_t_half_h", "ke_basis",
             "curated_half_life", "curated_half_life_phase",
             "alpha_t_half_reported", "alpha_t_half_h", "other_t_half_reported",
             "group_n", "n_records", "bioanalytical_method", "pk_quote",
             "figure_image", "figure_caption",
             "evidence", "extraction_status", "needs_curation"]
    return out[[c for c in order if c in out.columns]].sort_values(
        ["PMID", "agent"], key=lambda s: pd.to_numeric(s, errors="coerce")
        if s.name == "PMID" else s)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", default="mouse", choices=["mouse", "rat"])
    args = ap.parse_args()

    out = build(args.species)
    dest = ROOT / args.species / "pk_parameters.csv"
    out.to_csv(dest, index=False)

    def present(col: str) -> pd.Series:
        s = out[col]
        return s.notna() & s.astype("string").fillna("").str.strip().ne("")

    print(f"[write] {dest} - {len(out)} arms across {out['PMID'].nunique()} publications")
    for col, label in [("dose", "dose"), ("route", "route"), ("cmax", "Cmax"),
                       ("tmax", "tmax"), ("ka_t_half_reported", "ka (absorption t1/2)"),
                       ("ke_t_half_reported", "ke (elimination t1/2)"),
                       ("figure_image", "conc-time figure")]:
        have = present(col)
        print(f"[have] {label:<24} {int(have.sum()):>4} / {len(out)} arms "
              f"({out.loc[have, 'PMID'].nunique():>3} of {out['PMID'].nunique()} pubs)")
    print("[ke basis]")
    print(out["ke_basis"].replace("", "none").value_counts().to_string())


if __name__ == "__main__":
    main()
