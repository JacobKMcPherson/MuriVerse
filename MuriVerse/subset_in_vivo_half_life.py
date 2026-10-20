"""Subset <species>/records_provenance.csv to live-dosed-animal half-life values.

The analysis set for in vivo pharmacokinetics is narrower than either of the columns it is
built from:

  * `design == "in vivo PK"` is not by itself the question - `design_of` buckets every tissue
    matrix as ex vivo before it reads the Vivo/Vitro flag, so a dosed animal whose organ was
    harvested is filed under `ex vivo (tissue homogenate)`. `animal_dosed` is the column that
    answers "did a live animal receive the peptide?" (see build_study_dossier.dosed_of), and
    both conditions are required here: dosed *and* a plasma/serum/blood/whole-animal readout.
  * a dosed record does not necessarily carry a half-life; PEPlife2 leaves the value blank on
    some rows, so `half_life_h` must be non-null as well.

Two further cleanups, applied in this order:

  * absorption-phase half-lives are dropped. They describe entry from a subcutaneous depot
    into circulation, not persistence in it, so they are a different parameter from the rest
    of the set - and the two papers that report them (28989813, 36630826) also report the
    elimination half-life of the same animals, which is retained. NB PEPlife2 labels the
    36630826 rows `Units Half Life == "Hours-1"`, but they are not rate constants: the source
    paper's Table 2 header reads "t1/2 absorption (h-1)" for what its Methods and body text
    give in hours ("t1/2 absorption=9.05 +/- 0.9h", t1/2 = 0.693/k already applied), so the
    typo is the paper's and the values need no conversion.
  * PEPlife2 repeats an identical measurement as several entries; the repeats are collapsed
    to the first copy so an unweighted median does not count those studies many times over.
    The duplicate key ignores whitespace, because PEPlife2 transcribed some tables twice with
    inconsistent spacing ("23.14 +/- 3.2(T1/2 Elimination)" vs "23.14 +/- 3.2 (T1/2
    Elimination)") and a literal key treats those as distinct measurements.

The surviving records are then classified along four axes, added as columns:

  pk_phase        the set still mixes PK phases. PEPlife2 carries the phase only as free text
                  inside the value ("2.12 +/- 0.3 (T1/2 Elimination)"), parsed here into
                  `pk_label`. Two-compartment fits contribute BOTH phases as separate
                  records, so an unfiltered median mixes a distribution half-life of minutes
                  with a terminal half-life of days from the same animal. The alpha/beta
                  marker is read BEFORE the word "elimination", because several papers label
                  both phases "elimination" and distinguish them only as T1/2,a vs T1/2,b -
                  the subscript is the reliable signal, not the noun. Half-lives of a
                  non-concentration readout (activity, excretion, radioactivity) are split
                  off: they answer a different question from a plasma-concentration t1/2.

  route_class     normalises `Dose Route`. Retro-orbital and tail-vein are intravenous.
                  "Injection" with no site named stays its own bucket rather than being
                  guessed at, and "N.A." is PEPlife2's blank, not a reported absence.

  topology        linear / cyclic / other, from the `Linear/Cyclic` prefix - the bracketed
                  detail ("Cyclic (C6-C11 disulfide bond in A chain)") is per-peptide text.

  hl_engineering  the half-life-extension strategy, pattern-matched across Name / Origin /
                  Chem-MOD / N-ter MOD / C-ter MOD, because it is the dominant driver of the
                  spread here (one sequence spans three orders of magnitude between its
                  native and its depot form). Buckets are ordered by how strongly they
                  determine clearance and the first match wins.

  python subset_in_vivo_half_life.py --species mouse
"""
import argparse
import re
import unicodedata
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent

# Same key build_study_dossier.py uses, re-applied post-filter so `is_repeat` keeps the first
# copy that survives the filters rather than one that was already dropped.
DUP_KEY = ["PMID", "sequence", "hl_sample", "hl_value", "hl_units"]

ABSORPTION = r"ab[os]{2}rption"  # PEPlife2 also spells it "Abosrption"

# PEPlife2 filed Table 2's absorption value for Ex-POEGMAopt in immunized mice under the
# elimination label; the paper's elimination value for that arm is 105.7 h (ID 4270/4281).
CORRECTIONS = {4269: "12.1 ± 1.0 (T1/2 Absorption)"}

# PEPlife2's export is mojibake in places: beta arrives as "Î2"/"Î²"/"Β", alpha as "Î±".
BETA = "(?:\u03b2|\u0392|\u00ce\u00b2|\u00ce2)"
ALPHA = "(?:\u03b1|\u0391|\u00ce\u00b1)"

PK_PHASE = [
    ("non-concentration", r"activity|excretion|radioactivity|metabolic"),
    # Absorption is a rate of entry from the injection site, not an elimination half-life,
    # and PEPlife2 spells it three ways; it used to fall through to "unlabelled", which the
    # training set treats as concentration decay.
    ("absorption", ABSORPTION + r"|\babsorb|\brelease\b|\blag\b"),
    ("distribution / alpha", rf"{ALPHA}|\balpha\b|t\s*1/2\s*,?\s*a\b|\binitial\b|"
                             r"first phase|1st phase|\bfast\b|distribution"),
    ("terminal / beta", rf"{BETA}|t\s*1/2\s*,?\s*b\b|terminal|elmination|elimination|"
                        r"second phase|2nd phase|\bslow\b|clearance"),
]

ROUTE = [
    ("iv or oral (not resolved)", r"intravenous or oral"),
    ("iv", r"\biv\b|intraven|tail vein|retro-?orbital|\bro\b"),
    ("sc", r"\bsc\b|subcutan"),
    ("ip", r"\bip\b|intraperit"),
    ("im", r"\bim\b|intramusc"),
    ("oral", r"\bpo\b|\bop\b|oral|gavage|intragastr"),
    ("intranasal", r"intranas"),
    ("injection (site not stated)", r"inject"),
]

TOPOLOGY = [("cyclic", r"^cyclic|^bicyclic"), ("linear", r"^linear")]

# Ordered by dominance over clearance: a depot or fusion overrides a bare sequence annotation.
ENGINEERING = [
    ("depot / particle formulation", r"microspher|nanoparticle|nanoghost|liposom|hydrogel|"
                                     r"\bdepot\b|\bms-\b|micelle|prepolymer"),
    ("fc fusion", r"\bfc\b|igg1|antibody|scfv|repebody"),
    ("albumin-based", r"\bhsa\b|\balbumin|\babd\b|albumin-bind"),
    ("polymer conjugate", r"\belp\b|zipp|polysarcosine|chondroitin|heparosan|poegma|"
                          r"pas\d|pasylat|xten|polyurethane|poly\(|dextran"),
    ("pegylated", r"\bpeg\b|pegyl|peg\d|peg-|-peg|\bpeg\d*[a-z]*\b"),
    ("lipidated", r"palmitoyl|myrist|stearo|fatty acid|lipida|\bc16\b|\bc18\b|octanoyl|acyl|"
                  r"cholester|\bchol\b"),
    ("radio/fluorescent label only", r"\d{2,3}i\b|i-?125|i125|3h |99mtc|177lu|64cu|68ga|"
                                     r"fitc|cy5|cy7|alexa|atto|dota|nota|fluoresc|radiolab"),
    ("backbone modification only", r"staple|lactam|d-amino|n-me|amidation|acetylation|"
                                   r"cyclic|disulfide"),
]

AXES = ["pk_phase", "route_class", "topology", "hl_engineering"]

COLUMNS = [
    "ID", "PMID", "YEAR",
    "Sequence", "sequence", "seq_length", "Name", "Length", "N-ter MOD", "C-ter MOD",
    "Linear/Cyclic", "topology", "Chirality", "Chem-MOD", "hl_engineering", "Origin", "Nature",
    "Half Life", "half_life_norm_h", "Units Half Life", "half-life(in Sec)", "half_life_h",
    "pk_label", "pk_phase",
    "Dose Route", "route_class",
    "Test Sample", "Vivo/Vitro", "setting", "matrix", "design", "animal_dosed",
    "Assay", "Protease", "Concentration", "Incubation Time", "Activity",
    "api_query", "is_repeat", "identical_copies",
    "study_title", "study_journal", "pub_year", "study_doi", "pmcid",
    "study_first_author", "study_publication_types",
    "evidence", "extraction_status", "curated_half_lives", "curated_doses",
]


def squash(s: pd.Series) -> pd.Series:
    # Python's re, not .str.replace: pandas 3's Arrow regex engine treats `\s` as ASCII-only,
    # so the NBSPs PEPlife2 left in its duplicated table rows survive and defeat the dup key.
    return s.fillna("").map(
        lambda v: re.sub(r"\s+", "", unicodedata.normalize("NFKC", str(v))).lower())


def tidy(s: pd.Series) -> pd.Series:
    return s.fillna("").map(
        lambda v: re.sub(r"\s+", " ", unicodedata.normalize("NFKC", str(v))).strip())


def bucket(text: pd.Series, rules: list[tuple[str, str]], default: str) -> pd.Series:
    # re.search, not Series.str.contains: pandas 3 dispatches to Arrow's RE2, which rejects
    # \u escapes outright and matches \s / \b as ASCII-only against this mojibake text.
    compiled = [(label, re.compile(pattern)) for label, pattern in rules]

    def first_match(value: str) -> str:
        for label, rx in compiled:
            if rx.search(value):
                return label
        return default

    return text.str.lower().map(first_match)


def add_axes(out: pd.DataFrame) -> pd.DataFrame:
    # `Half Life` is free text carrying its own units, an error term, an approximation marker
    # and the phase label ("97.3 +/- 3.2(T1/2 Elimination)", "~4", "53.6 +/- 5.8" minutes).
    # Same value as `half_life_h`, in hours, cut to the 3 significant figures the sources
    # actually report - the raw column carries float-division noise (0.5833333333333334).
    out["half_life_norm_h"] = out["half_life_h"].map(
        lambda v: v if pd.isna(v) else float(f"{v:.3g}"))
    out["pk_phase"] = bucket(out["pk_label"], PK_PHASE, "unlabelled")
    route = tidy(out["Dose Route"]).replace({"N.A.": "", "NA": "", "Not Reported": ""})
    out["route_class"] = bucket(route, ROUTE, "not reported")
    out["topology"] = bucket(tidy(out["Linear/Cyclic"]), TOPOLOGY, "not reported")
    chem = tidy(out["Name"])
    for col in ["Origin", "Chem-MOD", "N-ter MOD", "C-ter MOD"]:
        chem = chem + " | " + tidy(out[col])
    out["hl_engineering"] = bucket(chem, ENGINEERING, "unmodified / none stated")
    # `Length` is empty for every record in this subset, so measure the sequence directly.
    out["seq_length"] = out["sequence"].fillna("").astype(str).str.len().replace(0, pd.NA)
    return out


def summarise(out: pd.DataFrame, axis: str) -> pd.DataFrame:
    g = out.groupby(axis, dropna=False)
    return pd.DataFrame({
        "records": g.size(),
        "studies": g["PMID"].nunique(),
        "sequences": g["sequence"].nunique(),
        "median_h": g["half_life_h"].median().round(3),
        "min_h": g["half_life_h"].min().round(3),
        "max_h": g["half_life_h"].max().round(3),
    }).sort_values("records", ascending=False)


def build(species: str) -> tuple[pd.DataFrame, dict[str, int]]:
    src = ROOT / species / "records_provenance.csv"
    rec = pd.read_csv(src, low_memory=False)
    rec["sample"] = rec["Test Sample"].fillna("").astype(str)
    rec["Half Life"] = rec["Half Life"].mask(rec["ID"].isin(CORRECTIONS),
                                             rec["ID"].map(CORRECTIONS))
    # PEPlife2 tacks the PK parameter onto the value as free text: "2.12 ± 0.3 (T1/2 Elimination)".
    rec["pk_label"] = tidy(rec["Half Life"]).str.extract(r"\(([^)]*)\)$")[0].fillna("")
    rec["hl_sample"] = squash(rec["sample"])
    rec["hl_value"] = squash(rec["Half Life"])
    rec["hl_units"] = squash(rec["Units Half Life"])

    steps: dict[str, int] = {}
    out = rec[rec["in_cohort"].astype(bool)]
    steps["cohort records"] = len(out)
    out = out[out["animal_dosed"].eq("dosed")]
    steps["live animal dosed"] = len(out)
    out = out[out["half_life_h"].notna()]
    steps["half-life reported"] = len(out)
    out = out[out["design"].eq("in vivo PK")]
    steps["design = in vivo PK"] = len(out)
    out = out[~out["hl_value"].str.contains(ABSORPTION)]
    steps["dropped absorption phase"] = len(out)

    out = out.copy()
    out["identical_copies"] = out.groupby(DUP_KEY, dropna=False)["ID"].transform("size")
    out["is_repeat"] = out.duplicated(DUP_KEY, keep="first")
    out = out[~out["is_repeat"]]
    steps["deduplicated"] = len(out)

    out = add_axes(out)
    out = out[[c for c in COLUMNS if c in out.columns]]
    return (out.sort_values(["PMID", "ID"], key=lambda s: pd.to_numeric(s, errors="coerce")),
            steps)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", default="mouse", choices=["mouse", "rat"])
    args = ap.parse_args()

    out, steps = build(args.species)
    dest = ROOT / args.species / "in_vivo_half_life.csv"
    out.to_csv(dest, index=False)

    prev = None
    for label, n in steps.items():
        delta = "" if prev is None else f"  ({n - prev:+d})"
        print(f"[filter] {label:<26} {n:>4}{delta}")
        prev = n
    print(f"[write] {dest} - {len(out)} records, "
          f"{out['PMID'].nunique()} studies, {out['sequence'].nunique()} sequences")
    print(f"[half-life h] median {out['half_life_h'].median():.3g}, "
          f"range {out['half_life_h'].min():.3g}-{out['half_life_h'].max():.3g}")

    for axis in AXES:
        print(f"\n[{axis}]")
        print(summarise(out, axis).to_string())

    core = out[out["pk_phase"].isin(["terminal / beta", "unlabelled"])]
    print(f"\n[analysis set] single-phase + terminal (drops distribution-phase and "
          f"non-concentration half-lives) = {len(core)} records, "
          f"{core['PMID'].nunique()} studies, {core['sequence'].nunique()} sequences, "
          f"median {core['half_life_h'].median():.3g} h")


if __name__ == "__main__":
    main()
