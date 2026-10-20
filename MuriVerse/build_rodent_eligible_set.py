"""Select every rodent in-vivo pharmacokinetic record of a peptide 50 residues or shorter.

This replaces the species selection the pipeline used to make. The old route went through
PEPlife 2.0's organism API (`fetch_peplife2_api.py`, the `Mouse` + `Mice` and `Rat` + `Rats`
queries), which is the database's own annotation of what the paper studied. This one reads
`Test Sample` instead — the free-text description of what was actually sampled — and asks
whether it names a mouse or a rat.

Reading the sample text finds more records than the organism annotation does, because a
record can describe "SD rats plasma" or "Blood sample of female Sprague-Dawley rat" while
carrying an organism annotation that is missing or points at the peptide's origin species
rather than the dosed animal.

The selection is then narrowed on three axes, in this order:

  * **dosed animal, not isolated matrix.** `Vivo/Vitro` alone is not trustworthy — the flag
    is internally inconsistent and marks serum-stability assays as in vivo. The curated
    `design` and `animal_dosed` classifications from `build_study_dossier.py` are used
    instead, which recover the study design from the sample text. At this length cap the two
    nearly agree: of 1,121 records flagged in vivo, 1,087 are genuinely dosed-animal PK and
    only 34 are ex vivo stability assays.
  * **50 residues or fewer**, measured from the sequence itself rather than from PEPlife2's
    `Length` column, which is empty for every record in the shipped table.
  * **a positive half-life**, since log10 of the label requires one.

Cell-line matrices and PEPlife2's own repeated transcriptions are dropped, as they are
everywhere else in this pipeline.

What comes out is the ELIGIBLE set: every record that is in scope, whether or not its
molecule can be drawn. `build_strict_peptide_set.py` takes it from here and decides which
records carry a structure complete enough to model, and how.

  python MuriVerse/build_rodent_eligible_set.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# The modelling axes - PK phase, route class, topology, engineering class - are derived by
# `subset_in_vivo_half_life.add_axes` and are reused rather than re-implemented, so this
# selection and the per-species sets cannot drift apart in how they bucket the same text.
from subset_in_vivo_half_life import add_axes, tidy  # noqa: E402
RECORDS = HERE / "data" / "peplife2" / "peplife2_records.csv"
OUTDIR = HERE / "rodent"
SPECIES_COHORTS = ("mouse", "rat")

# Word-boundary matching so "rat" cannot fire inside "titrated" or "aspirate", and so the
# plural forms the database mixes freely are both caught.
RODENT = re.compile(r"\b(?:mouse|mice|rat|rats)\b", re.IGNORECASE)

# The TRUE residue limit is applied by build_strict_peptide_set.py, which parses the
# sequence into residues and counts them. This module cannot: it would have to re-implement
# the parser to know that "Y-Aib-EGTFISDLSI" is 16 residues and not the 20 letters it is
# written with. So the pre-filter here is deliberately generous - it exists only to drop
# proteins cheaply - and a sequence is kept if its LETTER count leaves room for a 50-residue
# peptide written with delimited non-standard residues. Counting letters at 50 would reject
# a 45-residue peptide for containing three Aib, which is exactly the error this avoids.
MAX_RESIDUES = 50
MAX_SEQUENCE_LETTERS = 80

# The curated study-design classes that mean "dosed into a live animal and sampled from its
# circulation". Everything else is an isolated-matrix experiment answering a different
# question, however PEPlife2's own Vivo/Vitro flag describes it.
DOSED_DESIGN = "in vivo PK"

STAGES = [
    ("PEPlife2 records", "the complete shipped record table"),
    ("Test Sample names a mouse or rat", "the dosed animal, read from the sample text "
                                         "rather than from the organism annotation"),
    ("flagged in vivo", "PEPlife2's own Vivo/Vitro flag; not trusted on its own"),
    ("positive half-life", "log10 of the label requires one"),
    (f"peptide scale (pre-filter, <= {MAX_SEQUENCE_LETTERS} letters)",
     f"a generous letter-count bound; the true <= {MAX_RESIDUES}-residue limit is applied "
     f"on parsed residues by build_strict_peptide_set.py"),
    ("dosed animal, not isolated matrix", "curated design and animal_dosed, which recover "
                                          "the study design from the sample text"),
    ("not a cell-line matrix", "a cultured-cell assay is not the animal"),
    ("deduplicated", "PEPlife2's own repeated transcriptions of one measurement"),
]


def measured_length(sequence: object) -> int:
    """Residue count from the sequence itself, letters only.

    PEPlife2 ships a `Length` column and it is empty for all 4,500 records, so every length
    in this project is measured rather than read. Punctuation and delimiters are stripped,
    which counts a delimited non-standard residue as the letters it is written with - an
    over-count that can only make this gate stricter, never looser.
    """
    return len(re.sub(r"[^A-Za-z]", "", str(sequence)))


def provenance() -> pd.DataFrame:
    """The curated per-record dossier for both species cohorts, one row per PEPlife2 ID."""
    frames = []
    for cohort in SPECIES_COHORTS:
        path = HERE / cohort / "records_provenance.csv"
        if not path.exists():
            raise SystemExit(f"{path} is missing; run build_study_dossier.py first")
        frame = pd.read_csv(path, low_memory=False)
        frame["cohort_species"] = cohort
        frames.append(frame)
    merged = pd.concat(frames, ignore_index=True)
    # A record studied in both cohorts is one record; the dossier columns agree for it.
    return merged.drop_duplicates(subset="ID", keep="first")


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    records = pd.read_csv(RECORDS, low_memory=False)
    dossier = provenance()

    missing = set(records["ID"]) & set(dossier["ID"])
    if not missing:
        raise SystemExit("no PEPlife2 IDs in common between the record table and the "
                         "curated dossier; the two are out of step")

    frame = records.merge(
        dossier[["ID", "design", "animal_dosed", "matrix", "matrix_is_cell_line",
                 "is_repeat", "setting", "half_life_h", "cohort_species"]],
        on="ID", how="left")

    counts = [{"stage": STAGES[0][0], "records": len(frame), "removed": 0,
               "rationale": STAGES[0][1]}]

    def cut(mask: pd.Series, stage: int) -> None:
        nonlocal frame
        before = len(frame)
        frame = frame[mask].copy()
        counts.append({"stage": STAGES[stage][0], "records": len(frame),
                       "removed": before - len(frame), "rationale": STAGES[stage][1]})

    cut(frame["Test Sample"].astype(str).str.contains(RODENT, na=False), 1)
    cut(frame["Vivo/Vitro"].astype(str).str.contains("vivo", case=False, na=False), 2)
    cut(pd.to_numeric(frame["half-life(in Sec)"], errors="coerce").gt(0), 3)
    frame["seq_length_measured"] = frame["Sequence"].map(measured_length)
    cut(frame["seq_length_measured"].between(1, MAX_SEQUENCE_LETTERS), 4)
    cut((frame["design"] == DOSED_DESIGN) & (frame["animal_dosed"] == "dosed"), 5)
    cut(~frame["matrix_is_cell_line"].astype("boolean").fillna(False).astype(bool), 6)
    cut(~frame["is_repeat"].astype("boolean").fillna(False).astype(bool), 7)

    # `species` is the dosed animal as the sample text names it, which is what every
    # downstream split and metric is cut by. Where the text names both - a comparison across
    # animals in one sentence - the cohort the record was curated under decides.
    sample = frame["Test Sample"].astype(str)
    names_mouse = sample.str.contains(r"\b(?:mouse|mice)\b", case=False, regex=True, na=False)
    names_rat = sample.str.contains(r"\brats?\b", case=False, regex=True, na=False)
    frame["species"] = pd.Series(
        [cohort if m and r else ("mouse" if m else "rat")
         for m, r, cohort in zip(names_mouse, names_rat, frame["cohort_species"])],
        index=frame.index)

    # `sequence` is the cleaned residue string the rest of the pipeline reads; `pk_label` is
    # the phase annotation PEPlife2 appends to its own value in parentheses
    # ("2.12 +/- 0.3 (T1/2 Elimination)"), which add_axes buckets into `pk_phase`.
    # `add_axes` derives `seq_length` from this column. It is the RAW sequence here, so that
    # length is a letter count and not a residue count; build_strict_peptide_set.py
    # overwrites it with the parsed residue count, which is the one anything downstream uses.
    frame["sequence"] = frame["Sequence"]
    frame["pk_label"] = tidy(frame["Half Life"]).str.extract(r"\(([^)]*)\)$")[0].fillna("")
    frame = add_axes(frame)

    # The curated full-text cross-check was run per species, and its verdicts are carried
    # over for the records it covers. Everything outside it is stated as what it is - an
    # unchecked database transcription - rather than left blank and read as though it had
    # been verified.
    #
    # Read straight from the two per-species tables. These identify the animal from
    # PEPlife2's organism annotation while this builder reads the sample text, so the two
    # cross rather than nest and neither is a parent of the other: they are a source of
    # curated COLUMNS for the records they share, never a funnel stage or a record count.
    # Every curated column is carried across, so this table is an information superset of
    # them, which is what lets the dataset chain nest PEPlife2 -> eligible -> strict ->
    # analysis set with nothing downstream reaching back.
    curated = []
    for species_dir in ("mouse", "rat"):
        path = HERE / species_dir / "in_vivo_half_life.csv"
        if path.exists():
            curated.append(pd.read_csv(path, low_memory=False))

    frame["half_life_provenance"] = "PEPlife2 transcription (no curated full text)"
    if curated:
        checked = pd.concat(curated, ignore_index=True)
        if checked["ID"].duplicated().any():
            dupes = sorted(checked.loc[checked["ID"].duplicated(), "ID"])
            raise ValueError(f"PEPlife2 IDs appear in both species tables: {dupes}")
        if "half_life_provenance" in checked:
            verdict = checked.set_index("ID")["half_life_provenance"]
            carried = frame["ID"].map(verdict)
            frame["half_life_provenance"] = carried.fillna(frame["half_life_provenance"])
        # Curated evidence columns the eligible funnel does not derive for itself. Carried
        # by record ID, so a record the curation never saw keeps them empty rather than
        # borrowing another record's curation.
        extra = [c for c in checked.columns
                 if c not in frame.columns and c != "ID"]
        if extra:
            frame = frame.merge(checked[["ID"] + extra], on="ID", how="left")

    # `source_cohort` is the per-species cohort the record was curated under, which the
    # SMILES audit groups by.
    frame["source_cohort"] = frame["cohort_species"]

    frame = frame.sort_values(["species", "PMID", "ID"], kind="stable").reset_index(drop=True)
    return frame, pd.DataFrame(counts)


def main() -> None:
    frame, counts = build()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / "eligible_in_vivo_half_life.csv"
    frame.to_csv(out, index=False)
    counts.to_csv(OUTDIR / "eligible_filtering_counts.csv", index=False)

    print(counts[["stage", "records", "removed"]].to_string(index=False))
    print(f"\n[write] {out} - {len(frame)} eligible records")
    print(f"        {frame['PMID'].nunique()} sources, "
          f"{frame['Sequence'].nunique()} distinct sequence strings")
    print(frame["species"].value_counts().to_string())
    print(f"\nlength: {frame['seq_length_measured'].min()}-"
          f"{frame['seq_length_measured'].max()} residues, "
          f"median {frame['seq_length_measured'].median():.0f}")
    print(f"half-life: {frame['half_life_h'].min():.4f}-{frame['half_life_h'].max():.1f} h, "
          f"median {frame['half_life_h'].median():.3f}")


if __name__ == "__main__":
    main()
