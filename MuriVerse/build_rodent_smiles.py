"""Run the SMILES builder over EVERY record in the eligible rodent set.

`build_halflife_training_set.py` builds SMILES only for the modelling subset: it first drops
non-decay PK phases, sequence-less rows and duplicate feature/label pairs, so its SMILES
ledger describes a subset of the eligible records. This script answers the separate
question "what happens if we try to draw all 1,159?" and writes one row per record with the
outcome, so the figures can show the full funnel rather than the modelling subset's tail.

It does not filter. Every record in `rodent/eligible_in_vivo_half_life.csv` gets a row in
`rodent/smiles_all_records.csv`, including the ones whose sequence column cannot support a
structure at all — those are reported, not silently dropped.

Output columns beyond the identifiers:

  smiles                the string, or empty
  smiles_status         first blocker, or 'complete'
  smiles_coverage       complete | partial | none — the PEPTIDE's own axes only; an
                        attached partner never demotes it (see peptide_smiles)
  smiles_blocker_class  why a record is not `complete`, grouped for plotting
  sequence_problem      when coverage is 'none', what is wrong with the sequence column
  smiles_*_encoded      whether each terminal/topology axis is in the emitted string;
                        blank when no SMILES was emitted
  smiles_chirality_resolved  whether residue stereochemistry is specified in the string
  smiles_other_mod_review   source Chem-MOD text is present but is not parsed by the builder
  smiles_reconstructed_sequence  sequence read back OUT of `smiles` by `reconstruct_sequence`,
                        independent of the `sequence` column above; one rotation of the true
                        sequence for a head-to-tail cyclic peptide, since those have no fixed
                        starting residue
  smiles_reconstructed_cyclic  whether that reverse read treated the molecule as a macrocycle
  smiles_roundtrip_status  'no smiles' | 'exact match' | 'match (cyclic rotation)' |
                        'mismatch' | 'decode failed' — the forward sequence-to-SMILES step
                        checked against an independent SMILES-to-sequence read

Usage: conda run -n peptiverse-local python MuriVerse/build_rodent_smiles.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_halflife_training_set import bucket_mod  # noqa: E402
from peptide_smiles import (  # noqa: E402
    SIDE_CHAIN,
    build_smiles,
    describe_termini,
    looks_like_prose,
    reconstruct_sequence,
)

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "rodent" / "eligible_in_vivo_half_life.csv"
OUT = HERE / "rodent" / "smiles_all_records.csv"

STANDARD = set(SIDE_CHAIN)


def sequence_problem(sequence: object) -> str:
    """Why the sequence column cannot support a structure. Empty when it can."""
    seq = "" if pd.isna(sequence) else " ".join(str(sequence).split()).upper()
    if not seq:
        return "no sequence recorded"
    if looks_like_prose(seq):
        return "residue names written as prose"
    letters = set(seq)
    # B (Asx) and Z (Glx) name a PAIR of residues and X names none; no structure follows from
    # them. O and U are IUPAC codes for ornithine and selenocysteine, but in this cohort they
    # are almost always prose leakage or a transcription slip - see the audit in the caption.
    if letters & set("BZX"):
        return "ambiguity code (B/Z/X) in sequence"
    if letters & set("OU"):
        return "O/U in sequence (not a residue here)"
    if letters - STANDARD:
        return "non-residue character in sequence"
    return ""


def blocker_class(status: str, coverage: str) -> str:
    if coverage == "complete" and status == "complete":
        return "complete"
    if status.startswith("N-terminal partner"):
        return "N-terminal partner not drawable"
    if status.startswith("C-terminal partner"):
        return "C-terminal partner not drawable"
    if status.startswith("cyclization"):
        return "cyclization not drawable"
    if status.startswith("chirality"):
        return "chirality unspecified"
    return "sequence unusable"


def roundtrip_status(sequence: object, smiles: str, reconstructed: str | None,
                     cyclic: bool) -> str:
    """Compare the ORIGINAL sequence column to what `reconstruct_sequence` reads back
    out of `smiles`, independent of the forward builder's own reasoning."""
    if not smiles:
        return "no smiles"
    if reconstructed is None:
        return "decode failed"
    original = "" if pd.isna(sequence) else " ".join(str(sequence).split()).upper()
    if cyclic:
        # a macrocycle has no fixed starting residue: a rotation IS the same sequence
        doubled = original + original
        return ("match (cyclic rotation)"
                if len(reconstructed) == len(original) and reconstructed in doubled
                else "mismatch")
    return "exact match" if reconstructed == original else "mismatch"


def main() -> None:
    df = pd.read_csv(SOURCE, low_memory=False)
    print(f"[in]  {SOURCE.relative_to(HERE.parent)}: {len(df)} records")
    print(df.groupby("source_cohort").size().to_string())

    built = [
        build_smiles(
            "" if pd.isna(r["sequence"]) else r["sequence"],
            r["N-ter MOD"], r["C-ter MOD"], r["Linear/Cyclic"], r["Chirality"],
        )
        for _, r in df.iterrows()
    ]

    out = df[["species", "source_cohort", "ID", "PMID", "Name", "sequence",
              "N-ter MOD", "C-ter MOD", "Linear/Cyclic", "Chirality", "Chem-MOD",
              "topology", "hl_engineering", "half_life_h", "pk_phase"]].copy()
    out["n_ter_mod"] = df["N-ter MOD"].map(bucket_mod)
    out["c_ter_mod"] = df["C-ter MOD"].map(bucket_mod)
    out["chirality"] = df["Chirality"].map(
        lambda v: "not reported" if str(v).strip().lower() in ("", "nan", "n.a.", "na")
        else str(v).strip().lower()
    )
    out["smiles"] = [b.smiles for b in built]
    out["smiles_status"] = [b.status for b in built]
    out["smiles_coverage"] = [b.coverage for b in built]
    out["smiles_n_partner"] = [b.n_partner for b in built]
    out["smiles_c_partner"] = [b.c_partner for b in built]
    out["smiles_stereo_unspecified"] = [b.stereo_unspecified for b in built]
    out["smiles_n_terminal_encoded"] = [b.n_terminal_encoded for b in built]
    out["smiles_c_terminal_encoded"] = [b.c_terminal_encoded for b in built]
    out["smiles_topology_encoded"] = [b.topology_encoded for b in built]
    out["smiles_chirality_resolved"] = [b.chirality_resolved for b in built]
    # What is ATTACHED, read from the source text rather than from the emitted string, so
    # these are unaffected by whether a partner happens to be drawable. The published
    # ledger carries them and `describe_termini` reproduces it exactly on all 1,159 rows;
    # regenerating without them would have silently dropped five columns from an artefact
    # other work reads.
    termini = [describe_termini(r["N-ter MOD"], r["C-ter MOD"]) for _, r in df.iterrows()]
    out["n_ter_class"] = [t.n_class for t in termini]
    out["c_ter_class"] = [t.c_class for t in termini]
    out["partner_class"] = [t.partner_class for t in termini]
    out["peg_kda"] = [t.peg_kda for t in termini]
    out["partner_kda"] = [t.partner_kda for t in termini]
    out["smiles_other_mod_review"] = ~df["Chem-MOD"].fillna("").astype(str).str.strip().str.lower().isin(
        ("", "n.a.", "na", "none", "free", "not reported")
    )
    out["sequence_problem"] = df["sequence"].map(sequence_problem)
    out["smiles_blocker_class"] = [
        blocker_class(b.status, b.coverage) for b in built
    ]

    decoded = [reconstruct_sequence(b.smiles) if b.smiles else (None, False) for b in built]
    out["smiles_reconstructed_sequence"] = [d[0] for d in decoded]
    out["smiles_reconstructed_cyclic"] = [d[1] for d in decoded]
    out["smiles_roundtrip_status"] = [
        roundtrip_status(seq, b.smiles, recon, cyclic)
        for seq, b, (recon, cyclic) in zip(df["sequence"], built, decoded)
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)

    print(f"\n[out] {OUT.relative_to(HERE.parent)}: {len(out)} records")
    print("\ncoverage")
    print(out["smiles_coverage"].value_counts().to_string())
    print("\nblocker")
    print(out["smiles_blocker_class"].value_counts().to_string())
    unusable = out[out["smiles_coverage"] == "none"]
    print(f"\nsequence problems ({len(unusable)} records)")
    print(unusable["sequence_problem"].value_counts().to_string())
    print("\nSMILES -> sequence round trip")
    print(out["smiles_roundtrip_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
