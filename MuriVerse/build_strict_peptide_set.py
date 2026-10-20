"""Build the peptide-only rodent half-life cohort with a verified SMILES for every record.

The source is `rodent/eligible_in_vivo_half_life.csv`: live-dosed mouse and rat records
with positive labels. This stage keeps a record only when it is a single peptide chain of at
most 50 residues that can be written down faithfully, carries no attached half-life partner,
reports an elimination half-life, and has a SMILES that survives an independent decode.

  * The sequence must parse residue by residue (`peptide_smiles.parse_residues`); a string
    reconstructed from prose, a repeat or polymer notation, a staple placeholder or a missing
    sequence is excluded, never handed to ESM-2 or HMMER as though it were a sequence.
  * Pegylated, polymer, Fc, albumin and depot/particle records, and any record with a named
    or unresolved partner, are excluded: their clearance is the partner's, not the peptide's.
  * Absorption, distribution and activity half-lives are excluded. Valid elimination values
    are not trimmed by magnitude; approximate, range and one-sided values are flagged in
    `half_life_qualifier`, and unusually long values receive a QC flag.
  * The bare peptide is independently decoded and attached-fragment masses are checked. A
    record whose structure cannot be drawn and verified is excluded, so every kept record has
    a verified SMILES and the molecular encoder never receives an unverified string.
    Unresolved modifications remain explicit tiers (`partly resolved`,
    `modification not drawn`).

Every excluded record is listed with its reason in `rodent/strict/exclusions.csv`. Repeated
observations are preserved. Candidate repeats with the same structure, species and label are
flagged in `repeat_observation_candidates.csv` for evidence-based review; the pipeline does
not infer duplication from sequence and half-life alone. The stage counts are written to
`rodent/strict/filtering_counts.csv`.

Usage: python MuriVerse/build_strict_peptide_set.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import halflife_condition_curation as curation  # noqa: E402
import peptide_mods as mods  # noqa: E402
import subset_in_vivo_half_life as subset  # noqa: E402
from half_life_text import half_life_qualifier, normalise_unit  # noqa: E402
from peptide_smiles import (  # noqa: E402
    MACROMOLECULAR,
    UNSUPPORTED_RESIDUES,
    build_smiles,
    chirality_pattern,
    SARCOSINE,
    esm_sequence,
    looks_like_prose,
    parse_residues,
    reconstruct_residues,
    render_residues,
    residue_template,
)

SOURCE = HERE / "rodent" / "eligible_in_vivo_half_life.csv"
OUTDIR = HERE / "rodent" / "strict"

# The engineering classes whose clearance is dominated by a non-peptide partner. The gate
# is applied on top of the structural gates rather than instead of them: a partner is
# often named only in `Name` or `Origin`, which the SMILES builder never reads, so 164 of
# the 928 structurally "complete" pooled records carry one.
ENGINEERED_PARTNER_CLASSES = {
    "pegylated",
    "fc fusion",
    "albumin-based",
    "polymer conjugate",
    "depot / particle formulation",
}

# Repeat and polymer notation in `Sequence`: "(VPGAG)120", "(hpcdypewqwlcelGk)3", "pSar200-b-...",
# "Poly(L-glutamic acid)". The residue parser would read the repeat unit and drop the count,
# so a 600-residue elastin-like polypeptide passes the 50-residue cap as a pentapeptide.
MULTIMER_NOTATION = re.compile(r"\)\s*\d+\b|\)\s*[nN]\b|\bpoly|[A-Za-z]\d{3,}|\d{3,}", re.I)

# Phases that are a concentration-decay (elimination) half-life. PEPlife2 leaves most records
# unlabelled, which is read as elimination.
ELIMINATION_PHASES = ("unlabelled", "terminal / beta")

# Blocking reasons that describe something BOLTED ON to the peptide rather than a change to
# the peptide itself. A record whose only blocking clauses are these still has a truthful
# structure for its peptide - the partner simply is not in it - so it is kept in the
# `partner carried` tier, with the partner reaching the model through `partner_class`,
# `peg_kda` and `partner_kda` instead of through the string.
#
# Everything NOT listed here alters the chain: a sarcosine or an N-methylated residue is a
# different backbone, a staple is a new covalent bond inside the peptide, and a clause the
# classifier cannot read at all could be either. Those still exclude, because for them the
# string would not be the peptide, which is the one thing this set guarantees.
PARTNER_REASONS = {
    "PEG conjugate",
    "protein or peptide partner",
    "particle or polymer vehicle",
    "dye structure not public",
    "conjugate of unstated composition",
    "conjugate composition not resolved",
    "drug or cytotoxin payload",
    "radiometal without a named chelator",
    "fatty acyl named without a chain length",
    "glycosylation",
}

# Modification text the classifier could not read but which names a payload, carrier, dye or
# conjugation handle: a drug (gemcitabine), a particle (sHDL), a protein or peptide partner
# (ghrelin 12-28, HCG), a targeting ligand (folate), a maleimide handle for protein
# conjugation, a fluorocarbon or cholesterol anchor. These are attached partners.
PARTNER_TEXT = re.compile(
    r"gem[ci]+tabine|\bshdl\b|evans ?blue|folate|\bfa-pep|cholesterol|sialic|glucosyl|\bhcg\b|"
    r"maleimid|mpb-mal|\bmpa\b|fluorocarbon|boronic|lactic acid|\bes2\b|\brccms?\b|"
    r"ghrelin|isobutyl citrate|nanoparticle|\bnps?\b", re.I)

# Label reasons of the form "<dye> label: attachment site not stated" are partners too; the
# dye is known, only where it went is not.
PARTNER_REASON_SUFFIX = "label: attachment site not stated"

# How completely the emitted string describes the dosed molecule, worst last. Every kept
# record has a verified structure; the label is a feature, so a model can learn to discount
# the lower tiers rather than being handed them silently.
#
#   complete              the string is the whole dosed molecule
#   partly resolved       the peptide is drawn but one reported feature - topology,
#                         configuration, a terminus - is not specified in the string
#   modification not drawn a clause that alters the chain could not be drawn, so the string
#                         is the peptide WITHOUT it
TIERS = ["complete", "partly resolved", "modification not drawn"]


def worst_tier(*tiers: str) -> str:
    """The least complete of the tiers given, which is what a record is labelled with."""
    return max((t for t in tiers if t), key=TIERS.index, default="complete")


def is_partner_reason(reason: str) -> bool:
    return reason in PARTNER_REASONS or reason.endswith(PARTNER_REASON_SUFFIX)

# Residues longer than this are proteins, whose renal-threshold clearance is a different
# regime. 50 is the working boundary: a chain that long is still cleared largely by
# glomerular filtration and peptidase action, whereas past it the molecule begins to behave
# like a small protein and the half-life is set by size rather than by the pharmacophore.
# (It was 100 until this cap was tightened; 100 was chosen as the point where ESM-2's
# 1022-token window starts truncating, which is a property of the encoder rather than of the
# biology, and is no longer the binding constraint at 50.)
MAX_RESIDUES = 50

# Records excluded by hand, by PEPlife2 ID, each with the reason. A curated list sits
# uncomfortably beside the rule-based gates, so it stays short and every entry has to name a
# defect a gate cannot see. The repo already works this way upstream, in
# subset_in_vivo_half_life.py's CURATED_DUPLICATES and CURATED_NOT_A_HALF_LIFE.
CURATED_EXCLUSIONS = {
    5555: ("'Free-ADI' is an 83-residue fragment of arginine deiminase, a ~410-residue "
           "enzyme, and its 1445 h label - the largest in the cohort by a wide margin - is "
           "an 'Activity Half Life', a stability readout rather than plasma concentration "
           "decay. Neither defect is visible to the structural gates: the fragment is a "
           "clean L-linear sequence of templated residues, with both termini free and no "
           "Chem-MOD, so it draws and verifies perfectly. It is the sequence's PROVENANCE "
           "and the label's MEANING that disqualify it. At 83 residues it would now also "
           "fail the peptide-scale gate, but it is kept listed here because the curated "
           "gate fires first and because the reason it does not belong is not its length."),
    4547: ("The source sequence is only 'E' while the named dosed agent is cyclic "
           "3PRGD2 with a PEG4 radiolabel; E is not its peptide sequence."),
    4549: ("The source sequence is only 'E' while the named dosed agent is cyclic "
           "3PRGD2 with a PEG4 radiolabel; E is not its peptide sequence."),
    4551: ("The source sequence is only 'E' while the named dosed agent is cyclic "
           "3PRGD2 with a PEG4 radiolabel; E is not its peptide sequence."),
    4553: ("The source sequence is only 'E' while the named dosed agent is cyclic "
           "3PRGD2 with a PEG4 radiolabel; E is not its peptide sequence."),
    4405: ("[3H]BPC157's 102 h is total radioactivity, not intact peptide. The SAME study "
           "(PMID 36588717, doi:10.3389/fphar.2022.1026182) reports the prototype's "
           "elimination half-life as under 30 min and contributes five LC-MS records of the "
           "identical sequence at 0.13-0.50 h, which stay in the cohort; the paper states "
           "that [3H]BPC157 is 'rapidly metabolized into a variety of small peptide "
           "fragments, thus forming single amino acids that entered normal amino acid "
           "metabolism', so the 358-fold gap is the tritium recycling into the amino-acid "
           "pool. No structural gate can see this: the molecule, the species and the matrix "
           "are the same as the five records kept, and `radio/fluorescent label only` is a "
           "KEPT engineering class because labelled records are not systematically biased "
           "(their median half-life is 0.22 h against 0.35 h for the rest). It is the "
           "label's MEANING - what quantity decayed - that disqualifies it, exactly as for "
           "5555 above."),
}

CURATED_EXCLUSIONS[2183] = (
    "The 1.5 h recorded for this intermedin/adrenomedullin analogue is not a measurement. "
    "WO2012138867A2 contains no rat pharmacokinetic experiment: its only rodent dosing is "
    "intraperitoneal, for a blood-pressure readout, and the sole '1.5 hours' anywhere in "
    "the document is claim 1 - 'A biologically active intermedin or adrenomedullin peptide "
    "having a serum half-life of greater than 1.5 hours' - which is a patentable threshold "
    "the invention must clear, not a value observed in an animal. The record also carries "
    "'Rat plasma (Intravenous)' and an ELISA, neither of which the patent supports. No "
    "structural gate can see any of this: the peptide draws and verifies cleanly, and a "
    "label that is a legal claim rather than an observation is exactly the kind of defect "
    "this list exists for. The route curation reaches the same verdict from the other "
    "direction - no primary source establishes a route - but that gate should not be what "
    "this record depends on, because the label would be wrong at any route.")

CURATED_EXCLUSIONS[5654] = (
    "GsMTx4's 168 h ('1 Week') is the half-life of the D-enantiomer GsMTx4-D washing out of "
    "gastrocnemius, heart and diaphragm muscle (PMID 30174173: 'GsMTx4-D showed a 1-week "
    "half-life when administered by subcutaneous injection over four weeks'), recorded as "
    "'Mice plasma'. It is a tissue washout, not a circulating half-life, and it is the "
    "largest label in the set.")

# Rule-based exclusion read from the `Test Sample` text: the animal or the sampled matrix is
# not a normal circulating-PK setting. Organ ablation and receptor knockout change clearance
# directly; disease models and tissue or lymph matrices do not report plasma elimination.
NON_PHYSIOLOGICAL_ANIMAL = re.compile(
    r"nephrectom|hepatectom|infarct|tumou?r[- ]?(?:model|bearing)|bearing|xenograft|db/db|"
    r"ob ?/ ?ob|obese|diabetic|spontaneously hypertensive|\bSHR\b|transgenic|-/-|Brattleboro",
    re.IGNORECASE)
NON_CIRCULATING_MATRIX = re.compile(
    r"tissue|lysate|homogenate|pancrea|stomach|brain|liver|kidney|lymph", re.IGNORECASE)
CIRCULATING_MATRIX = re.compile(r"plasma|serum|blood", re.IGNORECASE)

# Engineering class corrected by hand, by PEPlife2 ID. `hl_engineering` is pattern-matched in
# subset_in_vivo_half_life.py over Name / Origin / Chem-MOD / the two termini, so a half-life
# partner that is named only in the paper's text and never in those fields is invisible to it
# and the record is scored `unmodified / none stated`. These overrides are applied BEFORE the
# partner gate, so a record whose real class turns out to be an engineered partner is then
# excluded by that existing rule rather than by a second hand-maintained list. Every entry
# cites the paper that establishes the class.
CURATED_ENGINEERING = {
    4958: ("fc fusion",
           "TB01-3 is an antibody, not a bare peptide. PMID 33706686 "
           "(doi:10.1080/19420862.2021.1893425, mAbs) panned a synthetic antibody library "
           "built by 'incorporating GPCR binding motifs into the heavy chain "
           "complementarity-determining region 3'; the recorded 39-residue sequence is that "
           "GLP-1-derived CDR3 graft, while the 168 h that was measured by ELISA is the "
           "IgG scaffold's half-life - a week, the textbook rodent IgG value. The name "
           "'TB01-3' matches no ENGINEERING pattern, so the record reached the cohort as "
           "`unmodified / none stated` with `partner_class = none`."),
    5045: ("albumin-based",
           "LM06 is an oxyntomodulin derivative fused to an albumin-binding helix domain. "
           "PMID 33336763 (doi:10.26355/eurrev_202012_24038) states that 'different helix "
           "albumin-binding domains were fused to the mutated OXM via a thrombin-cleavable "
           "linker', and selected LM06 for its 'higher albumin-binding affinity, "
           "sustained-release efficiency'. Its 53-74 h is therefore albumin's clearance, "
           "which is the stated reason `albumin-based` is an excluded partner class. The "
           "name 'LM06' matches no ENGINEERING pattern."),
    5046: ("albumin-based", "Second LM06 dose group; see 5045."),
    5047: ("albumin-based", "Third LM06 dose group; see 5045."),
    5330: ("lipidated",
           "PP18 is a fatty-acid-acylated oxyntomodulin derivative: PMID 31389463 "
           "(doi:10.1039/c9ob01132c) 'site-specifically modified [PP11] with three lengths "
           "of fatty acid chains to provide long-acting conjugates PP13-PP24', of which "
           "PP18 is one. `lipidated` is a KEPT class, so this correction relabels the "
           "record rather than removing it - 29 lipidated records are already in the "
           "cohort - but leaving a 38.5 h acylated peptide labelled `unmodified` "
           "mis-prices the unmodified class in Figure 2d and in the model's covariates."),
}

GATES = [
    ("eligible rodent records", "mouse or rat, dosed in vivo, <= 50 residues; see build_rodent_eligible_set.py"),
    ("not curated out", "a defect no structural gate can see; see CURATED_EXCLUSIONS"),
    ("physiological circulating half-life", "organ-ablated, knockout, transgenic and disease-model "
                                            "animals, and tissue, lysate or lymph samples, are not "
                                            "normal plasma elimination"),
    ("no engineered partner class", "pegylated, polymer, Fc, albumin and depot/particle records: "
                                    "clearance is the partner's, not the peptide's"),
    ("single peptide chain", "repeat and polymer notation ('(VPGAG)120', 'pSar200') is read as its "
                             "repeat unit, which is not the molecule dosed"),
    ("faithful residue sequence", "every character of `Sequence` must be a residue; prose, "
                                  "staple placeholders, unsupported residues and a missing "
                                  "sequence cannot give an ESM-2 or HMMER input"),
    (f"peptide scale (<= {MAX_RESIDUES} aa)", "longer chains are proteins in a different "
                                              "clearance regime"),
    ("no other attached partner", "any record whose name, termini or modifications name or "
                                  "leave unresolved an attached partner (carrier, dye, drug, "
                                  "protein or polymer)"),
    ("elimination half-life", "absorption, distribution and activity half-lives; no value "
                              "magnitude cutoff is applied"),
    ("structure verified", "the SMILES could not be built, or the emitted string contradicts the "
                           "sequence it was built from, or its conjugate mass delta is wrong"),
]


def sequence_problem(record: pd.Series) -> tuple[list[tuple[str, str]] | None, str]:
    """(parsed residues, problem) for one record. Residues are None iff there is a problem.

    Read from `Sequence`, PEPlife2's ORIGINAL column, not from `sequence`, the flattened
    derivative the rest of the pipeline uses. The two differ on 403 of 1,300 rodent rows and
    the difference is not cosmetic:

      * `Sequence` marks a D residue by lower-casing its letter ("rPKPfQwFwLL"), so a mixed
        peptide's D positions are known individually. `sequence` is upper-cased, so that
        information is simply gone - which is why 32 records were excluded for having
        chirality the pipeline could not resolve when in fact the source states it exactly.
      * `Sequence` delimits a non-natural residue ("pGlu-LYENKPRRPYIL"); `sequence` strips
        the punctuation to "PGLULYENKPRRPYIL", which is then detected as prose and rejected.

    Where the flattened column disagrees with an upper-cased reading of the original, the
    record is refused: two columns of the same database contradicting each other about which
    residues a peptide contains is not a basis for drawing it.
    """
    residues, problem = parse_residues(record.get("Sequence"))
    if residues is None:
        return None, problem
    if looks_like_prose(render_residues(residues).upper()):
        return None, "residue names written as prose"

    unsupported = sorted({t for t, _ in residues if t in UNSUPPORTED_RESIDUES})
    if unsupported:
        return None, (f"no residue template for {', '.join(unsupported)}: it alters the "
                      f"backbone rather than the side chain")
    unknown = sorted({t for t, _ in residues if residue_template(t) is None})
    if unknown:
        return None, f"no residue template for {', '.join(unknown)}"

    flattened = re.sub(r"[^A-Z]", "", str(record.get("sequence") or "").upper())
    plain = "".join(t for t, _ in residues if len(t) == 1).upper()
    if flattened and plain and flattened != plain and len(flattened) == len(plain):
        return None, (f"the two sequence columns disagree: `Sequence` reads {plain} and "
                      f"`sequence` reads {flattened}")
    return residues, ""


HYDROGEN = 1.00782503207


def fragment_mass(fragment: str) -> float:
    """Mass the fragment contributes once bonded, which is TWO hydrogens short of its own.

    Both losses are easy to miss, and leaving either one out makes the check fail on every
    conjugated record by a flat 1.008 Da:

      * `_attach` removes one hydrogen from the anchor atom - the lysine N-zeta, say - to
        free a valence for the new bond. That hydrogen was counted in the bare peptide.
      * deleting the '*' placeholder from the fragment leaves its attachment carbon one
        bond short, and RDKit satisfies that valence with an implicit hydrogen. In the
        assembled molecule that bond goes to the anchor instead, so this hydrogen is not
        there either.
    """
    mol = Chem.MolFromSmiles(fragment)
    dummy_free = Chem.DeleteSubstructs(mol, Chem.MolFromSmiles("*"))
    return ExactMolWt(dummy_free) - 2 * HYDROGEN


def verify(residues: list[tuple[str, str]], bare: str, full: str,
           fragments: list[str]) -> str:
    """'' when both checks pass, else the reason they did not.

    The round trip is over residue IDENTITY, ORDER and CONFIGURATION: `reconstruct_residues`
    reads each alpha carbon's CIP descriptor back out of the string and maps it to L or D, so
    a peptide drawn with D positions taken from the sequence's own case is confirmed to carry
    them. A constitution-only check would pass an all-L string for a mixed-chirality peptide,
    which is exactly the error this set exists to exclude.
    """
    decoded, cyclic = reconstruct_residues(bare)
    if decoded is None:
        return "bare peptide did not decode"

    # Glycine, sarcosine and Aib have no alpha stereocentre, so the decoder reports '' for them
    # and the forward side's setting for those positions is not a claim about anything. A
    # multi-character residue (Nle, Orn, pGlu...) is parsed without a configuration - the
    # sequence spelling carries none - so only a configuration the record actually states is
    # compared; the decoder's reading of the others is not a contradiction.
    def tokens_of(items):
        return [t for t, _ in items]

    def configs_agree(want_items, got_items):
        return all((not w) or t in ("G", "Aib", SARCOSINE) or w == g
                   for (t, w), (_, g) in zip(want_items, got_items))

    def matches(want_items, got_items):
        return tokens_of(want_items) == tokens_of(got_items) and configs_agree(want_items, got_items)

    if not matches(residues, decoded):
        rotations = [residues[i:] + residues[:i] for i in range(len(residues))]
        if not (cyclic and any(matches(r, decoded) for r in rotations)):
            return (f"bare peptide decoded to {render_residues(decoded)}, not the recorded "
                    f"{render_residues(residues)}")

    if not fragments:
        return "" if bare == full else "conjugate-free record changed on assembly"
    bare_mol, full_mol = Chem.MolFromSmiles(bare), Chem.MolFromSmiles(full)
    if bare_mol is None or full_mol is None:
        return "assembled structure did not parse"
    observed = ExactMolWt(full_mol) - ExactMolWt(bare_mol)
    expected = sum(fragment_mass(f) for f in fragments)
    if abs(observed - expected) > 0.01:
        return (f"conjugate mass delta {observed:.3f} does not match the expected "
                f"{expected:.3f}")
    return ""


def assess(row: pd.Series) -> dict:
    """Run every gate on one record, in funnel order, stopping at the first failure."""
    out = {
        "gate_failed": "",
        "gate_detail": "",
        "smiles": "",
        "smiles_bare": "",
        "conjugates": "",
        "tracer_label": "",
        "sequence_tokens": "",
        "sequence_esm": "",
        "nonstandard_positions": "",
        "n_nonstandard": 0,
        "chirality_from_sequence": "",
        "stereo_source": "",
        "n_conjugates": 0,
        "n_residues": 0,
        "structure_tier": "complete",
        "unresolved": "",
        "substitution_unconfirmed": "",
        "partner_text": "",
        "partner_source": "",
        "hl_engineering_corrected": "",
        "hl_engineering_correction_reason": "",
    }

    curated = CURATED_EXCLUSIONS.get(row.get("ID"))
    if curated is not None:
        return {**out, "gate_failed": "not curated out", "gate_detail": curated[0]}

    sample = str(row.get("Test Sample") or "")
    animal = NON_PHYSIOLOGICAL_ANIMAL.search(sample)
    if animal:
        return {**out, "gate_failed": "physiological circulating half-life",
                "gate_detail": f"non-physiological animal: {animal.group(0)!r} in {sample!r}"}
    if NON_CIRCULATING_MATRIX.search(sample) and not CIRCULATING_MATRIX.search(sample):
        return {**out, "gate_failed": "physiological circulating half-life",
                "gate_detail": f"non-circulating matrix: {sample!r}"}

    # Correct the engineering class BEFORE the partner gate reads it, so that a partner the
    # pattern matcher could not see is excluded by that gate and not by a separate list.
    correction = CURATED_ENGINEERING.get(row.get("ID"))
    engineering = correction[0] if correction else row.get("hl_engineering")
    if correction:
        out["hl_engineering_corrected"] = engineering
        out["hl_engineering_correction_reason"] = correction[1]

    if engineering in ENGINEERED_PARTNER_CLASSES:
        detail = f"engineering class: {engineering}"
        if correction:
            detail += (f" (corrected from {row.get('hl_engineering')!r}: {correction[1]})")
        return {**out, "gate_failed": "no engineered partner class", "gate_detail": detail}

    raw_sequence = row.get("Sequence")
    if isinstance(raw_sequence, str) and MULTIMER_NOTATION.search(raw_sequence):
        return {**out, "gate_failed": "single peptide chain",
                "gate_detail": f"repeat or polymer notation: {raw_sequence[:60]}"}

    residues, problem = sequence_problem(row)
    if residues is None:
        # A string that is not a residue-by-residue sequence is not an input for ESM-2 or
        # HMMER, and a letter-by-letter fallback would feed them fragments of residue names.
        return {**out, "gate_failed": "faithful residue sequence", "gate_detail": problem}
    tokens = [t for t, _ in residues]
    # The canonical spelling, carrying the D positions in its case, which is what the
    # builder reads and what the round trip is checked against.
    sequence = render_residues(residues)
    # What ESM-2 is actually given. A residue it was never pretrained on becomes X, and
    # which positions those are is recorded rather than left implicit: the molecular branch
    # sees the true structure while the sequence branch sees an unknown token, and that
    # asymmetry is a measurable quality gradient on these rows.
    esm, blurred = esm_sequence(tokens)
    out = {**out, "n_residues": len(tokens),
           "sequence_tokens": sequence, "sequence_esm": esm,
           "nonstandard_positions": " ".join(str(i) for i in blurred),
           "n_nonstandard": len(blurred),
           "chirality_from_sequence": chirality_pattern(residues)}

    if len(tokens) > MAX_RESIDUES:
        return {**out, "gate_failed": f"peptide scale (<= {MAX_RESIDUES} aa)",
                "gate_detail": f"{len(tokens)} residues"}

    # A partner no longer ends the assessment. The peptide is still drawn, truthfully and
    # in full; what is bolted onto it is recorded here and reaches the model as a covariate.
    partners: list[tuple[str, str]] = []

    # The terminal columns belong to the SMILES builder, which recognises caps that match a
    # whole field. A field it does not recognise is handed to `peptide_mods` instead, which
    # can resolve a site-specific conjugate written there ("Lys4 linked with Palm fatty
    # acid"); whichever module takes a field, something must account for it.
    from peptide_smiles import _c_terminal_kind, _n_terminal_fragment, _normalise
    n_text, c_text = _normalise(row.get("N-ter MOD")), _normalise(row.get("C-ter MOD"))
    n_raw, c_raw = row.get("N-ter MOD"), row.get("C-ter MOD")
    # A transposed pair of terminal fields, which the cohort carries at least once:
    # "N-ter MOD = Amidation, C-ter MOD = Acetylation". Neither cap belongs where it is
    # written. Amidation converts a CARBOXYL to an amide and has no meaning at an amino
    # terminus, and an acyl cap needs an amine to acylate, so the two fields can only have
    # been swapped. Each field is recognised at the OTHER terminus and at neither its own,
    # which is the test applied here - a field that fits where it stands is never moved.
    if (_c_terminal_kind(n_text)[1] and _n_terminal_fragment(c_text)[1]
            and not _n_terminal_fragment(n_text)[1] and not _c_terminal_kind(c_text)[1]):
        n_text, c_text = c_text, n_text
        n_raw, c_raw = c_raw, n_raw
    _, n_is_cap = _n_terminal_fragment(n_text)
    _, c_is_cap = _c_terminal_kind(c_text)

    extra_fields = [row.get("Chem-MOD")]
    n_arg, c_arg = n_raw, c_raw
    for text, is_cap in ((n_text, n_is_cap), (c_text, c_is_cap)):
        # peptide_smiles already knows which terminal descriptions name a partner no small
        # molecule can stand in for, and that knowledge must be consulted before the field
        # is handed to peptide_mods: "AFF" is an affibody, and a generic clause classifier
        # reading it as a residue name let a 15 kDa protein fusion through this gate.
        if not is_cap and MACROMOLECULAR.search(text):
            partners.append(("terminus", text))
    # A terminal field the builder does not recognise as a cap is handed to peptide_mods
    # below and the terminus itself drawn free, whether that field named a partner or a
    # site-specific conjugate. That is what makes the peptide drawable in both tiers.
    terminal_partner = {text for kind, text in partners if kind == "terminus"}
    if not n_is_cap:
        if n_text not in terminal_partner:
            extra_fields.append(n_text)
        n_arg = "Free"
    if not c_is_cap:
        if c_text not in terminal_partner:
            extra_fields.append(c_text)
        c_arg = "Free"

    verdict = mods.report(tokens, *extra_fields)
    out["tracer_label"] = verdict.tracer_label
    # Reported separately from the other blockers, because it is a different kind of
    # problem: not a modification this module cannot draw, but an annotation the sequence
    # column cannot be reconciled with. See the note above peptide_mods.SUBSTITUTION.
    unconfirmed = [c for c in verdict.blocking if c.unconfirmed]
    if unconfirmed:
        # The sequence can neither confirm nor refute the claim. The peptide still draws as
        # the sequence states it; what is uncertain is recorded rather than used to discard.
        out["substitution_unconfirmed"] = (f"{unconfirmed[0].reason}: "
                                           f"{unconfirmed[0].text}")
    altering = [c for c in verdict.blocking if not is_partner_reason(c.reason)]
    payload = [c for c in altering if PARTNER_TEXT.search(str(c.text or ""))]
    partners.extend(("modification", f"{c.reason}: {c.text}") for c in payload)
    altering = [c for c in altering if c not in payload]
    tier_from_mods = ""
    if altering:
        # The clause changes the chain and cannot be drawn, so the string will be the
        # peptide WITHOUT it. That is a real limitation and is labelled as one.
        first = altering[0]
        tier_from_mods = "modification not drawn"
        out["unresolved"] = f"{first.reason}: {first.text}"
    partners.extend(("modification", f"{c.reason}: {c.text}")
                    for c in verdict.blocking if is_partner_reason(c.reason))

    conjugates = [(site.position, c.fragment)
                  for c in verdict.drawable for site in c.sites]
    out["conjugates"] = "; ".join(
        f"{site.residue}{site.position} {c.reason}"
        for c in verdict.drawable for site in c.sites)
    out["n_conjugates"] = len(conjugates)

    # The `Name` column is the last place a partner hides, and nothing upstream reads it.
    # "Au-AR pep-PROTAC" is a gold-nanoparticle conjugate with a 26.3 h half-life - ~100x
    # this set's median - whose every annotation column says Free and unmodified, so it
    # passed every other gate. `hl_engineering` misses it too, because its ENGINEERING
    # patterns have no term for a gold carrier.
    #
    # Only Name is screened, not Origin: Origin is provenance prose about where the
    # SEQUENCE came from ("Derived from the tenth domain of fibronectin type III"), which
    # describes nothing attached to the molecule. Screening it would reject a plain peptide
    # for its ancestry.
    name = " ".join(str(row.get("Name") or "").split())
    for match in MACROMOLECULAR.finditer(name):
        matched = match.group(0).lower()
        # A PEG named here is accounted for when the builder actually drew a discrete one.
        if "peg" in matched and any(c.reason.startswith("discrete PEG")
                                    for c in verdict.drawable):
            continue
        partners.append(("name", f"{matched} in Name: {name}"))
        break

    if partners:
        return {**out, "gate_failed": "no other attached partner",
                "gate_detail": " | ".join(text for _, text in partners)[:300]}

    phase = row.get("pk_phase")
    if phase not in ELIMINATION_PHASES:
        return {**out, "gate_failed": "elimination half-life",
                "gate_detail": f"{phase}: {row.get('pk_label')}"}
    built = build_smiles(sequence, n_arg, c_arg, row.get("Linear/Cyclic"),
                         row.get("Chirality"), conjugates=conjugates)
    bare = build_smiles(sequence, n_arg, c_arg, row.get("Linear/Cyclic"),
                        row.get("Chirality"))

    unspecified = []
    if not built.n_terminal_encoded or not built.c_terminal_encoded:
        unspecified.append("terminus: " +
                           (n_text if not built.n_terminal_encoded else c_text))
    if not built.topology_encoded:
        unspecified.append("topology: " + _normalise(row.get("Linear/Cyclic")))
    if not built.chirality_resolved:
        unspecified.append("chirality: " + _normalise(row.get("Chirality")))
    if unspecified:
        out["unresolved"] = "; ".join(filter(None, [out.get("unresolved", "")] + unspecified))

    if built.smiles is None or bare.smiles is None:
        return {**out, "gate_failed": "structure verified",
                "gate_detail": built.status if built.smiles is None else bare.status}

    out["smiles"], out["smiles_bare"] = built.smiles, bare.smiles
    out["stereo_source"] = built.stereo_source
    failure = verify(residues, bare.smiles, built.smiles,
                     [f for _, f in conjugates])
    if failure:
        # A string that failed the independent check is never passed to PeptideCLM, and a
        # record without a verified structure is not kept.
        return {**out, "smiles": "", "smiles_bare": "",
                "gate_failed": "structure verified", "gate_detail": failure}

    # The tier is decided last, on a record whose PEPTIDE has already drawn and verified.
    # `complete` means the string is the whole dosed molecule. `partner carried` means the
    # string is the whole PEPTIDE and something else is attached to it, named here and
    # encoded as a covariate rather than silently dropped or silently drawn wrong.
    out["structure_tier"] = worst_tier(
        "partly resolved" if unspecified else "",
        tier_from_mods)
    return out


def build() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    source = pd.read_csv(SOURCE, low_memory=False)
    # The phase buckets were widened (absorption, "alpha" spelled out), so the phase is
    # re-derived from the label here rather than trusted from the eligible file.
    source["pk_phase"] = subset.bucket(source["pk_label"].fillna(""), subset.PK_PHASE,
                                       "unlabelled")
    source["half_life_qualifier"] = [
        half_life_qualifier(v, u) for v, u in zip(source["Half Life"], source["Units Half Life"])]
    source["half_life_unit"] = source["Units Half Life"].map(normalise_unit)
    assessed = pd.DataFrame([assess(row) for _, row in source.iterrows()],
                            index=source.index)
    frame = pd.concat([source, assessed], axis=1)
    # Carry the curated class onto the record itself, so a correction that does NOT trip the
    # partner gate (PP18, `lipidated`) travels with the kept row into the model's covariates
    # and into Figure 2d, instead of only explaining an exclusion.
    corrected = frame["hl_engineering_corrected"].astype(str) != ""
    frame.loc[corrected, "hl_engineering"] = frame.loc[corrected, "hl_engineering_corrected"]

    counts = [{"stage": GATES[0][0], "records": len(frame), "removed": 0,
               "rationale": GATES[0][1]}]
    kept = frame
    for stage, rationale in GATES[1:]:
        before = len(kept)
        kept = kept[kept["gate_failed"] != stage]
        counts.append({"stage": stage, "records": len(kept),
                       "removed": before - len(kept), "rationale": rationale})

    kept = kept[kept["gate_failed"] == ""].copy()

    # A sequence/species/label match cannot distinguish an independent study, a new route or
    # dose arm, or a different molecular topology. Preserve source records and expose likely
    # repeats as review candidates instead of deleting them automatically.
    kept["half_life_extreme_qc"] = kept["half_life_h"].astype(float) > 168.0
    # These verbatim PEPlife2 descriptors distinguish many treatment arms. Their absence
    # from the key would mislabel different dose/route arms as duplicates. They are a first-
    # pass review key; the full-text evidence table is authoritative when available.
    condition_key = ["PMID", "species", "smiles", "Dose Route", "Concentration",
                     "Test Sample", "Assay"]
    within_key = condition_key + ["half_life_h"]
    across_key = ["species", "smiles", "half_life_h"]
    kept["same_study_structure_condition_n"] = kept.groupby(
        condition_key, dropna=False)["ID"].transform("size").astype(int)
    kept["same_study_structure_condition_label_n"] = kept.groupby(
        within_key, dropna=False)["ID"].transform("size").astype(int)
    kept["same_structure_label_across_sources_n"] = kept.groupby(
        across_key, dropna=False)["ID"].transform("size").astype(int)
    kept["exact_duplicate_candidate"] = (
        kept["same_study_structure_condition_label_n"] > 1)
    kept["within_context_repeated_measurement"] = (
        kept["same_study_structure_condition_n"] > 1)
    kept["cross_source_repeated_measurement"] = (
        kept["same_structure_label_across_sources_n"] > 1)
    kept["repeat_observation_candidate"] = (
        kept["exact_duplicate_candidate"] | kept["cross_source_repeated_measurement"]
    )
    kept["repeat_candidate_scope"] = ""
    kept.loc[kept["exact_duplicate_candidate"], "repeat_candidate_scope"] = (
        "same publication, structure, source arm and label: review for duplicate entry")
    across = kept["cross_source_repeated_measurement"]
    kept.loc[across & ~kept["exact_duplicate_candidate"], "repeat_candidate_scope"] = (
        "same structure/species/label across sources: preserve as independent evidence; "
        "review citation overlap")
    counts.append({"stage": "observations retained; repeat candidates flagged",
                   "records": len(kept), "removed": 0,
                   "rationale": "preserve source IDs; review possible repeats using full "
                               "experimental context"})

    # The dose route, resolved against the primary sources wherever PEPlife2's `Dose Route`
    # column is blank. Done here rather than in the training builder so the published kept
    # file carries the route it was curated to, with the evidence level beside it: a table
    # built from this file and a figure built from the modelling file then say the same
    # thing about the same record. Nothing is dropped for its route here - the cohort-level
    # parenteral requirement belongs to build_halflife_training_set.py.
    kept = curation.apply_route(kept)
    curated = int((kept["route_evidence_level"] != "database_route_column").sum())
    counts.append({"stage": "dose route curated from the sources",
                   "records": len(kept), "removed": 0,
                   "rationale": f"{curated} records had a blank Dose Route and were "
                                f"resolved against the primary sources; see "
                                f"data/curation/route_curation.csv"})

    columns = ["species", "ID", "PMID", "Name", "sequence", "hl_engineering", "half_life_h",
               "N-ter MOD", "C-ter MOD", "Linear/Cyclic", "Chirality", "Chem-MOD",
               "gate_failed", "gate_detail"]
    excluded = frame[frame["gate_failed"] != ""][columns].copy()
    return kept, excluded, pd.DataFrame(counts)


def main() -> None:
    problems = mods.selfcheck()
    if problems:
        raise SystemExit("fragment library is broken:\n" + "\n".join(problems))

    kept, excluded, counts = build()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    kept.to_csv(OUTDIR / "in_vivo_half_life.csv", index=False)
    excluded.to_csv(OUTDIR / "exclusions.csv", index=False)
    counts.to_csv(OUTDIR / "filtering_counts.csv", index=False)
    kept[kept["repeat_observation_candidate"]].to_csv(
        OUTDIR / "repeat_observation_candidates.csv", index=False)

    print(counts[["stage", "records", "removed"]].to_string(index=False))
    print(f"\nkept {len(kept)} records over {kept['PMID'].nunique()} studies")
    print(kept["species"].value_counts().to_string())
    print("\nhalf-life extension classes retained:")
    print(kept["hl_engineering"].value_counts().to_string())
    print(f"\nrecords carrying a drawn side-chain conjugate: "
          f"{(kept['n_conjugates'] > 0).sum()}")
    print(f"records flagged as radiotracer-labelled: "
          f"{(kept['tracer_label'] != '').sum()}")
    print(f"records carrying any PEG: "
          f"{kept['Chem-MOD'].fillna('').str.contains('PEG', case=False).sum()}")
    print(f"records containing a non-standard residue (X to ESM-2): "
          f"{(kept['n_nonstandard'] > 0).sum()}")
    print(f"records whose D positions came from the sequence's own case: "
          f"{(kept['stereo_source'] == 'sequence case').sum()}")
    print("chirality as the sequence states it:")
    print(kept["chirality_from_sequence"].value_counts().to_string())


if __name__ == "__main__":
    main()
