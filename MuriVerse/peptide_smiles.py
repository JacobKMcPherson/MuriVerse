"""Build a SMILES for a modified peptide.

PeptiVerse encodes N-terminal modification, C-terminal modification, cyclization and
chirality in the SMILES string rather than as separate model inputs — its half-life
metadata carries only `sequence, SMILES, half_life_*`. This module reconstructs that
representation for the MuriVerse mouse cohort from PEPlife2's free-text modification
columns.

Roughly a third of the cohort is a peptide conjugated to a partner that has no honest
small-molecule SMILES at all: IgG Fc, serum albumin, GFP, His-tag, XTEN, DARPins,
chondroitin, heparosan, 20-55 kDa PEG, iron-oxide and gold particles. Emitting a bare
backbone for those would claim the peptide was unmodified, and dropping them discards
most of the long-half-life end of the label distribution — which is precisely the
signal the model is meant to learn.

So the builder reports *coverage* alongside the SMILES:

  complete      parsed constitution is represented, possibly with unspecified
                stereochemistry; this does not audit the separate Chem-MOD field
  partial       the peptide is drawn but one of its own axes is not pinned down: an
                unpositioned cyclization, or an L/D series the source leaves unresolved. A
                partner named in `n_partner`/
                `c_partner` is attached in reality and is NOT in the string
  none          no SMILES, because the sequence itself is not usable: it contains a
                residue with no determinate structure (B = Asx, Z = Glx, X = unknown,
                and O/U per the note under SIDE_CHAIN), or the sequence column holds
                prose rather than a sequence

Callers must propagate `coverage` into the model as a feature. A `partial` string
is not a claim about the intact molecule, and treating it as one would be a silent lie
of exactly the kind this module exists to avoid. `stereo_unspecified` is reported the same
way and for the same reason.

The backbone is assembled from residue templates rather than `Chem.MolFromSequence`,
which accepts only the 20 standard amino acids and always emits the L enantiomer. The
templates reproduce RDKit's output atom-for-atom on standard sequences (`selfcheck()`
asserts this, and also that a disulfide's ring label never collides with an aromatic side
chain's), and additionally carry all-D and stereochemically unspecified peptides.

`reconstruct_sequence` runs the builder in reverse: given only a generated SMILES string
(no sequence, no annotation columns), it walks the backbone and re-derives the one-letter
sequence, so the forward direction can be checked by a genuine independent read of the
string rather than by re-deriving the same thing twice. It is not a general SMILES parser —
it relies on the same residue templates and only has to undo what `build_smiles` does:
an optional N-terminal acyl cap, the residue chain (optionally disulfide-bridged or
head-to-tail cyclised), and an optional C-terminal amide or methyl ester. Chirality is
deliberately not decoded (that axis is already checked directly from the annotation
columns, in `smiles_chirality_resolved`); only residue identity is.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# B = Asx (Asn or Asp), Z = Glx (Gln or Glu), X = unknown: no determinate structure.
# O and U are here for a different reason — see the note under SIDE_CHAIN.
AMBIGUOUS_AA = set("BZXOU")

# Side chain hung off CA. Glycine and proline are special-cased in `_backbone`.
SIDE_CHAIN = {
    "A": "C",
    "C": "CS",
    "D": "CC(=O)O",
    "E": "CCC(=O)O",
    "F": "Cc1ccccc1",
    "G": None,
    "H": "Cc1c[nH]cn1",
    "I": "[C@@H](C)CC",
    "K": "CCCCN",
    "L": "CC(C)C",
    "M": "CCSC",
    "N": "CC(N)=O",
    "P": None,
    "Q": "CCC(N)=O",
    "R": "CCCNC(=N)N",
    "S": "CO",
    "T": "[C@@H](C)O",
    "V": "C(C)C",
    "W": "Cc1c[nH]c2ccccc12",
    "Y": "Cc1ccc(O)cc1",
}
# O and U deliberately have no template. PEPlife2 does define them in its own free text
# ("O denotes L-ornithine"), but in this cohort most O and U characters are not residue
# codes at all — they are fragments of residue names that were transcribed into the
# sequence column as prose: PROPARGYLGLYCINE, NMEVALETPRONMETHR, PGLULYENK, PSARBPCYSSOET,
# HOOFOPHLEPFFF. Mapping O to ornithine there would invent a structure for a row whose
# sequence is not a sequence. `AMBIGUOUS_AA` below therefore covers them, and the prose
# rows are reported separately by `looks_like_prose`.

# Side-chain templates carrying an atom map number on the atom that accepts a conjugate.
# The map number survives `MolFromSmiles`, so a residue's anchor can be found by label
# after RDKit renumbers the molecule, rather than by index arithmetic over the emitted
# token stream — which breaks the moment a ring closure or a cap shifts the order.
# `clear_tags` strips the labels before the final string is written.
#
# The anchor is the residue's own nucleophile: lysine's N-zeta, cysteine's sulfur, and the
# hydroxyl oxygen of serine, threonine and tyrosine. H counts are written explicitly
# because a bracket atom in SMILES carries no implicit hydrogens.
TAGGED_SIDE_CHAIN = {
    "K": "CCCC[NH2:{tag}]",
    "C": "C[SH:{tag}]",
    "S": "C[OH:{tag}]",
    "T": "{beta}(C)[OH:{tag}]",
    "Y": "Cc1ccc([OH:{tag}])cc1",
}

# Residues whose side chain carries its own stereocentre. `_backbone` writes these for the
# natural series; a D peptide is the enantiomer, so every centre inverts, not just CA.
BETA_STEREO = {"I": "[C@@H](C)CC", "T": "[C@@H](C)O"}
BETA_STEREO_D = {"I": "[C@H](C)CC", "T": "[C@H](C)O"}


# --- non-standard residues -------------------------------------------------------------
# PEPlife2 writes a non-natural residue as a hyphen-delimited three-letter code inside an
# otherwise one-letter sequence: "H-Aib-QGTFTSDYSKYL...". Upper-casing and stripping
# punctuation - which is what the pipeline used to do - turns ONE residue into THREE
# standard ones (Aib -> A, I, B), lengthening the sequence and inventing residues the
# peptide does not contain. Those rows were therefore excluded outright for having no
# residue template.
#
# The ones below are added because the cohort actually names them and because each is a
# plain side-chain substitution: the backbone is untouched, so one template entry makes the
# residue both buildable and - since `_side_chain_signatures` is derived from this table -
# decodable, which is what the round-trip check requires. Counts are legend mentions in the
# rodent cohort's Chem-MOD plus occurrences in PEPlife2's own structure archive.
#
# Aib is the exception and the biggest single win (48 legend mentions, 90 in the archive).
# It is quaternary - two methyls on the alpha carbon, no hydrogen - so it alters the
# backbone and is special-cased in `_backbone` and in `_identify_residue`. Being
# quaternary and symmetric it is also achiral, so it takes no stereo token.
NONSTANDARD_SIDE_CHAIN = {
    "Orn": "CCCN",                    # ornithine, lysine one methylene shorter
    "Nle": "CCCC",                    # norleucine
    "Nva": "CCC",                     # norvaline
    "Abu": "CC",                      # 2-aminobutyric acid
    "Dab": "CCN",                     # 2,4-diaminobutyric acid
    "Dpr": "CN",                      # 2,3-diaminopropionic acid
    "hArg": "CCCCNC(=N)N",            # homoarginine
    "Cha": "CC1CCCCC1",               # cyclohexylalanine
    "Nal": "Cc1cccc2ccccc12",         # 3-(1-naphthyl)alanine
    "Tys": "Cc1ccc(OS(=O)(=O)O)cc1",  # O-sulfotyrosine
}

# Quaternary alpha carbons: no hydrogen, two side chains, therefore no stereocentre.
QUATERNARY = {"Aib": ("C", "C")}

# Backbone-altering residues. Sarcosine and the N-methylated residues (one methyl on the amide
# nitrogen) and N-terminal pyroglutamate (a lactam between the N-terminus and the side chain)
# each have their own construction in `_backbone` and decoding in `_identify_residue`, so the
# round-trip check covers them like any other residue.
SARCOSINE = "Sar"
PYROGLUTAMATE = "pGlu"
# Residues that may carry an N-methyl, keyed by the suffix of the token ("NMeVal" -> "Val").
N_METHYLATABLE = {
    "Leu": "L", "Val": "V", "Thr": "T", "Phe": "F", "Nva": "Nva",
}
N_METHYL_TOKENS = tuple("NMe" + name for name in N_METHYLATABLE)
# Spelling variants in `Sequence` that denote the same residue.
RESIDUE_ALIASES = {
    "meleu": "NMeLeu", "meval": "NMeVal", "methr": "NMeThr",
    "mephe": "NMePhe", "nmeleu": "NMeLeu", "nmeval": "NMeVal", "nmethr": "NMeThr",
    "nmephe": "NMePhe", "nmenva": "NMeNva", "pglu": "pGlu", "sar": "Sar", "1nal": "Nal",
}

# Named in the cohort but deliberately NOT supported, because each alters the backbone or
# the ring system in a way that has no construction here: Hyp / Aze / Oic (fused or altered
# proline ring), Iva (quaternary AND chiral). Records naming these stay excluded, and say so.
UNSUPPORTED_RESIDUES = ("Hyp", "Aze", "Oic", "Iva", "Ibu", "Npa", "Agp",
                        # Spelled in upper case in `Sequence`, so they would otherwise be
                        # expanded as a run of one-letter codes: AMHA -> A, M, H, A.
                        "Pff", "AMHA", "EtPro", "MePro", "DBip", "pIPhe", "Tz4", "Pal",
                        "Cbm", "Dpro", "Bip", "Cit", "Phg", "Tic", "Thi", "Nma")

# The 20 standard residues also appear spelled out in `Sequence` ("...-Leu-MePro-..."),
# where they are unambiguous and map straight to their one-letter code.
STANDARD_THREE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# Longest code first, so "hArg" is not matched as "Arg" and "NMeLeu" not as "MeLeu".
_RESIDUE_CODES = sorted(
    list(NONSTANDARD_SIDE_CHAIN) + list(QUATERNARY) + list(UNSUPPORTED_RESIDUES)
    + [SARCOSINE, PYROGLUTAMATE] + list(N_METHYL_TOKENS),
    key=len, reverse=True)
_CODE_LOOKUP = {c.lower(): c for c in _RESIDUE_CODES}
_CODE_LOOKUP.update({alias: code for alias, code in RESIDUE_ALIASES.items() if " " not in alias})


def residue_template(token: str) -> tuple[str, str] | None:
    """(kind, side chain) for a token: kind is 'standard', 'nonstandard' or 'quaternary'."""
    if token in NONSTANDARD_SIDE_CHAIN:
        return "nonstandard", NONSTANDARD_SIDE_CHAIN[token]
    if token in QUATERNARY:
        return "quaternary", ""
    if token == SARCOSINE:
        return "sarcosine", ""
    if token == PYROGLUTAMATE:
        return "pyroglutamate", ""
    if token in N_METHYL_TOKENS:
        base = N_METHYLATABLE[token[3:]]
        return "n-methyl", (SIDE_CHAIN[base] if len(base) == 1 else NONSTANDARD_SIDE_CHAIN[base])
    if len(token) == 1 and token in SIDE_CHAIN:
        return "standard", SIDE_CHAIN[token] or ""
    return None


# The letter ESM-2 is given for a residue it was never pretrained on. Recorded per position
# rather than applied silently: the molecular branch gets the true structure while the
# sequence branch gets an unknown-residue token, and that asymmetry is a measurable quality
# gradient, not a detail.
ESM_UNKNOWN = "X"


def esm_sequence(tokens: list[str]) -> tuple[str, list[int]]:
    """(sequence as ESM-2 will see it, 1-based positions replaced by X)."""
    out, blurred = [], []
    for i, token in enumerate(tokens, start=1):
        if len(token) == 1 and token in SIDE_CHAIN:
            out.append(token)
        else:
            out.append(ESM_UNKNOWN)
            blurred.append(i)
    return "".join(out), blurred


# --- per-residue chirality, from the column that still carries it ----------------------
# PEPlife2 ships TWO sequence columns and they are not the same data. `Sequence` is the
# original: it keeps hyphen-delimited non-natural residue codes AND marks a D residue by
# writing its letter in lower case ("rPKPfQwFwLL", "SCSlPQT"). `sequence` is a flattened
# derivative - upper-cased and punctuation-stripped - which turns "pGlu-LYENKPRRPYIL" into
# "PGLULYENKPRRPYIL", i.e. into exactly the prose that `looks_like_prose` then rejects, and
# throws the D/L pattern away. 403 of 1,300 rodent rows differ between the two.
#
# Reading the original column is what makes a mixed-chirality peptide representable. The
# `Chirality` column can only say "Mix", which is not a structure and cost 32 records their
# place in the strict set; the lower-case letters say WHICH residues are D, which is a
# structure. 62 rows carry that information (52 of them marked `Mix`).
#
# Delimited sequences are handled conservatively, because case cannot be trusted to mean
# D there: "pIPhe" is one residue (p-iodo-phenylalanine) whose letters p, I, P, h, e are
# all individually valid residue codes, so expanding it would silently invent five residues
# from one. A delimited sequence is therefore accepted only when every chunk is either a
# supported residue code or a single letter, and refused otherwise.

# In this cohort `O` is a residue code only where the source states it. Elsewhere an `O` is a
# fragment of a residue name transcribed into the sequence column, or a typo (a motilin
# chimera whose other copy reads Q at the same position), so a blanket O -> ornithine would
# invent structures. The sequences below are the ones whose papers define O as L-ornithine:
# apidaecin analogues Api88/Api137, oncocin analogue Onc72, and the C5a antagonist PMX53.
ORNITHINE_SEQUENCES = {
    "ONNRPVYIPRPRPPHPRL", "VDKPPYLPRPRPPROIYNO", "FOP-Cha-WR", "OP-(d-Cha)-WR",
}

# Multi-letter chunks that are NOT strings of one-letter residue codes even though every one
# of their letters is a valid code ("Dmt" would otherwise expand to Asp, Met, Thr; "desH" to
# D-Asp, D-Glu, D-Ser, His). They are refused with their own reason rather than expanded.
NON_RESIDUE_WORDS = {
    "des", "desh", "gu", "pr", "hoo", "hle", "sip", "dmt", "tmsala", "ahx", "cbc", "dcm",
    "dallo", "oeg", "nh2", "oh", "ile", "diacid", "agp", "pal", "dab", "dpr", "cha",
    # Two residue codes run together with no delimiter: "SarLys-RPPGFSPf" is
    # sarcosine-lysine, not Ser-Ala-Arg-Leu-Tyr-Ser. Expanding it letter by letter turned
    # one dipeptide into six invented residues and lengthened a 10-mer to a 14-mer, and the
    # round-trip check could not catch it because the wrong sequence decodes consistently
    # back to itself. Refused with a reason rather than silently mis-expanded; splitting
    # undelimited codes would need the tokenizer to guess where one ends.
    "sarlys",
}

_PAREN_CODE = re.compile(r"\(\s*([A-Za-z][A-Za-z0-9]*)\s*\)")


def parse_residues(raw: object) -> tuple[list[tuple[str, str]] | None, str]:
    """([(residue token, 'L' | 'D' | '')], problem) read from PEPlife2's `Sequence`.

    Case carries configuration for one-letter codes: upper case is the natural series,
    lower case its enantiomer. A multi-character code takes '' - unknown from the spelling
    alone - and the caller supplies it from the `Chirality` column.

    A sequence is accepted only when every character belongs to a residue. Anything that is
    prose, a repeat or polymer notation, a staple placeholder or a residue this module has
    no construction for is refused with the reason, never expanded letter by letter.
    """
    text = " ".join(str(raw).split())
    if not text or text.upper() == "NAN":
        return None, "no sequence recorded"

    # The residue name is the whole chunk, so a parenthesised code ("(Nva)", "(NMeNva)") is
    # unwrapped before splitting, and "N-MeVal" - an N-methyl marker, not asparagine - is
    # joined to its residue. A repeat count or any other bracket is not a residue and stays,
    # so the chunk check below refuses it.
    text = _PAREN_CODE.sub(r"\1", text)
    text = re.sub(r"\bN-\s*Me(?=[A-Za-z])", "NMe", text, flags=re.I)
    text = re.sub(r"\bp-Gln\b", "pGlu", text)
    ornithine = text.replace(" ", "") in ORNITHINE_SEQUENCES or text in ORNITHINE_SEQUENCES

    def one_letter(chunk: str) -> list[tuple[str, str]] | str:
        out_: list[tuple[str, str]] = []
        for c in chunk:
            if c.upper() == "O" and ornithine:
                out_.append(("Orn", ""))
            elif c.upper() in SIDE_CHAIN:
                out_.append((c.upper(), "D" if c.islower() else "L"))
            else:
                return c
        return out_

    if re.fullmatch(r"[A-Za-z]+", text):
        residues_or_bad = one_letter(text)
        if isinstance(residues_or_bad, str):
            return None, f"no residue template for {residues_or_bad.upper()}"
        return residues_or_bad, ""

    out: list[tuple[str, str]] = []
    for chunk in re.split(r"[-\s]+", text):
        chunk = chunk.strip()
        if not chunk:
            continue
        if chunk in STANDARD_THREE or chunk.capitalize() in STANDARD_THREE:
            letter = STANDARD_THREE.get(chunk, STANDARD_THREE.get(chunk.capitalize()))
            out.append((letter, "D" if chunk[0].islower() else "L"))
            continue
        code = _CODE_LOOKUP.get(chunk.lower())
        if code is not None:
            out.append((code, ""))
            continue
        if len(chunk) == 1:
            got = one_letter(chunk)
            if isinstance(got, str):
                return None, f"no residue template for {chunk.upper()}"
            out.extend(got)
            continue
        # A multi-letter chunk is expanded one residue per letter only when it is a plain
        # alphabetic run of one-letter codes and not a known non-residue word. Case then
        # marks configuration ("GViT" is Gly-Val-D-Ile-Thr; "HWSYlLRPG" has a D-Leu).
        if (chunk.isalpha() and chunk.lower() not in NON_RESIDUE_WORDS
                and all(c.upper() in SIDE_CHAIN or (c.upper() == "O" and ornithine)
                        for c in chunk)):
            out.extend(one_letter(chunk))
            continue
        return None, (f"sequence contains {chunk!r}, which is neither a residue code nor a "
                      f"run of one-letter codes")
    if not out:
        return None, "no residues found"
    return out, ""


def chirality_pattern(residues: list[tuple[str, str]]) -> str:
    """'L', 'D', 'mix' or '' for the configurations a parsed sequence actually specifies."""
    given = {c for t, c in residues if c and t in SIDE_CHAIN and t != "G"}
    if not given:
        return ""
    if given == {"L"}:
        return "L"
    if given == {"D"}:
        return "D"
    return "mix"


BLANK = re.compile(r"^(|free|n\.?a\.?|none|not reported|not mentioned|nan)$", re.I)

# Acyl caps expressible as a plain carbon chain. Patterns must match the WHOLE field: a
# description like "Lys1 linked with (palmitoyl-Glu-OH)" names a side-chain conjugate, not
# N-terminal palmitoylation, and drawing it as the latter would be a silent lie.
N_TERM_ACYL = [
    (r"(n-?(terminal )?)?acetyl(ation|ated)?", "*C(C)=O"),
    (r"(n-?(terminal )?)?formyl(ation|ated)?", "*C=O"),
    (r"propionyl(ation|ated)?", "*C(=O)CC"),
    (r"butyryl(ation|ated)?", "*C(=O)CCC"),
    (r"hexanoyl(ation|ated)?|caproyl(ation|ated)?", "*C(=O)CCCCC"),
    (r"octanoyl(ation|ated)?|caprylyl(ation|ated)?", "*C(=O)CCCCCCC"),
    (r"decanoyl(ation|ated)?", "*C(=O)CCCCCCCCC"),
    (r"lauroyl(ation|ated)?|dodecanoyl(ation|ated)?", "*C(=O)CCCCCCCCCCC"),
    (r"myristoyl(ation|ated)?|tetradecanoyl(ation|ated)?", "*C(=O)CCCCCCCCCCCCC"),
    (r"palmitoyl(ation|ated)?|hexadecanoyl(ation|ated)?", "*C(=O)CCCCCCCCCCCCCCC"),
    (r"stearoyl(ation|ated)?|octadecanoyl(ation|ated)?", "*C(=O)CCCCCCCCCCCCCCCCC"),
    (r"hydrocinnamate|hydrocinnamoyl(ation|ated)?|hc ?= ?hydrocinnamate.*",
     "*C(=O)CCc1ccccc1"),
    (r"alpha-?hydroxybutyric acid|2-?hydroxybutyr(yl|ic acid)", "*C(=O)C(O)CC"),
]

# The N-terminal guanidinylation of the apidaecin analogues Api88 and Api137, written in
# the cohort as the paper's own legend "Gu: N,N,N',N' tetramethylguanidino", is NOT here
# and is deliberately left undrawn. Every other entry above is an acyl bonded to the
# terminal nitrogen through a single bond, which is what `_attach` builds. Gu is not an
# acyl: the terminal amine becomes the third nitrogen of a tetramethylguanidine, whose only
# neutral tautomer doubles the bond to that nitrogen, and whose physiological form is a
# guanidinium cation. Drawing it would mean either a single bond RDKit rejects as
# pentavalent nitrogen, or the cohort's only charged structure next to an arginine side
# chain this module draws neutral. Those records stay `modification not drawn`.

# Terminal descriptions that are chemically real but whose partner cannot be drawn as a
# small molecule. Matching one is reported through `n_partner`/`c_partner` and the
# terminal-encoded flags - never by demoting `coverage`, and never by dropping the record.
MACROMOLECULAR = re.compile(
    r"\balbumin\b|\bhsa\b|\bmsa\b|\bfc\b|\bigg\b|\bgfp\b|repebody|his-?tag|strep-?tag|"
    r"transferrin|\belp\b|\bzipp\b|nanobody|antibody|\babd\b|\baff\b|darpin|nanofitin|"
    r"\bxten\b|\bpas#?\d|\bebp\b|scfv|scfc|\btrail\b|fn3|lectin|interferon|rhifn|"
    r"chondroitin|heparosan|polysarcosine|\bhpn\b|\bpcl\b|dextran|hyaluron|polyglutam|"
    r"poegma|oegma|\bpvp\b|polyvinyl|phosphorylcholine|chitosan|chitooligosaccharide|"
    r"\bpeg\b|peg\d|peg-|pegylat|polyethylene glycol|\bmpeg\b|\bdspe\b|"
    r"microsphere|nanoparticle|nanoghost|liposom|hydrogel|implant|\bplga\b|micelle|depot|"
    r"encapsulat|iron oxide|\bionp\b|gold particle|\bau\b|\bsenp\b|\bpamam\b|streptavidin|"
    r"glycosylat|\bctp\b|genetic fusion|fusion|fused|recombinantly linked|leader sequence",
    re.I,
)

C_TERM_AMIDE = re.compile(r"(c-?(terminal )?)?(amidation|amidated|amide|-?nh2|conh2)", re.I)
C_TERM_METHYL_ESTER = re.compile(r"(methyl ?ester(ification)?|-?ome)", re.I)
# A C-terminal N-ethyl amide, which PEPlife2 writes as "NHEt (Ethylamine)". It is the
# C terminus of the whole GnRH superagonist family - leuprolide, buserelin, histrelin are
# all -Pro-NHEt - so it is a cap the cohort needs rather than an exotic one. "Ethylamine"
# standing alone names the same cap: an amine in the C-TERMINAL column is the amide's
# substituent, not a free base.
C_TERM_ETHYLAMIDE = re.compile(
    r"(c-?(terminal )?)?(nh\s*-?\s*et|n-?ethyl\s*amide|ethyl\s*amide|ethyl\s*amine)"
    r"(\s*\(\s*ethyl\s*amine\s*\))?", re.I)

# "Amidation, Orn at position 23" carries a representable C-terminal amide plus a note
# naming a residue the sequence column already spells with 'O'. Honour the amide.
C_TERM_AMIDE_WITH_NOTE = re.compile(r"^\s*amidation\s*,.*|.*,\s*amidation\s*$", re.I)

# PEG with its size stated, as a molecular weight or as an explicit number of oxyethylene
# units. Both forms appear: "Pegylated with 20,000 mw PEG", "PEG(20k)", "40 kDa PEG",
# "mPEG5000", "(PEG)36", "linker PEG24". The size reaches the model as a number on
# `Termini.peg_kda`; no PEG is ever drawn into the SMILES.
PEG_MENTION = re.compile(r"\bm?peg\b|peg\d|\bpegylat\w*|polyethylene glycol", re.I)
PEG_KDA = re.compile(r"(\d[\d,\.]*)\s*(k(?:da)?\b|kilodalton)", re.I)
PEG_DALTON = re.compile(r"(\d[\d,\.]*)\s*(?:mw\b|da\b|dalton)|\bm?peg[\s\-(]*(\d{3,})", re.I)
PEG_UNITS = re.compile(r"\(?\bpeg\)?[\s\-]*\(?(\d{1,3})\b", re.I)
PEG_UNIT_MASS = 44.05    # one oxyethylene unit
PEG_MAX_UNITS = 1200     # ceiling on the parsed size, so a mis-read number cannot explode

# Nominal masses for partners PEPlife2 names but never gives a sequence or a size for.
# These are textbook values for the named entity, not measurements of the specific
# construct, and they exist so that "an Fc is attached" and "a His-tag is attached" are not
# the same number to a model. Ordered longest-pattern-first is unnecessary here because the
# lookup is a search over an ordered list, as with every other bucket in this module.
PARTNER_MASS_KDA = [
    (r"\bigg\b|antibody|\bmab\b", 150.0),
    (r"transferrin", 80.0),
    (r"\balbumin\b|\bhsa\b|\bmsa\b|\bbsa\b", 66.5),
    (r"\bxten\b", 60.0),
    (r"\belp\b|\bzipp\b", 55.0),
    (r"\bfc\b|scfc", 50.0),
    (r"\bgfp\b", 27.0),
    (r"scfv", 27.0),
    (r"repebody|darpin|nanofitin|fn3|\babd\b|nanobody", 15.0),
    (r"\bctp\b", 4.0),
    (r"his-?tag|strep-?tag", 1.0),
]

# "C6-C22", and also "Cys2-Cys10 And Cys4-Cys8" - the three-letter spelling is the one
# PEPlife2 uses whenever it writes the bond out in words rather than in shorthand.
DISULFIDE_PAIR = re.compile(r"\b(?:cys|c)\s*(\d+)\s*[-–]\s*(?:cys|c)\s*(\d+)", re.I)
HEAD_TO_TAIL = re.compile(r"n\s*-\s*c\s*terminal|head[ -]to[ -]tail|lactam.*n.*c-?term", re.I)

# A topology field that claims ONE ring and gives no positions: "Cyclic", "Cyclic peptide",
# "Cyclic (Cys-Cys Disulfide Bond)". "Bicyclic" is deliberately excluded - it claims two
# rings, which a single inferred bond cannot supply.
SINGLE_RING = re.compile(
    r"cyclic(?:\s*peptide)?(?:\s*\((?:\s*cys\s*[-–]\s*cys)?\s*(?:disulphide|disulfide)?"
    r"\s*bond(?:ed)?\s*\))?", re.I)
# Says the ring is a disulfide without saying between which residues, so a cysteine-free
# sequence carrying this note is contradictory rather than head-to-tail.
DISULFIDE_NOTE = re.compile(r"disulphide|disulfide|\bs-?s\b|cys\s*[-–]\s*cys", re.I)


@dataclass(frozen=True)
class Build:
    """Outcome of one reconstruction."""
    smiles: str | None
    status: str
    coverage: str          # 'complete' | 'partial' | 'none'
    n_partner: str = ""    # undrawable partner at the N terminus, verbatim
    c_partner: str = ""    # undrawable partner at the C terminus, verbatim
    # True when PEPlife2 reports a mixed L/D peptide without saying which residues are D.
    # The string is then drawn with every alpha carbon unspecified, so it is complete in
    # constitution but not in configuration. That distinction matters here rather than being
    # a technicality: D residues are one of the standard ways to blunt proteolysis, so a
    # model given only the flat structure cannot see one of the reasons the label is high.
    stereo_unspecified: bool = False
    # Whether each recorded structural axis is represented in the emitted string. None
    # means that no SMILES was emitted, so the axis was not evaluated.
    # How many side-chain conjugates from `peptide_mods` were bonded into the string.
    conjugates_encoded: int = 0
    # 'sequence case' when the D positions were read off the sequence itself, which is the
    # only way a mixed-chirality peptide becomes a structure; 'Chirality column' otherwise.
    stereo_source: str = ""
    n_terminal_encoded: bool | None = None
    c_terminal_encoded: bool | None = None
    topology_encoded: bool | None = None
    chirality_resolved: bool | None = None

    def __iter__(self):
        """Back-compatible with the old `smiles, status = build_smiles(...)` unpacking."""
        return iter((self.smiles, self.status))


def _normalise(value: object) -> str:
    return " ".join(str(value).split())


# Residue and modifier names that were transcribed into PEPlife2's sequence column as
# prose instead of one-letter codes. A hit means the row's "sequence" is a description,
# so its length, its ESM-2 embedding and its homology cluster are all meaningless — not
# just its SMILES. Matched only against short sequences: real protein sequences contain
# several of these as substrings by chance (an antibody heavy chain hits NME and ENK).
PROSE_TOKENS = re.compile(
    r"CYSTEIN|PHENYLALAN|TRYPTOP|THREONY|GLYCINE|PROPARGYL|TYROSY|LEUCY|VALY|ALANY|"
    r"LYSYL|ARGINY|GLUTAM|ASPART|PROLY|HYDROXY|AMIDED|MEPRO|"
    r"PSAR|SSOET|HOOFOP",
    re.I,
)
PROSE_MAX_LENGTH = 120


def looks_like_prose(sequence: str) -> bool:
    """True if the sequence column holds a description rather than a sequence."""
    seq = _normalise(sequence).upper()
    return len(seq) <= PROSE_MAX_LENGTH and bool(PROSE_TOKENS.search(seq))


def _centre(chirality: str) -> str:
    """CA stereo token. 'l' is the natural series, 'd' its enantiomer, else unspecified."""
    return {"l": "[C@@H]", "d": "[C@H]"}.get(chirality, "C")


def _backbone(sequence, chirality, head_to_tail: bool,
              disulfides: tuple[tuple[int, int], ...],
              tags: dict[int, int] | None = None) -> str | None:
    """Assemble the peptide as a SMILES string, or None if a residue has no template.

    `tags` maps a 1-based residue position to an atom map number, which is written onto
    that residue's side-chain anchor atom so a conjugate can be bonded there later.

    `chirality` is either one setting for the whole chain or one per residue, so a peptide
    whose D positions are known individually is drawn with each centre as reported rather
    than with every centre unspecified.
    """
    tokens = list(sequence) if isinstance(sequence, str) else list(sequence)
    if not tokens or any(residue_template(t) is None for t in tokens):
        return None
    tags = tags or {}
    if any(not 1 <= pos <= len(tokens) for pos in tags):
        return None
    if any(tokens[pos - 1] not in TAGGED_SIDE_CHAIN for pos in tags):
        return None

    # Ring-closure labels. The aromatic side-chain templates above already use digits 1
    # and 2, and a digit is only free once its own ring has closed — so a disulfide, whose
    # closure stays open across many residues, must not reuse them or it will pair with an
    # aromatic ring opening instead. Two-digit %nn labels keep the two sets disjoint.
    closures: dict[int, str] = {}
    ring = 10
    for i, j in disulfides:
        closures[i] = closures[j] = f"%{ring}"
        ring += 1
    head_ring = f"%{ring}" if head_to_tail else None
    if head_to_tail:
        ring += 1
    proline_ring = f"%{ring}"
    if ring > 99:
        return None

    per_residue = (list(chirality) if not isinstance(chirality, str)
                   else [chirality] * len(tokens))
    if len(per_residue) != len(tokens):
        return None
    parts = []
    for pos, (aa, chir_here) in enumerate(zip(tokens, per_residue), start=1):
        ca = _centre(chir_here)
        chirality = chir_here
        if aa == "G":
            unit = "NCC(=O)"
        elif aa == SARCOSINE:
            if pos in tags:
                return None
            unit = "N(C)CC(=O)"
        elif aa == PYROGLUTAMATE:
            # Only an N-terminal residue can be pyroglutamate: its ring closes on the very
            # nitrogen that would otherwise carry the preceding residue.
            if pos != 1 or head_ring is not None or pos in tags:
                return None
            unit = "O=C1CC" + ({"l": "[C@H]", "d": "[C@@H]"}.get(chir_here, "C")) + "(N1)C(=O)"
        elif aa in N_METHYL_TOKENS:
            if pos in tags:
                return None
            base = N_METHYLATABLE[aa[3:]]
            if len(base) == 1:
                side = SIDE_CHAIN[base]
                if base in BETA_STEREO:
                    side = (BETA_STEREO_D if chirality == "d" else BETA_STEREO)[base]
            else:
                side = NONSTANDARD_SIDE_CHAIN[base]
            unit = f"N(C){ca}({side})C(=O)"
        elif aa == "P":
            unit = f"N{proline_ring}{ca}(CCC{proline_ring})C(=O)"
        elif aa in QUATERNARY:
            # Two side chains and no alpha hydrogen, so no stereo token: Aib is achiral
            # whatever the record's Chirality column says.
            first, second = QUATERNARY[aa]
            unit = f"NC({first})({second})C(=O)"
        elif aa in NONSTANDARD_SIDE_CHAIN:
            side = NONSTANDARD_SIDE_CHAIN[aa]
            if pos in tags:
                return None  # no tagged template for a non-standard side chain
            unit = f"N{ca}({side})C(=O)"
        else:
            side = SIDE_CHAIN[aa]
            if aa in BETA_STEREO:
                side = (BETA_STEREO_D if chirality == "d" else BETA_STEREO)[aa]
            if aa == "C" and pos in closures:
                if pos in tags:
                    # A cysteine sulfur cannot both close a disulfide and carry a
                    # conjugate. The record says it does, so the record is wrong.
                    return None
                side = f"CS{closures[pos]}"
            elif pos in tags:
                template = TAGGED_SIDE_CHAIN[aa]
                beta = ("[C@H]" if chirality == "d" else "[C@@H]")
                side = template.format(tag=tags[pos], beta=beta)
            unit = f"N{ca}({side})C(=O)"
        if pos == 1 and head_ring is not None:
            unit = unit.replace("N", f"N{head_ring}", 1)
        parts.append(unit)

    return "".join(parts) + (head_ring if head_ring is not None else "O")


def _attach(mol: Chem.Mol, target_idx: int, fragment_smiles: str) -> Chem.Mol:
    """Bond `fragment_smiles` (whose '*' marks the attachment point) onto `target_idx`."""
    fragment = Chem.MolFromSmiles(fragment_smiles)
    offset = mol.GetNumAtoms()
    combined = Chem.RWMol(Chem.CombineMols(mol, fragment))

    dummy = next(a for a in combined.GetAtoms()
                 if a.GetIdx() >= offset and a.GetAtomicNum() == 0)
    anchor = dummy.GetNeighbors()[0].GetIdx()
    target = combined.GetAtomWithIdx(target_idx)
    target.SetNumExplicitHs(max(target.GetTotalNumHs() - 1, 0))
    target.SetNoImplicit(True)
    combined.AddBond(target_idx, anchor, Chem.BondType.SINGLE)
    combined.RemoveAtom(dummy.GetIdx())
    return combined.GetMol()


def tagged_atom(mol: Chem.Mol, tag: int) -> int | None:
    """Index of the atom carrying atom map number `tag`, or None."""
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == tag:
            return atom.GetIdx()
    return None


def clear_tags(mol: Chem.Mol) -> Chem.Mol:
    """Drop every atom map number, so the emitted SMILES carries no bookkeeping."""
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol


def _terminal_indices(mol: Chem.Mol) -> tuple[int | None, int | None, int | None]:
    """(alpha-amino N, C-terminal carbonyl C, C-terminal OH oxygen) of the assembled chain.

    The builder emits the chain N-to-C, so the first atom is the alpha-amino nitrogen and
    the last is the terminal hydroxyl oxygen — unless the chain was closed head-to-tail,
    in which case there is no free terminus to modify.
    """
    n_idx = 0 if mol.GetAtomWithIdx(0).GetSymbol() == "N" else None
    last = mol.GetNumAtoms() - 1
    oxt_idx = last if mol.GetAtomWithIdx(last).GetSymbol() == "O" else None
    c_idx = None
    if oxt_idx is not None:
        neighbours = mol.GetAtomWithIdx(oxt_idx).GetNeighbors()
        c_idx = neighbours[0].GetIdx() if len(neighbours) == 1 else None
    return n_idx, c_idx, oxt_idx


def _n_terminal_fragment(text: str) -> tuple[str | None, bool]:
    """Returns (fragment SMILES or None, drawable)."""
    if BLANK.fullmatch(text):
        return None, True
    for pattern, fragment in N_TERM_ACYL:
        if re.fullmatch(pattern, text, re.I):
            return fragment, True
    return None, False


def _c_terminal_kind(text: str) -> tuple[str | None, bool]:
    """Returns (None | 'amide' | 'ethylamide' | 'ester', drawable)."""
    if BLANK.fullmatch(text):
        return None, True
    # Before the plain amide, whose pattern would match "N-ethylamide" on "amide" alone
    # and drop the ethyl.
    if C_TERM_ETHYLAMIDE.fullmatch(text):
        return "ethylamide", True
    if C_TERM_AMIDE.fullmatch(text) or C_TERM_AMIDE_WITH_NOTE.fullmatch(text):
        return "amide", True
    if C_TERM_METHYL_ESTER.fullmatch(text):
        return "ester", True
    return None, False




def build_smiles(sequence: str, n_ter: object, c_ter: object,
                 linear_cyclic: object, chirality: object,
                 conjugates: "list[tuple[int, str]] | None" = None) -> Build:
    """Reconstruct the molecule. Always returns a Build; see the module docstring.

    `conjugates` is an optional list of (1-based residue position, fragment SMILES) pairs
    resolved by `peptide_mods`, each bonded onto that residue's side-chain anchor. They are
    attached after the termini, so a record can carry a cap and a side-chain conjugate at
    once.
    """
    raw = _normalise(sequence)

    # A missing sequence must be refused BEFORE the residue check, because pandas hands a
    # blank cell over as the float nan, whose string form upper-cases to "NAN" - three
    # perfectly valid residue codes. Without this, 17 rodent records with no sequence at
    # all were drawn as the tripeptide Asn-Ala-Asn and reported `complete`.
    if BLANK.fullmatch(raw) or raw.upper() == "NAN":
        return Build(None, "no sequence recorded", "none")

    if looks_like_prose(raw.upper()):
        return Build(None, "sequence column holds prose", "none")

    # Residues, not characters: a hyphen-delimited non-natural code is ONE residue, and a
    # lower-case letter is a D residue rather than a different one. Accepts a plain
    # upper-case sequence unchanged, so every existing caller behaves as before.
    residues, problem = parse_residues(raw)
    if residues is None:
        return Build(None, problem, "none")
    tokens = [t for t, _ in residues]
    unsupported = sorted({t for t in tokens if t in UNSUPPORTED_RESIDUES})
    if unsupported:
        return Build(None, f"no residue template for {', '.join(unsupported)}", "none")
    if any(t in AMBIGUOUS_AA for t in tokens if len(t) == 1):
        return Build(None, "ambiguous residue", "none")
    if any(residue_template(t) is None for t in tokens):
        bad = sorted({t for t in tokens if residue_template(t) is None})
        return Build(None, f"no residue template for {', '.join(bad)}", "none")
    sequence = tokens

    chir = _normalise(chirality).lower()
    if BLANK.fullmatch(chir):
        chir = "l"

    # Configuration comes from the SEQUENCE's own case where the sequence carries case, and
    # from the `Chirality` column otherwise. The two paths matter in opposite directions:
    #
    #  * `Sequence` writes a D residue in lower case, so a mixed peptide's D positions are
    #    known individually. That is a structure, and it is drawn.
    #  * `sequence`, the flattened column, is entirely upper case, so it carries no case
    #    signal at all and the column has to govern. Every existing caller passing that
    #    column therefore behaves exactly as before.
    #
    # 'mix' on the COLUMN still means unresolved - the pattern is reported but not located -
    # and is drawn with every centre unspecified rather than discarded.
    from_case = chirality_pattern(residues)
    if from_case == "mix":
        per_residue = [(c or "l").lower() for _, c in residues]
        stereo_known = True
        stereo_source = "sequence case"
    elif from_case == "D":
        per_residue = ["d"] * len(residues)
        stereo_known = True
        stereo_source = "sequence case"
    else:
        stereo_known = chir in ("l", "d")
        setting = chir if stereo_known else "mix"
        per_residue = [setting] * len(residues)
        stereo_source = "Chirality column"

    n_text, c_text = _normalise(n_ter), _normalise(c_ter)
    n_fragment, n_drawable = _n_terminal_fragment(n_text)
    c_kind, c_drawable = _c_terminal_kind(c_text)

    n_partner = "" if n_drawable else n_text
    c_partner = "" if c_drawable else c_text

    disulfides, head_to_tail, topo_ok = _cyclization(_normalise(linear_cyclic), sequence)
    topology_note = "" if topo_ok else _normalise(linear_cyclic)

    conjugates = list(conjugates or [])
    tags = {position: 900 + i for i, (position, _) in enumerate(conjugates)}
    if len(tags) != len(conjugates):
        return Build(None, "two conjugates at one residue", "none")

    smiles = _backbone(sequence, per_residue, head_to_tail, tuple(disulfides), tags)
    if smiles is None:
        return Build(None, "backbone not assembled", "none")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return Build(None, "backbone did not parse", "none")

    n_idx, c_idx, oxt_idx = _terminal_indices(mol)

    # `_attach` appends the fragment and deletes only the dummy, which sits above every
    # backbone atom, so the terminal indices found above stay valid afterwards.
    if n_fragment is not None and n_idx is not None:
        mol = _attach(mol, n_idx, n_fragment)

    if c_kind is not None and oxt_idx is not None:
        rw = Chem.RWMol(mol)
        if c_kind in ("amide", "ethylamide"):
            atom = rw.GetAtomWithIdx(oxt_idx)
            atom.SetAtomicNum(7)
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
            if c_kind == "ethylamide":
                # The hydroxyl oxygen has just become the amide nitrogen, so the ethyl
                # hangs off the same atom index. Spell the primary amide's two hydrogens
                # out first: `_attach` takes one hydrogen off the atom it bonds to, and on
                # an unsanitised molecule the implicit count it would read there is zero,
                # which leaves a nitrogen radical rather than a secondary amide.
                atom = rw.GetAtomWithIdx(oxt_idx)
                atom.SetNumExplicitHs(2)
                atom.SetNoImplicit(True)
            mol = rw.GetMol()
            if c_kind == "ethylamide":
                mol = _attach(mol, oxt_idx, "*CC")
        else:
            mol = _attach(rw.GetMol(), oxt_idx, "*C")

    for position, fragment in conjugates:
        anchor = tagged_atom(mol, tags[position])
        if anchor is None:
            return Build(None, "side-chain anchor lost", "none")
        mol = _attach(mol, anchor, fragment)

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return Build(None, "sanitization failed", "none")

    out = Chem.MolToSmiles(clear_tags(mol))

    # What is missing from the string, in the order it matters for clearance.
    missing = []
    if not n_drawable:
        missing.append("N-terminal partner")
    if not c_drawable:
        missing.append("C-terminal partner")
    if not topo_ok:
        missing.append("cyclization")
    if not stereo_known:
        missing.append("chirality")

    axis_flags = dict(conjugates_encoded=len(conjugates), stereo_source=stereo_source,
                      n_terminal_encoded=n_drawable,
                      c_terminal_encoded=c_drawable,
                      topology_encoded=topo_ok,
                      chirality_resolved=stereo_known)
    if not missing:
        return Build(out, "complete", "complete", **axis_flags)

    status = {
        "N-terminal partner": "N-terminal partner not drawable",
        "C-terminal partner": "C-terminal partner not drawable",
        "cyclization": "cyclization not drawable",
        "chirality": "chirality unspecified",
    }[missing[0]]
    if len(missing) > 1:
        status += f" (+{len(missing) - 1} more)"

    # Coverage is a statement about the PEPTIDE and nothing else: whether its own
    # constitution, topology and configuration are in the string. An attached partner does
    # NOT demote it.
    #
    # The earlier `peptide_core` value broke that rule, and `build_halflife_training_set`'s
    # docstring records why it was abandoned: driven by partners, it measured what RDKit
    # could draw rather than what the molecule is, and it decayed every time the builder
    # improved. Partner presence reaches the caller through `n_partner`/`c_partner` and the
    # `n_terminal_encoded`/`c_terminal_encoded` flags, and reaches the model through the
    # termini feature block, so nothing is lost by keeping it out of this field.
    #
    # Both peptide axes count. The published column is 928 complete / 158 partial / 231
    # none, and that 158 is exactly 83 unpositioned cyclizations plus 75 unresolved L/D
    # series - so an unresolved configuration demotes the string as an unpositioned
    # macrocycle does. `stereo_unspecified` still distinguishes which of the two it was.
    coverage = "partial" if any(m in ("cyclization", "chirality") for m in missing) \
        else "complete"
    return Build(out, status, coverage, n_partner, c_partner or topology_note,
                 stereo_unspecified=not stereo_known, **axis_flags)


def _peg_units(text: str) -> int | None:
    """Number of oxyethylene units named in a PEG description, or None if unstated.

    A number is read as a molecular weight when it is large or carries a mass unit, and as
    a literal unit count only in the `PEG24` / `(PEG)36` form, where the discrete PEGs used
    as spacers are written that way. The distinction matters: `PEG24` is 24 units (~1 kDa)
    while `PEG 20000` is 20 kDa, i.e. 454 units, and reading one as the other misstates the
    size by two orders of magnitude - which, for PEG, is the whole point of the annotation.
    """
    if not PEG_MENTION.search(text):
        return None
    compact = text.replace(",", "")
    m = PEG_KDA.search(compact)
    if m:
        return min(int(round(float(m.group(1)) * 1000 / PEG_UNIT_MASS)), PEG_MAX_UNITS)
    m = PEG_DALTON.search(compact)
    if m:
        value = float(next(g for g in m.groups() if g))
        return min(int(round(value / PEG_UNIT_MASS)), PEG_MAX_UNITS)
    m = PEG_UNITS.search(compact)
    if m:
        return min(int(m.group(1)), PEG_MAX_UNITS)
    return None


def _cyclization(text: str, sequence: str) -> tuple[list[tuple[int, int]], bool, bool]:
    """Returns (disulfide residue pairs, head-to-tail, representable).

    `sequence` is a token list here, so an index into it is a RESIDUE offset rather than a
    character offset - a Cys after an Aib sits at a different character offset than residue
    offset.

    `Linear/Cyclic` often says only "Cyclic", with the ring left to the reader. Two of
    those cases are not guesses but the only ring the molecule can close, and they are
    inferred rather than abandoned:

      * exactly two cysteines and no other detail - the disulfide between them is the only
        disulfide the peptide has the thiols for, and "Cyclic (Cys-Cys Disulfide Bond)"
        says so in words without giving the positions
      * no cysteine at all - with no thiol there is no disulfide, so a cyclic peptide is
        closed through its backbone. Both records this applies to say as much in their own
        sources: c[wLwReQeR] is written in head-to-tail bracket notation (PMID 28691076)
        and the YSNSG peptide is named a cyclopeptide in its title (PMID 31866565).

    Anything else underdetermined is left unrepresentable. A bare "Bicyclic" needs two
    rings and four cysteines, and a lone cysteine, or a pair of positions that are not both
    cysteines, means the numbering belongs to a parent the sequence column does not hold.
    """
    if BLANK.fullmatch(text) or re.fullmatch(r"linear", text, re.I):
        return [], False, True
    if HEAD_TO_TAIL.search(text):
        return [], True, True

    pairs = [(int(a), int(b)) for a, b in DISULFIDE_PAIR.findall(text)]
    if not pairs:
        cys = [i + 1 for i, token in enumerate(sequence) if token == "C"]
        if not SINGLE_RING.fullmatch(text.strip()):
            return [], False, False
        if len(cys) == 2:
            return [(cys[0], cys[1])], False, True
        if not cys and not DISULFIDE_NOTE.search(text):
            return [], True, True
        return [], False, False
    seen: set[int] = set()
    for i, j in pairs:
        if not (1 <= i <= len(sequence) and 1 <= j <= len(sequence)):
            return [], False, False
        if sequence[i - 1] != "C" or sequence[j - 1] != "C":
            return [], False, False
        # A cysteine sulfur forms at most one disulfide, so a residue appearing in two
        # pairs means the field was mis-transcribed; don't guess which pair was meant.
        if i in seen or j in seen or i == j:
            return [], False, False
        seen.update((i, j))
    return pairs, False, True


@dataclass(frozen=True)
class Termini:
    """What is attached to the peptide, as features rather than as atoms."""
    n_class: str           # bucket for the N-terminal modification
    c_class: str           # bucket for the C-terminal modification
    partner_class: str     # the heaviest thing attached across both termini
    peg_kda: float         # stated PEG mass, 0.0 when none is named or no size is given
    partner_kda: float     # nominal mass of a named protein/polymer partner, else 0.0
    n_text: str = ""       # verbatim source text, for audit
    c_text: str = ""


# Buckets for a single terminus, ordered so the chemistry that dominates clearance wins.
# Deliberately finer at the small end than PARTNER_CLASSES, which answers a different
# question ("what is the heaviest thing on this molecule").
TERMINAL_CLASSES = [
    ("particle / depot", r"microspher|nanoparticle|nanoghost|liposom|hydrogel|implant|"
                         r"\bplga\b|micelle|depot|encapsulat|\bsenp\b|\bpamam\b|iron oxide"),
    ("protein fusion", r"\balbumin\b|\bhsa\b|\bmsa\b|\bfc\b|\bigg\b|\bgfp\b|repebody|"
                       r"transferrin|\belp\b|\bzipp\b|nanobody|antibody|\babd\b|darpin|"
                       r"\bxten\b|scfv|scfc|fn3|interferon|\bctp\b|fusion|fused|"
                       r"leader sequence|extension peptide|his-?tag|strep-?tag|\bhdl\b"),
    ("polymer", r"chondroitin|heparosan|polysarcosine|\bhpn\b|dextran|hyaluron|polyglutam|"
                r"poegma|oegma|\bpvp\b|polyvinyl|chitosan|\bpas#?\d|phosphorylcholine"),
    ("pegylation", r"\bm?peg\b|peg\d|pegylat|polyethylene glycol"),
    ("lipidation", r"palmito|myrist|stearo|lauro|octanoy|decanoy|fatty acid|cholester|"
                   r"\bchol\b|lipid|\bc16\b|\bc18\b"),
    ("label / imaging", r"fitc|cy\d|alexa|rhodamine|\bdota\b|nota|do3a|fluoresc|gadoteridol|"
                        r"dbco|biotin|radiolab|radioiodin|tritium|iodinat|"
                        r"\b\d{1,3}\s*-?\s*[a-z]{1,2}\b[\s-]*label"),
    ("linker / small conjugate", r"maleimid|succinim|\bmpa\b|\bahx\b|citrate|gemcitabine|"
                                 r"diflunisal|linker"),
    ("acetylation", r"acetylat|^ac$"),
    ("amidation", r"amidat|amide|\bnh2\b|\bconh2\b"),
    ("esterification", r"esterif|methyl ester|\bome\b"),
    ("formylation", r"formyl"),
]


def _bucket_terminus(text: str) -> str:
    if BLANK.fullmatch(text.strip()):
        return "free"
    lowered = " ".join(text.lower().split())
    for label, pattern in TERMINAL_CLASSES:
        if re.search(pattern, lowered):
            return label
    return "other / unparsed"


def _partner_kda(text: str) -> float:
    for pattern, mass in PARTNER_MASS_KDA:
        if re.search(pattern, text, re.I):
            return mass
    return 0.0


def describe_termini(n_ter: object, c_ter: object) -> Termini:
    """Everything attached to the peptide, read from the two terminal-modification fields.

    Read from the SOURCE text and never from the emitted SMILES, so this is unaffected by
    whether a partner happens to be drawable - which was the defect in the previous design,
    where `smiles_coverage` silently doubled as a \"has a big partner\" flag and therefore
    measured the builder's capability rather than the molecule.
    """
    n_text, c_text = _normalise(n_ter), _normalise(c_ter)
    joined = f"{n_text} | {c_text}"
    n_class, c_class = _bucket_terminus(n_text), _bucket_terminus(c_text)

    # The heaviest thing attached anywhere, on the same ordering as a single terminus, so
    # "PEG on one end, Fc on the other" reports the Fc.
    partner = "none"
    for label, pattern in TERMINAL_CLASSES:
        if re.search(pattern, " ".join(joined.lower().split())):
            partner = label
            break

    peg_units = _peg_units(n_text) or _peg_units(c_text) or 0
    return Termini(
        n_class=n_class, c_class=c_class, partner_class=partner,
        peg_kda=round(peg_units * PEG_UNIT_MASS / 1000.0, 3),
        partner_kda=_partner_kda(joined), n_text=n_text, c_text=c_text,
    )

# --- reverse direction: SMILES -> amino acid sequence --------------------------------

# A carbon bonded directly to both an amine-type nitrogen and a carbonyl carbon is, in
# every structure this builder emits, an alpha carbon — no side chain, cap or partner
# fragment it ever attaches has that combination locally. This lets the backbone be found
# by substructure search rather than by retracing the original construction order, which
# is lost once RDKit renumbers atoms on `MolToSmiles`.
_ALPHA_CARBON = Chem.MolFromSmarts("[CX4]([#7])[CX3]=O")


def _side_chain_signatures() -> dict[str, str]:
    """Constitutional (stereo-blind) canonical SMILES of each side chain, as it reads
    once cut off its alpha carbon. Ile and Thr carry their own side-chain stereocentre in
    `SIDE_CHAIN`, but it never varies with the residue's own L/D setting, so constitution
    alone already distinguishes all 20 residues - no need to match configuration too, which
    would require reproducing the extracted fragment's exact bond-insertion order below.
    """
    sigs = {}
    for table in (SIDE_CHAIN, NONSTANDARD_SIDE_CHAIN):
        for aa, frag in table.items():
            if frag is None:  # Gly, Pro: identified structurally, not by a side-chain match
                continue
            frag_mol = Chem.MolFromSmiles("*" + frag)
            signature = Chem.MolToSmiles(frag_mol, isomericSmiles=False)
            # A clash would make one side chain decode as two different residues, which
            # would silently corrupt every round-trip check that relies on this table.
            if signature in sigs and sigs[signature] != aa:
                raise AssertionError(
                    f"side chains of {sigs[signature]} and {aa} are indistinguishable")
            sigs[signature] = aa
    return sigs


_SIDE_CHAIN_SIGNATURES = _side_chain_signatures()


def _carbonyl_neighbour(mol: Chem.Mol, atom_idx: int) -> int | None:
    """The neighbouring carbon bearing a C=O, if any."""
    for nb in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        if nb.GetSymbol() != "C":
            continue
        for bond in nb.GetBonds():
            other = bond.GetOtherAtom(nb)
            if other.GetSymbol() == "O" and bond.GetBondTypeAsDouble() == 2.0:
                return nb.GetIdx()
    return None


def _amine_neighbour(mol: Chem.Mol, ca_idx: int) -> int | None:
    for nb in mol.GetAtomWithIdx(ca_idx).GetNeighbors():
        if nb.GetSymbol() == "N":
            return nb.GetIdx()
    return None


_N_METHYL_BASE_TO_NAME = {base: name for name, base in N_METHYLATABLE.items()}


def _is_pyroglutamate(mol: Chem.Mol, ca_idx: int, amine_idx: int) -> bool:
    """A five-membered N,C,C,C,C ring through CA and its amine with a ring carbonyl: the
    lactam of N-terminal pyroglutamate. Proline has the same ring without the carbonyl."""
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 5 or ca_idx not in ring or amine_idx not in ring:
            continue
        symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
        if symbols.count("N") != 1 or symbols.count("C") != 4:
            continue
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "C" and any(
                    b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetSymbol() == "O"
                    for b in atom.GetBonds()):
                return True
    return False


def _identify_residue(mol: Chem.Mol, ca_idx: int, amine_idx: int,
                      carbonylc_idx: int) -> str | None:
    """Residue token for the alpha carbon at `ca_idx`, or None if unrecognised.

    Backbone variants are read first: an N-terminal pyroglutamate lactam, and a methyl on
    the amide nitrogen (sarcosine, or N-methyl + the residue the side chain names).
    """
    if _is_pyroglutamate(mol, ca_idx, amine_idx):
        return PYROGLUTAMATE
    amine = mol.GetAtomWithIdx(amine_idx)
    methyl = [n for n in amine.GetNeighbors() if n.GetIdx() != ca_idx
              and n.GetSymbol() == "C" and n.GetDegree() == 1]
    base = _identify_core(mol, ca_idx, amine_idx, carbonylc_idx)
    if not methyl:
        return base
    if base == "G":
        return SARCOSINE
    name = _N_METHYL_BASE_TO_NAME.get(base)
    return f"NMe{name}" if name else None


def _identify_core(mol: Chem.Mol, ca_idx: int, amine_idx: int,
                   carbonylc_idx: int) -> str | None:
    """One-letter code for the alpha carbon at `ca_idx`, or None if unrecognised."""
    others = [n.GetIdx() for n in mol.GetAtomWithIdx(ca_idx).GetNeighbors()
             if n.GetIdx() not in (amine_idx, carbonylc_idx)]

    # Proline: a genuine pyrrolidine ring (N-CA-CB-CG-CD, size 5, all C but the one N)
    # containing both CA and its amine. A disulfide loop or a head-to-tail macrocycle can
    # also put CA and amine in a shared ring, but never one of this size and composition.
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 5 or ca_idx not in ring or amine_idx not in ring:
            continue
        symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
        if symbols.count("N") == 1 and symbols.count("C") == 4:
            return "P"

    if not others:
        return "G"
    if len(others) > 1:
        # A quaternary alpha carbon. Aib is the only one built, and its two side chains are
        # both methyls; anything else with two side chains is not a template this module
        # emits, so it stays unrecognised rather than being guessed at.
        if len(others) == 2 and all(
                mol.GetAtomWithIdx(i).GetSymbol() == "C"
                and mol.GetAtomWithIdx(i).GetDegree() == 1 for i in others):
            return "Aib"
        return None
    side0 = others[0]

    # Cystine: CB bonded to an S that itself bonds onward to a second S. Handled before
    # the generic fragment match below, because cutting across that second S-S bond would
    # add a dummy atom the (single-cut) side-chain signatures don't have.
    for nb in mol.GetAtomWithIdx(side0).GetNeighbors():
        if nb.GetSymbol() == "S":
            if any(nb2.GetSymbol() == "S" and nb2.GetIdx() != side0
                  for nb2 in nb.GetNeighbors()):
                return "C"

    # Every other side chain is a closed fragment: walk out from side0, away from CA.
    # The disulfide-S-S guard below never fires here (only Cys ever has one), but costs
    # nothing to keep for any future residue that might.
    visited = {ca_idx}
    atoms: set[int] = set()
    bonds: set[int] = set()
    frontier = [side0]
    while frontier:
        nxt = []
        for idx in frontier:
            if idx in visited:
                continue
            visited.add(idx)
            atoms.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                other = bond.GetOtherAtom(atom)
                if other.GetIdx() in visited:
                    continue
                if atom.GetSymbol() == "S" and other.GetSymbol() == "S":
                    continue
                bonds.add(bond.GetIdx())
                nxt.append(other.GetIdx())
        frontier = nxt

    # Copy just the (handful of) fragment atoms into a fresh, small molecule rather than
    # fragmenting the whole peptide: `Chem.FragmentOnBonds`/`GetMolFrags` on the full mol
    # scale with the WHOLE molecule's size, called once per residue, which is minutes slow
    # on a several-hundred-residue engineered fusion partner (XTEN, ELP repeats).
    fragment = Chem.RWMol()
    dummy_idx = fragment.AddAtom(Chem.Atom(0))
    old_to_new = {}
    for idx in atoms:
        old = mol.GetAtomWithIdx(idx)
        new = Chem.Atom(old.GetAtomicNum())
        new.SetIsAromatic(old.GetIsAromatic())
        new.SetFormalCharge(old.GetFormalCharge())
        new.SetNoImplicit(old.GetNoImplicit())
        new.SetNumExplicitHs(old.GetNumExplicitHs())
        old_to_new[idx] = fragment.AddAtom(new)
    fragment.AddBond(dummy_idx, old_to_new[side0], Chem.BondType.SINGLE)
    for bond_idx in bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        fragment.AddBond(old_to_new[bond.GetBeginAtomIdx()], old_to_new[bond.GetEndAtomIdx()],
                         bond.GetBondType())

    try:
        Chem.SanitizeMol(fragment)
        signature = Chem.MolToSmiles(fragment, isomericSmiles=False)
    except Exception:
        return None
    return _SIDE_CHAIN_SIGNATURES.get(signature)


# CIP code an alpha carbon takes in the natural series, established by building all 20
# residues in each series and reading RDKit's assignment back (see `selfcheck`). Every
# residue is S except cysteine, whose sulfur outranks its carbonyl and flips the
# descriptor; glycine has no centre at all. The D series is the mirror of this.
L_SERIES_CIP = {"C": "R"}
L_SERIES_DEFAULT = "S"


def _configuration(token: str, cip: str) -> str:
    """'L', 'D' or '' for a residue whose alpha carbon carries CIP code `cip`."""
    if not cip or token in QUATERNARY or token in ("G", SARCOSINE):
        return ""
    base = N_METHYLATABLE[token[3:]] if token in N_METHYL_TOKENS else token
    natural = L_SERIES_CIP.get(base, L_SERIES_DEFAULT)
    return "L" if cip == natural else "D"


def reconstruct_residues(smiles: str) -> tuple[list[tuple[str, str]] | None, bool]:
    """Read (residue token, configuration) pairs back out of a generated SMILES string.

    This is the reverse half of the interconversion check, and unlike the constitutional
    walk it also recovers CONFIGURATION, so a peptide whose D positions were read from the
    sequence's own case is verified to have been drawn with those positions D - rather than
    the round trip passing on constitution alone while the stereochemistry goes unchecked.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, False
    order, neighbours, cyclic = _walk_backbone(mol)
    if order is None:
        return None, bool(cyclic)
    # Stereochemistry is perceived once, on the whole molecule, only in this function;
    # the constitutional walk deliberately avoids it because legacy CIP assignment is
    # pathologically slow on peptides with many centres.
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    out: list[tuple[str, str]] = []
    for ca in order:
        token = _identify_residue(mol, ca, *neighbours[ca])
        if token is None:
            return None, cyclic
        cip = mol.GetAtomWithIdx(ca).GetPropsAsDict().get("_CIPCode", "")
        out.append((token, _configuration(token, cip)))
    return out, cyclic


def render_residues(residues: list[tuple[str, str]]) -> str:
    """The canonical `Sequence` spelling of parsed residues: lower case marks a D residue."""
    out: list[str] = []
    for token, config in residues:
        text = token.lower() if (config == "D" and len(token) == 1) else token
        if len(token) == 1:
            out.append(text)
        else:
            if out and not out[-1].endswith("-"):
                out.append("-")
            out.append(text + "-")
    return "".join(out).strip("-")


def join_tokens(tokens: list[str]) -> str:
    """Render residue tokens the way PEPlife2 writes them.

    A run of one-letter codes is concatenated; a multi-character code is set off by hyphens,
    so "H-Aib-QGTF" round-trips to itself rather than to an ambiguous "HAibQGTF" in which
    the Aib could be read as three residues.
    """
    out: list[str] = []
    for token in tokens:
        if len(token) == 1:
            out.append(token)
        else:
            if out and not out[-1].endswith("-"):
                out.append("-")
            out.append(token + "-")
    return "".join(out).strip("-") if len(tokens) != 1 else tokens[0]


def reconstruct_sequence(smiles: str) -> tuple[str | None, bool]:
    """Read the residue sequence back out of a generated SMILES string.

    Returns (sequence, is_cyclic). Sequence is None if the string doesn't parse, its
    backbone can't be walked as one single chain (or one ring, for a head-to-tail
    cyclic peptide), or any residue's side chain doesn't match a known template. For a
    cyclic peptide, `sequence` is ONE valid rotation of the true sequence: a macrocycle
    has no intrinsic starting residue, so equality must be checked up to rotation, not
    exactly, by the caller.
    """
    residues, cyclic = reconstruct_residues(smiles)
    if residues is None:
        return None, cyclic
    return join_tokens([t for t, _ in residues]), cyclic


def _walk_backbone(mol):
    """(alpha carbons in chain order, their amine/carbonyl neighbours, is_cyclic)."""
    # No whole-molecule stereochemistry perception here: it isn't needed for the
    # constitutional walk below, and legacy CIP assignment is pathologically slow on
    # peptides with many stereocentres (each `_identify_residue` call perceives
    # stereochemistry only on its own few-atom fragment, which stays cheap).
    # `maxMatches` defaults to 1000 and silently truncates anything past it - a real risk
    # here, since some engineered fusion partners (XTEN, ELP repeats) run past 1000 residues.
    alpha_carbons = sorted({m[0] for m in mol.GetSubstructMatches(
        _ALPHA_CARBON, maxMatches=mol.GetNumAtoms())})
    if not alpha_carbons:
        return None, None, False

    neighbours = {ca: (_amine_neighbour(mol, ca), _carbonyl_neighbour(mol, ca))
                 for ca in alpha_carbons}
    amine_to_ca = {amine: ca for ca, (amine, _) in neighbours.items()}

    # ca -> the next residue's CA, found via its own carbonyl carbon's OTHER nitrogen
    # neighbour (the next peptide bond) - not the amide/ester oxygen a C-terminal cap
    # may have replaced, and not present at all past a free or capped C-terminus.
    next_ca: dict[int, int] = {}
    has_predecessor: set[int] = set()
    for ca, (amine, carbonylc) in neighbours.items():
        if carbonylc is None:
            continue
        for nb in mol.GetAtomWithIdx(carbonylc).GetNeighbors():
            if nb.GetSymbol() == "N" and nb.GetIdx() != amine:
                downstream = amine_to_ca.get(nb.GetIdx())
                if downstream is not None:
                    next_ca[ca] = downstream
                    has_predecessor.add(downstream)

    starts = [ca for ca in alpha_carbons if ca not in has_predecessor]
    limit = len(alpha_carbons) + 1
    if len(starts) == 1:
        order, cyclic = [starts[0]], False
        while order[-1] in next_ca:
            order.append(next_ca[order[-1]])
            if len(order) > limit:
                return None, None, False
    elif not starts:
        # every CA has a predecessor: head-to-tail cyclic, walk once around the ring
        # starting anywhere (the caller compares up to rotation)
        order, cyclic = [alpha_carbons[0]], True
        while next_ca.get(order[-1]) != order[0]:
            if order[-1] not in next_ca or len(order) > limit:
                return None, None, False
            order.append(next_ca[order[-1]])
    else:
        return None, None, False  # more than one chain: not a shape this builder emits

    return order, neighbours, cyclic



def selfcheck() -> list[str]:
    """Residue templates must reproduce RDKit's own peptide builder atom-for-atom."""
    problems = []
    for seq in ["A", "G", "P", "AG", "AP", "PA", "GG", "TTT", "IIT",
                "ACDEFGHIKLMNPQRSTVWY", "SQDSAFRIQERLRNSKMAHSSSCFGQKIDRIGAVSRLGCDGLRLF",
                "P" * 12]:                      # more prolines than there are ring digits
        mine = _backbone(seq, "l", False, ())
        ref = Chem.MolFromSequence(seq)
        if mine is None or ref is None:
            problems.append(f"{seq[:24]}: not built")
            continue
        a, b = Chem.MolToSmiles(Chem.MolFromSmiles(mine)), Chem.MolToSmiles(ref)
        if a != b:
            problems.append(f"{seq[:24]}:\n  ours  {a}\n  rdkit {b}")

    # A disulfide closure stays open across the residues it spans. The aromatic side-chain
    # templates open their own rings in between, so the two label sets must not overlap —
    # if they do, the SMILES still parses and keeps its formula but wires the rings wrongly.
    for seq, pairs in [("CFC", ((1, 3),)), ("CWYHFC", ((1, 6),)),
                       ("CFCWC", ((1, 5),)), ("CFCWCAC", ((1, 3), (5, 7)))]:
        built = _backbone(seq, "l", False, pairs)
        mol = Chem.MolFromSmiles(built) if built else None
        if mol is None:
            problems.append(f"{seq} {pairs}: did not parse")
            continue
        aromatic_rings = sum(1 for r in mol.GetRingInfo().AtomRings()
                             if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r))
        expected = sum(seq.count(aa) * (2 if aa == "W" else 1) for aa in "FWYH")
        if aromatic_rings != expected:
            problems.append(f"{seq} {pairs}: {aromatic_rings} aromatic rings, expected "
                            f"{expected} — ring labels collided\n  {built}")
        for i, j in pairs:
            sg = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "S"]
            if not any(mol.GetBondBetweenAtoms(x, y) for x in sg for y in sg if x != y):
                problems.append(f"{seq} {pairs}: no S-S bond formed")
                break

    # Round trip: every sequence/modification combination above, run forward through
    # `build_smiles` and back through `reconstruct_sequence`, must return exactly (or, for
    # a cyclic peptide, up to rotation of) the sequence that went in.
    cases = [
        dict(sequence="ACDEFGHIKLMNPQRSTVWY", n_ter="free", c_ter="free",
            linear_cyclic="linear", chirality="l"),
        dict(sequence="GAG", n_ter="free", c_ter="free", linear_cyclic="linear", chirality="l"),
        dict(sequence="PAP", n_ter="free", c_ter="free", linear_cyclic="linear", chirality="l"),
        dict(sequence="CFC", n_ter="free", c_ter="free", linear_cyclic="C1-C3", chirality="l"),
        dict(sequence="CWYHFC", n_ter="free", c_ter="free", linear_cyclic="C1-C6", chirality="l"),
        dict(sequence="AAA", n_ter="acetylation", c_ter="amidation",
            linear_cyclic="linear", chirality="l"),
        dict(sequence="AAA", n_ter="free", c_ter="free",
            linear_cyclic="head-to-tail", chirality="l"),
        dict(sequence="KKKKK", n_ter="palmitoylation", c_ter="methyl ester",
            linear_cyclic="linear", chirality="l"),
        dict(sequence="AAA", n_ter="free", c_ter="free", linear_cyclic="linear", chirality="d"),
        dict(sequence="AAA", n_ter="free", c_ter="free", linear_cyclic="linear", chirality="mix"),
        dict(sequence="IIT", n_ter="free", c_ter="free", linear_cyclic="linear", chirality="l"),
        dict(sequence="SQDSAFRIQERLRNSKMAHSSSCFGQKIDRIGAVSRLGCDGLRLF", n_ter="free",
            c_ter="free", linear_cyclic="linear", chirality="l"),
    ]
    # Backbone-altering residues go through the same forward-and-back check as every other.
    for raw, chir in [("Sar-GViT-Nva-IRP", "l"), ("pGlu-HWSYlLRPG", "l"),
                      ("NMeVal-A-NMeLeu-NMeThr-NMePhe-G", "l"), ("Sar-Sar", "l"),
                      ("pGlu-RPRLSHKGPMPF", "l"), ("Sar-A", "d")]:
        built = build_smiles(raw, "free", "free", "linear", chir)
        want, _ = parse_residues(raw)
        if built.smiles is None or want is None:
            problems.append(f"round trip {raw}: no SMILES ({built.status})")
            continue
        got, _ = reconstruct_residues(built.smiles)
        def norm(items):
            return [(t, c if t not in ("G", SARCOSINE, "Aib") and c else "") for t, c in items]
        if got is None or [t for t, _ in got] != [t for t, _ in want]:
            problems.append(f"round trip {raw}: decoded {got}")
        elif chir == "l" and any(c != w and w for (t, c), (_, w) in zip(norm(got), norm(want))):
            problems.append(f"round trip {raw}: configuration {got} vs {want}")
    for case in cases:
        built = build_smiles(case["sequence"], case["n_ter"], case["c_ter"],
                             case["linear_cyclic"], case["chirality"])
        if built.smiles is None:
            problems.append(f"round trip {case['sequence']}: builder gave no SMILES")
            continue
        recon, cyclic = reconstruct_sequence(built.smiles)
        if recon is None:
            problems.append(f"round trip {case['sequence']} {case['linear_cyclic']}: "
                            f"reconstruction failed\n  {built.smiles}")
            continue
        ok = (recon in case["sequence"] * 2 and len(recon) == len(case["sequence"])) \
            if cyclic else recon == case["sequence"]
        if not ok:
            problems.append(f"round trip {case['sequence']} {case['linear_cyclic']}: "
                            f"got {recon!r} (cyclic={cyclic})\n  {built.smiles}")
    return problems


if __name__ == "__main__":
    issues = selfcheck()
    print("\n".join(issues) if issues else "selfcheck: residue templates match RDKit")
