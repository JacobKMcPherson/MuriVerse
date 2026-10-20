"""Classify PEPlife2 modification free text, and turn the drawable part into SMILES.

`peptide_smiles.py` reconstructs the peptide itself: residues, topology, chirality and the
two termini. It reads `N-ter MOD`, `C-ter MOD` and `Linear/Cyclic` and nothing else, which
leaves `Chem-MOD` — 226 distinct strings over 483 of the 983 rodent peptide-scale records —
entirely unexamined. That gap is why `smiles_coverage == 'complete'` has been overstating
completeness: a record reading "Lys24 conjugated with (-(gE)3-C18)" got a bare-peptide
string and was called complete, and a record reading "Glycosylation at 84th position" got
the same string and the same verdict.

This module closes that gap by forcing every clause of every modification field into one of
three verdicts, with no silent fourth option:

  inert     the clause adds no atoms to the molecule, so the peptide-only string is already
            the whole structure. Two sub-cases, both real rather than convenient:
            * sequence-descriptive text — legends ("Aib = a-aminoisobutyric acid"),
              substitution notes ("V24L modifications", "Substituition of G for A at
              position 2"), chirality notes ("All D-amino acids except Gly"). These
              describe the sequence column rather than adding to it, and the strict
              sequence gate has already rejected the rows whose sequence cannot carry them.
            * isotopic substitution — 3H, 14C, 13C, 15N. An isotope changes mass, not
              constitution or connectivity.
  drawable  the clause names a real fragment AND a resolvable attachment site, so it is
            built and bonded on. Returns the fragment SMILES and the site.
  blocking  the clause is chemically real but cannot be drawn truthfully — a polydisperse
            kDa-scale PEG, a glycan, a protein partner, a proprietary dye, a drug payload,
            or a structure the source only gives as "structure given in paper". A record
            with any blocking clause is excluded from the strict set, by design.

Radiotracers are the one judgement call worth stating outright, because they decide 87 of
the 121 label-only records. A radioiodination ("I125 labeling", "Radioiodinated peptide at
Tyr") IS constitutional — it hangs an iodine off a tyrosyl ring, unlike a 3H or 14C
substitution. It is nonetheless classed `inert` here, for a reason about the dosed material
rather than about drawing difficulty: tracer radioiodination runs at specific activities far
below one iodine per peptide, so the overwhelming majority of the administered molecules
carry no iodine at all. Drawing the iodine would describe the minority species and
misstate the bulk of what was injected. The same argument covers a radiometal, which in
addition is a coordination complex rather than a covalent adduct. Every such record is
flagged in `tracer_label` so the choice can be inspected, and reversed, without re-deriving
it.

A radiometal named WITHOUT its chelator ("Cu-labelled", "Mn and 89Zr labeling") is blocking
instead: the chelator is a real covalent adduct of a few hundred daltons, and its absence
from the text means the structure is simply unknown.

Conjugate recipes are assembled from the Lys/Cys/Ser anchor outward. Spacer order inside a
recipe is NOT inferred when two different spacer types appear: PEPlife2 writes
"N(eps)hexadecanedioyl-gGlu-2xOEG" distal-first and "(ADO-linker-Glu-C18)" proximal-first,
and nothing in the text distinguishes them, so a multi-spacer recipe is blocking rather
than a guess. Single-spacer recipes are unambiguous — a fatty acyl has no amine, so it can
only be the distal cap — and are built.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from rdkit import Chem

from peptide_smiles import (N_METHYL_TOKENS, NONSTANDARD_SIDE_CHAIN, PYROGLUTAMATE,
                            QUATERNARY, SARCOSINE)

# Residues the SMILES builder now has a template for. A clause that merely names one of
# these is a LEGEND explaining a code the sequence already uses, not an unrepresentable
# modification, so it must not block the record: "Aib = a-aminoisobutyric acid" alongside a
# sequence reading "H-Aib-QGTF" adds no atoms the structure does not already have.
#
# Sarcosine and pyroglutamate belong here with the rest. They alter the backbone, so they
# are constructed in `peptide_smiles._backbone` rather than through the side-chain table,
# but they are templated and decodable all the same, and a legend naming one is as inert as
# a legend naming Aib.
SUPPORTED_RESIDUES = (set(NONSTANDARD_SIDE_CHAIN) | set(QUATERNARY)
                      | {SARCOSINE, PYROGLUTAMATE} | set(N_METHYL_TOKENS))
_SUPPORTED_LOWER = {r.lower() for r in SUPPORTED_RESIDUES}

# A legend spells the residue out - "O denotes L-ornithine", "pGlu = Pyro-glutamate" - and
# the code it defines ("O", "pGlu") is not always a word the clause repeats. Matching the
# NAME back to the code is what makes those legends recognisable; matching only the code
# left the broad non-standard-residue blocking rule to reject records whose residue the
# builder templates and whose sequence already carries it.
RESIDUE_NAME_SPELLINGS = {
    "Orn": (r"ornithine",),
    "Nle": (r"norleucine",),
    "Nva": (r"norvaline",),
    "Abu": (r"aminobutyric acid", r"aminobutanoic acid"),
    "Dab": (r"diaminobutyric acid",),
    "Dpr": (r"diaminopropionic acid",),
    "hArg": (r"homoarginine", r"homo-arginine"),
    "Cha": (r"cyclohexylalanine",),
    "Nal": (r"naphthylalanine", r"naphthyl-alanine"),
    "Tys": (r"sulfotyrosine", r"sulpho-?tyrosine"),
    "Aib": (r"aminoisobutyric acid", r"aminoisobutyrate"),
    SARCOSINE: (r"sarcosine", r"n-methyl ?glycine"),
    PYROGLUTAMATE: (r"pyroglutamate", r"pyroglutamic acid", r"pyro-?glutam\w*",
                    r"5-oxo-?prolin\w*", r"5-oxoproline"),
    # An N-methylated residue is templated, so a legend naming one is inert - but ONLY
    # where the sequence carries that exact token, which the `named <= present` test below
    # enforces. The distinction matters: "Methylation of alpha-C of Pro" is alpha-methyl
    # proline, a quaternary centre with no construction here, and its sequence carries no
    # NMePro token, so it keeps blocking.
    "NMePhe": (r"n-?methyl[- ]?phenylalanine", r"methylation at phenylalanine"),
    "NMeLeu": (r"n-?methyl[- ]?leucine", r"methylation at leucine"),
    "NMeVal": (r"n-?methyl[- ]?valine", r"methylation at valine"),
    "NMeThr": (r"n-?methyl[- ]?threonine", r"methylation at threonine"),
    "NMeNva": (r"n-?methyl[- ]?norvaline",),
}
_RESIDUE_NAME_RX = tuple(
    (code, re.compile(r"\b(?:" + "|".join(spellings) + r")\b"))
    for code, spellings in RESIDUE_NAME_SPELLINGS.items())


def named_templated_residues(text: str) -> set[str]:
    """Lower-cased codes of templated residues this clause names, by code or by name."""
    found = {w for w in re.findall(r"[a-z]+", text) if w in _SUPPORTED_LOWER}
    for code, rx in _RESIDUE_NAME_RX:
        if rx.search(text):
            found.add(code.lower())
    return found

# --- fragment library -----------------------------------------------------------------
# Every fragment is written with '*' marking the atom that bonds to the peptide. Fragments
# attaching to a nucleophile (Lys N-zeta, Cys S, Ser/Thr/Tyr O) lead with their acyl
# carbon, so '*' sits on that carbon and the bond formed is an amide, thioester or ester
# as the anchor dictates.

# Terminal acyl caps: a fatty acid or fatty diacid, which has no amine and therefore can
# only ever be the distal end of a conjugate.
ACYL = {
    "c8": "*C(=O)CCCCCCC",
    "c10": "*C(=O)CCCCCCCCC",
    "c12": "*C(=O)CCCCCCCCCCC",
    "c14": "*C(=O)CCCCCCCCCCCCC",
    "c16": "*C(=O)CCCCCCCCCCCCCCC",
    "c18": "*C(=O)CCCCCCCCCCCCCCCCC",
    "c18:1": "*C(=O)CCCCCCC/C=C/CCCCCCCC",
    # Diacids keep a free distal carboxylate, which is the point of them: it is what binds
    # albumin. C-number counts the whole chain including both carbonyls.
    "c14-diacid": "*C(=O)CCCCCCCCCCCC(=O)O",
    "c16-diacid": "*C(=O)CCCCCCCCCCCCCC(=O)O",
    "c18-diacid": "*C(=O)CCCCCCCCCCCCCCCC(=O)O",
}

ACYL_SYNONYM = {
    "octanoyl": "c8", "octanoic": "c8", "caprylyl": "c8",
    "decanoyl": "c10", "lauroyl": "c12", "dodecanoyl": "c12",
    "myristoyl": "c14", "myristic": "c14", "tetradecanoyl": "c14",
    "palmitoyl": "c16", "palmitic": "c16", "palm": "c16", "hexadecanoyl": "c16",
    "stearoyl": "c18", "stearic": "c18", "octadecanoyl": "c18",
    "oleic": "c18:1", "ole": "c18:1",
    "tetradecanedioyl": "c14-diacid", "tetradecanedioic": "c14-diacid",
    "hexadecanedioyl": "c16-diacid", "hexadecanedioic": "c16-diacid",
    "octadecanedioyl": "c18-diacid", "octadecanedioic": "c18-diacid",
}

# Spacers: each carries an acyl carbon at its proximal end and an amine at its distal end,
# so units concatenate directly into an amide chain. Written without a leading '*' because
# they are composed, not attached, by `assemble`; the '*' is added to the finished chain.
SPACER = {
    # gamma-glutamyl: linked through its side-chain carboxyl, alpha-amine left free.
    "gglu": "C(=O)CC[C@@H](C(=O)O)N",
    # 8-amino-3,6-dioxaoctanoyl, written AEEA / ADO / OEG / Novo's "PEG2".
    "oeg": "C(=O)COCCOCCN",
    # 6-aminohexanoyl.
    "ahx": "C(=O)CCCCCN",
    # beta-alanyl.
    "bala": "C(=O)CCN",
}

SPACER_SYNONYM = {
    "ggl": "gglu", "ggla": "gglu", "gglu": "gglu", "ge": "gglu", "glu": "gglu",
    "gammaglu": "gglu", "gammaglutamyl": "gglu", "glutamyl": "gglu",
    "ye": "gglu", "yglu": "gglu", "ygl": "gglu",
    "oeg": "oeg", "ado": "oeg", "aeea": "oeg", "peg2": "oeg",
    "ahx": "ahx", "aminohexanoyl": "ahx",
    "bala": "bala", "betaala": "bala", "betaalanine": "bala",
}

# PEG is excluded outright, in every form. The earlier rule admitted a discrete, countable
# oligomer - a defined molecule with an exact structure, which `PEG12` genuinely is - and
# excluded only the polydisperse kDa-scale material. In this cohort that distinction bought
# exactly ONE record out of 134 PEGylated ones, because 42 state no size at all and the rest
# are 5-46 kDa chains, frequently at an unstated or side-chain site. Keeping a single
# PEGylated record meant the set's PEG policy could not be stated in one sentence and that
# record was the only thing standing between "no polymer conjugates" and a caveat. The
# oligomer builder and the DEG unit table are gone with it, rather than left unreachable.

# Covalent organic labels whose structure is exact and public. A dye that is proprietary or
# named only by trade name is NOT here, and is blocking rather than approximated.
LABEL = {
    # Fluorescein isothiocyanate, isomer I: the thiocarbamoyl bonds to the anchor amine.
    "fitc": "*C(=S)Nc1ccc2c(c1)C1(OC2=O)c2ccc(O)cc2Oc2cc(O)ccc21",
    # 5-carboxytetramethylrhodamine, amide-linked.
    "tamra": "*C(=O)c1ccc2c(c1)C(=O)O[C@]21c2ccc(N(C)C)cc2Oc2cc(N(C)C)ccc21",
    # d-Biotin, amide-linked.
    "biotin": "*C(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12",
    # 7-hydroxycoumarin-3-carboxamide.
    "coumarin": "*C(=O)C1=Cc2ccc(O)cc2OC1=O",
    # DOTA and NOTA, amide-linked through one pendant arm; the radiometal itself is a
    # coordination complex and is recorded in `tracer_label`, not drawn.
    "dota": "*C(=O)CN1CCN(CC(=O)O)CCN(CC(=O)O)CCN(CC(=O)O)CC1",
    "nota": "*C(=O)CN1CCN(CC(=O)O)CCN(CC(=O)O)CC1",
    # Dibenzocyclooctyne, amide-linked.
    "dbco": "*C(=O)CCC(=O)N1CC2=CC=CC=C2C#CC2=CC=CC=C2C1",
    # Aryl sulfonate, as "SO3H group at Tyr-4": bonds to a ring carbon, not a nucleophile.
    "so3h": "*S(=O)(=O)O",
}

LABEL_SYNONYM = {
    "fitc": "fitc", "fluorescein": "fitc", "fluoresceinisothiocyanate": "fitc",
    "tamra": "tamra", "carboxytetramethylrhodamine": "tamra",
    "biotin": "biotin", "biotinylation": "biotin",
    "coumarin": "coumarin", "hydroxycoumarin": "coumarin",
    "dota": "dota", "nota": "nota", "dbco": "dbco",
    # NB no "so3h" entry - see the note on LABEL["so3h"].
}


# --- anchors ---------------------------------------------------------------------------
# Which atom of a residue accepts a conjugate, and which residues can accept one at all.
# Keyed by the one-letter code; the value names the anchor for `peptide_smiles.tag_site`.
ANCHOR = {"K": "amine", "C": "thiol", "S": "hydroxyl", "T": "hydroxyl", "Y": "hydroxyl"}

RESIDUE_NAME = {
    "lys": "K", "lysine": "K", "k": "K",
    "cys": "C", "cysteine": "C", "c": "C",
    "ser": "S", "serine": "S", "s": "S",
    "thr": "T", "threonine": "T", "t": "T",
    "tyr": "Y", "tyrosine": "Y", "y": "Y",
    "his": "H", "histidine": "H", "h": "H",
    "asn": "N", "glu": "E", "asp": "D", "gly": "G", "pro": "P", "arg": "R",
}


@dataclass(frozen=True)
class Site:
    """A resolved attachment point: 1-based residue position and the anchor atom kind."""
    position: int
    anchor: str
    residue: str


@dataclass(frozen=True)
class Clause:
    """One clause of a modification field, with its verdict.

    `sites` holds every attachment point the clause names: a dye written as "at Lysine 15,
    23 and 27" is one clause with three sites, and the fragment goes on at each.
    """
    text: str
    verdict: str                      # 'inert' | 'drawable' | 'blocking'
    reason: str
    sites: tuple[Site, ...] = ()
    fragment: str | None = None
    tracer: str = ""                  # 'isotope' | 'radiohalogen' | 'radiometal' | ''
    # A positional substitution claim that is consistent with the sequence but that the
    # sequence cannot positively confirm. See the note above SUBSTITUTION.
    unconfirmed: bool = False

    def __post_init__(self):
        object.__setattr__(self, "sites", tuple(self.sites))


@dataclass
class ModReport:
    """The verdict over a whole record's modification text."""
    clauses: list[Clause] = field(default_factory=list)
    tracer_label: str = ""

    @property
    def substitution_unconfirmed(self) -> bool:
        return any(c.unconfirmed for c in self.clauses)

    @property
    def blocking(self) -> list[Clause]:
        return [c for c in self.clauses if c.verdict == "blocking"]

    @property
    def drawable(self) -> list[Clause]:
        return [c for c in self.clauses if c.verdict == "drawable"]

    @property
    def status(self) -> str:
        if self.blocking:
            return self.blocking[0].reason
        return "all modifications accounted for"


# --- clause classification -------------------------------------------------------------

BLANK = re.compile(r"^(|free|n\.?a\.?|none|not reported|not mentioned|nan|-)$", re.I)

# The absence of a gamma-carboxyl, which is the plain glutamate the sequence already spells.
UNCARBOXYLATED = re.compile(r"\buncarboxylat\w*|\bun-?gamma-?carboxylat\w*", re.I)

# Blocking is tested FIRST, so that a clause naming both a drawable fragment and an
# undrawable partner ("Lys27 conjugaged with (PEG(10KDa)-HGEG...)") is blocked rather than
# half-drawn.
BLOCKING = [
    # Any PEG, of any size, discrete or polydisperse. `peg2`/OEG as a SPACER UNIT inside a
    # fatty-acyl recipe is not caught here: that is an 8-amino-3,6-dioxaoctanoyl residue
    # with an exact structure, reached through SPACER_SYNONYM during recipe assembly, and
    # it never appears as a polymer conjugate on its own.
    ("PEG conjugate",
     r"\bm?peg\b|peg\s*\d|peg\s*\(|pegylat|polyethylene glycol|"
     r"diethylene glycol|triethylene glycol|tetraethylene glycol"),
    ("glycosylation", r"glycosylat|lactose|\bglycan\b"),
    ("protein or peptide partner",
     r"\balbumin\b|\bhsa\b|\babd\b|abdcon|\bfc\b|\bigg\b|scfv|streptavidin|\bstrep\b|"
     r"transferrin|interferon|\bifn\b|gp41|apoc|\bga3\b|\btt1\b|\bip protein|nanobody|"
     r"repebody|\belp\b|xten|\bebp\b|chain linked with disulfide"),
    ("drug or cytotoxin payload",
     r"\bsn38\b|\bmmae\b|\bdtx\b|thapsigargin|indomethacin|diflunisal|\bdif\b|auristatin|"
     r"doxorubicin|\bdox\b|vitamin d3|cotinine|\bcot\b|l12adt"),
    ("particle or polymer vehicle",
     r"nanoparticle|microspher|liposom|polyplex|\bmnp|micelle|hydrogel|\bpcl\b|dextran|"
     r"polysarcosine|chondroitin|heparosan|poegma|chitosan|\bpamam\b|polyurethane|poly\("),
    ("dye structure not public",
     r"alexa|\bdir\b|flamma|zw800|fluorescence dye|fluorescent dye|fluorescently labeled|"
     r"fluorescently labelled|\bcy5\b|cy5\.5|\bcy7\b|indocarbocyanine|cyanine"),
    ("structure withheld by source",
     r"structure (is )?given in paper|structure givven in paper|given in paper|"
     r"structure given|\br[1-4]\b\s*,|=\s*zr\b"),
    ("radiometal without a named chelator",
     r"cu-?label|\bmn and\b|89zr|\b99mtc\b(?!.*hynic)|177lu(?!.*(peg|dota|nota))|"
     r"\b64cu\b(?!.*(dota|nota))|\b68ga\b(?!.*(dota|nota))|\b18f\b|\bal18f\b|noda"),
    ("backbone surrogate not representable",
     r"ch2nh|reduced amine bond|silaprolin|\bsip\b|tmsala|methylation of alpha|"
     r"methylation at|\bn-?me|n-methyl|peptoid"),
    ("non-standard residue not in the template set",
     r"\bagp\b|\bdpr\b|propargylglycine|naphthyl|napthyl|cyclohexylalanin|\bcha\b|"
     r"\boic\b|isovalin|norvalin|\bnva\b|homoleucin|\bhle\b|\bpff\b|phe\(4-f\)|"
     r"hydroxyprolin|\bhyp\b|\b1nal\b|\bnal\b|\bamha\b|ethylprolin|\betpro\b|"
     r"aminoisobutyric|\baib\b|\bibu\b|\biva\b|sarcosin|\bsar\b|ornithin|\born\b|"
     r"homo-?arginin|\bharg\b|diaminobutyric|\bdab\b|\bnle\b|\babu\b|pyroglutam|"
     r"pyro-?glutam|\bpglu\b|aminoisobutyr|beta alanine|triazole|\btz\d|\bdeap\b|"
     r"hynic|\bzipp\b|dehydro-?pro|cystine|\bmodified central region\b"),
    ("multimerisation", r"tetrameri|dimeri|interconnected through"),
    ("conjugate of unstated composition",
     r"\bx\d?\s*=\s*linker|\blinker\b(?!.*\b(ggl|oeg|ado|aeea|gly)\b)|\bl3\b|\bs4\b|"
     r"\bc2\s*=|\bc1\s*=|eliminative linker|\bmal-|maleimide|vinylsulfone|"
     r"isothiocyanate(?! )|conjugated with drug"),
]

# Isotopic substitution and tracer labelling. Tested before the sequence-descriptive and
# drawable rules so that the `tracer` flag is always set.
TRACER = [
    ("isotope", r"\b3h\b|\[3h\]|tritiat|tritium|\bc14\b|\b14c\b|\bc-14\b|\b13c\b|"
                r"\bn15\b|\b15n\b"),
    ("radiohalogen", r"\b12[35]i\b|\bi-?12[35]\b|\b131i\b|\(12[35]\)i|iodinat|radioiod|"
                     r"\biodo\b"),
    ("radiometal", r"\b177lu\b|\b99mtc\b|\b64cu\b|\b68ga\b|\bdota\b|\bnota\b"),
    # "Radiolabelling of exendin 3 at position 1 (histidine)" names no isotope but is
    # unambiguously a tracer preparation, and is treated as one.
    ("radiolabel", r"radiolab\w*|radioactive label"),
]

# Text that describes the sequence column rather than adding atoms to the molecule.
# Residue tokens a bare clause may legitimately consist of: the 20 three-letter codes, the
# full names, and the non-standard codes this cohort writes out. Anything else that is three
# letters long is not a residue and must not be waved through as one.
_RESIDUE_TOKEN = "(?:" + "|".join([
    "ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", "his", "ile", "leu", "lys",
    "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val",
    "alanine", "arginine", "asparagine", "aspartate", "aspartic", "cysteine", "glutamine",
    "glutamate", "glutamic", "glycine", "histidine", "isoleucine", "leucine", "lysine",
    "methionine", "phenylalanine", "proline", "serine", "threonine", "tryptophan",
    "tryptophane", "tyrosine", "valine",
    "aib", "nle", "orn", "dab", "sar", "cha", "abu", "nal", "iva", "nva", "hyp", "harg",
    "pglu", "aze", "oic", "dpr", "agp",
]) + ")"

SEQUENCE_DESCRIPTIVE = [
    # A clause that is nothing but a residue token is the tail of a substitution list:
    # "Substituition of amino acids G,L,K,Q,M,E,..." splits into one clause per letter, and
    # "Gly2", "Nle10", "D-Phe11", "desB27", "D-serine" are the same thing spelled out.
    # The token must be an actual residue: an unrestricted three-letter rule here read
    # "AFF" - an affibody fused to the C terminus - as a residue name and called a 15 kDa
    # protein partner atom-free.
    r"^(?:d-)?(?:[acdefghiklmnpqrstvwy]|" + _RESIDUE_TOKEN + r")\d*$",
    r"^des-?(?:" + _RESIDUE_TOKEN + r"|[a-z])?\d+$",
    # Legends: "Aib = a-aminoisobutyric acid", "O denotes L-ornithine", "y = D-Tyr".
    r"^[^=]{1,30}=[^=]+$",
    r"\bdenotes\b",
    # Point-substitution shorthand and prose: "V24L modifications", "T639I mutation",
    # "Substituition of G for A at position 2", "Ala2 to Gly2", "des-Gly".
    r"^[a-z]\d+[a-z]( (substitu\w*|modification|mutation)s?)?$",
    r"substi\w*|replac\w*|\bmutation\b|^introduction of|\bintroduced\b|\bmodifications?$|\bmodificat\w*$|"
    r"\bdes-?[a-z]{3}\b|^removing |^removal|^modif\w*|^addition at|"
    r"^insertion of|^incorporation of|^inserting |^post translational",
    # Chirality notes, already carried by the Chirality column.
    r"\bd-amino acids?\b|\ball d-amino\b|\bd-isomers?\b|\bd isoform\b|all d-",
    # Positional notes naming a residue the sequence already spells.
    r"^[a-z]{3}\d+(,\s*[a-z]{3}\d+)* (modification|substitu\w*)",
    r"^(two )?amphipathic helices|^modified by adding|^at \d+ position",
    r"\baddition at\b|\baddition of\b|^[a-z]{3} addition\b|^two amphipathic",
]

ACYL_SITE_VERB = (r"(acylat\w*|conjugat\w*|link\w*|modif\w*|appended|attach\w*|"
                  r"present at|lipidat\w*|\bat\b|\bof\b)")


# " and " joins two independent statements only when the right-hand side carries its own
# verb: "Lys-34 is replaced by Arg and Lys-26 is acylated with a C16" is two clauses, while
# "Lysine 15,23 and 27" is one site list. Splitting the latter produced clause fragments
# like "23 and 27", which were the single largest source of unrecognised text.
CLAUSE_VERB = re.compile(
    r"\b(acylat\w*|conjugat\w*|link\w*|label\w*|replac\w*|substitut\w*|modif\w*|"
    r"attach\w*|append\w*|pegylat\w*|is |was |were )", re.I)


# --- substitution claims, checked rather than assumed -----------------------------------
# A clause like "V24L modifications" or "Lys-34 is replaced by Arg" adds no atoms, so the
# obvious verdict is `inert`. That verdict carries a hidden assumption: that the sequence
# column already HAS the substitution applied. For part of this cohort it does not.
#
# Liraglutide is the clearest case. PEPlife2 stores the sequence
# HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG - native GLP-1(7-37), with lysine at the position the
# drug has arginine - and records the difference only as free text, "K34R substituition at
# GLP-1". Treating that text as inert emits native GLP-1 and labels it liraglutide:
# C172H265N41O51 against the real C172H265N43O51, two nitrogens and a guanidine short of
# the molecule whose half-life is on the row.
#
# So every substitution claim that names a position is checked against the sequence. The
# product residue must be the one actually sitting there. Where the numbering is relative
# to a parent peptide - liraglutide's "34" indexes GLP-1 from residue 1, not the stored
# 31-mer - the position falls outside the sequence and the claim cannot be checked at all,
# which is also a refusal: an unverifiable claim about which residues are present is not a
# basis for drawing a structure.
#
# One class of error survives this check and cannot be caught from a single row. A claim
# whose stated position ALREADY holds the product residue is self-consistent whether or not
# the sequence column is right. Three liraglutide rows read "Lys-30 is replaced by Arg",
# position 30 does hold Arg, and the check passes - yet the stored sequence is still native
# GLP-1(7-37), carrying Lys where the drug has Arg two residues earlier, and the built
# molecule comes out as C172H265N41O51 against liraglutide's C172H265N43O51. Detecting that
# needs the parent sequence, which the row does not carry. Such claims are therefore
# flagged `unconfirmed` rather than silently trusted, and the flag reaches the built
# dataset as `substitution_unconfirmed` so these rows can be found in one filter.
AA3 = {"ala": "A", "arg": "R", "asn": "N", "asp": "D", "cys": "C", "gln": "Q",
       "glu": "E", "gly": "G", "his": "H", "ile": "I", "leu": "L", "lys": "K",
       "met": "M", "phe": "F", "pro": "P", "ser": "S", "thr": "T", "trp": "W",
       "tyr": "Y", "val": "V"}
AA_FULL = {"alanine": "A", "arginine": "R", "asparagine": "N", "aspartate": "D",
           "aspartic acid": "D", "cysteine": "C", "glutamine": "Q", "glutamate": "E",
           "glutamic acid": "E", "glycine": "G", "histidine": "H", "isoleucine": "I",
           "leucine": "L", "lysine": "K", "methionine": "M", "phenylalanine": "F",
           "proline": "P", "serine": "S", "threonine": "T", "tryptophan": "W",
           "tyrosine": "Y", "valine": "V"}
_NAME3 = "|".join(AA3)
_NAMEF = "|".join(sorted(AA_FULL, key=len, reverse=True))

SUBSTITUTION = [
    # "K34R", "V24L modifications", "T639I mutation", "A14E"
    (r"\b([acdefghiklmnpqrstvwy])(\d{1,4})([acdefghiklmnpqrstvwy])\b", "src_pos_dst"),
    # "Lys-34 is replaced by Arg", "Lys-30 is replaced by Arg"
    (rf"\b({_NAME3})\w*[\s\-]*(\d{{1,4}})\s+(?:is|was|were)?\s*replaced?\s+"
     rf"(?:by|with)\s+({_NAME3})", "src_pos_dst_named"),
    # "Ala2 to Gly2", "Trp3 -> Phe3", "Ala2 was substituted with Gly"
    (rf"\b({_NAME3})\w*\s*(\d{{1,4}})\s*(?:->|to|was substituted with)\s*"
     rf"({_NAME3})", "src_pos_dst_named"),
    # "Substituting arginine 14 with alanine"
    (rf"\b({_NAMEF})\s+(\d{{1,4}})\s+with\s+({_NAMEF})", "src_pos_dst_full"),
    # "Ala5 substituition in GpTx-1" - product unknown, but the source residue must be gone.
    (rf"\b({_NAME3})\w*\s*(\d{{1,4}})\s+substi\w*", "src_pos_only"),
]


def _substitution_claims(text: str) -> list[tuple[int, str | None, str | None]]:
    """Every (position, expected residue, displaced residue) claim in a clause."""
    claims = []
    for pattern, kind in SUBSTITUTION:
        for match in re.finditer(pattern, text, re.I):
            first, position, last = match.group(1), int(match.group(2)), match.group(3) \
                if match.lastindex and match.lastindex >= 3 else None
            if kind == "src_pos_dst":
                src, dst = first.upper(), (last or "").upper()
            elif kind == "src_pos_dst_named":
                src = AA3.get(first.lower()[:3])
                dst = AA3.get((last or "").lower()[:3])
            elif kind == "src_pos_dst_full":
                src, dst = AA_FULL.get(first.lower()), AA_FULL.get((last or "").lower())
            else:
                src, dst = AA3.get(first.lower()[:3]), None
            if src in STANDARD or dst in STANDARD:
                claims.append((position, dst, src))
    return claims


STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def _check_substitutions(text: str, sequence: str) -> str:
    """'' when every positional substitution claim holds, else why it does not."""
    for position, expected, displaced in _substitution_claims(text):
        if not 1 <= position <= len(sequence):
            return (f"substitution at position {position} cannot be checked against a "
                    f"{len(sequence)}-residue sequence, so the numbering is relative to a "
                    f"parent the sequence column does not hold")
        actual = sequence[position - 1]
        if expected is not None and expected in STANDARD:
            if actual != expected:
                return (f"source reports {expected} at position {position} but the "
                        f"sequence column has {actual}: the substitution is not applied")
        elif displaced is not None and actual == displaced:
            return (f"source reports {displaced}{position} was substituted but the "
                    f"sequence column still has {displaced} there")
    return ""


def _split_clauses(text: str) -> list[str]:
    """Split on commas outside brackets, then on ' and ' only between real clauses."""
    out, depth, current = [], 0, []
    for i, ch in enumerate(text):
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth = max(depth - 1, 0)
        # A comma directly in front of a bare number continues a position list
        # ("Lysine 15,23 and 27") and must not break the clause.
        if ch == "," and depth == 0 and not re.match(r"\s*\d", text[i + 1:]):
            out.append("".join(current))
            current = []
        else:
            current.append(ch)
    out.append("".join(current))

    clauses = []
    for chunk in out:
        parts, buffer = [], chunk
        while True:
            match = re.search(r"\s+and\s+", buffer)
            if not match:
                break
            head, tail = buffer[:match.start()], buffer[match.end():]
            if CLAUSE_VERB.search(tail) and CLAUSE_VERB.search(head):
                parts.append(head)
                buffer = tail
            else:
                break
        parts.append(buffer)
        clauses.extend(parts)
    return [c.strip() for c in clauses if c.strip()]


def _normalise(value: object) -> str:
    """Lower-case, collapse whitespace, and fold the Greek and mojibake PEPlife2 carries."""
    text = " ".join(str(value).split()).lower()
    for src, dst in (("γ", "g"), ("α", "a"), ("β", "b"), ("μ", "u"),
                     ("Î³", "g"), ("Î±", "a"), ("Î²", "b"),
                     ("Ã®", ""), ("×", "x"), ("→", "->"),
                     ("Ψ", "psi"), ("–", "-"), ("—", "-")):
        text = text.replace(src, dst)
    return text


# Fatty-acyl and PEG tokens look exactly like residue-with-position tokens: "C16" is a
# palmitoyl chain, not cysteine 16, and "PEG4" is not proline 4. They are masked out before
# the site search so a chain length can never be read as an attachment position.
CHEMICAL_TOKEN = re.compile(
    r"\bc[\s\-]?\d{1,2}\b(?![a-z])|\bpeg[\s\-]?\d{1,2}\b|\(?-?\(?[gG]e\)?\d\b|"
    r"\b\d+\s*x\s*(?:gglu|ggl|oeg|ado|aeea|peg2)\b|\b(?:gglu|ggl|oeg|ado|aeea)\b|"
    r"\b\d{1,3}\s*kda\b|\b1[0-9]{2}[a-z]{1,2}\b", re.I)

# "Lys24", "Lysine16", "Lys-26", "K-26", "Cys14". No trailing \b after the residue name:
# in "Lysine16" there is no word boundary between the 'e' and the '1', which silently cost
# every record written that way.
SITE_NAMED = re.compile(
    r"\b(lys(?:ine)?|cys(?:teine)?|ser(?:ine)?|thr(?:eonine)?|tyr(?:osine)?|"
    r"his(?:tidine)?)[\s\-]*(\d+)", re.I)
# A bare one-letter form, which only ever appears hyphenated or glued: "K-26", "K12".
SITE_LETTER = re.compile(r"\b([kcsty])[\s\-]?(\d{1,3})\b", re.I)
# Trailing members of a position list: "Lysine 15,23 and 27" -> 23, 27.
SITE_EXTRA = re.compile(r"(?:,|\band\b)\s*(\d{1,3})\b")


def _resolve(code: str, position: int, sequence: str) -> tuple[Site | None, str]:
    """Check one (residue, position) claim against the sequence."""
    if not 1 <= position <= len(sequence):
        return None, f"{code}{position} lies outside the sequence"
    actual = sequence[position - 1]
    if actual != code:
        return None, f"text names {code}{position} but the sequence has {actual} there"
    if code not in ANCHOR:
        return None, f"residue {code}{position} cannot carry a conjugate"
    return Site(position, ANCHOR[code], code), ""


def _find_sites(clause: str, sequence: str) -> tuple[list[Site], str]:
    """Resolve every attachment residue named in `clause`, or say why it cannot be.

    A clause may legitimately name several ("acetylation at Lys10 and Lys23",
    "Cy5.5 at Lysine 15,23 and 27"), in which case the fragment goes on at each.
    """
    masked = CHEMICAL_TOKEN.sub(" ", clause)

    named = SITE_NAMED.findall(masked) or [
        (letter, number) for letter, number in SITE_LETTER.findall(masked)]
    if named:
        codes = {RESIDUE_NAME[name.lower()] for name, _ in named}
        if len(codes) > 1:
            return [], "conjugate stated across two different residue types"
        code = codes.pop()
        positions = [int(number) for _, number in named]
        # Trailing bare numbers after the last named residue belong to the same list.
        tail = masked[masked.rfind(str(positions[-1])) + len(str(positions[-1])):]
        positions += [int(n) for n in SITE_EXTRA.findall(tail)]

        sites = []
        for position in sorted(set(positions)):
            site, problem = _resolve(code, position, sequence)
            if site is None:
                return [], problem
            sites.append(site)
        return sites, ""

    at_position = re.search(
        r"(?:at|of)?\s*position\s*(?:at\s*)?(\d+)|\bat (\d+)\s*(?:st|nd|rd|th)?\s*"
        r"position|\b(\d+)\s*(?:st|nd|rd|th) pos", masked, re.I)
    if at_position:
        position = int(next(g for g in at_position.groups() if g))
        if not 1 <= position <= len(sequence):
            return [], "named position lies outside the sequence"
        code = sequence[position - 1]
        if code not in ANCHOR:
            return [], f"residue {code}{position} cannot carry a conjugate"
        return [Site(position, ANCHOR[code], code)], ""

    residue_only = {RESIDUE_NAME[m.lower()] for m in re.findall(
        r"\b(lys(?:ine)?|cys(?:teine)?|ser(?:ine)?|thr(?:eonine)?|tyr(?:osine)?)\b",
        masked, re.I)}
    if len(residue_only) == 1:
        code = residue_only.pop()
        hits = [i + 1 for i, aa in enumerate(sequence) if aa == code]
        if len(hits) == 1:
            return [Site(hits[0], ANCHOR[code], code)], ""
        return [], f"site named only as {code} and the sequence has {len(hits)} of them"
    return [], "attachment site not stated"


def _recipe_units(clause: str) -> tuple[list[tuple[str, int]], str | None, str]:
    """Parse a conjugate recipe into (spacer units, terminal acyl key, problem)."""
    bracketed = re.findall(r"\(([^()]*)\)", clause)
    body = bracketed[-1] if bracketed else clause

    acyl_key = None
    for word, key in ACYL_SYNONYM.items():
        if re.search(rf"\b{word}", body) or re.search(rf"\b{word}", clause):
            acyl_key = key
            break
    if acyl_key is None:
        chain = re.search(r"\bc[\s\-]?(\d{1,2})\b[\s\-]*(di-?acid|dioic|monoacid)?", body) \
            or re.search(r"\bc[\s\-]?(\d{1,2})\b[\s\-]*(di-?acid|dioic|monoacid)?", clause)
        if chain:
            key = f"c{int(chain.group(1))}"
            if chain.group(2) and "mono" not in chain.group(2):
                key += "-diacid"
            if key in ACYL:
                acyl_key = key
    if acyl_key is None and re.search(r"fatty acid|fatty chain|lipidat|acylat", clause):
        return [], None, "fatty acyl named without a chain length"

    # "2xOEG" glues the multiplier to the token, and the spacer pattern below needs a word
    # boundary in front of the token - between "x" and "oeg" there is none, so the whole
    # spacer went unmatched and was silently dropped. A dropped spacer is the one failure
    # mode this function must never have: it does not refuse the recipe, it assembles a
    # SHORTER conjugate and reports success, so "C18 diacid-gGlu-2xOEG" came out as
    # C18-diacid-gGlu with both OEG units missing. Separating the multiplier first makes
    # the count reachable, and the audit below refuses anything still unaccounted for.
    body = re.sub(r"(\d+)\s*[x×]\s*", r"\1x ", body)

    units: list[tuple[str, int]] = []
    consumed = 0
    for count_prefix, token, count_suffix in re.findall(
            r"(\d+)?\s*x?\s*\(?\b(ggl?u?|[gy]e|[gy]glu|oeg|ado|aeea|peg2|ahx|bala)\b"
            r"\)?\s*x?(\d+)?",
            body):
        key = SPACER_SYNONYM.get(token)
        if key is None:
            continue
        consumed += 1
        count = int(count_prefix or count_suffix or 1)
        units.append((key, count))
    # Every spacer name the recipe mentions must have become a unit. Where one did not, the
    # recipe is not understood, and saying so is the only safe answer - assembling what was
    # understood would emit a conjugate missing a piece of itself.
    mentioned = len(re.findall(
        r"\b(?:ggl?u?|[gy]e|[gy]glu|oeg|ado|aeea|peg2|ahx|bala)\b", body))
    if mentioned != consumed:
        return [], acyl_key, "conjugate names a spacer this recipe cannot place"
    if len({key for key, _ in units}) > 1:
        return [], acyl_key, "conjugate has two spacer types in an unstated order"
    return units, acyl_key, ""


def assemble(units: list[tuple[str, int]], acyl_key: str | None) -> str | None:
    """Build a conjugate SMILES from the anchor outward. Spacers proximal, acyl distal."""
    if acyl_key is None and not units:
        return None
    chain = "".join(SPACER[key] * count for key, count in units)
    if acyl_key is None:
        return None
    return "*" + chain + ACYL[acyl_key].lstrip("*")


def classify(clause: str, sequence: str) -> Clause:
    """One clause in, one verdict out."""
    text = _normalise(clause)
    if BLANK.fullmatch(text):
        return Clause(clause, "inert", "blank")

    tracer = ""
    for kind, pattern in TRACER:
        if re.search(pattern, text):
            tracer = kind
            break

    # A discrete PEG oligomer must be recognised before the polydisperse-PEG block, which
    # is deliberately broad.
    # A residue legend for a residue that IS templated, and that the sequence actually
    # contains, is descriptive. Checked before the blocking table, whose non-standard
    # residue rule is deliberately broad and would otherwise reject the very records the
    # new templates exist to admit.
    present = {t.lower() for t in sequence} if not isinstance(sequence, str) else set()
    named = named_templated_residues(text)
    if named and named <= present and not re.search(
            r"conjugat|acylat|lipidat|pegylat|link|label|fusion|glycosylat", text):
        return Clause(clause, "inert",
                      f"legend for {', '.join(sorted(named))}, which the sequence carries "
                      f"and the builder templates")

    # "Uncarboxylation at all glutamic acid" takes nothing off the structure the sequence
    # drew. The modified residue in this family is gamma-carboxyglutamate; its
    # UNcarboxylated form is ordinary glutamate, which is exactly the E the sequence column
    # already carries, so the note describes the sequence rather than a change to it.
    # Settled here rather than among the sequence-descriptive rules at the end, because the
    # acylation branch between the two reads "...at all glutamic acid" as a site-named
    # gamma-glutamyl recipe and blocks on a composition it cannot assemble.
    if UNCARBOXYLATED.search(text):
        return Clause(clause, "inert", "uncarboxylated glutamate is the glutamate the "
                                       "sequence already spells; adds no atoms")

    for reason, pattern in BLOCKING:
        if re.search(pattern, text):
            return Clause(clause, "blocking", reason, tracer=tracer)

    if tracer in ("isotope", "radiohalogen", "radiometal", "radiolabel"):
        # A chelator named alongside the metal is a real covalent adduct and is drawn; the
        # metal is not. DOTA/NOTA reach here only when the earlier blocking rule for a
        # chelator-less metal did not fire.
        for word, key in LABEL_SYNONYM.items():
            if key in ("dota", "nota") and re.search(rf"\b{word}\b", text):
                sites, problem = _find_sites(text, sequence)
                if sites:
                    return Clause(clause, "drawable", f"{key} chelator", sites,
                                  LABEL[key], tracer)
        return Clause(clause, "inert", f"{tracer} tracer label", tracer=tracer)

    for word, key in LABEL_SYNONYM.items():
        if re.search(rf"\b{word}\b", text):
            sites, problem = _find_sites(text, sequence)
            if not sites:
                return Clause(clause, "blocking", f"{key} label: {problem}", tracer=tracer)
            return Clause(clause, "drawable", f"{key} label", sites, LABEL[key], tracer)

    # Side-chain acylation, lipidation and spacer conjugates.
    if re.search(ACYL_SITE_VERB, text) and re.search(
            r"\bc[\s\-]?\d{1,2}\b|fatty|acyl|lipid|" + "|".join(ACYL_SYNONYM) + "|"
            + "|".join(SPACER_SYNONYM), text):
        units, acyl_key, problem = _recipe_units(text)
        if problem:
            return Clause(clause, "blocking", problem, tracer=tracer)
        fragment = assemble(units, acyl_key)
        if fragment is None:
            return Clause(clause, "blocking", "conjugate composition not resolved",
                          tracer=tracer)
        sites, problem = _find_sites(text, sequence)
        if not sites:
            return Clause(clause, "blocking", f"side-chain conjugate: {problem}",
                          tracer=tracer)
        return Clause(clause, "drawable", "side-chain conjugate", sites, fragment,
                      tracer)

    # Side-chain acetylation, which the terminal columns do not cover.
    if re.search(r"acetylation at (lys|k)", text):
        sites, problem = _find_sites(text, sequence)
        if not sites:
            return Clause(clause, "blocking", f"side-chain acetylation: {problem}")
        return Clause(clause, "drawable", "side-chain acetylation", sites, "*C(C)=O")

    for pattern in SEQUENCE_DESCRIPTIVE:
        if re.search(pattern, text):
            # "Adds no atoms" is only true if the sequence column already carries the
            # change the clause describes. Where the clause says which residue sits where,
            # that is checkable, so it is checked.
            mismatch = _check_substitutions(text, sequence)
            if mismatch:
                return Clause(clause, "blocking", mismatch, tracer=tracer)
            claims = _substitution_claims(text)
            vacuous = any(expected is not None and 1 <= position <= len(sequence)
                          and sequence[position - 1] == expected
                          and displaced is not None and displaced != expected
                          for position, expected, displaced in claims)
            if vacuous:
                # A claim whose stated position already holds the product residue is
                # consistent with the sequence whether or not the sequence is right, so the
                # sequence cannot confirm it. Liraglutide is the proof that "consistent" is
                # not good enough: three rows read "Lys-30 is replaced by Arg", position 30
                # does hold Arg, and the stored sequence is still native GLP-1(7-37) - two
                # nitrogens and a guanidine short of the drug named on the row. Such a
                # record is excluded rather than flagged and kept, because a structure that
                # may silently be the parent molecule is not a structure for the label.
                return Clause(clause, "blocking",
                              "substitution claim cannot be confirmed by the sequence: the "
                              "stated position already holds the product residue, so the "
                              "sequence is consistent with the change either way",
                              unconfirmed=True)
            return Clause(clause, "inert", "describes the sequence, adds no atoms")

    return Clause(clause, "blocking", "modification text not recognised", tracer=tracer)


def report(sequence: str, *fields: object) -> ModReport:
    """Classify every clause of every modification field for one record."""
    out = ModReport()
    tracers: list[str] = []
    for value in fields:
        text = " ".join(str(value).split())
        if BLANK.fullmatch(_normalise(text)):
            continue
        for clause in _split_clauses(text):
            verdict = classify(clause, sequence)
            out.clauses.append(verdict)
            if verdict.tracer:
                tracers.append(verdict.tracer)
    for kind in ("radiometal", "radiohalogen", "radiolabel", "isotope"):
        if kind in tracers:
            out.tracer_label = kind
            break
    return out


def selfcheck() -> list[str]:
    """Every fragment in the library must parse, and carry exactly one attachment dummy."""
    problems = []
    library = {**ACYL, **LABEL}
    for name, smiles in library.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            problems.append(f"{name}: does not parse ({smiles})")
            continue
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummies) != 1:
            problems.append(f"{name}: {len(dummies)} attachment dummies, expected 1")
    # A spacer is only ever used mid-chain, so it is checked the way it is used: capped at
    # the proximal end and with a second copy of itself hung off its distal amine.
    for name, smiles in SPACER.items():
        if Chem.MolFromSmiles("C" + smiles) is None:
            problems.append(f"spacer {name}: does not parse ({smiles})")
        elif Chem.MolFromSmiles("C" + smiles + smiles) is None:
            problems.append(f"spacer {name}: does not chain to itself ({smiles})")
    # `so3h` is in LABEL with no synonym on purpose; every other entry must be reachable.
    unreachable = (set(ACYL) | set(SPACER) | (set(LABEL) - {"so3h"})) - (
        set(ACYL_SYNONYM.values()) | set(SPACER_SYNONYM.values())
        | set(LABEL_SYNONYM.values()))
    problems += [f"{name}: defined but no synonym reaches it" for name in sorted(unreachable)]
    for name, key in {**ACYL_SYNONYM, **SPACER_SYNONYM, **LABEL_SYNONYM}.items():
        if key not in ACYL and key not in SPACER and key not in LABEL:
            problems.append(f"synonym {name} points at missing fragment {key}")
    return problems


if __name__ == "__main__":
    issues = selfcheck()
    print("\n".join(issues) if issues else "fragment library OK")
