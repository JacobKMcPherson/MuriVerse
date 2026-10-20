"""Physicochemical descriptors of the dosed molecule, for the half-life feature matrix.

This block exists because the cohort had no feature that describes the MOLECULE in terms a
clearance mechanism is written in. The molecule-annotation block - `hl_engineering`, `topology`
and `chirality` (`smiles_coverage` is constant on the strict cohort and is dropped; `species`
now lives in the study-context block) - holds bucketed CLASS labels, so two different peptides
of the same class receive the identical vector and therefore the identical prediction, and on
its own it scores rho = 0.104 on held-out clusters and studies. The two learned
representations that can see the molecule (ESM-2 over the sequence, PeptideCLM over the SMILES)
score no better on their own (0.145 and 0.177), which leaves the molecule side of the model
carrying almost nothing.

The descriptors here are deliberately NOT another learned embedding. They are the handful of
quantities that renal filtration, protease exposure and plasma-protein binding actually
depend on - size, lipophilicity, polar surface, hydrogen bonding, flexibility, and above all
net charge, since cationicity drives both glomerular handling and tubular reuptake for
peptides of this size. Being interpretable and low-dimensional (16 columns against
PeptideCLM's 768) matters on a few hundred rows, and being computed from the verified SMILES
rather than from the sequence means they are defined for molecules that have no faithful
single-letter sequence at all - which is the case that Figure 4 turns on, because six of the
ten residues of a polymyxin are 2,4-diaminobutyric acid and have no code.

THE BLOCK IS PART OF THE SHIPPED MODEL. The shipped configuration is molecule annotation +
study context + terminal classes + these 16 descriptors (71 columns), fitted on the 213
intravenous records. Measured on that cohort under the protocol every other number in Figure 3
uses (5-fold cross-validation holding out whole groups under the union of homology cluster and
source, three seeds, 95% group bootstrap; `table31_muriverse_halflife_feature_ablation_strict_iv.csv`):

    molecule annotation only                      9 cols   rho = 0.104 [-0.37, 0.37]
    physchem only                                16 cols   rho = 0.168 [ 0.01, 0.47]
    study context only                           33 cols   rho = 0.347 [ 0.12, 0.49]
    physchem + study context                     49 cols   rho = 0.522 [ 0.41, 0.63]
    physchem + study context + termini           62 cols   rho = 0.587 [ 0.40, 0.68]
    annotation + context + termini + physchem    71 cols   rho = 0.574 [ 0.35, 0.69]  (shipped)
    SMILES (PeptideCLM) only                    768 cols   rho = 0.177 [ 0.02, 0.44]
    SMILES + physchem + study context           817 cols   rho = 0.247 [ 0.12, 0.48]

Two things follow. Adding these descriptors to study context takes it from 0.347 to 0.522,
which is the largest single gain in the ablation. And putting the learned structural embedding
back on top of them makes it much WORSE (0.247 against 0.522), so on a cohort this size
PeptideCLM is not adding structural information to the descriptors, it is adding variance.
The best-scoring configuration overall drops the annotation block (0.587); the shipped one
keeps it, at a cost of 0.013, so that a prediction always rests on the engineering class,
topology and chirality a curator recorded (see REQUIRED_SHIPPED_BLOCKS in
evaluate_halflife_model.py).

Earlier versions of this docstring quoted measurements on a 488-row route-pooled cohort and
stated that the block "does NOT improve the shipped configuration (annotation + termini)".
Both were true of that cohort and that configuration and are superseded: the shipped
configuration changed (it now carries study context and these descriptors) and so did the
cohort.

n_basic and n_acidic count SMARTS MATCHES, not charge centres: the amine pattern and the
guanidine pattern overlap on guanidine nitrogens, so an arginine contributes more than one
to n_basic. net_charge_ph74 is n_basic - n_acidic, a structural count with no pKa in it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem import rdMolDescriptors as rdMD

RDLogger.DisableLog("rdApp.*")

# Kyte-Doolittle hydropathy, for the sequence-side columns. Only the 20 standard residues
# carry a value; a non-standard residue is skipped rather than guessed, so GRAVY on a
# peptide full of Dab is the mean over the residues that do have one.
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# A primary/secondary amine that is not an amide nitrogen, and a guanidine. Counted as
# SMARTS rather than read off `formal_charge`, which is zero for every record: the builders
# emit neutral SMILES, so the ionisation state at pH 7.4 has to be inferred from the groups
# that carry it.
BASIC_AMINE = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
GUANIDINE = Chem.MolFromSmarts("[NX3][CX3](=[NX2,NX3+])[NX3]")
CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")


def structure_descriptors(smiles: str) -> dict[str, float]:
    """RDKit descriptors of one molecule. Empty dict if the SMILES does not parse."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return {}
    n_basic = (len(mol.GetSubstructMatches(BASIC_AMINE))
               + len(mol.GetSubstructMatches(GUANIDINE)))
    n_acidic = len(mol.GetSubstructMatches(CARBOXYLIC_ACID))
    return {
        "mw": Descriptors.MolWt(mol),
        "clogp": Crippen.MolLogP(mol),
        "molmr": Crippen.MolMR(mol),
        "tpsa": rdMD.CalcTPSA(mol),
        "labute_asa": rdMD.CalcLabuteASA(mol),
        "hbd": rdMD.CalcNumHBD(mol),
        "hba": rdMD.CalcNumHBA(mol),
        "rotatable_bonds": rdMD.CalcNumRotatableBonds(mol),
        "rings": rdMD.CalcNumRings(mol),
        "aromatic_rings": rdMD.CalcNumAromaticRings(mol),
        "frac_csp3": rdMD.CalcFractionCSP3(mol),
        "heavy_atoms": float(mol.GetNumHeavyAtoms()),
        "amide_bonds": rdMD.CalcNumAmideBonds(mol),
        "n_basic": float(n_basic),
        "n_acidic": float(n_acidic),
        "net_charge_ph74": float(n_basic - n_acidic),
    }


def sequence_descriptors(sequence: str) -> dict[str, float]:
    """Composition of one sequence. Non-standard residues are excluded from the fractions."""
    residues = [c for c in str(sequence).upper() if c in KYTE_DOOLITTLE]
    denom = len(residues) or 1
    positive = sum(c in "KR" for c in residues)
    negative = sum(c in "DE" for c in residues)
    return {
        "seq_length": float(len(str(sequence))),
        "gravy": sum(KYTE_DOOLITTLE[c] for c in residues) / denom,
        "frac_positive": positive / denom,
        "frac_negative": negative / denom,
        "net_sequence_charge": float(positive - negative),
        "frac_aromatic": sum(c in "FWY" for c in residues) / denom,
        "frac_aliphatic": sum(c in "AVLIM" for c in residues) / denom,
        "frac_polar": sum(c in "STNQ" for c in residues) / denom,
        "frac_gly_pro": sum(c in "GP" for c in residues) / denom,
        "frac_cys": sum(c == "C" for c in residues) / denom,
    }


def build(meta: pd.DataFrame, smiles_column: str = "smiles",
          sequence_column: str | None = None) -> pd.DataFrame:
    """One row of descriptors per metadata row, in metadata order.

    STRUCTURE ONLY by default. `sequence_descriptors` is kept in this module and is tested,
    but the block does not include it, for two independent reasons that happen to agree:

    * It scores WORSE. Under the same grouped CV, on the earlier 488-row cohort, adding the ten
      residue-composition columns to the sixteen structural ones costs the configuration
      Figure 4 applies about two points of Spearman (structure only + annotation + termini
      0.517 [0.24, 0.63]; with composition 0.500 [0.22, 0.62]) and the block on its own
      about two as well (0.178 against 0.161). On a few hundred rows ten extra columns of a
      quantity the ESM-2 block already describes is dilution, not information.

    * It is not computable for the molecules the model is FOR. A composition fraction needs
      a faithful single-letter sequence, and the application case does not have one: six of
      the ten residues of a polymyxin are 2,4-diaminobutyric acid, which has no code, so
      GRAVY would be averaged over the four residues that happen to have one. A
      configuration offered for application has to be one that can actually be computed on
      the molecule it is applied to, otherwise the model shipped is not the model scored.

    A SMILES that does not parse yields NaN rather than a zero row, because zero is a
    meaningful value for most of these columns and XGBoost handles NaN as missing natively.
    Constant columns are dropped: the builders emit neutral SMILES, so `formal_charge` is
    zero for every record and a constant column is a split the tree can never use.
    """
    sequences = (meta[sequence_column] if sequence_column
                 else [None] * len(meta))
    rows = []
    for smiles, sequence in zip(meta[smiles_column], sequences):
        row = structure_descriptors(smiles)
        if sequence_column:
            row.update(sequence_descriptors(sequence))
        rows.append(row)
    frame = pd.DataFrame(rows, index=meta.index).astype(float)
    constant = [c for c in frame.columns if frame[c].nunique(dropna=True) <= 1]
    return frame.drop(columns=constant)


def matrix(meta: pd.DataFrame, **kwargs) -> tuple[np.ndarray, list[str]]:
    """The block as a float32 matrix plus its column names, for the feature layout."""
    frame = build(meta, **kwargs)
    return frame.to_numpy(np.float32), list(frame.columns)
