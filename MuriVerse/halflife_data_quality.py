"""Data-quality metadata for the MuriVerse rodent half-life model.

This module makes label evidence, structure fidelity, and observation identity explicit.
A/B labels require an exact PEPlife2-ID link to curated primary-source PK evidence. C is
unverified database evidence and D is conflicting or unresolved primary-source evidence.
Only an exact whole-molecule structure may feed whole-molecule model feature blocks.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem


STRUCTURE_FIDELITY = (
    "exact_whole_molecule",
    "exact_peptide_partial_modifier",
    "peptide_core_only",
    "unresolved",
)

LEGACY_STRUCTURE_MAP = {
    "complete": "exact_whole_molecule",
    "whole dosed molecule": "exact_whole_molecule",
    "partly resolved": "exact_peptide_partial_modifier",
    "partial": "exact_peptide_partial_modifier",
    "modification not drawn": "peptide_core_only",
    "partner carried": "peptide_core_only",
    "peptide core only": "peptide_core_only",
    "sequence only": "unresolved",
    "none": "unresolved",
    "": "unresolved",
}

LABEL_QUALITY = ("A", "B", "C", "D")


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return " ".join(str(value).strip().split())


def _stable_id(prefix: str, *parts: object) -> str:
    payload = "\x1f".join(_text(p) for p in parts)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]
    return f"{prefix}_{digest}"


def _canonical_smiles(value: object) -> str:
    text = _text(value)
    if not text:
        return ""
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def add_structure_fidelity(meta: pd.DataFrame) -> pd.DataFrame:
    """Add explicit structure fidelity and a model-safe whole-molecule SMILES column."""
    out = meta.copy()
    legacy = out.get("structure_tier", pd.Series("", index=out.index)).map(_text).str.lower()
    fidelity = legacy.map(LEGACY_STRUCTURE_MAP).fillna("unresolved")

    smiles = out.get("smiles", pd.Series("", index=out.index)).map(_text)
    status = out.get("smiles_status", pd.Series("", index=out.index)).map(_text).str.lower()
    missing = smiles.eq("") | status.isin({"none", "failed", "invalid"})
    fidelity = fidelity.mask(missing, "unresolved")

    out["structure_fidelity"] = pd.Categorical(
        fidelity, categories=STRUCTURE_FIDELITY, ordered=True
    ).astype(str)
    out["smiles_whole_molecule_verified"] = out["structure_fidelity"].eq(
        "exact_whole_molecule"
    )
    # Preserve smiles as provenance. Only smiles_for_model may feed whole-molecule blocks.
    out["smiles_for_model"] = out["smiles"].where(
        out["smiles_whole_molecule_verified"], np.nan
    )
    return out


def add_label_quality(meta: pd.DataFrame) -> pd.DataFrame:
    """Assign A/B/C/D evidence tiers without inferring primary-source verification."""
    out = meta.copy()
    evidence = out.get(
        "label_evidence_tier", pd.Series("database_only", index=out.index)
    ).fillna("database_only").astype(str)
    qualifier = out.get(
        "half_life_qualifier", pd.Series("", index=out.index)
    ).fillna("").astype(str).str.lower()

    tier = pd.Series("C", index=out.index, dtype=object)
    reason = pd.Series(
        "database transcription; no exact arm-level primary-source match",
        index=out.index, dtype=object
    )

    matched = evidence.eq("full_text_match")
    exact = qualifier.eq("exact")
    tier.loc[matched & exact] = "A"
    reason.loc[matched & exact] = (
        "exact-ID primary-source arm; numeric half-life agrees; exact label"
    )
    tier.loc[matched & ~exact] = "B"
    reason.loc[matched & ~exact] = (
        "exact-ID primary-source arm; numeric half-life agrees; non-exact label"
    )

    conflict = evidence.isin({"full_text_conflict", "full_text_value_unresolved"})
    tier.loc[conflict] = "D"
    reason.loc[evidence.eq("full_text_conflict")] = (
        "primary-source arm conflicts with database half-life"
    )
    reason.loc[evidence.eq("full_text_value_unresolved")] = (
        "primary-source arm identified but half-life value could not be resolved"
    )

    out["label_quality_tier"] = pd.Categorical(
        tier, categories=LABEL_QUALITY, ordered=True
    ).astype(str)
    out["label_quality_reason"] = reason
    out["gold_training_label"] = out["label_quality_tier"].isin(["A", "B"])
    out["exact_label"] = qualifier.eq("exact")
    return out


def add_identifiers(meta: pd.DataFrame) -> pd.DataFrame:
    """Add molecule, experiment, study, arm and homology-cluster identifiers."""
    out = meta.copy()
    source = out.get("PMID", pd.Series("", index=out.index)).map(_text)
    out["study_id"] = [_stable_id("study", x) for x in source]

    molecule_ids = []
    for _, row in out.iterrows():
        if _text(row.get("structure_fidelity")) == "exact_whole_molecule":
            canonical = _canonical_smiles(row.get("smiles_for_model"))
            if canonical:
                molecule_ids.append(_stable_id("mol", "smiles", canonical))
                continue
        # Conservative identity for incomplete structures: never claim peptide-core SMILES
        # is the administered molecule.
        molecule_ids.append(_stable_id(
            "mol_partial",
            row.get("sequence"),
            row.get("n_ter_mod"),
            row.get("c_ter_mod"),
            row.get("topology"),
            row.get("chirality"),
            row.get("hl_engineering"),
            row.get("partner_text"),
            row.get("unresolved"),
            row.get("Name"),
        ))
    out["molecule_id"] = molecule_ids

    record_source = out.get("ID", out.get("peplife_id", out.get("id")))
    record_id = record_source.map(_text)
    out["arm_id"] = [
        _stable_id("arm", study, rid) for study, rid in zip(out["study_id"], record_id)
    ]

    def first_value(row: pd.Series, primary: str, fallback: str = "") -> str:
        value = _text(row.get(primary))
        if value:
            return value
        return _text(row.get(fallback)) if fallback else ""

    experiment_ids = []
    for _, row in out.iterrows():
        experiment_ids.append(_stable_id(
            "exp",
            row.get("study_id"),
            row.get("species"),
            first_value(row, "evidence_route", "route_model"),
            first_value(row, "evidence_matrix", "sample_matrix"),
            first_value(row, "evidence_bioanalytical_method", "assay_method"),
            first_value(row, "evidence_strain", "strain"),
            first_value(row, "evidence_sex", "sex"),
            first_value(row, "evidence_dose", "pk_dose_text"),
            first_value(row, "evidence_formulation", "formulation_class"),
            row.get("evidence_sampling_design"),
        ))
    out["experiment_id"] = experiment_ids
    out["homology_cluster"] = out.get(
        "cluster", pd.Series("", index=out.index)
    ).map(_text)
    return out


def add_priority_source_rank(meta: pd.DataFrame, n: int = 10) -> pd.DataFrame:
    """Flag the largest sources in the current model table."""
    out = meta.copy()
    counts = (
        out.assign(_source=out["PMID"].astype(str))
        .groupby("_source", dropna=False)
        .size()
        .rename("records")
        .reset_index()
        .sort_values(["records", "_source"], ascending=[False, True], kind="stable")
        .reset_index(drop=True)
    )
    counts["priority_source_rank"] = np.arange(1, len(counts) + 1)
    rank = counts.set_index("_source")["priority_source_rank"]
    out["priority_source_rank"] = out["PMID"].astype(str).map(rank).astype(int)
    out["priority_source_top10"] = out["priority_source_rank"].le(n)
    return out


def prepare_training_metadata(meta: pd.DataFrame, priority_n: int = 10) -> pd.DataFrame:
    out = add_structure_fidelity(meta)
    out = add_label_quality(out)
    out = add_identifiers(out)
    out = add_priority_source_rank(out, priority_n)
    return out


def priority_source_summary(meta: pd.DataFrame, n: int = 10) -> list[dict]:
    """JSON-serialisable curation status for the current largest sources."""
    rows = meta[meta["priority_source_rank"].le(n)].copy()
    result = []
    for source, group in rows.groupby("PMID", sort=False):
        id_column = "peplife_id" if "peplife_id" in group else "ID"
        matched = group.get(
            "pk_observation_matched", pd.Series(False, index=group.index)
        ).fillna(False).astype(bool)
        result.append({
            "source": str(source),
            "rank": int(group["priority_source_rank"].iloc[0]),
            "records": int(len(group)),
            "label_quality": {
                str(k): int(v)
                for k, v in group["label_quality_tier"].value_counts().sort_index().items()
            },
            "exact_labels": int(group["exact_label"].sum()),
            "gold_labels": int(group["gold_training_label"].sum()),
            "arm_level_full_text_links": int(matched.sum()),
            "structure_fidelity": {
                str(k): int(v)
                for k, v in group["structure_fidelity"].value_counts().sort_index().items()
            },
            "remaining_curation_ids": [
                str(x) for x in group.loc[~group["gold_training_label"], id_column].tolist()
            ],
        })
    return sorted(result, key=lambda x: x["rank"])


def write_quality_manifest(meta: pd.DataFrame, path: Path, priority_n: int = 10) -> None:
    """Write a small non-LFS JSON status file for review and reproducibility."""
    payload = {
        "rows": int(len(meta)),
        "label_quality_counts": {
            str(k): int(v)
            for k, v in meta["label_quality_tier"].value_counts().sort_index().items()
        },
        "exact_label_counts": {
            str(k).lower(): int(v)
            for k, v in meta["exact_label"].value_counts().sort_index().items()
        },
        "structure_fidelity_counts": {
            str(k): int(v)
            for k, v in meta["structure_fidelity"].value_counts().sort_index().items()
        },
        "priority_sources": priority_source_summary(meta, priority_n),
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
