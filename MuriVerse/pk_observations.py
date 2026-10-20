"""Load arm-level PK evidence and join only by explicit PEPlife2 record ID.

The publication-level `agents` list summarizes molecules. `pk_observations` carries the
species-by-treatment-arm outcome and conditions that belong to individual source records.
No sequence/name fuzzy join is performed here.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pandas as pd


TABLE_COLUMNS = [
    "peplife2_id", "PMID", "species", "agent_name", "strain", "sex", "age",
    "body_weight", "group_n", "dose", "route", "formulation", "formulation_class",
    "matrix", "sampling_design", "bioanalytical_method", "half_life",
    "half_life_phase", "cmax", "tmax", "auc", "clearance",
    "volume_of_distribution", "mrt", "bioavailability", "evidence",
]


def _read_extractions(root: Path) -> list[dict]:
    rows: dict[str, dict] = {}
    for species in ("mouse", "rat"):
        folder = root / species / "extractions"
        if not folder.exists():
            continue
        for path in sorted(folder.glob("*.json")):
            extraction = json.loads(path.read_text())
            pmid = str(extraction.get("pmid") or path.stem)
            for observation in extraction.get("pk_observations", []):
                record_id = str(observation.get("peplife2_id", "")).strip()
                if not record_id:
                    raise ValueError(f"{path}: every pk_observation needs a PEPlife2 ID")
                row = {column: observation.get(column, "") for column in TABLE_COLUMNS}
                row["PMID"] = pmid
                row["evidence"] = json.dumps(observation.get("evidence", []), ensure_ascii=False)
                if record_id in rows and rows[record_id] != row:
                    raise ValueError(f"Conflicting PK evidence for PEPlife2 ID {record_id}")
                rows[record_id] = row
    return list(rows.values())


def write_table(root: Path, path: Path | None = None) -> pd.DataFrame:
    """Write and return the long-form PK evidence table."""
    rows = _read_extractions(root)
    table = pd.DataFrame(rows, columns=TABLE_COLUMNS)
    if path is None:
        path = root / "rodent" / "pk_observations.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    return table


def _parse_dose(text: object) -> tuple[float, str]:
    value = str(text or "").strip().lower().replace("μ", "µ")
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*(mg|µg|ug|mcg|ng|g)\s*/\s*kg\s*", value)
    if match:
        number, unit = float(match.group(1)), match.group(2)
        factor = {"g": 1000.0, "mg": 1.0, "µg": 0.001, "ug": 0.001,
                  "mcg": 0.001, "ng": 0.000001}[unit]
        return number * factor, "mg/kg"
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*(nmol|µmol|umol|pmol)\s*/\s*kg\s*", value)
    if match:
        return float(match.group(1)), match.group(2) + "/kg"
    return math.nan, "not reported"


def _parse_scalar(text: object, pattern: str, scale: float = 1.0) -> float:
    match = re.fullmatch(pattern, str(text or "").strip(), re.I)
    return float(match.group(1)) * scale if match else math.nan


def _dose_features(text: object) -> tuple[float, str]:
    amount, unit = _parse_dose(text)
    return (math.log10(amount) if unit == "mg/kg" and amount > 0 else math.nan, unit)


def _age_years(text: object) -> float:
    value = str(text or "").strip().lower()
    match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*(days?|d|weeks?|wk|wks|months?|years?|yrs?)", value)
    if not match:
        return math.nan
    number, unit = float(match.group(1)), match.group(2)
    if unit.startswith("d"):
        return number / 365.25
    if unit.startswith(("w", "wk")):
        return number / 52.18
    if unit.startswith("month"):
        return number / 12
    return number


def _half_life_hours(text: object) -> float:
    value = str(text or "").strip().lower().replace("μ", "µ")
    match = re.search(r"(?<![\w.])(\d+(?:\.\d+)?)\s*(seconds?|secs?|s|min(?:ute)?s?|h(?:ours?)?|days?)\b", value)
    if not match:
        return math.nan
    number, unit = float(match.group(1)), match.group(2)
    if unit.startswith("s"):
        return number / 3600
    if unit.startswith("min"):
        return number / 60
    if unit.startswith("d"):
        return number * 24
    return number


def _route(value: object) -> str:
    text = str(value or "").strip().lower()
    for result, pattern in (("iv", r"intraven|\biv\b"), ("sc", r"subcutan|\bsc\b"),
                            ("ip", r"intraperiton|\bip\b"), ("oral", r"\boral|gavage"),
                            ("intranasal", r"intranas"), ("im", r"intramuscular|\bim\b")):
        if re.search(pattern, text):
            return result
    return ""


def _assay(value: object) -> str:
    text = str(value or "").lower()
    if re.search(r"lc.?ms|mass spect", text):
        return "LC-MS"
    if re.search(r"radioimmuno|\bria\b", text):
        return "radioimmunoassay"
    if re.search(r"elisa|immunoassay", text):
        return "ELISA / immunoassay"
    if "hplc" in text:
        return "HPLC"
    if re.search(r"radio|gamma|scintill|pet|fluoresc|imaging", text):
        return "radiometric / imaging"
    return ""


def enrich_training_rows(meta: pd.DataFrame, table: pd.DataFrame) -> pd.DataFrame:
    """Add exact-ID evidence, QC status and explicit arm conditions to model metadata."""
    out = meta.copy()
    if table.empty:
        merged = out
    else:
        table = table.copy()
        if table["peplife2_id"].duplicated().any():
            raise ValueError("A PEPlife2 record ID maps to multiple PK observations")
        table["peplife2_id"] = table["peplife2_id"].astype(str)
        table = table.rename(columns={
            column: f"evidence_{column}" for column in table.columns
            if column != "peplife2_id"
        })
        out["ID"] = out["ID"].astype(str)
        merged = out.merge(table, how="left", left_on="ID", right_on="peplife2_id",
                           validate="one_to_one", suffixes=("", "_evidence"))

    # Preserve database-derived context beside any full-text value. Body evidence is used
    # when an exact source ID has been curated; disagreements are surfaced, never hidden.
    for field in ("route_model", "sample_matrix", "assay_method", "strain", "sex"):
        if field not in merged:
            merged[field] = "unknown"
        merged[f"database_{field}"] = merged[field]

    if not table.empty:
        for index, row in merged.iterrows():
            if pd.isna(row.get("peplife2_id")):
                continue
            conflicts = []
            if str(row.get("species")) != str(row.get("evidence_species", row.get("species"))):
                conflicts.append("species")
            explicit = {
                "route_model": _route(row.get("evidence_route")),
                "sample_matrix": str(row.get("evidence_matrix") or "").lower(),
                "assay_method": _assay(row.get("evidence_bioanalytical_method")),
                "strain": str(row.get("evidence_strain") or "").strip(),
                "sex": str(row.get("evidence_sex") or "").strip().lower(),
            }
            for field, value in explicit.items():
                if value:
                    old = str(row.get(f"database_{field}", "")).strip()
                    if old and old.lower() not in ("unknown", "not reported", value.lower()):
                        conflicts.append(field)
                    merged.at[index, field] = value
            merged.at[index, "pk_condition_conflicts"] = ";".join(sorted(set(conflicts)))

    if "formulation_class" not in merged:
        merged["formulation_class"] = "not reported"
    if "evidence_formulation_class" in merged:
        has_evidence_form = merged["evidence_formulation_class"].fillna("").ne("")
        merged.loc[has_evidence_form, "formulation_class"] = merged.loc[
            has_evidence_form, "evidence_formulation_class"]
    merged["formulation_class"] = merged["formulation_class"].fillna("not reported").replace("", "not reported")
    merged["pk_dose_text"] = merged.get("evidence_dose", pd.Series("", index=merged.index)).fillna("")
    parsed = merged["pk_dose_text"].map(_dose_features)
    merged["pk_dose_log10_mg_kg"] = [x[0] for x in parsed]
    merged["pk_dose_unit"] = [x[1] for x in parsed]
    merged["pk_age_years"] = merged.get("evidence_age", pd.Series("", index=merged.index)).map(_age_years)
    merged["pk_body_weight_g"] = merged.get("evidence_body_weight", pd.Series("", index=merged.index)).map(
        lambda s: _parse_scalar(s, r"(\d+(?:\.\d+)?)\s*g", 1.0))
    merged["pk_group_n"] = merged.get("evidence_group_n", pd.Series("", index=merged.index)).map(
        lambda s: _parse_scalar(s, r"n\s*=\s*(\d+)(?:\s+.*)?"))
    merged["pk_observation_matched"] = merged.get(
        "peplife2_id", pd.Series(index=merged.index, dtype=object)).notna()
    merged["pk_label_reported"] = merged.get(
        "evidence_half_life", pd.Series("", index=merged.index)).fillna("")
    merged["pk_evidence_json"] = merged.get(
        "evidence_evidence", pd.Series("", index=merged.index)).fillna("")
    merged["pk_label_reported_h"] = merged["pk_label_reported"].map(_half_life_hours)
    numeric = pd.to_numeric(merged["half_life_h"], errors="coerce")
    delta = (numeric - merged["pk_label_reported_h"]).abs()
    comparable = merged["pk_label_reported_h"].notna() & numeric.notna()
    merged["label_evidence_tier"] = "database_only"
    merged.loc[comparable & (delta <= numeric.abs().clip(lower=0.01) * 0.02),
               "label_evidence_tier"] = "full_text_match"
    merged.loc[comparable & (delta > numeric.abs().clip(lower=0.01) * 0.02),
               "label_evidence_tier"] = "full_text_conflict"
    merged.loc[merged["pk_observation_matched"] & ~comparable,
               "label_evidence_tier"] = "full_text_value_unresolved"
    if "half_life_provenance" in merged:
        merged.loc[merged["label_evidence_tier"] == "full_text_match",
                   "half_life_provenance"] = "matches curated full text"
        merged.loc[merged["label_evidence_tier"] == "full_text_conflict",
                   "half_life_provenance"] = "conflicts with curated full text"
    return merged
