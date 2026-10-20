"""Build lightweight mouse and rat indexes for the pooled half-life training table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FIELDS = ("metadata_row", "id", "species", "PMID", "cluster", "split",
          "half_life_h", "structure_tier")


def write_index(meta_path: Path, out_dir: Path) -> dict:
    with meta_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = (set(FIELDS) - {"metadata_row"}) | {
            "smiles", "smiles_coverage", "smiles_status"
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing metadata columns: {sorted(required - set(reader.fieldnames or []))}")
        rows = list(reader)

    ids = [row["id"] for row in rows]
    if not rows or any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("Training metadata must have nonempty, unique IDs")
    invalid_labels = [row["id"] for row in rows
                      if not row.get("half_life_h") or
                      not _is_positive_finite(row["half_life_h"])]
    if invalid_labels:
        raise ValueError(f"Training metadata has invalid half-life labels for {invalid_labels[:5]}")
    species = {row["species"] for row in rows}
    if species != {"mouse", "rat"}:
        raise ValueError(f"Expected a pooled mouse and rat set, found {sorted(species)}")
    splits = {row["split"] for row in rows}
    if splits != {"train", "val"}:
        raise ValueError(f"Expected train and val splits, found {sorted(splits)}")
    invalid_structure = [row["id"] for row in rows
                         if row["structure_tier"] not in
                         {"complete", "partly resolved", "modification not drawn"}
                         or not row["smiles"]
                         or row["smiles_status"] != "verified by build_strict_peptide_set"
                         or (row["structure_tier"] == "complete") !=
                         (row["smiles_coverage"] == "complete")]
    if invalid_structure:
        raise ValueError("Strict training rows need a recognized structure tier and a "
                         f"verified peptide SMILES: {invalid_structure[:5]}")
    cluster_splits = defaultdict(set)
    for row in rows:
        cluster_splits[row["cluster"]].add(row["split"])
    leaking = [cluster for cluster, assigned in cluster_splits.items() if len(assigned) > 1]
    if leaking:
        raise ValueError(f"Homology clusters cross train and val: {leaking[:5]}")

    out_dir.mkdir(parents=True, exist_ok=True)
    counts = {}
    for animal in ("mouse", "rat"):
        selected = [(position, row) for position, row in enumerate(rows) if row["species"] == animal]
        with (out_dir / f"{animal}_index.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for position, row in selected:
                writer.writerow({field: position if field == "metadata_row" else row[field]
                                 for field in FIELDS})
        counts[animal] = dict(sorted(Counter(row["split"] for _, row in selected).items()))

    quality_counts = {}
    for field in ("half_life_qualifier", "label_evidence_tier", "repeat_candidate_scope"):
        if field in (reader.fieldnames or []):
            quality_counts[field] = dict(sorted(Counter(
                row.get(field, "") or "not reported" for row in rows).items()))
    for field in ("pk_observation_matched", "half_life_extreme_qc",
                  "exact_duplicate_candidate", "cross_source_repeated_measurement"):
        if field in (reader.fieldnames or []):
            values = Counter(str(row.get(field, "False")).lower() == "true" for row in rows)
            quality_counts[field] = {"true": values[True], "false": values[False]}

    manifest = {
        "canonical_metadata": str(meta_path.relative_to(ROOT.parent)),
        "canonical_metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
        "rows": len(rows),
        "species_split_counts": counts,
        "homology_clusters": len(cluster_splits),
        "quality_counts": quality_counts,
        "index_columns": list(FIELDS),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def _is_positive_finite(value: str) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return number > 0 and number < float("inf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path,
                        default=ROOT / "rodent/strict/halflife_meta_with_split.csv")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "rodent/strict/training_data")
    args = parser.parse_args()
    print(json.dumps(write_index(args.metadata.resolve(), args.output.resolve()), indent=2))
