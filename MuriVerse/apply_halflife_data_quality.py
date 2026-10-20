"""Materialize MuriVerse half-life quality metadata after the normal training-data build."""

from pathlib import Path

import pandas as pd

from halflife_data_quality import prepare_training_metadata, write_quality_manifest

ROOT = Path(__file__).resolve().parent
STRICT = ROOT / "rodent" / "strict"


def main() -> None:
    path = STRICT / "halflife_meta_with_split.csv"
    meta = pd.read_csv(path, low_memory=False)
    meta = prepare_training_metadata(meta)
    write_quality_manifest(meta, STRICT / "data_quality_manifest.json")
    columns = [
        "id", "peplife_id", "PMID", "species",
        "molecule_id", "experiment_id", "study_id", "arm_id", "homology_cluster",
        "label_quality_tier", "label_quality_reason", "gold_training_label", "exact_label",
        "structure_fidelity", "smiles_whole_molecule_verified",
        "priority_source_rank", "priority_source_top10",
    ]
    meta[columns].to_json(
        STRICT / "halflife_quality_metadata.jsonl",
        orient="records", lines=True,
    )


if __name__ == "__main__":
    main()
