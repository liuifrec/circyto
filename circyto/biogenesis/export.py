from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import pandas as pd

from circyto import __version__
from circyto.biogenesis.schema import (
    BIOGENESIS_SCHEMA_VERSION,
    validate_biogenesis_bundle,
)


BIOGENESIS_OUTPUT_FILENAMES = {
    "candidates": "circRNA_candidates.v1.tsv",
    "cell_contexts": "cell_contexts.v1.tsv",
    "observations": "cell_circ_observations.v1.tsv",
    "provenance": "biogenesis_provenance.v1.json",
}


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Biogenesis source table not found: {path}")
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def _unique_values(frame: pd.DataFrame, column: str) -> list[str]:
    if column not in frame.columns:
        return []
    return sorted(
        {
            str(value).strip()
            for value in frame[column].tolist()
            if str(value).strip()
        }
    )


def export_biogenesis_bundle(
    *,
    candidates_path: Path,
    cell_contexts_path: Path,
    observations_path: Path,
    outdir: Path,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Validate three source TSVs and export a canonical schema-v1 bundle."""
    bundle = validate_biogenesis_bundle(
        _read_tsv(candidates_path),
        _read_tsv(cell_contexts_path),
        _read_tsv(observations_path),
    )

    outputs = {
        key: outdir / filename
        for key, filename in BIOGENESIS_OUTPUT_FILENAMES.items()
    }
    existing = [path for path in outputs.values() if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Biogenesis export would replace existing files: "
            + ", ".join(str(path) for path in existing)
            + ". Use --overwrite to replace them."
        )

    outdir.mkdir(parents=True, exist_ok=True)
    candidates = bundle.candidates.sort_values("circ_id").reset_index(drop=True)
    contexts = bundle.cell_contexts.sort_values("cell_id").reset_index(drop=True)
    observations = bundle.observations.sort_values(
        ["cell_id", "circ_id"]
    ).reset_index(drop=True)
    candidates.to_csv(outputs["candidates"], sep="\t", index=False)
    contexts.to_csv(outputs["cell_contexts"], sep="\t", index=False)
    observations.to_csv(outputs["observations"], sep="\t", index=False)

    provenance = {
        "command_name": "circyto export-biogenesis",
        "circyto_version": __version__,
        "schema_family": "circyto-biogenesis",
        "schema_versions": {
            "candidates": BIOGENESIS_SCHEMA_VERSION,
            "cell_contexts": BIOGENESIS_SCHEMA_VERSION,
            "observations": BIOGENESIS_SCHEMA_VERSION,
        },
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_files": {
            "candidates": str(candidates_path.resolve()),
            "cell_contexts": str(cell_contexts_path.resolve()),
            "observations": str(observations_path.resolve()),
        },
        "output_files": {
            key: str(path.resolve())
            for key, path in outputs.items()
            if key != "provenance"
        },
        "record_counts": {
            "candidates": int(len(candidates)),
            "cell_contexts": int(len(contexts)),
            "observations": int(len(observations)),
        },
        "provenance_values": {
            "datasets": _unique_values(observations, "dataset_id"),
            "donors": _unique_values(observations, "donor_id"),
            "protocols": _unique_values(observations, "protocol"),
            "detectors": _unique_values(observations, "detector"),
            "reference_genomes": _unique_values(observations, "reference_genome"),
            "reference_annotations": _unique_values(
                observations, "reference_annotation"
            ),
            "workflow_uuids": _unique_values(observations, "workflow_uuid"),
        },
        "label_semantics": {
            "allowed": ["positive", "unlabelled"],
            "non_detection_is_negative": False,
        },
    }
    outputs["provenance"].write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return {
        "outdir": str(outdir.resolve()),
        "schema_version": BIOGENESIS_SCHEMA_VERSION,
        "record_counts": provenance["record_counts"],
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
    }
