from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from circyto import __version__
from circyto.pipeline.scrr_cnv_import import canonical_scrr_cell_id
from circyto.pipeline.workflow_reporting import sanitize_frame_for_anndata, write_json

try:
    import anndata as ad  # type: ignore

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False

try:
    import mudata as mu  # type: ignore

    HAS_MUDATA = True
except Exception:
    HAS_MUDATA = False


SCRR_CELL_MAP_COLUMNS = [
    "gsm_id",
    "rna_cell_id",
    "dna_cell_id",
    "canonical_cell_id",
    "sample_title",
    "molecule",
    "treatment",
    "source_name",
]


def _open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _first(values: Iterable[str] | None) -> str:
    if values is None:
        return ""
    for value in values:
        return str(value).strip()
    return ""


def _soft_value(fields: dict[str, list[str]], key: str) -> str:
    return _first(fields.get(key))


def _extract_treatment(fields: dict[str, list[str]]) -> str:
    for key in ("Sample_treatment_protocol_ch1", "Sample_treatment_protocol"):
        value = _soft_value(fields, key)
        if value:
            return value
    for value in fields.get("Sample_characteristics_ch1", []):
        text = str(value).strip()
        lowered = text.lower()
        if lowered.startswith("treatment:"):
            return text.split(":", 1)[1].strip()
        if lowered.startswith("treatment ="):
            return text.split("=", 1)[1].strip()
    return ""


def _infer_molecule(sample_title: str, fields: dict[str, list[str]]) -> str:
    title = str(sample_title).strip()
    if title.startswith("RNA_"):
        return "RNA"
    if title.startswith("DNA_"):
        return "DNA"
    value = _soft_value(fields, "Sample_molecule_ch1") or _soft_value(fields, "Sample_molecule")
    return value.upper() if value.upper() in {"RNA", "DNA"} else value


def _iter_soft_samples(path: Path) -> list[dict[str, list[str]]]:
    samples: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line.startswith("^SAMPLE"):
                if current is not None:
                    samples.append(current)
                current = {"SAMPLE": [line.split("=", 1)[1].strip() if "=" in line else ""]}
                continue
            if current is None or not line.startswith("!"):
                continue
            if "=" not in line:
                continue
            key, value = line[1:].split("=", 1)
            current.setdefault(key.strip(), []).append(value.strip())
    if current is not None:
        samples.append(current)
    return samples


def _validate_no_duplicates(values: list[str], *, label: str) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if not value:
            continue
        if value in seen:
            duplicates.append(value)
        seen.add(value)
    if duplicates:
        preview = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(f"Duplicate {label}: {preview}")


def _require_columns(df: pd.DataFrame, required: list[str], *, source: Path) -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")


def build_scrr_cell_map(*, soft_path: Path, output_path: Path) -> dict[str, Any]:
    raw_samples = _iter_soft_samples(soft_path)
    rows: list[dict[str, str]] = []
    for fields in raw_samples:
        gsm_id = _soft_value(fields, "Sample_geo_accession") or _first(fields.get("SAMPLE"))
        sample_title = _soft_value(fields, "Sample_title")
        molecule = _infer_molecule(sample_title, fields)
        if molecule not in {"RNA", "DNA"}:
            continue
        canonical_cell_id = canonical_scrr_cell_id(sample_title)
        rows.append(
            {
                "gsm_id": gsm_id,
                "rna_cell_id": sample_title if molecule == "RNA" else "",
                "dna_cell_id": sample_title if molecule == "DNA" else "",
                "canonical_cell_id": canonical_cell_id,
                "sample_title": sample_title,
                "molecule": molecule,
                "treatment": _extract_treatment(fields),
                "source_name": _soft_value(fields, "Sample_source_name_ch1")
                or _soft_value(fields, "Sample_source_name"),
            }
        )

    if not rows:
        raise ValueError(f"No RNA_* or DNA_* scRR sample titles found in {soft_path}.")

    df = pd.DataFrame(rows, columns=SCRR_CELL_MAP_COLUMNS)
    _validate_no_duplicates(df["gsm_id"].astype(str).tolist(), label="GSM IDs")
    _validate_no_duplicates(df["sample_title"].astype(str).tolist(), label="Sample_title values")
    for molecule in ("RNA", "DNA"):
        molecule_df = df[df["molecule"] == molecule]
        _validate_no_duplicates(
            molecule_df["canonical_cell_id"].astype(str).tolist(),
            label=f"{molecule} canonical cell IDs",
        )

    rna_by_canonical = {
        str(row["canonical_cell_id"]): str(row["sample_title"])
        for _, row in df[df["molecule"] == "RNA"].iterrows()
    }
    dna_by_canonical = {
        str(row["canonical_cell_id"]): str(row["sample_title"])
        for _, row in df[df["molecule"] == "DNA"].iterrows()
    }
    df["rna_cell_id"] = [
        str(row["sample_title"])
        if row["molecule"] == "RNA"
        else rna_by_canonical.get(str(row["canonical_cell_id"]), "")
        for _, row in df.iterrows()
    ]
    df["dna_cell_id"] = [
        str(row["sample_title"])
        if row["molecule"] == "DNA"
        else dna_by_canonical.get(str(row["canonical_cell_id"]), "")
        for _, row in df.iterrows()
    ]

    all_canonical = sorted(set(df["canonical_cell_id"].astype(str)))
    rna_only = [cell_id for cell_id in all_canonical if cell_id in rna_by_canonical and cell_id not in dna_by_canonical]
    dna_only = [cell_id for cell_id in all_canonical if cell_id in dna_by_canonical and cell_id not in rna_by_canonical]
    paired = [cell_id for cell_id in all_canonical if cell_id in rna_by_canonical and cell_id in dna_by_canonical]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return {
        "description": "Built scRR GSM-to-biological-cell mapping from GEO SOFT metadata.",
        "source_soft": str(soft_path.resolve()),
        "output_cell_map": str(output_path.resolve()),
        "columns": SCRR_CELL_MAP_COLUMNS,
        "n_rows": int(df.shape[0]),
        "n_rna_rows": int((df["molecule"] == "RNA").sum()),
        "n_dna_rows": int((df["molecule"] == "DNA").sum()),
        "n_paired_canonical_cells": int(len(paired)),
        "n_rna_only_canonical_cells": int(len(rna_only)),
        "n_dna_only_canonical_cells": int(len(dna_only)),
        "paired_canonical_cells": paired,
        "rna_only_canonical_cells": rna_only,
        "dna_only_canonical_cells": dna_only,
    }


def _read_cell_map(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    _require_columns(df, SCRR_CELL_MAP_COLUMNS, source=path)
    _validate_no_duplicates(df["gsm_id"].astype(str).tolist(), label="mapping GSM IDs")
    return df


def _target_column(target_id: str) -> str:
    normalized = str(target_id).strip().lower()
    if normalized in {"canonical", "canonical_cell_id"}:
        return "canonical_cell_id"
    if normalized in {"rna", "rna_cell_id", "rna_title"}:
        return "rna_cell_id"
    raise ValueError("target_id must be canonical_cell_id or rna_cell_id")


def _safe_uns(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _safe_uns(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_safe_uns(item) for item in value]
    if isinstance(value, tuple):
        return [_safe_uns(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _remap_obs_frame(
    obs: pd.DataFrame,
    *,
    mapping_by_gsm: dict[str, dict[str, str]],
    target_column: str,
    allow_partial: bool,
    axis_label: str,
) -> tuple[pd.DataFrame, list[str]]:
    original_ids = [str(value) for value in obs.index.astype(str).tolist()]
    target_ids: list[str] = []
    missing: list[str] = []
    metadata_rows: list[dict[str, str]] = []
    for original_id in original_ids:
        mapped = mapping_by_gsm.get(original_id)
        if mapped is None:
            missing.append(original_id)
            if not allow_partial:
                continue
            target_id = original_id
            metadata_rows.append(
                {
                    "gsm_id": original_id if original_id.startswith("GSM") else "",
                    "canonical_cell_id": "",
                    "rna_cell_id": "",
                    "dna_cell_id": "",
                    "sample_title": "",
                    "molecule": "",
                }
            )
        else:
            target_id = mapped[target_column]
            if not target_id:
                missing.append(original_id)
                if not allow_partial:
                    continue
                target_id = original_id
            metadata_rows.append(
                {
                    "gsm_id": mapped["gsm_id"],
                    "canonical_cell_id": mapped["canonical_cell_id"],
                    "rna_cell_id": mapped["rna_cell_id"],
                    "dna_cell_id": mapped["dna_cell_id"],
                    "sample_title": mapped["sample_title"],
                    "molecule": mapped["molecule"],
                }
            )
        target_ids.append(target_id)

    if missing and not allow_partial:
        preview = ", ".join(missing[:5])
        raise ValueError(f"{axis_label} has obs IDs missing from scRR cell map: {preview}")
    _validate_no_duplicates(target_ids, label=f"{axis_label} remapped obs IDs")

    out = obs.copy()
    out["original_obs_id"] = original_ids
    for column in ("gsm_id", "canonical_cell_id", "rna_cell_id", "dna_cell_id", "sample_title", "molecule"):
        out[column] = [row[column] for row in metadata_rows]
    out.index = target_ids
    out.index.name = None
    return out, missing


def remap_scrr_mudata_obs(
    *,
    input_h5mu: Path,
    cell_map_path: Path,
    output_h5mu: Path,
    target_id: str = "canonical_cell_id",
    allow_partial: bool = False,
    overwrite: bool = False,
) -> dict[str, Any]:
    if not HAS_MUDATA:
        raise RuntimeError("mudata is not installed; install circyto[mudata] or pip install mudata")
    if output_h5mu.exists() and not overwrite:
        raise ValueError(f"Output already exists: {output_h5mu}. Use --overwrite to replace it.")

    target_column = _target_column(target_id)
    mapping_df = _read_cell_map(cell_map_path)
    mapping_by_gsm = {
        str(row["gsm_id"]): {column: str(row[column]) for column in SCRR_CELL_MAP_COLUMNS}
        for _, row in mapping_df.iterrows()
    }

    source = mu.read_h5mu(str(input_h5mu))
    remapped_modalities = {}
    missing_by_axis: dict[str, list[str]] = {}
    for modality, modality_adata in source.mod.items():
        copied = modality_adata.copy()
        remapped_obs, missing = _remap_obs_frame(
            copied.obs,
            mapping_by_gsm=mapping_by_gsm,
            target_column=target_column,
            allow_partial=allow_partial,
            axis_label=f"{modality} modality",
        )
        remapped_obs[f"original_{modality}_obs_id"] = remapped_obs["original_obs_id"].astype(str)
        copied.obs = sanitize_frame_for_anndata(remapped_obs)
        copied.obs_names = list(remapped_obs.index.astype(str))
        remapped_modalities[str(modality)] = copied
        missing_by_axis[str(modality)] = missing

    top_obs, top_missing = _remap_obs_frame(
        source.obs,
        mapping_by_gsm=mapping_by_gsm,
        target_column=target_column,
        allow_partial=allow_partial,
        axis_label="MuData",
    )
    missing_by_axis["mudata"] = top_missing

    remapped = mu.MuData(remapped_modalities)
    remapped.obs = sanitize_frame_for_anndata(top_obs.reindex(remapped.obs_names))
    remapped.uns.update(source.uns)
    remapped.uns["circyto_scrr_remap"] = _safe_uns(
        {
            "command_name": "circyto remap-scrr-mudata-obs",
            "circyto_version": __version__,
            "source_h5mu": str(input_h5mu.resolve()),
            "cell_map": str(cell_map_path.resolve()),
            "target_id": target_column,
            "allow_partial": allow_partial,
            "missing_by_axis": missing_by_axis,
        }
    )

    output_h5mu.parent.mkdir(parents=True, exist_ok=True)
    remapped.write_h5mu(str(output_h5mu))
    modality_shapes = {
        modality: [int(adata.n_obs), int(adata.n_vars)]
        for modality, adata in remapped.mod.items()
    }
    return {
        "description": "Remapped scRR RNA/circ MuData obs IDs through GEO GSM metadata.",
        "source_h5mu": str(input_h5mu.resolve()),
        "cell_map": str(cell_map_path.resolve()),
        "output_h5mu": str(output_h5mu.resolve()),
        "target_id": target_column,
        "allow_partial": allow_partial,
        "n_obs": int(remapped.n_obs),
        "modalities": list(remapped.mod.keys()),
        "modality_shapes": modality_shapes,
        "missing_by_axis": missing_by_axis,
    }


def _obs_set(adata_or_mdata: Any) -> set[str]:
    return set(str(value) for value in adata_or_mdata.obs_names.astype(str).tolist())


def _overlap_counts(modality_obs: dict[str, set[str]]) -> dict[str, int]:
    rna = modality_obs.get("rna", set())
    circ = modality_obs.get("circ", set())
    cnv = modality_obs.get("cnv", set())
    return {
        "n_rna_obs": int(len(rna)),
        "n_circ_obs": int(len(circ)),
        "n_cnv_obs": int(len(cnv)),
        "n_rna_circ_overlap": int(len(rna & circ)),
        "n_rna_cnv_overlap": int(len(rna & cnv)),
        "n_circ_cnv_overlap": int(len(circ & cnv)),
        "n_trimodal_overlap": int(len(rna & circ & cnv)),
        "n_union_obs": int(len(rna | circ | cnv)),
    }


def merge_scrr_cnv(
    *,
    input_h5mu: Path,
    cnv_h5ad: Path,
    output_h5mu: Path,
    summary_json: Path | None = None,
    allow_partial: bool = False,
    overwrite: bool = False,
) -> dict[str, Any]:
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for merge-scrr-cnv")
    if not HAS_MUDATA:
        raise RuntimeError("mudata is not installed; install circyto[mudata] or pip install mudata")
    if output_h5mu.exists() and not overwrite:
        raise ValueError(f"Output already exists: {output_h5mu}. Use --overwrite to replace it.")
    if summary_json is None:
        summary_json = output_h5mu.with_suffix(".summary.json")
    if summary_json.exists() and not overwrite:
        raise ValueError(f"Summary already exists: {summary_json}. Use --overwrite to replace it.")

    rna_circ = mu.read_h5mu(str(input_h5mu))
    cnv = ad.read_h5ad(str(cnv_h5ad))
    required = {"rna", "circ"}
    missing_modalities = sorted(required - set(rna_circ.mod.keys()))
    if missing_modalities:
        raise ValueError(f"{input_h5mu} is missing required modalities: {', '.join(missing_modalities)}")

    modality_obs = {
        "mudata": _obs_set(rna_circ),
        "rna": _obs_set(rna_circ.mod["rna"]),
        "circ": _obs_set(rna_circ.mod["circ"]),
        "cnv": _obs_set(cnv),
    }
    counts = _overlap_counts(modality_obs)
    all_equal = (
        modality_obs["mudata"]
        == modality_obs["rna"]
        == modality_obs["circ"]
        == modality_obs["cnv"]
    )
    if not all_equal and not allow_partial:
        raise ValueError(
            "RNA, circ, and CNV obs IDs do not match exactly. "
            f"overlap_counts={counts}. Use --allow-partial to permit partial overlap."
        )

    modalities = {modality: rna_circ.mod[modality].copy() for modality in rna_circ.mod.keys()}
    cnv_copy = cnv.copy()
    top_obs = rna_circ.obs.copy()
    if all_equal:
        target_order = [str(value) for value in rna_circ.obs_names.astype(str).tolist()]
        for modality, modality_adata in list(modalities.items()):
            modalities[modality] = modality_adata[target_order, :].copy()
        cnv_copy = cnv_copy[target_order, :].copy()
        cnv_obs = cnv_copy.obs.copy()
        top_obs = top_obs.reindex(target_order)
        for column in cnv_obs.columns:
            if column not in top_obs.columns:
                top_obs[column] = cnv_obs[column]
    modalities["cnv"] = cnv_copy

    merged = mu.MuData(modalities)
    if all_equal:
        merged.obs = sanitize_frame_for_anndata(top_obs)
    merged.uns.update(rna_circ.uns)
    merged.uns["circyto_scrr_trimodal"] = _safe_uns(
        {
            "command_name": "circyto merge-scrr-cnv",
            "circyto_version": __version__,
            "source_rna_circ_h5mu": str(input_h5mu.resolve()),
            "source_cnv_h5ad": str(cnv_h5ad.resolve()),
            "allow_partial": allow_partial,
            "overlap_counts": counts,
        }
    )

    output_h5mu.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5mu(str(output_h5mu))
    modality_shapes = {
        modality: [int(adata.n_obs), int(adata.n_vars)]
        for modality, adata in merged.mod.items()
    }
    summary = {
        "description": "Merged remapped scRR RNA/circ MuData with CNV AnnData into tri-modal MuData.",
        "source_rna_circ_h5mu": str(input_h5mu.resolve()),
        "source_cnv_h5ad": str(cnv_h5ad.resolve()),
        "output_h5mu": str(output_h5mu.resolve()),
        "output_summary_json": str(summary_json.resolve()),
        "allow_partial": allow_partial,
        "n_obs": int(merged.n_obs),
        "modalities": list(merged.mod.keys()),
        "modality_shapes": modality_shapes,
        "overlap_counts": counts,
    }
    write_json(summary_json, summary)
    return summary
