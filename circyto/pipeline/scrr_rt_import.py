from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from circyto import __version__
from circyto.pipeline.scrr_cnv_import import (
    canonical_scrr_cell_id,
    infer_paired_rna_cell_id,
)
from circyto.pipeline.workflow_reporting import (
    apply_standard_provenance,
    sanitize_frame_for_anndata,
    utc_now_iso,
    write_json,
)

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


CELL_PREFIXES = ("DNA_", "RNA_")
MAPPING_COLUMNS = (
    "gsm_id",
    "rna_cell_id",
    "dna_cell_id",
    "canonical_cell_id",
    "sample_title",
    "molecule",
    "treatment",
    "source_name",
)
GENE_ID_CANDIDATES = (
    "feature_id",
    "bin_id",
    "gene_id",
    "Geneid",
    "gene",
    "gene_name",
    "name",
)
CHROM_CANDIDATES = ("seqname", "chrom", "chromosome", "chr")


def infer_paired_dna_cell_id(cell_id: str) -> str:
    value = str(cell_id).strip()
    if value.startswith("RNA_"):
        return "DNA_" + value[len("RNA_") :]
    return value


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


def _validate_unique(values: list[str], *, label: str) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value in seen:
            duplicates.append(value)
        seen.add(value)
    if duplicates:
        preview = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(f"Duplicate {label} identifiers: {preview}")


def _read_wide_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", keep_default_na=False, compression="infer")
    df = df.rename(columns={column: str(column).strip() for column in df.columns})
    if df.empty:
        raise ValueError(f"{path} is empty.")
    return df


def _is_cell_column(column: object) -> bool:
    value = str(column).strip()
    return any(value.startswith(prefix) for prefix in CELL_PREFIXES)


def _detect_cell_columns(df: pd.DataFrame) -> list[str]:
    cell_columns = [str(column) for column in df.columns if _is_cell_column(column)]
    if not cell_columns:
        raise ValueError(
            "Could not identify scRR RT cell columns. Expected sample columns beginning with DNA_ or RNA_."
        )
    _validate_unique(cell_columns, label="RT cell")
    return cell_columns


def _chrom_column(columns: list[str]) -> str | None:
    lookup = {column.lower(): column for column in columns}
    for candidate in CHROM_CANDIDATES:
        match = lookup.get(candidate.lower())
        if match is not None:
            return match
    return None


def _feature_type_from_metadata(metadata: pd.DataFrame) -> str:
    columns = [str(column) for column in metadata.columns]
    has_coords = _chrom_column(columns) is not None and "start" in columns and "end" in columns
    has_gene = any(column.lower() in {"gene_id", "geneid", "gene", "gene_name"} for column in columns)
    if has_coords and has_gene:
        return "gene_intersect_genomic_feature"
    if has_coords:
        return "genomic_bin"
    if has_gene:
        return "gene_intersect_feature"
    return "feature"


def _make_unique_feature_ids(feature_ids: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for idx, feature_id in enumerate(feature_ids):
        value = feature_id or f"feature_{idx + 1}"
        count = seen.get(value, 0)
        seen[value] = count + 1
        if count:
            out.append(f"{value}_{count + 1}")
        else:
            out.append(value)
    return out


def _build_feature_frame(df: pd.DataFrame, cell_columns: list[str]) -> pd.DataFrame:
    metadata_columns = [str(column) for column in df.columns if str(column) not in set(cell_columns)]
    metadata = df.loc[:, metadata_columns].copy()
    metadata = metadata.rename(columns={column: str(column).strip() for column in metadata.columns})
    feature_type = _feature_type_from_metadata(metadata)
    chrom_col = _chrom_column([str(column) for column in metadata.columns])
    if chrom_col is not None and "start" in metadata.columns and "end" in metadata.columns:
        metadata["seqname"] = metadata[chrom_col].astype(str)
        metadata["start"] = pd.to_numeric(metadata["start"], errors="raise").astype(np.int64)
        metadata["end"] = pd.to_numeric(metadata["end"], errors="raise").astype(np.int64)
        if (metadata["end"] <= metadata["start"]).any():
            raise ValueError("RT features with genomic coordinates must have end > start.")
        metadata["bin_size"] = metadata["end"] - metadata["start"]
        raw_feature_ids = (
            metadata["seqname"]
            + ":"
            + metadata["start"].astype(str)
            + "-"
            + metadata["end"].astype(str)
        ).astype(str).tolist()
    else:
        chosen = None
        for candidate in GENE_ID_CANDIDATES:
            if candidate in metadata.columns:
                chosen = candidate
                break
        if chosen is None and metadata_columns:
            chosen = metadata_columns[0]
        if chosen is None:
            raw_feature_ids = [f"feature_{idx + 1}" for idx in range(df.shape[0])]
        else:
            raw_feature_ids = metadata[chosen].astype(str).str.strip().tolist()

    feature_ids = _make_unique_feature_ids(raw_feature_ids)
    var = metadata.copy()
    var["feature_id"] = feature_ids
    var["feature_type"] = feature_type
    var = var.set_index("feature_id", drop=False)
    var.index.name = None
    return var


def _matrix_from_cells(df: pd.DataFrame, cell_columns: list[str]) -> tuple[np.ndarray, dict[str, Any]]:
    raw = df.loc[:, cell_columns]
    try:
        numeric = raw.apply(pd.to_numeric, errors="raise")
        values = numeric.to_numpy(dtype=np.float32).T
        unique_values = sorted({float(value) for value in values.ravel().tolist()})
        integral = all(float(value).is_integer() for value in unique_values)
        if set(unique_values).issubset({0.0, 1.0}):
            return values.astype(np.int8), {
                "value_semantics": "binary replication state",
                "X_dtype": "int8",
                "observed_values": [int(value) for value in unique_values],
                "category_encoding": None,
            }
        if integral:
            return values.astype(np.int16), {
                "value_semantics": "integer encoded replication state",
                "X_dtype": "int16",
                "observed_values": [int(value) for value in unique_values[:50]],
                "category_encoding": None,
            }
        return values.astype(np.float32), {
            "value_semantics": "numeric replication timing/state signal",
            "X_dtype": "float32",
            "observed_values": [float(value) for value in unique_values[:50]],
            "category_encoding": None,
        }
    except Exception:
        categorical = raw.astype(str)
        categories = sorted({value for value in categorical.to_numpy().ravel().tolist() if value != ""})
        if not categories:
            raise ValueError("RT table has no numeric or categorical replication-state values.")
        encoding = {category: idx for idx, category in enumerate(categories)}
        encoded = categorical.apply(lambda column: column.map(encoding))
        if encoded.isna().any().any():
            raise ValueError("RT table contains empty categorical replication-state values.")
        return encoded.to_numpy(dtype=np.int16).T, {
            "value_semantics": "categorical replication state encoded as integer labels",
            "X_dtype": "int16",
            "observed_values": categories[:50],
            "category_encoding": encoding,
        }


def _load_cell_mapping(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    missing = [column for column in MAPPING_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required mapping columns: {', '.join(missing)}")
    mapping: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        payload = {column: str(row[column]).strip() for column in MAPPING_COLUMNS}
        for key in ("dna_cell_id", "rna_cell_id", "sample_title"):
            value = payload.get(key, "")
            if value:
                if value in mapping and mapping[value] != payload:
                    raise ValueError(f"Duplicate conflicting mapping entry for {value}.")
                mapping[value] = payload
    return mapping


def _source_molecule(cell_id: str) -> str:
    value = str(cell_id).strip()
    if value.startswith("DNA_"):
        return "DNA"
    if value.startswith("RNA_"):
        return "RNA"
    return ""


def _build_obs(
    cell_columns: list[str],
    *,
    cell_mapping_path: Path | None,
    obs_id_strategy: str,
) -> tuple[pd.DataFrame, int]:
    mapping = _load_cell_mapping(cell_mapping_path)
    rows: list[dict[str, str]] = []
    mapped_count = 0
    for source_cell_id in cell_columns:
        mapped = mapping.get(source_cell_id)
        source_molecule = _source_molecule(source_cell_id)
        if mapped is not None:
            mapped_count += 1
            canonical_cell_id = mapped["canonical_cell_id"] or canonical_scrr_cell_id(source_cell_id)
            dna_cell_id = mapped["dna_cell_id"] or infer_paired_dna_cell_id(source_cell_id)
            rna_cell_id = mapped["rna_cell_id"] or infer_paired_rna_cell_id(source_cell_id)
            pairing_strategy = "cell_mapping"
        else:
            canonical_cell_id = canonical_scrr_cell_id(source_cell_id)
            dna_cell_id = source_cell_id if source_cell_id.startswith("DNA_") else infer_paired_dna_cell_id(source_cell_id)
            rna_cell_id = source_cell_id if source_cell_id.startswith("RNA_") else infer_paired_rna_cell_id(source_cell_id)
            pairing_strategy = "strip_modality_prefix"

        if obs_id_strategy == "canonical":
            cell_id = canonical_cell_id
        elif obs_id_strategy == "dna":
            cell_id = dna_cell_id
        elif obs_id_strategy == "rna":
            cell_id = rna_cell_id
        else:
            raise ValueError(f"Unsupported obs_id_strategy: {obs_id_strategy}")

        rows.append(
            {
                "cell_id": cell_id,
                "canonical_cell_id": canonical_cell_id,
                "dna_cell_id": dna_cell_id,
                "rna_cell_id": rna_cell_id,
                "source_cell_id": source_cell_id,
                "source_molecule": source_molecule,
                "pairing_strategy": pairing_strategy,
            }
        )

    obs = pd.DataFrame(rows)
    _validate_unique(obs["cell_id"].astype(str).tolist(), label="RT obs")
    obs = obs.set_index("cell_id", drop=False)
    obs.index.name = None
    return obs, mapped_count


def _read_avg_rt_bedgraph(path: Path) -> pd.DataFrame:
    bedgraph = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["seqname", "start", "end", "avg_rt"],
        comment="#",
        compression="infer",
        keep_default_na=False,
    )
    bedgraph = bedgraph[~bedgraph["seqname"].astype(str).str.startswith(("track", "browser"))].copy()
    if bedgraph.empty:
        raise ValueError(f"{path} has no bedGraph rows.")
    bedgraph["seqname"] = bedgraph["seqname"].astype(str)
    bedgraph["start"] = pd.to_numeric(bedgraph["start"], errors="raise").astype(np.int64)
    bedgraph["end"] = pd.to_numeric(bedgraph["end"], errors="raise").astype(np.int64)
    bedgraph["avg_rt"] = pd.to_numeric(bedgraph["avg_rt"], errors="raise").astype(np.float32)
    return bedgraph


def _attach_avg_rt_if_matching(var: pd.DataFrame, avg_rt_bedgraph_path: Path | None) -> tuple[pd.DataFrame, dict[str, Any]]:
    summary: dict[str, Any] = {
        "source_avg_rt_bedgraph": str(avg_rt_bedgraph_path.resolve()) if avg_rt_bedgraph_path is not None else None,
        "avg_rt_storage": None,
        "n_avg_rt_rows": 0,
    }
    if avg_rt_bedgraph_path is None:
        return var, summary
    bedgraph = _read_avg_rt_bedgraph(avg_rt_bedgraph_path)
    summary["n_avg_rt_rows"] = int(bedgraph.shape[0])
    if {"seqname", "start", "end"}.issubset(var.columns) and int(var.shape[0]) == int(bedgraph.shape[0]):
        left = var.loc[:, ["seqname", "start", "end"]].reset_index(drop=True)
        right = bedgraph.loc[:, ["seqname", "start", "end"]].reset_index(drop=True)
        if left.equals(right):
            out = var.copy()
            out["avg_rt"] = bedgraph["avg_rt"].to_numpy(dtype=np.float32)
            summary["avg_rt_storage"] = "var['avg_rt']"
            return out, summary
    summary["avg_rt_storage"] = "uns['circyto']['avg_rt_bedgraph'] metadata only"
    return var, summary


def _overlap_counts(modality_obs: dict[str, set[str]]) -> dict[str, int]:
    rna = modality_obs.get("rna", set())
    circ = modality_obs.get("circ", set())
    rt = modality_obs.get("rt", set())
    return {
        "n_rna_obs": int(len(rna)),
        "n_circ_obs": int(len(circ)),
        "n_rt_obs": int(len(rt)),
        "n_rna_circ_overlap": int(len(rna & circ)),
        "n_rna_rt_overlap": int(len(rna & rt)),
        "n_circ_rt_overlap": int(len(circ & rt)),
        "n_trimodal_overlap": int(len(rna & circ & rt)),
        "n_union_obs": int(len(rna | circ | rt)),
    }


def _obs_set(adata_or_mdata: Any) -> set[str]:
    return set(str(value) for value in adata_or_mdata.obs_names.astype(str).tolist())


def import_scrr_rt(
    *,
    rt_table_path: Path,
    outdir: Path,
    avg_rt_bedgraph_path: Path | None = None,
    cell_mapping_path: Path | None = None,
    obs_id_strategy: str = "canonical",
    write_h5ad: bool = True,
) -> dict[str, Any]:
    df = _read_wide_table(rt_table_path)
    cell_columns = _detect_cell_columns(df)
    var_df = _build_feature_frame(df, cell_columns)
    matrix, semantics = _matrix_from_cells(df, cell_columns)
    var_df, avg_rt_summary = _attach_avg_rt_if_matching(var_df, avg_rt_bedgraph_path)
    obs_df, mapped_count = _build_obs(
        cell_columns,
        cell_mapping_path=cell_mapping_path,
        obs_id_strategy=obs_id_strategy,
    )

    outdir.mkdir(parents=True, exist_ok=True)
    cells_out = outdir / "rt_cells.tsv"
    features_out = outdir / "rt_features.tsv"
    obs_df.reset_index(drop=True).to_csv(cells_out, sep="\t", index=False)
    var_df.to_csv(features_out, sep="\t", index=False)

    h5ad_out = outdir / "rt.h5ad"
    h5ad_written = False
    if write_h5ad:
        if not HAS_ANNDATA:
            raise RuntimeError("anndata is required for import-scrr-rt h5ad output")
        adata = ad.AnnData(
            X=matrix,
            obs=sanitize_frame_for_anndata(obs_df),
            var=sanitize_frame_for_anndata(var_df),
        )
        adata.uns["circyto"] = _safe_uns(
            {
                "command_name": "circyto import-scrr-rt",
                "circyto_version": __version__,
                "source_rt_table": str(rt_table_path.resolve()),
                "source_avg_rt_bedgraph": str(avg_rt_bedgraph_path.resolve())
                if avg_rt_bedgraph_path is not None
                else "",
                "cell_mapping": str(cell_mapping_path.resolve()) if cell_mapping_path is not None else "",
                "obs_id_strategy": obs_id_strategy,
                "modality_name": "rt",
                "workflow_type": "scrr-replication-timing-import",
                "dna_modality_interpretation": "scRR-seq DNA replication timing/state profile",
                "X_semantics": semantics["value_semantics"],
                "category_encoding": semantics["category_encoding"],
                "avg_rt_bedgraph": avg_rt_summary,
            }
        )
        adata.write_h5ad(str(h5ad_out))
        h5ad_written = True

    feature_type_counts = {
        str(key): int(value)
        for key, value in var_df["feature_type"].astype(str).value_counts().sort_index().items()
    }
    summary: dict[str, Any] = {
        "description": "Imported processed scRR-seq DNA replication timing/state table into an RT modality.",
        "source_rt_table": str(rt_table_path.resolve()),
        "source_avg_rt_bedgraph": str(avg_rt_bedgraph_path.resolve()) if avg_rt_bedgraph_path is not None else None,
        "cell_mapping": str(cell_mapping_path.resolve()) if cell_mapping_path is not None else None,
        "obs_id_strategy": obs_id_strategy,
        "n_cells": int(len(cell_columns)),
        "n_features": int(var_df.shape[0]),
        "n_cells_with_mapping": int(mapped_count),
        "source_cell_id_prefixes": sorted({str(cell_id).split("_", 1)[0] for cell_id in cell_columns}),
        "feature_type_counts": feature_type_counts,
        "output_rt_cells": str(cells_out.resolve()),
        "output_rt_features": str(features_out.resolve()),
        "output_rt_h5ad": str(h5ad_out.resolve()) if h5ad_written else None,
        "output_summary": str((outdir / "scrr_rt_import_summary.json").resolve()),
        "modality_name": "rt",
        "X_semantics": f"cell x feature {semantics['value_semantics']}",
        "value_semantics": semantics["value_semantics"],
        "X_dtype": semantics["X_dtype"],
        "observed_values": semantics["observed_values"],
        "category_encoding": semantics["category_encoding"],
        "avg_rt_bedgraph": avg_rt_summary,
    }
    summary_path = outdir / "scrr_rt_import_summary.json"
    write_json(summary_path, summary)
    return apply_standard_provenance(
        summary,
        command_name="circyto import-scrr-rt",
        workflow_type="scrr-replication-timing-import",
        protocol="scRR-seq",
        read_layout="processed-geo-summary",
        genome_fasta=None,
        gtf=None,
        detector_backend="scRepli-seq-replication-timing-summary",
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )


def merge_scrr_rt(
    *,
    input_h5mu: Path,
    rt_h5ad: Path,
    output_h5mu: Path,
    summary_json: Path | None = None,
    allow_partial: bool = False,
    overwrite: bool = False,
) -> dict[str, Any]:
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for merge-scrr-rt")
    if not HAS_MUDATA:
        raise RuntimeError("mudata is not installed; install circyto[mudata] or pip install mudata")
    if output_h5mu.exists() and not overwrite:
        raise ValueError(f"Output already exists: {output_h5mu}. Use --overwrite to replace it.")
    if summary_json is None:
        summary_json = output_h5mu.with_suffix(".summary.json")
    if summary_json.exists() and not overwrite:
        raise ValueError(f"Summary already exists: {summary_json}. Use --overwrite to replace it.")

    rna_circ = mu.read_h5mu(str(input_h5mu))
    rt = ad.read_h5ad(str(rt_h5ad))
    required = {"rna", "circ"}
    missing_modalities = sorted(required - set(rna_circ.mod.keys()))
    if missing_modalities:
        raise ValueError(f"{input_h5mu} is missing required modalities: {', '.join(missing_modalities)}")

    modality_obs = {
        "mudata": _obs_set(rna_circ),
        "rna": _obs_set(rna_circ.mod["rna"]),
        "circ": _obs_set(rna_circ.mod["circ"]),
        "rt": _obs_set(rt),
    }
    counts = _overlap_counts(modality_obs)
    all_equal = (
        modality_obs["mudata"]
        == modality_obs["rna"]
        == modality_obs["circ"]
        == modality_obs["rt"]
    )
    if not all_equal and not allow_partial:
        raise ValueError(
            "RNA, circ, and RT obs IDs do not match exactly. "
            f"overlap_counts={counts}. Use --allow-partial to permit partial overlap."
        )

    modalities = {modality: rna_circ.mod[modality].copy() for modality in rna_circ.mod.keys()}
    rt_copy = rt.copy()
    top_obs = rna_circ.obs.copy()
    if all_equal:
        target_order = [str(value) for value in rna_circ.obs_names.astype(str).tolist()]
        for modality, modality_adata in list(modalities.items()):
            modalities[modality] = modality_adata[target_order, :].copy()
        rt_copy = rt_copy[target_order, :].copy()
        rt_obs = rt_copy.obs.copy()
        top_obs = top_obs.reindex(target_order)
        for column in rt_obs.columns:
            if column not in top_obs.columns:
                top_obs[column] = rt_obs[column]
    modalities["rt"] = rt_copy

    merged = mu.MuData(modalities)
    if all_equal:
        merged.obs = sanitize_frame_for_anndata(top_obs)
    merged.uns.update(rna_circ.uns)
    merged.uns["circyto_scrr_rt_trimodal"] = _safe_uns(
        {
            "command_name": "circyto merge-scrr-rt",
            "circyto_version": __version__,
            "source_rna_circ_h5mu": str(input_h5mu.resolve()),
            "source_rt_h5ad": str(rt_h5ad.resolve()),
            "allow_partial": allow_partial,
            "overlap_counts": counts,
            "modality_name": "rt",
        }
    )

    output_h5mu.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5mu(str(output_h5mu))
    modality_shapes = {
        modality: [int(adata.n_obs), int(adata.n_vars)]
        for modality, adata in merged.mod.items()
    }
    summary = {
        "description": "Merged remapped scRR RNA/circ MuData with replication timing/state AnnData into tri-modal MuData.",
        "source_rna_circ_h5mu": str(input_h5mu.resolve()),
        "source_rt_h5ad": str(rt_h5ad.resolve()),
        "output_h5mu": str(output_h5mu.resolve()),
        "output_summary_json": str(summary_json.resolve()),
        "allow_partial": allow_partial,
        "n_obs": int(merged.n_obs),
        "modalities": list(merged.mod.keys()),
        "modality_shapes": modality_shapes,
        "overlap_counts": counts,
        "modality_name": "rt",
    }
    write_json(summary_json, summary)
    return summary
