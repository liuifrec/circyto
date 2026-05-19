from __future__ import annotations

from pathlib import Path
from statistics import median
from typing import Any, Optional
import json
import os
import socket
import sys
import uuid
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import io as scio
from scipy import sparse as sp

from circyto import __version__


try:
    import anndata as ad

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False

try:
    import mudata as mu

    HAS_MUDATA = True
except Exception:
    HAS_MUDATA = False


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def generate_workflow_uuid() -> str:
    return str(uuid.uuid4())


def summarize_read_layouts(read_layouts: list[str]) -> str:
    normalized = sorted({str(item).strip() for item in read_layouts if str(item).strip()})
    if not normalized:
        return "unknown"
    if len(normalized) == 1:
        return normalized[0]
    return "mixed"


def stage_status_lists(stage_graph: list[dict[str, Any]]) -> dict[str, list[str]]:
    completed: list[str] = []
    skipped: list[str] = []
    failed: list[str] = []
    for item in stage_graph:
        stage = str(item.get("stage", "")).strip()
        status = str(item.get("status", "")).strip().lower()
        if not stage:
            continue
        if status == "completed":
            completed.append(stage)
        elif status in {"skipped", "disabled", "delegated"}:
            skipped.append(stage)
        elif status == "failed":
            failed.append(stage)
    return {
        "completed_stages": completed,
        "skipped_stages": skipped,
        "failed_stages": failed,
    }


def detect_partial_outputs(paths: dict[str, Any]) -> list[str]:
    missing: list[str] = []
    for label, value in dict(paths).items():
        if value in (None, ""):
            continue
        path = Path(str(value))
        if not path.exists():
            missing.append(label)
    return sorted(missing)


def cleanup_summary_block(
    *,
    cleanup: dict[str, Any] | None,
    cleanup_scope: str | None = None,
    cleanup_performed: bool = False,
    cleanup_deleted_paths: list[str] | None = None,
    cleanup_reclaimed_bytes: int = 0,
) -> dict[str, Any]:
    deleted = list(cleanup_deleted_paths or [])
    if cleanup is None:
        return {
            "enabled": cleanup_scope is not None,
            "scope": cleanup_scope,
            "performed": cleanup_performed,
            "deleted_paths": deleted,
            "reclaimed_bytes": int(cleanup_reclaimed_bytes),
        }
    return {
        "enabled": bool(cleanup.get("enabled", cleanup_scope is not None)),
        "scope": cleanup.get("cleanup_scope", cleanup.get("planned_scope", cleanup_scope)),
        "performed": bool(cleanup.get("cleanup_performed", cleanup_performed)),
        "deleted_paths": list(cleanup.get("cleanup_deleted_paths", deleted)),
        "reclaimed_bytes": int(cleanup.get("cleanup_reclaimed_bytes", cleanup_reclaimed_bytes) or 0),
        "candidate_count": int(cleanup.get("candidate_count", 0) or 0),
        "candidate_bytes": int(cleanup.get("candidate_bytes", 0) or 0),
    }


def apply_standard_provenance(
    payload: dict[str, Any],
    *,
    command_name: str,
    workflow_type: str,
    protocol: str,
    read_layout: str,
    genome_fasta: str | None,
    gtf: str | None,
    detector_backend: str | None,
    started_at: str,
    completed_at: str,
    elapsed_seconds: float,
    cleanup: dict[str, Any] | None = None,
    cleanup_scope: str | None = None,
    cleanup_performed: bool = False,
    cleanup_deleted_paths: list[str] | None = None,
    cleanup_reclaimed_bytes: int = 0,
    workflow_uuid: str | None = None,
) -> dict[str, Any]:
    payload["command_name"] = command_name
    payload["circyto_version"] = __version__
    payload["workflow_type"] = workflow_type
    payload["protocol"] = protocol
    payload["read_layout"] = read_layout
    payload["genome_fasta"] = genome_fasta
    payload["gtf"] = gtf
    payload["detector_backend"] = detector_backend
    payload["started_at"] = started_at
    payload["completed_at"] = completed_at
    payload["elapsed_seconds"] = round(float(elapsed_seconds), 3)
    payload["hostname"] = socket.gethostname()
    payload["python_version"] = sys.version.split()[0]
    payload["workflow_uuid"] = workflow_uuid or payload.get("workflow_uuid") or generate_workflow_uuid()
    payload["cleanup_summary"] = cleanup_summary_block(
        cleanup=cleanup,
        cleanup_scope=cleanup_scope,
        cleanup_performed=cleanup_performed,
        cleanup_deleted_paths=cleanup_deleted_paths,
        cleanup_reclaimed_bytes=cleanup_reclaimed_bytes,
    )
    stage_lists = stage_status_lists(list(payload.get("stage_graph", [])))
    payload["completed_stages"] = stage_lists["completed_stages"]
    payload["skipped_stages"] = stage_lists["skipped_stages"]
    payload["failed_stages"] = stage_lists["failed_stages"]
    payload["partial_outputs_detected"] = detect_partial_outputs(dict(payload.get("paths", {})))
    return payload


def read_index_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def numeric_summary(values: list[int | float]) -> dict[str, float | int | None]:
    if not values:
        return {"min": None, "median": None, "mean": None, "max": None}
    vals = [float(value) for value in values]
    return {
        "min": min(vals),
        "median": float(median(vals)),
        "mean": float(sum(vals) / len(vals)),
        "max": max(vals),
    }


def top_mapping_items(
    mapping: dict[str, int],
    *,
    count_key: str,
    limit: int = 10,
) -> list[dict[str, Any]]:
    ranked = sorted(mapping.items(), key=lambda item: (-int(item[1]), item[0]))
    return [{"cell_id": key, count_key: int(value)} for key, value in ranked[:limit]]


def load_circ_matrix(
    *,
    matrix_path: Path,
    circ_index_path: Path,
    cell_index_path: Path,
) -> tuple[sp.csr_matrix, list[str], list[str]]:
    X = scio.mmread(str(matrix_path)).tocsr()
    circ_ids = read_index_lines(circ_index_path)
    cell_ids = read_index_lines(cell_index_path)
    return X, circ_ids, cell_ids


def expand_cells(
    X_cells_by_circ: sp.csr_matrix,
    source_cell_ids: list[str],
    target_cell_ids: list[str],
) -> sp.csr_matrix:
    if not target_cell_ids:
        return sp.csr_matrix((0, X_cells_by_circ.shape[1]), dtype=X_cells_by_circ.dtype)
    source_pos = {cell_id: idx for idx, cell_id in enumerate(source_cell_ids)}
    zero_row = sp.csr_matrix((1, X_cells_by_circ.shape[1]), dtype=X_cells_by_circ.dtype)
    rows: list[sp.csr_matrix] = []
    for cell_id in target_cell_ids:
        idx = source_pos.get(cell_id)
        rows.append(X_cells_by_circ[idx] if idx is not None else zero_row)
    return sp.vstack(rows, format="csr")


def align_matrix_to_cells(
    X_cells_by_feature: sp.csr_matrix,
    source_cell_ids: list[str],
    target_cell_ids: list[str],
    *,
    join: str,
) -> tuple[sp.csr_matrix, list[str], set[str]]:
    source_pos = {cell_id: idx for idx, cell_id in enumerate(source_cell_ids)}
    if join not in {"inner", "outer"}:
        raise ValueError(f"Unsupported cell join: {join}")
    if join == "inner":
        matched = [cell_id for cell_id in target_cell_ids if cell_id in source_pos]
        X_aligned = X_cells_by_feature[[source_pos[cell_id] for cell_id in matched], :] if matched else sp.csr_matrix((0, X_cells_by_feature.shape[1]), dtype=X_cells_by_feature.dtype)
        return X_aligned.tocsr(), matched, set(matched)

    zero_row = sp.csr_matrix((1, X_cells_by_feature.shape[1]), dtype=X_cells_by_feature.dtype)
    rows: list[sp.csr_matrix] = []
    matched_set: set[str] = set()
    for cell_id in target_cell_ids:
        idx = source_pos.get(cell_id)
        if idx is None:
            rows.append(zero_row)
        else:
            rows.append(X_cells_by_feature[idx])
            matched_set.add(cell_id)
    return sp.vstack(rows, format="csr"), list(target_cell_ids), matched_set


def load_circ_feature_table(circ_ids: list[str], feature_path: Path) -> pd.DataFrame:
    default_df = pd.DataFrame(
        {
            "circ_id": circ_ids,
            "chrom": [""] * len(circ_ids),
            "start": [pd.NA] * len(circ_ids),
            "end": [pd.NA] * len(circ_ids),
            "strand": [""] * len(circ_ids),
            "host_gene": [""] * len(circ_ids),
        }
    ).set_index("circ_id")
    if not feature_path.exists():
        return default_df
    df = pd.read_csv(feature_path, sep="\t", keep_default_na=False)
    if "circ_id" not in df.columns:
        return default_df
    if "host_gene" not in df.columns and "gene_name" in df.columns:
        df["host_gene"] = df["gene_name"]
    df = df.set_index("circ_id")
    for column in ("chrom", "start", "end", "strand", "host_gene"):
        if column not in df.columns:
            df[column] = pd.NA if column in {"start", "end"} else ""
    df = df.reindex(circ_ids)[["chrom", "start", "end", "strand", "host_gene"]]
    for column in ("chrom", "strand", "host_gene"):
        df[column] = df[column].fillna("").astype(str)
    for column in ("start", "end"):
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def build_cell_qc_table(
    *,
    selected_cell_ids: list[str],
    assigned_reads: dict[str, int],
    alignment_cells: dict[str, dict[str, Any]],
    detector_cells: dict[str, dict[str, Any]],
    X_cells_by_circ: sp.csr_matrix,
) -> pd.DataFrame:
    circ_counts = np.asarray(X_cells_by_circ.getnnz(axis=1)).ravel().tolist() if X_cells_by_circ.shape[0] else [0] * len(selected_cell_ids)
    total_support = np.asarray(X_cells_by_circ.sum(axis=1)).ravel().tolist() if X_cells_by_circ.shape[0] else [0] * len(selected_cell_ids)
    records: list[dict[str, Any]] = []
    for idx, cell_id in enumerate(selected_cell_ids):
        align_record = alignment_cells.get(cell_id, {})
        detector_record = detector_cells.get(cell_id, {})
        detector_status = str(detector_record.get("status", "missing"))
        total_support_value = int(total_support[idx]) if idx < len(total_support) else 0
        records.append(
            {
                "cell_id": cell_id,
                "assigned_reads": int(assigned_reads.get(cell_id, 0)),
                "alignment_status": str(align_record.get("status", "missing")),
                "detector_status": detector_status,
                "circRNA_count": int(circ_counts[idx]) if idx < len(circ_counts) else 0,
                "total_circRNA_support": total_support_value,
                "empty_flag": bool(total_support_value == 0 or detector_status == "empty"),
                "detector_seconds": detector_record.get("seconds"),
                "alignment_seconds": align_record.get("seconds"),
            }
        )
    return pd.DataFrame.from_records(records).set_index("cell_id")


def build_circ_qc_table(
    *,
    circ_ids: list[str],
    feature_df: pd.DataFrame,
    X_cells_by_circ: sp.csr_matrix,
) -> pd.DataFrame:
    output_columns = [
        "chrom",
        "start",
        "end",
        "strand",
        "host_gene",
        "n_cells_detected",
        "total_support",
        "max_support",
        "mean_support_detected_cells",
    ]
    if not circ_ids:
        return pd.DataFrame(columns=output_columns).rename_axis("circ_id")

    total_support = np.asarray(X_cells_by_circ.sum(axis=0)).ravel()
    n_cells = np.asarray(X_cells_by_circ.getnnz(axis=0)).ravel()
    max_support = np.asarray(X_cells_by_circ.max(axis=0).toarray()).ravel() if X_cells_by_circ.shape[0] else np.zeros(len(circ_ids), dtype=int)
    mean_support = np.divide(total_support, n_cells, out=np.zeros_like(total_support, dtype=float), where=n_cells > 0)
    df = pd.DataFrame(index=pd.Index(circ_ids, name="circ_id"))
    df["chrom"] = feature_df.reindex(circ_ids)["chrom"].fillna("").astype(str)
    df["start"] = pd.to_numeric(feature_df.reindex(circ_ids)["start"], errors="coerce")
    df["end"] = pd.to_numeric(feature_df.reindex(circ_ids)["end"], errors="coerce")
    df["strand"] = feature_df.reindex(circ_ids)["strand"].fillna("").astype(str)
    df["host_gene"] = feature_df.reindex(circ_ids)["host_gene"].fillna("").astype(str)
    df["n_cells_detected"] = n_cells.astype(int)
    df["total_support"] = total_support.astype(int)
    df["max_support"] = max_support.astype(int)
    df["mean_support_detected_cells"] = mean_support.astype(float)
    df.index.name = "circ_id"
    if len(n_cells) and int(n_cells.max()) > X_cells_by_circ.shape[0]:
        raise AssertionError("n_cells_detected cannot exceed the number of cells in the matrix")
    return df[output_columns]


def matrix_section(X_cells_by_circ: sp.csr_matrix, circ_qc: pd.DataFrame) -> dict[str, Any]:
    n_cells, n_circs = X_cells_by_circ.shape
    nnz = int(X_cells_by_circ.nnz)
    cell_nnz = np.asarray(X_cells_by_circ.getnnz(axis=1)).ravel().astype(int).tolist() if n_cells else []
    circ_ncells = np.asarray(X_cells_by_circ.getnnz(axis=0)).ravel().astype(int).tolist() if n_circs else []
    recurrent = circ_qc.reset_index()[["circ_id", "n_cells_detected", "total_support"]].sort_values(
        by=["n_cells_detected", "total_support", "circ_id"],
        ascending=[False, False, True],
    )
    return {
        "n_circRNAs": int(n_circs),
        "n_cells": int(n_cells),
        "nnz": nnz,
        "sparsity": float(1.0 - (nnz / (n_cells * n_circs))) if n_cells and n_circs else 1.0,
        "circRNAs_per_cell_summary": numeric_summary(cell_nnz),
        "cells_per_circRNA_summary": numeric_summary(circ_ncells),
        "top_recurrent_circRNAs": recurrent.head(10).to_dict("records"),
    }


def directory_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for root, _, filenames in os.walk(path):
        for filename in filenames:
            file_path = Path(root) / filename
            try:
                total += file_path.stat().st_size
            except OSError:
                continue
    return int(total)


def largest_files_under(root: Path, *, limit: int = 10) -> list[dict[str, Any]]:
    if not root.exists():
        return []
    ranked: list[tuple[Path, int]] = []
    for walk_root, _, filenames in os.walk(root):
        for filename in filenames:
            file_path = Path(walk_root) / filename
            try:
                size = file_path.stat().st_size
            except OSError:
                continue
            ranked.append((file_path, int(size)))
    ranked.sort(key=lambda item: (-item[1], str(item[0])))
    return [
        {
            "path": str(path.resolve()),
            "bytes": size,
        }
        for path, size in ranked[: max(1, limit)]
    ]


def load_gene_counts(
    *,
    gene_counts: Path,
    gene_counts_format: str,
) -> tuple[sp.csr_matrix, list[str], pd.DataFrame]:
    format_name = gene_counts_format.lower()
    if format_name == "tsv":
        df = pd.read_csv(gene_counts, sep="\t")
        if df.shape[1] < 2:
            raise ValueError("Gene-count TSV must have at least one identifier column and one cell column")
        feature_col = df.columns[0]
        gene_ids = df[feature_col].astype(str).tolist()
        cell_ids = [str(column) for column in df.columns[1:]]
        matrix = sp.csr_matrix(df.iloc[:, 1:].to_numpy(dtype=np.int64).T)
        var = pd.DataFrame({feature_col: gene_ids}, index=gene_ids)
        var.index.name = "gene_id"
        return matrix, cell_ids, var

    if format_name != "mtx-dir":
        raise ValueError(f"Unsupported gene-count format: {gene_counts_format}")
    matrix_path = gene_counts / "matrix.mtx"
    barcodes_path = gene_counts / "barcodes.tsv"
    features_path = gene_counts / "features.tsv"
    genes_path = gene_counts / "genes.tsv"
    if not matrix_path.exists():
        raise FileNotFoundError(f"Gene-count MatrixMarket directory is missing {matrix_path}")
    feature_source = features_path if features_path.exists() else genes_path
    if not feature_source.exists():
        raise FileNotFoundError(f"Gene-count MatrixMarket directory is missing {features_path.name} or {genes_path.name}")
    if not barcodes_path.exists():
        raise FileNotFoundError(f"Gene-count MatrixMarket directory is missing {barcodes_path}")

    X = scio.mmread(str(matrix_path)).tocsr()
    cell_ids = read_index_lines(barcodes_path)
    feature_rows = pd.read_csv(feature_source, sep="\t", header=None)
    gene_ids = feature_rows.iloc[:, 0].astype(str).tolist()
    if X.shape[0] == len(gene_ids) and X.shape[1] == len(cell_ids):
        matrix = X.T.tocsr()
    elif X.shape[0] == len(cell_ids) and X.shape[1] == len(gene_ids):
        matrix = X.tocsr()
    else:
        raise ValueError(
            f"Gene-count shape mismatch: matrix={X.shape} cells={len(cell_ids)} genes={len(gene_ids)}"
        )
    var = pd.DataFrame(index=gene_ids)
    var.index.name = "gene_id"
    if feature_rows.shape[1] > 1:
        var["gene_name"] = feature_rows.iloc[:, 1].astype(str).tolist()
    return matrix, cell_ids, var


def build_shared_obs(
    *,
    base_cell_qc: pd.DataFrame,
    shared_cell_ids: list[str],
    circ_present: set[str],
    rna_present: set[str],
) -> pd.DataFrame:
    obs = base_cell_qc.reindex(shared_cell_ids).copy()
    obs["circ_present"] = [cell_id in circ_present for cell_id in shared_cell_ids]
    obs["rna_present"] = [cell_id in rna_present for cell_id in shared_cell_ids]
    for column in ("assigned_reads", "circRNA_count", "total_circRNA_support"):
        if column in obs.columns:
            obs[column] = obs[column].fillna(0).astype(int)
    if "empty_flag" in obs.columns:
        obs["empty_flag"] = obs["empty_flag"].fillna(True).astype(bool)
    for column in ("alignment_status", "detector_status"):
        if column in obs.columns:
            obs[column] = obs[column].fillna("missing")
    return obs


def sanitize_frame_for_anndata(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for column in out.columns:
        series = out[column]
        if pd.api.types.is_bool_dtype(series):
            out[column] = series.fillna(False).astype(bool)
            continue
        parsed = pd.to_numeric(series, errors="coerce")
        if int(parsed.notna().sum()) == int(series.notna().sum()) and int(series.notna().sum()) > 0:
            out[column] = parsed
            continue
        if pd.api.types.is_numeric_dtype(series):
            out[column] = pd.to_numeric(series, errors="coerce")
            continue
        out[column] = series.fillna("").astype(str)
    return out


def export_circ_h5ad(
    *,
    out_path: Path,
    X_cells_by_circ: sp.csr_matrix,
    cell_qc: pd.DataFrame,
    circ_qc: pd.DataFrame,
    uns_payload: dict[str, Any],
) -> Optional[Path]:
    if not HAS_ANNDATA:
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.AnnData(
        X=X_cells_by_circ.tocsr().astype(np.int32),
        obs=sanitize_frame_for_anndata(cell_qc),
        var=sanitize_frame_for_anndata(circ_qc),
    )
    adata.uns["circyto"] = uns_payload
    adata.write_h5ad(str(out_path))
    return out_path


def export_mudata_bundle(
    *,
    out_path: Path,
    circ_X: sp.csr_matrix,
    circ_qc: pd.DataFrame,
    cell_qc: pd.DataFrame,
    uns_payload: dict[str, Any],
    gene_counts: Path,
    gene_counts_format: str,
    cell_join: str,
) -> dict[str, Any]:
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for MuData export")
    if not HAS_MUDATA:
        raise RuntimeError("mudata is not installed; install mudata or use --no-export-mudata")

    rna_X, rna_cell_ids, rna_var = load_gene_counts(gene_counts=gene_counts, gene_counts_format=gene_counts_format)
    circ_cell_ids = list(cell_qc.index)
    circ_set = set(circ_cell_ids)
    rna_set = set(rna_cell_ids)
    shared = sorted(circ_set & rna_set)
    if cell_join == "inner":
        joint_cell_ids = shared
    elif cell_join == "outer":
        joint_cell_ids = sorted(circ_set | rna_set)
    else:
        raise ValueError(f"Unsupported cell join: {cell_join}")

    circ_aligned, circ_obs_ids, circ_present = align_matrix_to_cells(
        circ_X,
        circ_cell_ids,
        joint_cell_ids,
        join=cell_join,
    )
    rna_aligned, rna_obs_ids, rna_present = align_matrix_to_cells(
        rna_X,
        rna_cell_ids,
        joint_cell_ids,
        join=cell_join,
    )
    if circ_obs_ids != rna_obs_ids:
        raise AssertionError("circ and rna cell alignment produced different obs orders")
    shared_obs = build_shared_obs(
        base_cell_qc=cell_qc,
        shared_cell_ids=circ_obs_ids,
        circ_present=circ_present,
        rna_present=rna_present,
    )

    shared_obs_clean = sanitize_frame_for_anndata(shared_obs)
    circ_adata = ad.AnnData(
        X=circ_aligned.tocsr().astype(np.int32),
        obs=shared_obs_clean.copy(),
        var=sanitize_frame_for_anndata(circ_qc),
    )
    rna_adata = ad.AnnData(
        X=rna_aligned.tocsr().astype(np.int32),
        obs=shared_obs_clean.copy(),
        var=sanitize_frame_for_anndata(rna_var),
    )
    mdata = mu.MuData({"rna": rna_adata, "circ": circ_adata})
    mdata.obs = shared_obs_clean.copy()
    mdata.uns["circyto"] = uns_payload
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mdata.write_h5mu(str(out_path))
    return {
        "path": str(out_path.resolve()),
        "cell_join": cell_join,
        "n_cells_circ": len(circ_cell_ids),
        "n_cells_rna": len(rna_cell_ids),
        "n_cells_shared": len(shared),
        "dropped_circ_only_cells": sorted(circ_set - rna_set) if cell_join == "inner" else [],
        "dropped_rna_only_cells": sorted(rna_set - circ_set) if cell_join == "inner" else [],
    }


def summarize_h5ad(input_path: Path, outdir: Path) -> dict[str, Any]:
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for summarize-h5ad")
    adata = ad.read_h5ad(str(input_path))
    X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    circ_per_cell = np.asarray(X.getnnz(axis=1)).ravel().astype(int)
    support_per_cell = np.asarray(X.sum(axis=1)).ravel().astype(int)
    cells_per_circ = np.asarray(X.getnnz(axis=0)).ravel().astype(int)
    support_per_circ = np.asarray(X.sum(axis=0)).ravel().astype(int)
    outdir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame({"cell_id": list(adata.obs_names), "circRNA_count": circ_per_cell, "total_support": support_per_cell}).to_csv(
        outdir / "circRNAs_per_cell.tsv", sep="\t", index=False
    )
    pd.DataFrame({"circ_id": list(adata.var_names), "n_cells_detected": cells_per_circ, "total_support": support_per_circ}).to_csv(
        outdir / "cells_per_circRNA.tsv", sep="\t", index=False
    )
    pd.DataFrame({"circ_id": list(adata.var_names), "n_cells_detected": cells_per_circ, "total_support": support_per_circ}).sort_values(
        by=["n_cells_detected", "total_support", "circ_id"], ascending=[False, False, True]
    ).head(20).to_csv(outdir / "top_recurrent_circRNAs.tsv", sep="\t", index=False)
    pd.DataFrame({"cell_id": list(adata.obs_names), "total_support": support_per_cell}).to_csv(
        outdir / "total_support_distribution.tsv", sep="\t", index=False
    )

    payload = {
        "input": str(input_path.resolve()),
        "n_cells": int(adata.n_obs),
        "n_circRNAs": int(adata.n_vars),
        "nnz": int(X.nnz),
        "circRNAs_per_cell_summary": numeric_summary(circ_per_cell.tolist()),
        "cells_per_circRNA_summary": numeric_summary(cells_per_circ.tolist()),
        "total_support_per_cell_summary": numeric_summary(support_per_cell.tolist()),
    }
    write_json(outdir / "summary.json", payload)
    return payload
