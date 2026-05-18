from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional
import csv
import time

from circyto.pipeline.align_manifest import read_source_manifest
from circyto.pipeline.run_ciri3 import (
    RunCiri3Params,
    _resolve_row_protocol,
    _resolve_row_strandedness,
    run_ciri3_workflow,
)
from circyto.pipeline.workflow_reporting import (
    build_cell_qc_table,
    build_circ_qc_table,
    expand_cells,
    export_circ_h5ad,
    load_circ_feature_table,
    load_circ_matrix,
    load_json,
    matrix_section,
    numeric_summary,
    top_mapping_items,
    write_json,
)


ProgressFn = Callable[[str], None]


@dataclass(frozen=True)
class FullLengthCircRNAWorkflowParams:
    manifest: Path
    outdir: Path
    protocol: str
    genome_fasta: Path
    gtf: Path
    detector: str = "ciri3"
    star_index: Optional[Path] = None
    threads: int = 8
    parallel: int = 1
    chunk_size: int = 1
    skip_demux: bool = False
    export_h5ad: bool = True
    dry_run: bool = False
    fail_fast: bool = False
    command_template: Optional[str] = None
    experimental_paired_ramda: bool = False


def _workflow_paths(outdir: Path) -> dict[str, Path]:
    return {
        "root": outdir,
        "align": outdir / "align",
        "ciri3": outdir / "ciri3",
        "matrix": outdir / "matrix",
        "qc": outdir / "qc",
        "anndata": outdir / "anndata",
        "logs": outdir / "logs",
        "manifests": outdir / "manifests",
        "workflow_summary": outdir / "workflow_summary.json",
    }


def _emit(progress: ProgressFn, message: str) -> None:
    progress(f"[workflow full-length-circrna] {message}")


def _read_manifest_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = [{key: "" if value is None else str(value) for key, value in row.items()} for row in reader]
    return header, rows


def _n_input_reads_by_cell(raw_rows: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in raw_rows:
        cell_id = (row.get("cell_id") or row.get("sample_id") or "").strip()
        if not cell_id:
            continue
        raw = (row.get("n_input_reads") or "").strip()
        try:
            counts[cell_id] = int(raw) if raw else 0
        except ValueError:
            counts[cell_id] = 0
    return counts


def _per_cell_record_map(summary: dict[str, Any]) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for record in summary.get("cells", []):
        cell_id = str(record.get("cell_id", "")).strip()
        if cell_id:
            records[cell_id] = dict(record)
    return records


def _detector_status_counts(path: Path) -> dict[str, int]:
    summary = load_json(path)
    raw = summary.get("status_counts", {})
    return {str(key): int(value) for key, value in dict(raw).items()}


def _matrix_stats_header(matrix_path: Path) -> dict[str, int]:
    with matrix_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("%"):
                continue
            n_rows, n_cols, nnz = (int(part) for part in text.split())
            return {"n_rows": n_rows, "n_cols": n_cols, "nnz": nnz}
    raise ValueError(f"Could not read MatrixMarket header from {matrix_path}")


def _stage_graph(*, skip_demux_effective: bool, export_h5ad: bool, dry_run: bool) -> list[dict[str, Any]]:
    graph = [
        {
            "stage": "demux",
            "status": "skipped" if skip_demux_effective else "delegated",
            "detail": (
                "Manifest-based full-length workflow skips demux."
                if skip_demux_effective
                else "Pooled SMART-Seq3 demux remains the existing workflow smartseq3-ciri3 path."
            ),
        },
        {"stage": "alignment", "status": "planned" if dry_run else "completed"},
        {"stage": "detector", "status": "planned" if dry_run else "completed"},
        {"stage": "matrix", "status": "planned" if dry_run else "completed"},
        {"stage": "h5ad_export", "status": "planned" if dry_run and export_h5ad else ("completed" if export_h5ad else "disabled")},
        {"stage": "summary_qc", "status": "planned" if dry_run else "completed"},
    ]
    return graph


def _build_dry_run_summary(
    *,
    params: FullLengthCircRNAWorkflowParams,
    paths: dict[str, Path],
    source_rows,
    resolved_protocols: dict[str, str],
    resolved_strandedness: dict[str, str],
    skip_demux_effective: bool,
    warnings: list[str],
) -> dict[str, Any]:
    read_layout_counts = {"single-end": 0, "paired-end": 0}
    for row in source_rows:
        read_layout_counts[row.read_layout] = read_layout_counts.get(row.read_layout, 0) + 1
    return {
        "workflow": "full-length-circrna",
        "command_name": "circyto workflow full-length-circrna",
        "dry_run": True,
        "detector": params.detector,
        "protocol": params.protocol,
        "skip_demux_effective": skip_demux_effective,
        "stage_graph": _stage_graph(skip_demux_effective=skip_demux_effective, export_h5ad=params.export_h5ad, dry_run=True),
        "planned_cells": len(source_rows),
        "read_layout_counts": read_layout_counts,
        "protocols": sorted(set(resolved_protocols.values())),
        "strandedness": sorted({value for value in resolved_strandedness.values() if value}),
        "allow_paired_ramda": params.experimental_paired_ramda,
        "experimental_paired_ramda": params.experimental_paired_ramda,
        "warnings": warnings,
        "paths": {
            "manifest": str(params.manifest.resolve()),
            "align_dir": str(paths["align"].resolve()),
            "ciri3_dir": str(paths["ciri3"].resolve()),
            "matrix_dir": str(paths["matrix"].resolve()),
            "anndata_dir": str(paths["anndata"].resolve()) if params.export_h5ad else None,
            "run_ciri3_summary": str((paths["logs"] / "run_ciri3_summary.json").resolve()),
        },
    }


def _build_execution_summary(
    *,
    params: FullLengthCircRNAWorkflowParams,
    paths: dict[str, Path],
    source_rows,
    raw_rows: list[dict[str, str]],
    resolved_protocols: dict[str, str],
    resolved_strandedness: dict[str, str],
    skip_demux_effective: bool,
    elapsed_seconds: float,
    warnings: list[str],
) -> dict[str, Any]:
    selected_cell_ids = [row.cell_id for row in source_rows]
    assigned_reads = _n_input_reads_by_cell(raw_rows)
    alignment_summary_path = paths["align"] / "alignment_prepare_summary.json"
    detector_summary_path = paths["ciri3"] / "detector_run_summary.json"
    matrix_path = paths["matrix"] / "circ_counts.mtx"
    circ_index_path = paths["matrix"] / "circ_index.txt"
    cell_index_path = paths["matrix"] / "cell_index.txt"
    circ_feature_table_path = paths["matrix"] / "circ_feature_table.tsv"

    alignment_summary = load_json(alignment_summary_path)
    detector_summary = load_json(detector_summary_path)
    X_circ_by_cell, circ_ids, matrix_cell_ids = load_circ_matrix(
        matrix_path=matrix_path,
        circ_index_path=circ_index_path,
        cell_index_path=cell_index_path,
    )
    X_cells_by_circ = X_circ_by_cell.T.tocsr()
    X_selected_cells_by_circ = expand_cells(X_cells_by_circ, matrix_cell_ids, selected_cell_ids)
    feature_df = load_circ_feature_table(circ_ids, circ_feature_table_path)
    alignment_cells = _per_cell_record_map(alignment_summary)
    detector_cells = _per_cell_record_map(detector_summary)

    cell_qc = build_cell_qc_table(
        selected_cell_ids=selected_cell_ids,
        assigned_reads=assigned_reads,
        alignment_cells=alignment_cells,
        detector_cells=detector_cells,
        X_cells_by_circ=X_selected_cells_by_circ,
    )
    cell_qc["sample_id"] = selected_cell_ids
    cell_qc["protocol"] = [resolved_protocols[cell_id] for cell_id in selected_cell_ids]
    cell_qc["read_layout"] = [row.read_layout for row in source_rows]
    cell_qc["strandedness"] = [resolved_strandedness[cell_id] for cell_id in selected_cell_ids]

    circ_qc = build_circ_qc_table(
        circ_ids=circ_ids,
        feature_df=feature_df,
        X_cells_by_circ=X_selected_cells_by_circ,
    )

    cell_qc_path = paths["qc"] / "cell_qc.tsv"
    circ_qc_path = paths["qc"] / "circ_qc.tsv"
    cell_qc.reset_index().to_csv(cell_qc_path, sep="\t", index=False)
    circ_qc.reset_index().to_csv(circ_qc_path, sep="\t", index=False)

    h5ad_path = paths["anndata"] / "circ_counts.h5ad"
    h5ad_status = None
    if params.export_h5ad:
        exported_h5ad = export_circ_h5ad(
            out_path=h5ad_path,
            X_cells_by_circ=X_selected_cells_by_circ,
            cell_qc=cell_qc,
            circ_qc=circ_qc,
            uns_payload={
                "workflow_name": "full-length-circrna",
                "detector": params.detector,
                "reference": str(params.genome_fasta),
                "gtf": str(params.gtf),
                "protocol": params.protocol,
                "manifest": str(params.manifest.resolve()),
                "skip_demux_effective": skip_demux_effective,
                "dry_run": False,
                "allow_paired_ramda": params.experimental_paired_ramda,
                "experimental_paired_ramda": params.experimental_paired_ramda,
            },
        )
        h5ad_status = str(exported_h5ad.resolve()) if exported_h5ad is not None else None

    circ_counts_by_cell = {cell_id: int(value) for cell_id, value in zip(selected_cell_ids, cell_qc["circRNA_count"].tolist())}
    detector_status_counts = _detector_status_counts(detector_summary_path)

    return {
        "workflow": "full-length-circrna",
        "command_name": "circyto workflow full-length-circrna",
        "dry_run": False,
        "detector_backend": params.detector,
        "protocol": params.protocol,
        "skip_demux_effective": skip_demux_effective,
        "stage_graph": _stage_graph(skip_demux_effective=skip_demux_effective, export_h5ad=params.export_h5ad, dry_run=False),
        "planned_cells": len(source_rows),
        "completed_cells": int(detector_summary.get("completed_cells", 0)),
        "warnings": warnings,
        "protocols": sorted(set(resolved_protocols.values())),
        "strandedness": sorted({value for value in resolved_strandedness.values() if value}),
        "allow_paired_ramda": params.experimental_paired_ramda,
        "experimental_paired_ramda": params.experimental_paired_ramda,
        "alignment_status_counts": {
            str(key): int(value)
            for key, value in dict(alignment_summary.get("status_counts", {})).items()
        },
        "detector_status_counts": detector_status_counts,
        "detector": {
            "name": params.detector,
            "success": int(detector_status_counts.get("success", 0)),
            "empty": int(detector_status_counts.get("empty", 0)),
            "failed": int(detector_status_counts.get("failed", 0)),
            "circRNAs_per_cell_summary": numeric_summary(cell_qc["circRNA_count"].astype(int).tolist()),
            "top_cells_by_circRNA_count": top_mapping_items(circ_counts_by_cell, count_key="circRNA_count"),
            "elapsed_seconds": detector_summary.get("elapsed_seconds"),
        },
        "matrix": {
            **_matrix_stats_header(matrix_path),
            **matrix_section(X_selected_cells_by_circ, circ_qc),
        },
        "paths": {
            "manifest": str(params.manifest.resolve()),
            "alignment_manifest": str((paths["align"] / "alignment_manifest.tsv").resolve()),
            "alignment_summary": str(alignment_summary_path.resolve()),
            "detector_summary": str(detector_summary_path.resolve()),
            "matrix": str(matrix_path.resolve()),
            "circ_index": str(circ_index_path.resolve()),
            "cell_index": str(cell_index_path.resolve()),
            "circ_feature_table": str(circ_feature_table_path.resolve()),
            "cell_qc": str(cell_qc_path.resolve()),
            "circ_qc": str(circ_qc_path.resolve()),
            "h5ad": h5ad_status,
        },
        "workflow_timing": {
            "total_elapsed_seconds": round(elapsed_seconds, 3),
        },
    }


def run_full_length_circrna_workflow(
    params: FullLengthCircRNAWorkflowParams,
    *,
    progress: ProgressFn = print,
) -> dict[str, Any]:
    if params.detector != "ciri3":
        raise ValueError("full-length-circrna currently supports only --detector ciri3")
    if params.threads <= 0:
        raise ValueError("threads must be > 0")
    if params.parallel <= 0:
        raise ValueError("parallel must be > 0")
    if params.chunk_size <= 0:
        raise ValueError("chunk_size must be > 0")

    paths = _workflow_paths(params.outdir)
    for path in paths.values():
        if path.suffix:
            path.parent.mkdir(parents=True, exist_ok=True)
        else:
            path.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    source_rows = read_source_manifest(params.manifest, validate_files=not params.dry_run)
    _, raw_rows = _read_manifest_rows(params.manifest)
    if len(source_rows) != len(raw_rows):
        raise ValueError("Manifest parsing mismatch between structured and raw readers")

    resolved_protocols: dict[str, str] = {}
    resolved_strandedness: dict[str, str] = {}
    for row in source_rows:
        resolved_protocols[row.cell_id] = _resolve_row_protocol(row, params.protocol)
        resolved_strandedness[row.cell_id] = _resolve_row_strandedness(row, None)

    unique_protocols = sorted(set(resolved_protocols.values()))
    if len(unique_protocols) != 1:
        raise ValueError(f"full-length-circrna requires a single protocol per run; found: {', '.join(unique_protocols)}")

    effective_protocol = unique_protocols[0]
    skip_demux_effective = params.skip_demux or effective_protocol in {"ramda", "shin-ramda"}
    warnings: list[str] = []

    if effective_protocol == "smartseq3" and not skip_demux_effective:
        raise ValueError(
            "SMART-Seq3 pooled demux remains the existing `circyto workflow smartseq3-ciri3` path. "
            "Use --skip-demux only for already demultiplexed SMART-Seq3 per-cell FASTQ manifests."
        )

    paired_end_rows = [row for row in source_rows if row.read_layout == "paired-end"]
    if paired_end_rows and effective_protocol in {"ramda", "shin-ramda"}:
        message = (
            f"Protocol {effective_protocol} has paired-end rows. Paired-end rows use the validated STAR+CIRI3 paired-end "
            "route, which was locally validated on a real GSE278952 chr21 subset. Use --allow-paired-ramda for explicit "
            "opt-in before full hg38-scale biological validation."
        )
        if params.dry_run:
            warnings.append(message)
            _emit(progress, f"warning: {message}")
        elif not params.experimental_paired_ramda:
            raise ValueError(
                message
                + " Re-run with --allow-paired-ramda for real execution, or use --dry-run to inspect the planned commands. "
                + "--experimental-paired-ramda remains accepted as a deprecated alias."
            )
        else:
            warnings.append(message)
            _emit(progress, f"warning: {message}")

    _emit(progress, f"protocol={effective_protocol} cells={len(source_rows)} skip_demux={skip_demux_effective} dry_run={params.dry_run}")
    run_ciri3_workflow(
        RunCiri3Params(
            manifest=params.manifest,
            outdir=params.outdir,
            genome_fasta=params.genome_fasta,
            gtf=params.gtf,
            star_index=params.star_index,
            protocol=effective_protocol,
            strandedness=None,
            threads=params.threads,
            parallel=params.parallel,
            chunk_size=params.chunk_size,
            dry_run=params.dry_run,
            fail_fast=params.fail_fast,
            command_template=params.command_template,
        ),
        progress=progress,
    )

    if params.dry_run:
        summary = _build_dry_run_summary(
            params=params,
            paths=paths,
            source_rows=source_rows,
            resolved_protocols=resolved_protocols,
            resolved_strandedness=resolved_strandedness,
            skip_demux_effective=skip_demux_effective,
            warnings=warnings,
        )
    else:
        summary = _build_execution_summary(
            params=params,
            paths=paths,
            source_rows=source_rows,
            raw_rows=raw_rows,
            resolved_protocols=resolved_protocols,
            resolved_strandedness=resolved_strandedness,
            skip_demux_effective=skip_demux_effective,
            elapsed_seconds=time.perf_counter() - started,
            warnings=warnings,
        )

    write_json(paths["workflow_summary"], summary)
    _emit(progress, f"workflow summary: {paths['workflow_summary']}")
    return summary
