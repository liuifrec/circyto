from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional
import csv
import json

from circyto.demux.smartseq3 import SmartSeq3DemuxParams, demux_smartseq3
from circyto.detectors.ciri3 import Ciri3Detector
from circyto.manifest.alignment import read_alignment_manifest_tsv
from circyto.pipeline.align_manifest import (
    prepare_alignment_cache,
    read_source_manifest,
    run_detector_alignment_manifest,
)
from circyto.pipeline.collect import collect_matrix


ProgressFn = Callable[[str], None]


@dataclass(frozen=True)
class SmartSeq3Ciri3WorkflowParams:
    read1: Path
    read2: Path
    index1: Path
    index2: Path
    annotation: Path
    cell_id_column: str
    index1_column: str
    index2_column: str
    ref_fa: Path
    star_index: Path
    outdir: Path
    top_n: Optional[int] = None
    max_records: Optional[int] = None
    threads: int = 8
    parallel: int = 1
    chunk_size: int = 1
    max_mismatch: int = 0
    write_sink: bool = True
    resume: bool = False


def _workflow_paths(outdir: Path) -> dict[str, Path]:
    return {
        "root": outdir,
        "demux": outdir / "demux",
        "manifests": outdir / "manifests",
        "align": outdir / "align",
        "ciri3": outdir / "ciri3",
        "matrix": outdir / "matrix",
        "logs": outdir / "logs",
        "workflow_summary": outdir / "workflow_summary.json",
    }


def _emit(progress: ProgressFn, message: str) -> None:
    progress(f"[workflow smartseq3-ciri3] {message}")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _validate_demux_complete(demux_dir: Path) -> bool:
    summary_path = demux_dir / "demux_summary.json"
    manifest_path = demux_dir / "manifest.tsv"
    if not summary_path.exists() or not manifest_path.exists():
        return False
    _load_json(summary_path)
    read_source_manifest(manifest_path, validate_files=True)
    return True


def _validate_alignment_complete(align_dir: Path) -> bool:
    manifest_path = align_dir / "alignment_manifest.tsv"
    summary_path = align_dir / "alignment_prepare_summary.json"
    if not manifest_path.exists() or not summary_path.exists():
        return False
    _load_json(summary_path)
    read_alignment_manifest_tsv(manifest_path, validate_files=True)
    return True


def _validate_detector_complete(detector_dir: Path) -> bool:
    summary_path = detector_dir / "detector_run_summary.json"
    if not summary_path.exists():
        return False
    _load_json(summary_path)
    return True


def _validate_matrix_complete(matrix_dir: Path) -> bool:
    required = [
        matrix_dir / "circ_counts.mtx",
        matrix_dir / "circ_index.txt",
        matrix_dir / "cell_index.txt",
    ]
    return all(path.exists() for path in required)


def _read_manifest_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = [{key: "" if value is None else str(value) for key, value in row.items()} for row in reader]
    return header, rows


def _write_manifest_rows(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def _selected_manifest_path(manifests_dir: Path, top_n: Optional[int]) -> Path:
    if top_n is None:
        return manifests_dir / "all_manifest.tsv"
    return manifests_dir / f"top{top_n}_manifest.tsv"


def _compute_selected_cell_ids(
    *,
    demux_summary: dict[str, Any],
    manifest_rows: list[dict[str, str]],
    top_n: Optional[int],
) -> list[str]:
    manifest_cell_ids = [(row.get("cell_id") or "").strip() for row in manifest_rows if (row.get("cell_id") or "").strip()]
    reads_per_cell_raw = demux_summary.get("reads_per_cell", {})
    reads_per_cell = {str(cell_id): int(count) for cell_id, count in dict(reads_per_cell_raw).items()}
    if top_n is None:
        return sorted(manifest_cell_ids)
    ranked = sorted(reads_per_cell.items(), key=lambda item: (-item[1], item[0]))
    selected = {cell_id for cell_id, _ in ranked[:top_n]}
    return sorted(cell_id for cell_id in manifest_cell_ids if cell_id in selected)


def create_selected_manifest(
    *,
    demux_manifest: Path,
    demux_summary_path: Path,
    output_path: Path,
    top_n: Optional[int],
) -> dict[str, Any]:
    if top_n is not None and top_n <= 0:
        raise ValueError("top_n must be > 0 when provided")

    header, rows = _read_manifest_rows(demux_manifest)
    if not rows:
        raise ValueError(f"Demux manifest contains 0 rows: {demux_manifest}")
    demux_summary = _load_json(demux_summary_path)
    selected_cell_ids = _compute_selected_cell_ids(
        demux_summary=demux_summary,
        manifest_rows=rows,
        top_n=top_n,
    )
    selected_rows = [row for row in rows if (row.get("cell_id") or "").strip() in set(selected_cell_ids)]
    if not selected_rows:
        raise ValueError("Selected manifest would be empty")
    _write_manifest_rows(output_path, header, selected_rows)
    read_source_manifest(output_path, validate_files=True)
    return {
        "manifest_path": str(output_path.resolve()),
        "selected_cell_ids": selected_cell_ids,
        "selected_cell_count": len(selected_rows),
    }


def _selected_manifest_matches(
    *,
    manifest_path: Path,
    demux_manifest: Path,
    demux_summary_path: Path,
    top_n: Optional[int],
) -> bool:
    if not manifest_path.exists():
        return False
    read_source_manifest(manifest_path, validate_files=True)
    _, all_rows = _read_manifest_rows(demux_manifest)
    expected_ids = _compute_selected_cell_ids(
        demux_summary=_load_json(demux_summary_path),
        manifest_rows=all_rows,
        top_n=top_n,
    )
    _, selected_rows = _read_manifest_rows(manifest_path)
    actual_ids = sorted((row.get("cell_id") or "").strip() for row in selected_rows if (row.get("cell_id") or "").strip())
    return actual_ids == expected_ids


def _status_counts_from_summary(path: Path) -> dict[str, int]:
    summary = _load_json(path)
    raw = summary.get("status_counts", {})
    return {str(key): int(value) for key, value in dict(raw).items()}


def _matrix_stats(matrix_path: Path) -> dict[str, int]:
    with matrix_path.open("r", encoding="utf-8") as handle:
        size_line_found = False
        for line in handle:
            text = line.strip()
            if not text or text.startswith("%"):
                continue
            parts = text.split()
            if len(parts) != 3:
                raise ValueError(f"Malformed MatrixMarket size line in {matrix_path}: {text}")
            n_rows, n_cols, nnz = (int(part) for part in parts)
            size_line_found = True
            break
        if not size_line_found:
            raise ValueError(f"Could not read MatrixMarket header from {matrix_path}")
    return {"n_rows": n_rows, "n_cols": n_cols, "nnz": nnz}


def _build_workflow_summary(
    *,
    params: SmartSeq3Ciri3WorkflowParams,
    paths: dict[str, Path],
    selected_manifest_path: Path,
    completed_stages: list[str],
) -> dict[str, Any]:
    demux_summary = _load_json(paths["demux"] / "demux_summary.json")
    alignment_summary_path = paths["align"] / "alignment_prepare_summary.json"
    detector_summary_path = paths["ciri3"] / "detector_run_summary.json"
    matrix_path = paths["matrix"] / "circ_counts.mtx"
    circ_index_path = paths["matrix"] / "circ_index.txt"
    cell_index_path = paths["matrix"] / "cell_index.txt"
    selected_rows = read_source_manifest(selected_manifest_path, validate_files=True)
    matrix_stats = _matrix_stats(matrix_path)
    return {
        "workflow": "smartseq3-ciri3",
        "experimental": True,
        "command_options": {
            "read1": str(params.read1),
            "read2": str(params.read2),
            "index1": str(params.index1),
            "index2": str(params.index2),
            "annotation": str(params.annotation),
            "cell_id_column": params.cell_id_column,
            "index1_column": params.index1_column,
            "index2_column": params.index2_column,
            "ref_fa": str(params.ref_fa),
            "star_index": str(params.star_index),
            "outdir": str(params.outdir),
            "top_n": params.top_n,
            "max_records": params.max_records,
            "threads": params.threads,
            "parallel": params.parallel,
            "chunk_size": params.chunk_size,
            "max_mismatch": params.max_mismatch,
            "write_sink": params.write_sink,
            "resume": params.resume,
        },
        "demux": {
            "assigned_records": int(demux_summary.get("assigned_records", 0)),
            "total_records": int(demux_summary.get("total_records", 0)),
            "number_of_cells_detected": int(demux_summary.get("number_of_cells_detected", 0)),
        },
        "selected_cell_count": len(selected_rows),
        "alignment_status_counts": _status_counts_from_summary(alignment_summary_path),
        "detector_status_counts": _status_counts_from_summary(detector_summary_path),
        "matrix": matrix_stats,
        "completed_stages": completed_stages,
        "paths": {
            "demux_summary": str((paths["demux"] / "demux_summary.json").resolve()),
            "demux_manifest": str((paths["demux"] / "manifest.tsv").resolve()),
            "selected_manifest": str(selected_manifest_path.resolve()),
            "alignment_manifest": str((paths["align"] / "alignment_manifest.tsv").resolve()),
            "alignment_summary": str(alignment_summary_path.resolve()),
            "detector_summary": str(detector_summary_path.resolve()),
            "matrix": str(matrix_path.resolve()),
            "circ_index": str(circ_index_path.resolve()),
            "cell_index": str(cell_index_path.resolve()),
        },
    }


def run_smartseq3_ciri3_workflow(
    params: SmartSeq3Ciri3WorkflowParams,
    *,
    progress: ProgressFn = print,
) -> dict[str, Any]:
    if params.top_n is not None and params.top_n <= 0:
        raise ValueError("top_n must be > 0 when provided")
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

    completed_stages: list[str] = []

    demux_summary_path = paths["demux"] / "demux_summary.json"
    demux_manifest = paths["demux"] / "manifest.tsv"
    if params.resume and _validate_demux_complete(paths["demux"]):
        _emit(progress, f"resume skip: demux ({demux_summary_path})")
        completed_stages.append("demux_skipped")
    else:
        _emit(progress, "stage 1/5: SMART-Seq3 demux")
        demux_summary = demux_smartseq3(
            SmartSeq3DemuxParams(
                read1=params.read1,
                read2=params.read2,
                index1=params.index1,
                index2=params.index2,
                annotation=params.annotation,
                outdir=paths["demux"],
                cell_id_column=params.cell_id_column,
                index1_column=params.index1_column,
                index2_column=params.index2_column,
                max_mismatch=params.max_mismatch,
                max_records=params.max_records,
                write_sink=params.write_sink,
                compress_output=True,
                emit_manifest=True,
            )
        )
        _emit(
            progress,
            "demux complete: "
            f"assigned={demux_summary['assigned_records']} total={demux_summary['total_records']} "
            f"cells={demux_summary['number_of_cells_detected']}",
        )
        completed_stages.append("demux")

    selected_manifest = _selected_manifest_path(paths["manifests"], params.top_n)
    if params.resume and _selected_manifest_matches(
        manifest_path=selected_manifest,
        demux_manifest=demux_manifest,
        demux_summary_path=demux_summary_path,
        top_n=params.top_n,
    ):
        _emit(progress, f"resume skip: manifest selection ({selected_manifest})")
        completed_stages.append("manifest_selection_skipped")
    else:
        label = f"top {params.top_n}" if params.top_n is not None else "all cells"
        _emit(progress, f"stage 2/5: manifest selection ({label})")
        selection = create_selected_manifest(
            demux_manifest=demux_manifest,
            demux_summary_path=demux_summary_path,
            output_path=selected_manifest,
            top_n=params.top_n,
        )
        _emit(progress, f"selected {selection['selected_cell_count']} cells into {selected_manifest}")
        completed_stages.append("manifest_selection")

    alignment_manifest_path = paths["align"] / "alignment_manifest.tsv"
    if params.resume and _validate_alignment_complete(paths["align"]):
        _emit(progress, f"resume skip: alignment prep ({alignment_manifest_path})")
        completed_stages.append("alignment_skipped")
    else:
        _emit(progress, "stage 3/5: prepare STAR alignment cache for CIRI3")
        prepare_alignment_cache(
            manifest=selected_manifest,
            outdir=paths["align"],
            aligner="star",
            ref_fa=params.ref_fa,
            detector_hint="ciri3",
            threads=params.threads,
            parallel=params.parallel,
            chunk_size=params.chunk_size,
            extra_flags=f"--genomeDir {params.star_index} --readFilesCommand zcat",
            output_format="bam",
        )
        _emit(progress, f"alignment cache ready: {alignment_manifest_path}")
        completed_stages.append("alignment")

    detector_summary_path = paths["ciri3"] / "detector_run_summary.json"
    if params.resume and _validate_detector_complete(paths["ciri3"]):
        _emit(progress, f"resume skip: detector ({detector_summary_path})")
        completed_stages.append("detector_skipped")
    else:
        _emit(progress, "stage 4/5: run CIRI3 from alignments")
        run_detector_alignment_manifest(
            detector=Ciri3Detector(),
            manifest=alignment_manifest_path,
            outdir=paths["ciri3"],
            ref_fa=params.ref_fa,
            threads=params.threads,
            parallel=params.parallel,
            chunk_size=params.chunk_size,
        )
        _emit(progress, f"detector complete: {detector_summary_path}")
        completed_stages.append("detector")

    if params.resume and _validate_matrix_complete(paths["matrix"]):
        _emit(progress, f"resume skip: matrix ({paths['matrix'] / 'circ_counts.mtx'})")
        completed_stages.append("matrix_skipped")
    else:
        _emit(progress, "stage 5/5: collect matrix")
        collect_matrix(
            cirifull_dir=str(paths["ciri3"]),
            matrix_path=str(paths["matrix"] / "circ_counts.mtx"),
            circ_index_path=str(paths["matrix"] / "circ_index.txt"),
            cell_index_path=str(paths["matrix"] / "cell_index.txt"),
            min_count_per_cell=1,
        )
        _emit(progress, f"matrix complete: {paths['matrix'] / 'circ_counts.mtx'}")
        completed_stages.append("matrix")

    summary = _build_workflow_summary(
        params=params,
        paths=paths,
        selected_manifest_path=selected_manifest,
        completed_stages=completed_stages,
    )
    _write_json(paths["workflow_summary"], summary)
    _emit(progress, f"workflow summary: {paths['workflow_summary']}")
    return summary
