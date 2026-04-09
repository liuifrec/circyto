# circyto/pipeline/run_detector.py
from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from circyto.detectors import DetectorBase, DetectorRunInputs, DetectorResult
from circyto.detectors.base import get_detector_capabilities
from circyto.paths import resolve_manifest_path

def ensure_dir(path: Path) -> None:
    """
    Create the directory if it doesn't exist.
    Small local helper to avoid importing non-existent utils modules.
    """
    path.mkdir(parents=True, exist_ok=True)

def _pick_col(row: dict, keys: Tuple[str, ...]) -> Optional[str]:
    for k in keys:
        if k in row and row.get(k) not in (None, ""):
            return k
    return None


def _resolve_manifest_path(manifest_path: Path, value: str) -> Path:
    return resolve_manifest_path(manifest_path, value.strip())


def _detector_output_has_calls(path: Path, detector_name: str) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False

    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if detector_name in {"ciri-full", "ciri2"}:
        data_lines = [line for line in text[1:] if line.strip()]
        return len(data_lines) > 0

    data_lines = [line for line in text if line.strip() and not line.startswith("#")]
    return len(data_lines) > 0


def _count_output_rows(path: Path, detector_name: str) -> int | None:
    if not path.exists() or path.stat().st_size == 0:
        return 0

    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if detector_name in {"ciri-full", "ciri2", "ciri3"}:
        return len([line for line in text[1:] if line.strip()])
    return len([line for line in text if line.strip() and not line.startswith("#")])


def _outcome_category(status: str, *, raw_rows: int | None, normalized_rows: int | None) -> str:
    if status == "failed":
        return "failed"
    if status == "skipped_existing":
        return "skipped-existing"
    if normalized_rows and normalized_rows > 0:
        return "success-non-empty"
    if raw_rows is not None and raw_rows > 0 and (normalized_rows or 0) == 0:
        return "success-normalized-empty"
    return "success-empty"


def _flatten_detector_result(result: Any) -> List[DetectorResult]:
    if isinstance(result, DetectorResult):
        return [result]
    if isinstance(result, (list, tuple)):
        flat: List[DetectorResult] = []
        for item in result:
            flat.extend(_flatten_detector_result(item))
        return flat
    raise TypeError(f"Detector returned unsupported result type: {type(result)!r}")


def _summary_path(outdir: Path) -> Path:
    return outdir / "detector_run_summary.json"


def _write_summary(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _provenance_path(output_path: Path) -> Path:
    return Path(str(output_path) + ".provenance.json")


def _normalize_stamp_path(path: Path | None) -> str | None:
    if path is None:
        return None
    return str(path.resolve())


def _build_provenance_stamp(
    detector: DetectorBase,
    *,
    cell_id: str,
    r1: Path,
    r2: Path | None,
    ref_fa: Path | None,
    gtf: Path | None,
    threads: int,
    extra: dict[str, Any],
) -> dict[str, Any]:
    read_layout = "paired-end" if r2 is not None else "single-end"
    return {
        "detector": detector.name,
        "detector_class": detector.__class__.__name__,
        "detector_backend": detector.name,
        "cell_id": cell_id,
        "read1": _normalize_stamp_path(r1),
        "read2": _normalize_stamp_path(r2),
        "read_layout": read_layout,
        "execution_mode": "per-cell-fastq",
        "input_mode": "fastq",
        "reused_alignment": False,
        "ref_fa": _normalize_stamp_path(ref_fa),
        "gtf": _normalize_stamp_path(gtf),
        "threads": threads,
        "extra": extra,
    }


def _load_provenance_stamp(output_path: Path) -> dict[str, Any] | None:
    provenance_path = _provenance_path(output_path)
    if not provenance_path.exists():
        return None
    try:
        return json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _write_provenance_stamp(output_path: Path, stamp: dict[str, Any]) -> None:
    _provenance_path(output_path).write_text(
        json.dumps(stamp, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _provenance_matches(existing: dict[str, Any] | None, expected: dict[str, Any]) -> bool:
    if existing is None:
        return False
    if existing == expected:
        return True
    return all(expected.get(key) == value for key, value in existing.items())


def read_manifest(path: Path, *, validate_files: bool = True) -> List[Tuple[str, Path, Optional[Path]]]:
    """
    Read manifest rows as: (cell_id, r1_path, r2_path_or_None).

    Supports:
      - legacy:  cell_id, r1, r2
      - v1:      cell_id, read1, read2, (plus extra columns)

    Notes:
      - We intentionally do NOT require r2 (some platforms may be single-end).
      - Relative paths are resolved relative to the manifest location.
    """
    import csv

    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")

    rows: List[Tuple[str, Path, Optional[Path]]] = []
    seen_cells: set[str] = set()
    with path.open("r", newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        if rd.fieldnames is None:
            raise ValueError(f"Manifest has no header row: {path}")

        # Normalize header expectation
        if "cell_id" not in rd.fieldnames:
            raise KeyError(f"Manifest missing required column 'cell_id': {path}")

        for i, r in enumerate(rd, start=2):  # header is line 1
            cell = (r.get("cell_id") or "").strip()
            if not cell:
                raise ValueError(f"Empty cell_id at {path}:{i}")
            if cell in seen_cells:
                raise ValueError(f"Duplicate cell_id '{cell}' at {path}:{i}")
            seen_cells.add(cell)

            r1_key = _pick_col(r, ("r1", "read1"))
            r2_key = _pick_col(r, ("r2", "read2"))

            if r1_key is None:
                raise KeyError(f"Manifest missing r1/read1 for cell_id={cell} at {path}:{i}")

            r1_value = str(r[r1_key]).strip()
            if not r1_value:
                raise ValueError(f"Empty r1/read1 for cell_id={cell} at {path}:{i}")
            r1 = _resolve_manifest_path(path, r1_value)
            r2 = _resolve_manifest_path(path, str(r[r2_key]).strip()) if r2_key else None

            if validate_files and not r1.exists():
                raise FileNotFoundError(f"Manifest read1 not found for cell_id={cell}: {r1}")
            if validate_files and r2 is not None and not r2.exists():
                raise FileNotFoundError(f"Manifest read2 not found for cell_id={cell}: {r2}")

            rows.append((cell, r1, r2))

    if not rows:
        raise ValueError(f"Manifest contains 0 data rows: {path}")
    return rows
 


def run_detector_manifest(
    detector: DetectorBase,
    manifest: Path,
    outdir: Path,
    ref_fa: Path | None = None,
    gtf: Path | None = None,
    threads: int = 8,
    parallel: int = 4,
) -> list[DetectorResult]:
    """
    Run a single detector across all rows in a manifest.

    For detectors that are not process-parallel safe (e.g. CIRI-full),
    we respect detector.max_parallel if present and override the user-
    supplied `parallel` accordingly.
    """
    rows = read_manifest(manifest, validate_files=True)
    outdir.mkdir(parents=True, exist_ok=True)

    # NEW: limit effective parallelism based on detector capability
    det_max_parallel = get_detector_capabilities(detector).max_parallel
    effective_parallel = max(1, min(parallel, det_max_parallel, len(rows)))
    if effective_parallel < parallel:
        print(
            f"[circyto] Detector '{detector.name}' only supports parallel={effective_parallel}; "
            f"overriding requested parallel={parallel}."
        )

    started_at = time.perf_counter()
    summary_path = _summary_path(outdir)
    per_cell_records: list[dict[str, Any]] = []

    def _existing_output_for(cell_id: str) -> Path:
        if detector.name == "find-circ3":
            return outdir / cell_id / f"{cell_id}_splice_sites.bed"
        return outdir / f"{cell_id}.tsv"

    def _run_one(row: tuple[str, Path, Path | None]) -> tuple[list[DetectorResult], dict[str, Any]]:
        cell_id, r1, r2 = row
        read_layout = "paired-end" if r2 is not None else "single-end"
        existing_path = _existing_output_for(cell_id)
        expected_stamp = _build_provenance_stamp(
            detector,
            cell_id=cell_id,
            r1=r1,
            r2=r2,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
            extra={},
        )
        existing_stamp = _load_provenance_stamp(existing_path)
        if _detector_output_has_calls(existing_path, detector.name) and _provenance_matches(
            existing_stamp, expected_stamp
        ):
            result = DetectorResult(
                detector=detector.name,
                cell_id=cell_id,
                outdir=outdir,
                tsv_path=existing_path,
                meta={"skipped_existing": True},
            )
            return [result], {
                "cell_id": cell_id,
                "read_layout": read_layout,
                "status": "skipped_existing",
                "seconds": 0.0,
                "tsv_path": str(existing_path),
                "execution_mode": "per-cell-fastq",
                "input_mode": "fastq",
                "reused_alignment": False,
                "detector_backend": detector.name,
            }

        cell_started = time.perf_counter()
        inputs = DetectorRunInputs(
            cell_id=cell_id,
            r1=r1,
            r2=r2,
            outdir=outdir,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
            extra={},
        )
        raw_result = detector.run(inputs)
        flat_results = _flatten_detector_result(raw_result)
        if not flat_results:
            raise RuntimeError(f"Detector returned no results for cell_id={cell_id}")

        elapsed = round(time.perf_counter() - cell_started, 3)
        primary_path = flat_results[0].tsv_path
        normalized_rows = _count_output_rows(primary_path, detector.name)
        raw_output_path = flat_results[0].meta.get("raw_output_path") if flat_results[0].meta else None
        raw_rows = _count_output_rows(Path(raw_output_path), detector.name) if raw_output_path else None
        status = "success" if (normalized_rows or 0) > 0 else "empty"
        return flat_results, {
            "cell_id": cell_id,
            "read_layout": read_layout,
            "status": status,
            "outcome_category": _outcome_category(status, raw_rows=raw_rows, normalized_rows=normalized_rows),
            "seconds": elapsed,
            "tsv_path": str(primary_path),
            "log_path": str(flat_results[0].log_path) if flat_results[0].log_path else None,
            "execution_mode": "per-cell-fastq",
            "input_mode": "fastq",
            "reused_alignment": False,
            "detector_backend": detector.name,
            "reference_used": str(ref_fa) if ref_fa else None,
            "raw_output_path": raw_output_path,
            "raw_row_count": raw_rows,
            "normalized_row_count": normalized_rows,
        }

    results: list[DetectorResult] = []
    failures: list[dict[str, Any]] = []

    with ThreadPoolExecutor(max_workers=effective_parallel) as ex:
        future_map = {ex.submit(_run_one, row): row for row in rows}
        for fut in as_completed(future_map):
            cell_id, row_r1, row_r2 = future_map[fut]
            try:
                cell_results, record = fut.result()
                for cell_result in cell_results:
                    if _detector_output_has_calls(cell_result.tsv_path, detector.name):
                        _write_provenance_stamp(
                            cell_result.tsv_path,
                            _build_provenance_stamp(
                                detector,
                                cell_id=cell_result.cell_id,
                                r1=row_r1,
                                r2=row_r2,
                                ref_fa=ref_fa,
                                gtf=gtf,
                                threads=threads,
                                extra={},
                            ),
                        )
                results.extend(cell_results)
                per_cell_records.append(record)
            except Exception as exc:
                read_layout = "paired-end" if row_r2 is not None else "single-end"
                failures.append({"cell_id": cell_id, "read_layout": read_layout, "error": str(exc)})
                per_cell_records.append(
                    {
                        "cell_id": cell_id,
                        "read_layout": read_layout,
                        "status": "failed",
                        "outcome_category": "failed",
                        "seconds": None,
                        "error": str(exc),
                        "execution_mode": "per-cell-fastq",
                        "input_mode": "fastq",
                        "reused_alignment": False,
                        "detector_backend": detector.name,
                        "reference_used": str(ref_fa) if ref_fa else None,
                        "raw_output_path": None,
                        "raw_row_count": None,
                        "normalized_row_count": None,
                    }
                )

    per_cell_records.sort(key=lambda item: item["cell_id"])
    status_counts: dict[str, int] = {}
    for record in per_cell_records:
        status = str(record.get("status", "unknown"))
        status_counts[status] = status_counts.get(status, 0) + 1

    payload = {
        "detector": detector.name,
        "manifest": str(manifest.resolve()),
        "outdir": str(outdir.resolve()),
        "threads": threads,
        "parallel_requested": parallel,
        "parallel_effective": effective_parallel,
        "n_manifest_rows": len(rows),
        "status_counts": status_counts,
        "elapsed_seconds": round(time.perf_counter() - started_at, 3),
        "execution_mode": "per-cell-fastq",
        "input_mode": "fastq",
        "cells": per_cell_records,
    }
    if failures:
        payload["failures"] = failures
    _write_summary(summary_path, payload)

    success = status_counts.get("success", 0)
    skipped = status_counts.get("skipped_existing", 0)
    empty = status_counts.get("empty", 0)
    failed = status_counts.get("failed", 0)
    print(
        f"[circyto] Detector '{detector.name}' summary: "
        f"success={success} skipped={skipped} empty={empty} failed={failed} "
        f"(summary: {summary_path})"
    )

    if failures:
        failed_cells = ", ".join(item["cell_id"] for item in failures[:5])
        suffix = "" if len(failures) <= 5 else f" (+{len(failures) - 5} more)"
        raise RuntimeError(
            f"Detector '{detector.name}' failed for {len(failures)}/{len(rows)} cells: "
            f"{failed_cells}{suffix}. Summary: {summary_path}"
        )

    return results

def run_multidetector(
    detectors: Dict[str, DetectorBase],
    manifest: Path,
    root_outdir: Path,
    ref_fa: Optional[Path],
    gtf: Optional[Path],
    threads: int = 8,
    parallel: int = 4,
) -> Dict[str, List[DetectorResult]]:
    """
    Run multiple detectors over the same manifest.

    Layout:

      root_outdir/
        <detector_name>/
          <cell>.tsv

    Returns a dict: detector_name -> list[DetectorResult]
    """
    ensure_dir(root_outdir)

    results: Dict[str, List[DetectorResult]] = {}

    for name, det in detectors.items():
        det_out = root_outdir / name
        print(f"[circyto] Running detector {name} into {det_out}")
        res = run_detector_manifest(
            detector=det,
            manifest=manifest,
            outdir=det_out,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
            parallel=parallel,
        )
        results[name] = res

    return results
