# circyto/pipeline/run_detector.py
from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Optional, Dict

from circyto.detectors import DetectorBase, DetectorRunInputs, DetectorResult

def ensure_dir(path: Path) -> None:
    """
    Create the directory if it doesn't exist.
    Small local helper to avoid importing non-existent utils modules.
    """
    path.mkdir(parents=True, exist_ok=True)


def read_manifest(path: Path) -> List[tuple[str, Path, Optional[Path]]]:
    """
    Simple manifest reader.
    Expected columns: cell_id, r1, r2 (optional)
    """
    rows = []
    with path.open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            cell = r["cell_id"]
            r1 = Path(r["r1"])
            r2 = Path(r["r2"]) if "r2" in r and r["r2"] else None
            rows.append((cell, r1, r2))
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
    rows = read_manifest(manifest)
    outdir.mkdir(parents=True, exist_ok=True)

    # NEW: limit effective parallelism based on detector capability
    det_max_parallel = getattr(detector, "max_parallel", parallel)
    effective_parallel = min(parallel, det_max_parallel)
    if effective_parallel < parallel:
        print(
            f"[circyto] Detector '{detector.name}' only supports parallel={effective_parallel}; "
            f"overriding requested parallel={parallel}."
        )

    def _run_one(row: tuple[str, Path, Path | None]) -> DetectorResult:
        cell_id, r1, r2 = row
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
        return detector.run(inputs)

    results: list[DetectorResult] = []

    # Use the effective parallelism
    with ThreadPoolExecutor(max_workers=effective_parallel) as ex:
        futures = [ex.submit(_run_one, r) for r in rows]
        for fut in futures:
            results.append(fut.result())

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
