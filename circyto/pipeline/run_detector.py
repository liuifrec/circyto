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
    ref_fa: Optional[Path],
    gtf: Optional[Path],
    threads: int = 8,
    parallel: int = 4,
) -> List[DetectorResult]:
    """
    Dispatch a single detector across all rows in a manifest file.
    """

    ensure_dir(outdir)

    rows = read_manifest(manifest)
    results: List[DetectorResult] = []

    def _run_one(row):
        cell_id, r1, r2 = row
        inputs = DetectorRunInputs(
            cell_id=cell_id,
            r1=r1,
            r2=r2,
            outdir=outdir,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
        )
        return detector.run(inputs)

    with ThreadPoolExecutor(max_workers=parallel) as ex:
        futures = [ex.submit(_run_one, r) for r in rows]
        for fut in futures:
            results.append(fut.result())

    return results
