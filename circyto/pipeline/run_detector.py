# circyto/pipeline/run_detector.py
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from circyto.detectors import DetectorBase, DetectorRunInputs, DetectorResult

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


def read_manifest(path: Path) -> List[Tuple[str, Path, Optional[Path]]]:
    """
    Read manifest rows as: (cell_id, r1_path, r2_path_or_None).

    Supports:
      - legacy:  cell_id, r1, r2
      - v1:      cell_id, read1, read2, (plus extra columns)

    Notes:
      - We intentionally do NOT require r2 (some platforms may be single-end).
      - We do NOT force paths to exist here (HPC/remote mounts can be validated elsewhere).
    """
    import csv

    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")

    rows: List[Tuple[str, Path, Optional[Path]]] = []
    with path.open("r", newline="") as f:
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

            r1_key = _pick_col(r, ("r1", "read1"))
            r2_key = _pick_col(r, ("r2", "read2"))

            if r1_key is None:
                raise KeyError(f"Manifest missing r1/read1 for cell_id={cell} at {path}:{i}")

            r1 = Path(str(r[r1_key]).strip())
            r2 = Path(str(r[r2_key]).strip()) if r2_key else None

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
