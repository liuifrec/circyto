# circyto/pipeline/run_multidetector.py

from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

from circyto.detectors import available_detectors
from circyto.pipeline.run_detector import run_detector_manifest


def _flatten_results(detector_results):
    """
    Normalize detector_results into a flat List[DetectorResult].

    Some engines may (by mistake) return a List[List[DetectorResult]].
    Others correctly return List[DetectorResult].

    This helper makes both cases consistent.
    """
    flat = []
    for r in detector_results:
        if isinstance(r, (list, tuple)):
            flat.extend(r)
        else:
            flat.append(r)
    return flat


def run_multidetector_pipeline(
    detectors: Sequence[str],
    manifest: Path,
    outdir: Path,
    ref_fa: Path | None = None,
    gtf: Path | None = None,
    threads: int = 4,
    parallel: int = 1,
) -> dict[str, list]:
    """
    Run multiple detectors over the same manifest.

    Parameters
    ----------
    detectors
        Iterable of detector names (e.g. ["ciri-full", "find-circ3", "circexplorer2"]).
    manifest
        Path to a manifest TSV with columns like: cell_id, r1, r2.
    outdir
        Root output directory. Per-detector subdirectories are created here.
    ref_fa
        Reference FASTA passed through to run_detector_manifest.
    gtf
        Optional GTF for detectors that need it.
    threads
        Number of threads per detector run (forwarded to the detector).
    parallel
        Max number of cells to process in parallel for each detector.

    Returns
    -------
    dict[str, list]
        Mapping from detector name to the list of DetectorResult objects
        returned by run_detector_manifest for that detector.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    engines = available_detectors()

    # Validate requested detectors.
    missing = [d for d in detectors if d not in engines]
    if missing:
        available = ", ".join(sorted(engines.keys()))
        raise ValueError(
            f"Detector(s) not available: {', '.join(missing)}. "
            f"Available: {available}"
        )

    # What we return to callers/tests
    result_map: dict[str, list] = {}

    # What we serialize into summary.json (pure JSON-serializable)
    summary_json: dict[str, dict] = {}
    results_json: dict[str, list] = {}

    for det_name in detectors:
        det_engine = engines[det_name]
        det_outdir = outdir / det_name
        det_outdir.mkdir(parents=True, exist_ok=True)

        # Use the standard manifest runner (handles threading/parallel)
        raw_results = run_detector_manifest(
            detector=det_engine,
            manifest=manifest,
            outdir=det_outdir,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
            parallel=parallel,
        )

        # Normalize to a flat list of DetectorResult
        detector_results = _flatten_results(raw_results)

        # Populate return object
        result_map[det_name] = detector_results

        # Build lightweight JSON view
        results_json[det_name] = [
            {
                "cell_id": r.cell_id,
                "tsv_path": str(getattr(r, "tsv_path", "")) if getattr(r, "tsv_path", None) else None,
                "outdir": str(det_outdir),
            }
            for r in detector_results
        ]

        summary_json[det_name] = {
            "detector": det_name,
            "n_cells": len(detector_results),
            "outdir": str(det_outdir),
        }

    payload = {
        "summary": summary_json,
        "results": results_json,
    }

    summary_path = outdir / "summary.json"
    with summary_path.open("w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)

    return result_map
