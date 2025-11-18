from __future__ import annotations
from pathlib import Path
from typing import List, Dict

from circyto.detectors import build_default_engines
from circyto.pipeline.run_detector import run_detector_manifest
from circyto.utils import ensure_dir


def run_multidetector_pipeline(
    detectors: List[str],
    manifest: Path,
    outdir: Path,
    ref_fa: Path,
    gtf: Path,
    threads: int = 8,
    parallel: int = 1,
) -> Dict[str, dict]:
    """
    Internal orchestrator used by CLI.
    Returns dictionary of per-detector metadata.
    """
    ensure_dir(outdir)
    engines = build_default_engines()

    meta = {}

    for det in detectors:
        if det not in engines:
            raise ValueError(f"Detector '{det}' is not available")

        det_engine = engines[det]
        det_outdir = outdir / det
        ensure_dir(det_outdir)

        print(f"[multidetector] >>> Running {det}")
        results = run_detector_manifest(
            detector=det_engine,
            manifest=manifest,
            outdir=det_outdir,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
            parallel=parallel,
        )

        meta[det] = {
            "detector": det,
            "outdir": str(det_outdir),
            "n_cells": len(results),
        }

    return meta
