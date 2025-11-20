# circyto/pipeline/run_multidetector.py
from __future__ import annotations
import json
from pathlib import Path
from typing import List, Dict

from .run_detector import run_detector_manifest
from circyto.detectors import build_default_engines
from circyto.utils import ensure_dir


def run_multidetector_pipeline(
    detectors: List[str],
    manifest: Path,
    outdir: Path,
    ref_fa: Path | None,
    gtf: Path | None,
    threads: int = 8,
    parallel: int = 1,
) -> Dict[str, Dict]:

    outdir = Path(outdir)
    ensure_dir(outdir)

    engines = build_default_engines()

    summary: Dict[str, Dict] = {}

    for det_name in detectors:
        if det_name not in engines:
            raise ValueError(f"Detector '{det_name}' not available. Available: {list(engines.keys())}")

        det = engines[det_name]
        print(f"[multidetector] >>> Running {det_name}")

        det_out = outdir / det_name
        ensure_dir(det_out)

        results = run_detector_manifest(
            detector=det,
            manifest=manifest,
            outdir=det_out,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
            parallel=parallel,
        )

        summary[det_name] = {
            "detector": det_name,
            "n_cells": len(results),
            "outdir": str(det_out),
        }

    # ALWAYS write summary.json
    summary_path = outdir / "summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)

    print(f"[multidetector] Summary written to {summary_path}")

    return summary
