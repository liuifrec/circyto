from __future__ import annotations

import os
from pathlib import Path

import pytest

from circyto.detectors import build_default_engines
from circyto.pipeline.run_detector import run_detector_manifest
from circyto.pipeline.run_detector import read_manifest


@pytest.mark.slow
@pytest.mark.integration
def test_ciri2_single_end_chr21_known_positive() -> None:
    if os.environ.get("CIRCYTO_SKIP_INTEGRATION", "").lower() in {"1", "true", "yes"}:
        pytest.skip("CIRCYTO_SKIP_INTEGRATION is set, skipping integration test")

    repo_root = Path(__file__).resolve().parents[1]
    manifest = repo_root / "work" / "err1041421_se" / "manifest.tsv"
    ref_fa = repo_root / "ref" / "chr21.fa"
    gtf = repo_root / "ref" / "chr21.gtf"

    if not manifest.exists() or not ref_fa.exists() or not gtf.exists():
        pytest.skip("ERR1041421 SE manifest or chr21 reference assets are missing")
    try:
        read_manifest(manifest, validate_files=True)
    except FileNotFoundError as exc:
        pytest.skip(f"ERR1041421 SE manifest references missing local FASTQ assets: {exc}")

    engines = build_default_engines()
    det = engines["ciri2"]

    outdir = repo_root / "tests" / "ciri2_err1041421_integration"
    outdir.mkdir(parents=True, exist_ok=True)

    results = run_detector_manifest(
        detector=det,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=2,
        parallel=1,
    )

    assert len(results) == 1

    tsv_path = results[0].tsv_path
    lines = [line.strip() for line in tsv_path.read_text().splitlines() if line.strip()]
    assert len(lines) > 1, "CIRI2 SE known-positive regression returned header-only TSV"
    assert any("chr21\t39224408\t39238573\t-\t1" in line for line in lines[1:])
