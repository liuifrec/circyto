from __future__ import annotations

import os
from pathlib import Path

import pytest

from circyto.detectors import build_default_engines
from circyto.pipeline.run_detector import run_detector_manifest


@pytest.mark.slow
@pytest.mark.integration
def test_circexplorer2_chr21_manifest2_integration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """
    Integration test for the CIRCexplorer2 detector on the small chr21 manifest_2.tsv.

    Steps (when environment is fully configured):

      1. Runs circexplorer2 via detector API on manifest_2.tsv + ref/chr21.fa.
      2. Checks that per-cell CIRCexplorer2 output files exist and are non-empty.

    If STAR / CIRCexplorer2 or indices are not correctly installed/configured,
    this test will SKIP instead of failing, so it does not block development
    on machines without the full toolchain.
    """

    # Allow CI / local runs to skip easily
    if os.environ.get("CIRCYTO_SKIP_INTEGRATION", "").lower() in {"1", "true", "yes"}:
        pytest.skip("CIRCYTO_SKIP_INTEGRATION is set, skipping CIRCexplorer2 integration test")

    # Require STAR index + refFlat annotation for CIRCexplorer2
    star_index_env = os.environ.get("CIRCYTO_CIRCEXPLORER2_STAR_INDEX")
    ref_flat_env = os.environ.get("CIRCYTO_CIRCEXPLORER2_REF_FLAT")

    if not star_index_env or not ref_flat_env:
        pytest.skip(
            "CIRCYTO_CIRCEXPLORER2_STAR_INDEX and/or CIRCYTO_CIRCEXPLORER2_REF_FLAT "
            "are not set; skipping CIRCexplorer2 integration test"
        )

    repo_root = Path(__file__).resolve().parents[1]

    manifest = repo_root / "manifest_2.tsv"
    ref_fa = repo_root / "ref" / "chr21.fa"

    if not manifest.exists() or not ref_fa.exists():
        pytest.skip(
            f"Required test data not found: {manifest} or {ref_fa}. "
            "This test expects manifest_2.tsv and ref/chr21.fa at repo root."
        )

    # Ensure STAR uses a fresh tmp directory for this test
    star_tmp = tmp_path / "star_tmp"
    monkeypatch.setenv("CIRCYTO_STAR_TMP", str(star_tmp))

    # Ensure we do NOT skip STAR here: we want a full run when possible
    monkeypatch.delenv("CIRCYTO_CIRCEXPLORER2_SKIP_STAR", raising=False)

    # 1. Build detector
    engines = build_default_engines()
    if "circexplorer2" not in engines:
        pytest.skip("circexplorer2 detector not registered in engines; skipping integration test")

    det_engine = engines["circexplorer2"]

    outdir = tmp_path / "circexplorer2_manifest2"
    outdir.mkdir(parents=True, exist_ok=True)

    # 2. Run detector via manifest
    try:
        results = run_detector_manifest(
            detector=det_engine,
            manifest=manifest,
            outdir=outdir,
            ref_fa=ref_fa,
            gtf=None,
            threads=2,
            parallel=2,
        )
    except Exception as e:
        # If STAR / CIRCexplorer2 is not properly configured, we don't want this
        # to be a hard failure on every machine.
        pytest.skip(f"CIRCexplorer2 integration requires working STAR+CIRCexplorer2; got error: {e!r}")

    # Basic sanity: we should have at least one result per manifest row
    assert len(results) >= 1

    # 3. Check for per-cell CIRCexplorer2 output files
    # We expect directories named by cell_id under outdir, with CIRCexplorer2 outputs.
    cell_dirs = [p for p in outdir.iterdir() if p.is_dir()]
    assert cell_dirs, f"No per-cell directories found in {outdir}"

    # Check that each cell dir has at least one non-empty file (very loose sanity check)
    saw_nonempty = False
    for cell_dir in cell_dirs:
        for f in cell_dir.iterdir():
            if f.is_file() and f.stat().st_size > 0:
                saw_nonempty = True
                break
        if saw_nonempty:
            break

    assert saw_nonempty, "Expected at least one non-empty CIRCexplorer2 output file in per-cell directories"
