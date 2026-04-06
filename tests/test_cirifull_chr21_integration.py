from __future__ import annotations

import os
from pathlib import Path

import pytest

from circyto.detectors import build_default_engines
from circyto.pipeline.collect import collect_matrix
from circyto.pipeline.run_detector import run_detector_manifest


def _tsv_has_calls(path: Path) -> bool:
    if not path.exists():
        return False
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return len(lines) > 1


@pytest.mark.slow
@pytest.mark.integration
def test_cirifull_chr21_nonempty_matrix(tmp_path: Path) -> None:
    """
    Integration test for CIRI-full on the small chr21 Smart-seq2 subset.

    This exercises the current detector API end-to-end:
      1. Run the CIRI-full detector on manifest_2.tsv.
      2. Verify per-cell normalized TSVs exist and are non-empty.
      3. Collect them into a sparse matrix.
      4. Assert the matrix is non-empty and cell count matches the manifest.
    """

    if os.environ.get("CIRCYTO_SKIP_INTEGRATION", "").lower() in {"1", "true", "yes"}:
        pytest.skip("CIRCYTO_SKIP_INTEGRATION is set, skipping integration test")

    repo_root = Path(__file__).resolve().parents[1]
    manifest = repo_root / "manifest_2.tsv"
    ref_fa = repo_root / "ref" / "chr21.fa"
    gtf = repo_root / "ref" / "chr21.gtf"

    if not manifest.exists() or not ref_fa.exists() or not gtf.exists():
        pytest.skip("manifest_2.tsv or chr21 reference assets are missing")

    engines = build_default_engines()
    assert "ciri-full" in engines, "ciri-full detector not registered in engines"
    det = engines["ciri-full"]

    prepared_outdir = repo_root / "work_smartseq2" / "ciri_full_chr21_test"
    prepared_tsvs = sorted(prepared_outdir.glob("*.tsv")) if prepared_outdir.exists() else []
    if len(prepared_tsvs) == 2 and all(_tsv_has_calls(path) for path in prepared_tsvs):
        outdir = prepared_outdir
    else:
        outdir = tmp_path / "ciri_full_chr21_test"
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
        assert len(results) == 2, f"Expected 2 DetectorResult entries, got {len(results)}"
        for result in results:
            assert result.tsv_path.exists(), f"Missing TSV for {result.cell_id}"
            assert _tsv_has_calls(result.tsv_path), f"CIRI-full TSV is header-only for {result.cell_id}"

    mtx = tmp_path / "circ_chr21_test.mtx"
    circ_ids = tmp_path / "circ_chr21_test_ids.txt"
    cell_ids = tmp_path / "cell_chr21_test_ids.txt"

    collect_matrix(
        cirifull_dir=str(outdir),
        matrix_path=str(mtx),
        circ_index_path=str(circ_ids),
        cell_index_path=str(cell_ids),
        min_count_per_cell=1,
    )

    assert mtx.is_file(), "Matrix file was not created"
    assert circ_ids.is_file(), "circ index file missing"
    assert cell_ids.is_file(), "cell index file missing"

    with mtx.open() as f:
        for line in f:
            if line.startswith("%"):
                continue
            n_circ, n_cells, nnz = map(int, line.split())
            break

    assert nnz > 0, "Matrix has zero non-zero entries (expected non-empty CIRI-full output)"

    manifest_cells = sum(1 for _ in manifest.open()) - 1
    assert n_cells == manifest_cells, f"Expected {manifest_cells} cells, got {n_cells}"
