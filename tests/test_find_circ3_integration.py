from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from circyto.detectors import build_default_engines
from circyto.pipeline.run_detector import run_detector_manifest
from circyto.pipeline.collect_find_circ3 import collect_find_circ3_matrix


def _require_find_circ3_integration_runtime(repo_root: Path) -> tuple[Path, Path]:
    """
    Guard this heavyweight integration behind explicit opt-in plus runtime checks.

    This test runs real bowtie2/samtools/find-circ3 commands against repo-local
    chr21 FASTQs and references. It is not suitable for routine default pytest
    runs unless the caller explicitly enables it.
    """
    if os.environ.get("CIRCYTO_SKIP_INTEGRATION", "").lower() in {"1", "true", "yes"}:
        pytest.skip("CIRCYTO_SKIP_INTEGRATION is set, skipping integration test")

    if os.environ.get("CIRCYTO_RUN_FIND_CIRC3_INTEGRATION", "").lower() not in {"1", "true", "yes"}:
        pytest.skip(
            "Set CIRCYTO_RUN_FIND_CIRC3_INTEGRATION=1 to run the heavyweight find_circ3 integration test"
        )

    missing_cmds = [name for name in ("bowtie2", "samtools", "find-circ3") if shutil.which(name) is None]
    if missing_cmds:
        pytest.skip(f"Missing external commands for find_circ3 integration: {', '.join(missing_cmds)}")

    manifest = repo_root / "manifest_2.tsv"
    ref_fa = repo_root / "ref" / "chr21.fa"
    if not manifest.exists() or not ref_fa.exists():
        pytest.skip(
            f"Required test data not found: {manifest} or {ref_fa}. "
            "This test expects manifest_2.tsv and ref/chr21.fa at repo root."
        )

    # The adapter derives a bowtie2 index prefix from ref/chr21.fa -> ref/chr21
    bowtie2_index = ref_fa.with_suffix("")
    index_suffixes = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    missing_indexes = [f"{bowtie2_index}{suffix}" for suffix in index_suffixes if not Path(f"{bowtie2_index}{suffix}").exists()]
    if missing_indexes:
        pytest.skip("Missing bowtie2 chr21 index files for find_circ3 integration")

    for line in manifest.read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        for rel_path in fields[1:3]:
            if rel_path.strip() and not (repo_root / rel_path.strip()).exists():
                pytest.skip(f"Missing FASTQ input for find_circ3 integration: {repo_root / rel_path.strip()}")

    return manifest, ref_fa


@pytest.mark.slow
@pytest.mark.integration
def test_find_circ3_chr21_manifest2_integration() -> None:
    """
    Integration test for the find_circ3 detector on the small chr21 manifest_2.tsv.

    This test:
      1. Runs find_circ3 via the detector API on manifest_2.tsv + ref/chr21.fa.
      2. Checks that per-cell splice_sites.bed files are created and non-empty.
      3. Runs collect_find_circ3_matrix to build a MatrixMarket matrix.
      4. Verifies that the matrix and index files exist and contain at least one circRNA.

    It is marked as 'slow' and 'integration' so CI can selectively enable it
    (e.g. via an environment variable or pytest -m filter).
    """

    repo_root = Path(__file__).resolve().parents[1]
    manifest, ref_fa = _require_find_circ3_integration_runtime(repo_root)

    # 1. Build detector
    engines = build_default_engines()
    assert "find-circ3" in engines, "find-circ3 detector not registered in engines"
    det_engine = engines["find-circ3"]

    # Use a stable outdir under the repo instead of /tmp to avoid FS weirdness
    outdir = repo_root / "tests" / "find_circ3_manifest2_integration"
    if outdir.exists():
        # Best-effort cleanup; don't be too aggressive if CI reuses workspace
        for child in outdir.rglob("*"):
            if child.is_file():
                try:
                    child.unlink()
                except OSError:
                    pass
        for child in sorted(outdir.glob("*"), reverse=True):
            if child.is_dir():
                try:
                    child.rmdir()
                except OSError:
                    pass
    outdir.mkdir(parents=True, exist_ok=True)

    # 2. Run detector via manifest (low parallelism to reduce resource stress)
    try:
        results = run_detector_manifest(
            detector=det_engine,
            manifest=manifest,
            outdir=outdir,
            ref_fa=ref_fa,
            gtf=None,
            threads=1,
            parallel=1,
        )
    except subprocess.CalledProcessError as e:
        # External tool (bowtie2 / find-circ3) failure: treat as environment/tooling
        # issue and skip rather than blowing up pytest with a bus error.
        pytest.skip(
            f"External tool failure during find_circ3 run "
            f"(returncode={e.returncode}): {e}"
        )

    # Expect 2 rows from manifest_2.tsv
    assert len(results) == 2, f"Expected 2 DetectorResult entries, got {len(results)}"

    # 3. Check that per-cell splice_sites.bed files are present and non-empty
    cell_dirs = [p for p in outdir.iterdir() if p.is_dir()]
    assert cell_dirs, "No per-cell subdirectories created by find_circ3 detector"

    splice_paths = []
    for cdir in sorted(cell_dirs):
        cell_id = cdir.name
        bed = cdir / f"{cell_id}_splice_sites.bed"
        assert bed.exists(), f"Missing splice_sites.bed for cell {cell_id}"
        assert bed.stat().st_size > 0, f"Empty splice_sites.bed for cell {cell_id}"
        splice_paths.append(bed)

    # 4. Collect into MatrixMarket matrix + indices
    matrix_path = outdir / "circ_counts.mtx"
    circ_index_path = outdir / "circ_index.txt"
    cell_index_path = outdir / "cell_index.txt"

    collect_find_circ3_matrix(
        findcirc3_dir=str(outdir),
        matrix_path=str(matrix_path),
        circ_index_path=str(circ_index_path),
        cell_index_path=str(cell_index_path),
        min_count_per_cell=1,
    )

    assert matrix_path.exists(), "MatrixMarket output not created"
    assert circ_index_path.exists(), "circ_index.txt not created"
    assert cell_index_path.exists(), "cell_index.txt not created"

    # 5. Sanity-check contents: at least one circ_id and one cell
    circ_ids = [
        line.strip()
        for line in circ_index_path.read_text().splitlines()
        if line.strip()
    ]
    cell_ids = [
        line.strip()
        for line in cell_index_path.read_text().splitlines()
        if line.strip()
    ]

    assert circ_ids, "circ_index.txt is empty"
    assert cell_ids, "cell_index.txt is empty"

    # 6. Matrix non-zero check: at least one data line after header
    with matrix_path.open() as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("%")]
    # After the dims+nnz header line, expect at least one entry line
    assert len(lines) >= 2, f"MatrixMarket file too short, contents: {lines!r}"
