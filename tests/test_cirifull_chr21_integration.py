import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("[TEST] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def test_cirifull_chr21_nonempty_matrix():
    """
    Integration test:
    - run-manifest on a small chr21 Smart-seq2 subset (manifest_2.tsv)
    - collect CIRI-full TSVs into a sparse matrix
    - assert matrix is non-empty and cell count matches manifest
    """
    work_dir = ROOT
    manifest = work_dir / "manifest_2.tsv"
    outdir = work_dir / "work_smartseq2" / "ciri_full_chr21_test"
    ref_fa = work_dir / "ref" / "chr21.fa"
    gtf = work_dir / "ref" / "chr21.gtf"

    # 1) run-manifest (CIRI-full)
    run(
        [
            "circyto",
            "run-manifest",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--ref-fa",
            str(ref_fa),
            "--gtf",
            str(gtf),
            "--cmd-template",
            'bash -lc "'
            "export R1={r1} R2={r2} REF_FA={ref_fa} GTF={gtf} "
            "OUT_TSV={out_tsv} THREADS=4 ; "
            "tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh"
            '"',
        ],
        cwd=work_dir,
    )

    # 2) collect into matrix
    mtx = work_dir / "work_smartseq2" / "circ_chr21_test.mtx"
    circ_ids = work_dir / "work_smartseq2" / "circ_chr21_test_ids.txt"
    cell_ids = work_dir / "work_smartseq2" / "cell_chr21_test_ids.txt"

    run(
        [
            "circyto",
            "collect",
            "--cirifull-dir",
            str(outdir),
            "--matrix",
            str(mtx),
            "--circ-index",
            str(circ_ids),
            "--cell-index",
            str(cell_ids),
        ],
        cwd=work_dir,
    )

    # 3) assertions: files exist
    assert mtx.is_file(), "Matrix file was not created"
    assert circ_ids.is_file(), "circ index file missing"
    assert cell_ids.is_file(), "cell index file missing"

    # parse MatrixMarket header (first non-comment line)
    with mtx.open() as f:
        for line in f:
            if line.startswith("%"):
                continue
            n_circ, n_cells, nnz = map(int, line.split())
            break

    assert nnz > 0, "Matrix has zero non-zero entries (expected non-empty CIRI-full output)"

    # cells in matrix vs manifest_2.tsv (minus header)
    manifest_cells = sum(1 for _ in open(manifest)) - 1
    assert n_cells == manifest_cells, f"Expected {manifest_cells} cells, got {n_cells}"
