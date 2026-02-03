# tests/test_smoke_smartseq2_cli.py
from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_min_fastq_gz(path: Path, n_reads: int = 2) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        for i in range(n_reads):
            f.write(f"@read{i}\n")
            f.write("A" * 50 + "\n")
            f.write("+\n")
            f.write("I" * 50 + "\n")


def _write_barcodes_tsv(path: Path) -> None:
    path.write_text("sc01\tAACGTGAT\nsc02\tAAACATCG\n", encoding="utf-8")


def test_read_manifest_accepts_r1_r2(tmp_path: Path):
    m = tmp_path / "manifest_legacy.tsv"
    m.write_text(
        "cell_id\tr1\tr2\n"
        "sc01\t/a/b/c_R1.fastq.gz\t/a/b/c_R2.fastq.gz\n",
        encoding="utf-8",
    )

    from circyto.pipeline.run_detector import read_manifest

    rows = read_manifest(m)
    assert rows[0][0] == "sc01"
    assert rows[0][1].name.endswith("_R1.fastq.gz")
    assert rows[0][2].name.endswith("_R2.fastq.gz")


def test_read_manifest_accepts_read1_read2(tmp_path: Path):
    m = tmp_path / "manifest_v1.tsv"
    m.write_text(
        "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\n"
        "sc01\tsmartseq2\t/a/b/sc01_R1.fastq.gz\t/a/b/sc01_R2.fastq.gz\t\tLIB\t100\n",
        encoding="utf-8",
    )

    from circyto.pipeline.run_detector import read_manifest

    rows = read_manifest(m)
    assert rows[0][0] == "sc01"
    assert rows[0][1].name == "sc01_R1.fastq.gz"
    assert rows[0][2].name == "sc01_R2.fastq.gz"


@pytest.mark.parametrize("barcode_read", ["R2"])
def test_smoke_smartseq2_cli_wiring(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, barcode_read: str):
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    bar = tmp_path / "barcodes.tsv"
    outdir = tmp_path / "out"
    ref_fa = tmp_path / "chr21.fa"
    gtf = tmp_path / "genes.gtf"

    _write_min_fastq_gz(r1)
    _write_min_fastq_gz(r2)
    _write_barcodes_tsv(bar)

    ref_fa.write_text(">chr21\n" + "A" * 1000 + "\n", encoding="utf-8")
    gtf.write_text(
        'chr21\ttest\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        encoding="utf-8",
    )

    # --- monkeypatch demux (IMPORTANT: patch where smoke imports it) ---
    def fake_demux(params):
        out_fastq = params.outdir / "fastq"
        out_fastq.mkdir(parents=True, exist_ok=True)

        for cid in ["sc01", "sc02"]:
            for mate in ["R1", "R2"]:
                _write_min_fastq_gz(out_fastq / f"{cid}_{mate}.fastq.gz", n_reads=1)

        # write manifest exactly where smoke expects it
        params.manifest_path.write_text(
            "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\n"
            f"sc01\tsmartseq2\t{out_fastq/'sc01_R1.fastq.gz'}\t{out_fastq/'sc01_R2.fastq.gz'}\t\tLIB\t10\n"
            f"sc02\tsmartseq2\t{out_fastq/'sc02_R1.fastq.gz'}\t{out_fastq/'sc02_R2.fastq.gz'}\t\tLIB\t10\n",
            encoding="utf-8",
        )
        return {"sc01": {"assigned": 10}, "sc02": {"assigned": 10}}

    monkeypatch.setattr("circyto.cli.smoke.demux_smartseq2_pooled", fake_demux, raising=True)

    # --- monkeypatch run so we don't execute real tools ---
    def fake_run_detector_manifest(*, detector, manifest, outdir, ref_fa, gtf, threads, parallel):
        assert Path(manifest).exists()
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        # create dummy per-cell outputs (collector will be patched anyway)
        (out / "sc01").mkdir(parents=True, exist_ok=True)
        (out / "sc02").mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr("circyto.cli.smoke.run_detector_manifest", fake_run_detector_manifest, raising=True)

    # --- monkeypatch collect (signature in smoke is legacy collect_matrix(cirifull_dir, mtx, idx, idx)) ---
    def fake_collect(indir, matrix_path, circ_index_path, cell_index_path, *args, **kwargs):
        Path(matrix_path).write_text("%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n", encoding="utf-8")
        Path(circ_index_path).write_text("chr21:1-10:+\n", encoding="utf-8")
        Path(cell_index_path).write_text("sc01\n", encoding="utf-8")

    monkeypatch.setattr("circyto.cli.smoke.collect_matrix", fake_collect, raising=True)
    monkeypatch.setattr("circyto.cli.smoke.collect_find_circ3_matrix", fake_collect, raising=True)

    # --- monkeypatch convert ---
    def fake_convert_matrix_files(*, matrix_path, circ_index_path, cell_index_path, loom=None, h5ad=None):
        assert Path(matrix_path).exists()
        assert Path(circ_index_path).exists()
        assert Path(cell_index_path).exists()
        assert h5ad is not None
        Path(h5ad).write_text("dummy_h5ad\n", encoding="utf-8")

    monkeypatch.setattr("circyto.cli.smoke.convert_matrix_files", fake_convert_matrix_files, raising=True)

    res = runner.invoke(
        app,
        [
            "smoke",
            "smartseq2",
            "--r1", str(r1),
            "--r2", str(r2),
            "--barcodes", str(bar),
            "--outdir", str(outdir),
            "--library-id", "LIB",
            "--barcode-read", barcode_read,
            "--barcode-start", "0",
            "--barcode-length", "8",
            "--ref-fa", str(ref_fa),
            "--gtf", str(gtf),
            "--detector", "ciri2",
            "--limit-reads", "10",
            "--overwrite",
        ],
    )

    assert res.exit_code == 0, res.stdout + "\n" + (res.stderr or "")

    # sanity: output file exists
    assert (outdir / "h5ad" / "ciri2.h5ad").exists()
