from __future__ import annotations

from pathlib import Path

import pytest

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_manifest_v1(path: Path, r1: Path, r2: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\t".join(["cell_id", "platform", "read1", "read2", "bam", "library_id", "n_input_reads"]) + "\n"
        + "\t".join(["sc01", "smartseq2", str(r1), str(r2), "", "LIB", "10"]) + "\n"
    )


@pytest.mark.parametrize("detector", ["find-circ3", "ciri2"])
def test_smoke_smartseq2_wiring(monkeypatch, tmp_path: Path, detector: str):
    # Inputs (dummy)
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    barcodes = tmp_path / "barcodes.tsv"
    ref_fa = tmp_path / "ref.fa"
    r1.write_text("x")
    r2.write_text("y")
    barcodes.write_text("sc01\tAACGTGAT\n")
    ref_fa.write_text(">chr21\nACGT\n")

    outdir = tmp_path / "work"

    # --- monkeypatch build_default_engines to avoid pulling real engines ---
    class _FakeDet:
        def __init__(self, name: str):
            self.name = name
            self.max_parallel = 1

        def run(self, inputs):
            return {"cell_id": inputs.cell_id}

    def fake_build_default_engines():
        return {"find-circ3": _FakeDet("find-circ3"), "ciri2": _FakeDet("ciri2")}

    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", fake_build_default_engines)

    # --- monkeypatch demux: write manifest v1 and return fake stats ---
    def fake_demux(params):
        # create per-cell fastq files referenced by manifest
        fqdir = params.outdir / "fastq"
        fqdir.mkdir(parents=True, exist_ok=True)
        cell_r1 = fqdir / "sc01_R1.fastq.gz"
        cell_r2 = fqdir / "sc01_R2.fastq.gz"
        cell_r1.write_text("r1")
        cell_r2.write_text("r2")
        _write_manifest_v1(params.manifest_path, cell_r1, cell_r2)
        return {"sc01": {"assigned": 10}}

    monkeypatch.setattr("circyto.cli.smoke.demux_smartseq2_pooled", fake_demux)

    # --- monkeypatch run: just assert it received DetectorBase-like object and manifest exists ---
    def fake_run_detector_manifest(*, detector, manifest, outdir, ref_fa, gtf, threads, parallel):
        assert hasattr(detector, "name")
        assert Path(manifest).exists()
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        if detector.name == "find-circ3":
            det_dir = out / "sc01"
            det_dir.mkdir(parents=True, exist_ok=True)
            (det_dir / "sc01_splice_sites.bed").write_text("chr21\t1\t10\n")
        else:
            (out / "sc01.tsv").write_text(
                "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
                "circ1\tchr21\t1\t10\t+\t1\n"
            )

    monkeypatch.setattr("circyto.cli.smoke.run_detector_manifest", fake_run_detector_manifest)

    # --- monkeypatch collect: write minimal mtx + indexes ---
    def fake_collect(indir, matrix_path, circ_index_path, cell_index_path, *args, **kwargs):
        Path(matrix_path).write_text("%%MatrixMarket matrix coordinate integer general\n")
        Path(circ_index_path).write_text("circ1\n")
        Path(cell_index_path).write_text("sc01\n")

    monkeypatch.setattr("circyto.cli.smoke.collect_matrix", fake_collect)
    monkeypatch.setattr("circyto.cli.smoke.collect_find_circ3_matrix", fake_collect)

    # --- monkeypatch convert: write h5ad placeholder ---
    def fake_convert_matrix_files(*, matrix_path, circ_index_path, cell_index_path, loom=None, h5ad=None):
        assert Path(matrix_path).exists()
        assert Path(circ_index_path).exists()
        assert Path(cell_index_path).exists()
        assert h5ad is not None
        Path(h5ad).write_text("h5ad")

    monkeypatch.setattr("circyto.cli.smoke.convert_matrix_files", fake_convert_matrix_files)

    # Run CLI
    res = runner.invoke(
        app,
        [
            "smoke",
            "smartseq2",
            "--r1",
            str(r1),
            "--r2",
            str(r2),
            "--barcodes",
            str(barcodes),
            "--ref-fa",
            str(ref_fa),
            "--detector",
            detector,
            "--outdir",
            str(outdir),
            "--overwrite",
        ],
    )
    assert res.exit_code == 0, res.output

    # expected outputs exist
    assert (outdir / "demux" / "manifest.tsv").exists()
    assert (outdir / "matrix" / detector / "circ_counts.mtx").exists()
    assert (outdir / "h5ad" / f"{detector}.h5ad").exists()
