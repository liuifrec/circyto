from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def test_smoke_fails_when_detector_outputs_are_all_empty(
    monkeypatch,
    tmp_path: Path,
) -> None:
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    barcodes = tmp_path / "barcodes.tsv"
    ref_fa = tmp_path / "ref.fa"
    gtf = tmp_path / "genes.gtf"
    outdir = tmp_path / "work"

    r1.write_text("x", encoding="utf-8")
    r2.write_text("y", encoding="utf-8")
    barcodes.write_text("sc01\tAACGTGAT\n", encoding="utf-8")
    ref_fa.write_text(">chr21\nACGT\n", encoding="utf-8")
    gtf.write_text('chr21\ttest\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")

    class _FakeDet:
        def __init__(self, name: str):
            self.name = name
            self.max_parallel = 1

        def version(self):
            return "fake"

        def run(self, inputs):
            raise AssertionError("run() should be monkeypatched at the pipeline layer")

    monkeypatch.setattr(
        "circyto.cli.smoke.build_default_engines",
        lambda: {"ciri2": _FakeDet("ciri2")},
    )

    def fake_demux(params):
        fqdir = params.outdir / "fastq"
        fqdir.mkdir(parents=True, exist_ok=True)
        cell_r1 = fqdir / "sc01_R1.fastq.gz"
        cell_r2 = fqdir / "sc01_R2.fastq.gz"
        cell_r1.write_text("r1", encoding="utf-8")
        cell_r2.write_text("r2", encoding="utf-8")
        params.manifest_path.write_text(
            "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\n"
            f"sc01\tsmartseq2\t{cell_r1}\t{cell_r2}\t\tLIB\t10\n",
            encoding="utf-8",
        )
        return {"sc01": {"assigned": 10}}

    monkeypatch.setattr("circyto.cli.smoke.demux_smartseq2_pooled", fake_demux)

    def fake_run_detector_manifest(*, detector, manifest, outdir, ref_fa, gtf, threads, parallel):
        Path(outdir).mkdir(parents=True, exist_ok=True)
        # Header-only normalized TSV, which would previously flow through to an all-zero matrix.
        (Path(outdir) / "sc01.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\n",
            encoding="utf-8",
        )

    monkeypatch.setattr("circyto.cli.smoke.run_detector_manifest", fake_run_detector_manifest)

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
            "--gtf",
            str(gtf),
            "--detector",
            "ciri2",
            "--outdir",
            str(outdir),
            "--overwrite",
        ],
    )

    assert res.exit_code != 0
    assert res.exception is not None
    assert "only empty outputs" in str(res.exception)
