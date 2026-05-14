from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.align_manifest import read_source_manifest


runner = CliRunner()


def test_read_source_manifest_supports_protocol_alias_columns() -> None:
    manifest = Path("testdata/shin_ramda/manifest_chr21_smoke.tsv")
    rows = read_source_manifest(manifest, validate_files=True)
    assert len(rows) == 1
    row = rows[0]
    assert row.cell_id == "shin_ramda_chr21_se"
    assert row.protocol == "shin-ramda"
    assert row.strandedness == "forward"
    assert row.read_layout == "single-end"
    assert row.read2 is None


def test_run_ciri3_dry_run_prints_alignment_and_detector_commands(tmp_path: Path, monkeypatch) -> None:
    seen: dict[str, object] = {}

    def fake_plan_alignment_cache(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        return {
            "command_preview": [
                {
                    "cell_id": "ramda_chr21_pe",
                    "command": "STAR --genomeDir /ref/star --readFilesIn r1 r2",
                }
            ]
        }

    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "detector_alignment_plan.json").write_text(
            json.dumps(
                {
                    "command_preview": [
                        {
                            "cell_id": "ramda_chr21_pe",
                            "command": "java -jar CIRI3.jar -I planned.sam -O ciri3_raw.tsv",
                        }
                    ]
                }
            )
            + "\n",
            encoding="utf-8",
        )
        seen["manifest"] = kwargs["manifest"]
        seen["dry_run"] = kwargs["dry_run"]
        return []

    monkeypatch.setattr("circyto.pipeline.run_ciri3.plan_alignment_cache", fake_plan_alignment_cache)
    monkeypatch.setattr(
        "circyto.pipeline.run_ciri3.run_detector_alignment_manifest",
        fake_run_detector_alignment_manifest,
    )

    ref = tmp_path / "chr21.fa"
    gtf = tmp_path / "chr21.gtf"
    star_index = tmp_path / "star_index"
    ref.write_text(">chr21\nACGT\n", encoding="utf-8")
    gtf.write_text('chr21\ttest\texon\t1\t4\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")
    star_index.mkdir()

    result = runner.invoke(
        app,
        [
            "run-ciri3",
            "--manifest",
            "testdata/ramda/manifest_chr21_smoke.tsv",
            "--outdir",
            str(tmp_path / "run"),
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--star-index",
            str(star_index),
            "--protocol",
            "ramda",
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert "STAR --genomeDir /ref/star --readFilesIn r1 r2" in result.stdout
    assert "java -jar CIRI3.jar -I planned.sam -O ciri3_raw.tsv" in result.stdout
    assert seen["dry_run"] is True
    assert Path(seen["manifest"]).exists()


def test_run_ciri3_dry_run_accepts_single_end_shin_ramda(tmp_path: Path, monkeypatch) -> None:
    def fake_plan_alignment_cache(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        return {"command_preview": [{"cell_id": "shin_ramda_chr21_se", "command": "bwa mem ref.fa r1.fastq"}]}

    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "detector_alignment_plan.json").write_text(
            json.dumps(
                {
                    "command_preview": [
                        {"cell_id": "shin_ramda_chr21_se", "command": "java -jar CIRI3.jar -I planned.sam -O ciri3_raw.tsv"}
                    ]
                }
            )
            + "\n",
            encoding="utf-8",
        )
        return []

    monkeypatch.setattr("circyto.pipeline.run_ciri3.plan_alignment_cache", fake_plan_alignment_cache)
    monkeypatch.setattr(
        "circyto.pipeline.run_ciri3.run_detector_alignment_manifest",
        fake_run_detector_alignment_manifest,
    )

    ref = tmp_path / "chr21.fa"
    gtf = tmp_path / "chr21.gtf"
    star_index = tmp_path / "star_index"
    ref.write_text(">chr21\nACGT\n", encoding="utf-8")
    gtf.write_text('chr21\ttest\texon\t1\t4\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")
    star_index.mkdir()

    result = runner.invoke(
        app,
        [
            "run-ciri3",
            "--manifest",
            "testdata/shin_ramda/manifest_chr21_smoke.tsv",
            "--outdir",
            str(tmp_path / "run"),
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--star-index",
            str(star_index),
            "--protocol",
            "shin-ramda",
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert "bwa mem ref.fa r1.fastq" in result.stdout
