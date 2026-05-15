from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.align_manifest import read_source_manifest
from circyto.pipeline.run_ciri3 import _build_subset_rows


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


def test_read_source_manifest_normalizes_short_read_layout_aliases(tmp_path: Path) -> None:
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "sample_id\tfastq_1\tfastq_2\tprotocol\tstrandedness\tread_layout",
                f"cell1\t{fastq}\t\tramda\tunstranded\tsingle",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = read_source_manifest(manifest, validate_files=True)
    assert rows[0].read_layout == "single-end"


def test_build_subset_rows_absolutizes_manifest_relative_fastqs(tmp_path: Path) -> None:
    manifest = tmp_path / "ramda" / "manifest.tsv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "sample_id": "cell1",
            "fastq_1": "subset_100k/SRR5516584.100k.fastq.gz",
            "fastq_2": "",
            "protocol": "ramda",
            "strandedness": "unstranded",
            "read_layout": "single",
        }
    ]

    subset = _build_subset_rows(
        rows,
        source_manifest=manifest,
        keep_cell_ids={"cell1"},
        protocol_by_cell={"cell1": "ramda"},
        strandedness_by_cell={"cell1": "unstranded"},
    )

    assert subset[0]["fastq_1"] == str((manifest.parent / "subset_100k/SRR5516584.100k.fastq.gz").resolve())


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


def test_run_ciri3_dry_run_does_not_require_real_fastq(tmp_path: Path, monkeypatch) -> None:
    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "detector_alignment_plan.json").write_text(
            json.dumps(
                {
                    "command_preview": [
                        {
                            "cell_id": "ramda_se_1",
                            "command": "java -jar CIRI3.jar -I planned.sam -O ciri3_raw.tsv",
                        }
                    ]
                }
            )
            + "\n",
            encoding="utf-8",
        )
        return []

    monkeypatch.setattr(
        "circyto.pipeline.run_ciri3.run_detector_alignment_manifest",
        fake_run_detector_alignment_manifest,
    )
    monkeypatch.setattr("circyto.pipeline.align_manifest._tool_exists", lambda name: True)

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "sample_id\tfastq_1\tfastq_2\tprotocol\tstrandedness\tread_layout",
                "ramda_se_1\t/path/that/does/not/exist.fastq.gz\t\tramda\tunstranded\tsingle",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    ref = tmp_path / "chr21.fa"
    gtf = tmp_path / "chr21.gtf"
    ref.write_text(">chr21\nACGT\n", encoding="utf-8")
    gtf.write_text('chr21\ttest\texon\t1\t4\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "run-ciri3",
            "--manifest",
            str(manifest),
            "--outdir",
            str(tmp_path / "run"),
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--protocol",
            "ramda",
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert "/path/that/does/not/exist.fastq.gz" in result.stdout


def test_run_ciri3_single_end_real_execution_without_star_index(tmp_path: Path, monkeypatch) -> None:
    seen: dict[str, object] = {}

    def fake_prepare_alignment_cache(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / "alignment_manifest.tsv"
        manifest_path.write_text(
            "\n".join(
                [
                    "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness",
                    "shin_ramda_chr21_se\t\t/tmp/shin.sam\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tcache1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        seen["aligner"] = kwargs["aligner"]
        return manifest_path

    def fake_run_detector_alignment_manifest(**kwargs):
        seen["detector_manifest"] = kwargs["manifest"]
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "shin_ramda_chr21_se.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\n",
            encoding="utf-8",
        )
        return []

    def fake_collect_matrix(*args, **kwargs):
        matrix_dir = Path(args[1]).parent
        matrix_dir.mkdir(parents=True, exist_ok=True)
        for name in ("circ_counts.mtx", "circ_index.txt", "cell_index.txt"):
            (matrix_dir / name).write_text("", encoding="utf-8")

    monkeypatch.setattr("circyto.pipeline.run_ciri3.prepare_alignment_cache", fake_prepare_alignment_cache)
    monkeypatch.setattr(
        "circyto.pipeline.run_ciri3.run_detector_alignment_manifest",
        fake_run_detector_alignment_manifest,
    )
    monkeypatch.setattr("circyto.pipeline.run_ciri3.collect_matrix", fake_collect_matrix)

    ref = tmp_path / "chr21.fa"
    gtf = tmp_path / "chr21.gtf"
    ref.write_text(">chr21\nACGT\n", encoding="utf-8")
    gtf.write_text('chr21\ttest\texon\t1\t4\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")

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
            "--protocol",
            "shin-ramda",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert seen["aligner"] == "bwa-mem"
    assert Path(seen["detector_manifest"]).exists()
