from __future__ import annotations

import json
import os
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.detectors.base import DetectorCapabilities


runner = CliRunner()


class _FakeCiri3:
    name = "ciri3"
    input_type = "alignment"
    max_parallel = 1
    capabilities = DetectorCapabilities(
        accepts_fastq=False,
        accepts_alignment=True,
        prefers_paired=True,
        supports_single_end=True,
        supports_multisample_alignment=True,
        requires_unsorted_sam=True,
        supports_star=True,
        supports_bwa=True,
        max_parallel=1,
        recommended_execution_mode="alignment-first",
    )

    def validate_runtime(self):
        return True, [], {"preferred_mode": "direct-jar"}


def _make_packaged_demo(root: Path) -> Path:
    demo = root / "smoke_demo"
    demo.mkdir(parents=True, exist_ok=True)
    (demo / "reference.fa").write_text(">smokechr1\nACGTACGTACGT\n", encoding="utf-8")
    (demo / "reference.gtf").write_text(
        'smokechr1\ttest\texon\t1\t12\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        encoding="utf-8",
    )
    (demo / "single_end.fastq").write_text(
        "@se1\nACGTACGT\n+\nIIIIIIII\n",
        encoding="utf-8",
    )
    (demo / "paired_end_R1.fastq").write_text(
        "@pe1/1\nACGTACGT\n+\nIIIIIIII\n",
        encoding="utf-8",
    )
    (demo / "paired_end_R2.fastq").write_text(
        "@pe1/2\nTGCATGCA\n+\nIIIIIIII\n",
        encoding="utf-8",
    )
    return demo


def test_smoke_root_help_mentions_packaged_demo_options() -> None:
    result = runner.invoke(app, ["smoke", "--help"])
    assert result.exit_code == 0
    assert "--aligner" in result.stdout
    assert "--allow-empty" in result.stdout
    assert "--read-layout" in result.stdout


def test_smoke_bwa_ciri3_writes_summary_from_packaged_demo(tmp_path: Path, monkeypatch) -> None:
    demo_dir = _make_packaged_demo(tmp_path / "pkg")
    monkeypatch.setattr("circyto.cli.smoke.get_packaged_smoke_demo_dir", lambda: demo_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})
    monkeypatch.setattr("circyto.cli.smoke._ensure_bwa_index", lambda ref_fa: [str(ref_fa.with_suffix(".fa.bwt"))])

    def fake_prepare_alignment_cache(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / "alignment_manifest.tsv"
        manifest_path.write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\t"
            "chimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\tmapper_mode\tartifact_bucket\tsortedness\n"
            f"smoke_pe\t\t{outdir/'cache'/'alignment.sam'}\tsmoke\tpaired-end\t{kwargs['aligner']}\t{kwargs['ref_fa']}\t"
            "cache1\tsmoke_manifest.tsv\t\t\t\t\t0\tbwa_mem\tunsorted\n",
            encoding="utf-8",
        )
        (outdir / "alignment_prepare_summary.json").write_text(
            json.dumps({"aligner": kwargs["aligner"], "cells": [{"cell_id": "smoke_pe", "status": "aligned"}]}),
            encoding="utf-8",
        )
        return manifest_path

    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        tsv = outdir / "smoke_pe.tsv"
        tsv.write_text("circ_id\tchr\tstart\tend\tstrand\tsupport\n", encoding="utf-8")
        (outdir / "detector_run_summary.json").write_text(
            json.dumps(
                {"cells": [{"cell_id": "smoke_pe", "status": "empty", "normalized_row_count": 0, "input_file_type": "SAM", "mapper_mode": "0", "tsv_path": str(tsv)}]}
            ),
            encoding="utf-8",
        )
        return []

    def fake_collect_matrix(cirifull_dir, matrix_path, circ_index_path, cell_index_path, min_count_per_cell=1):
        Path(matrix_path).parent.mkdir(parents=True, exist_ok=True)
        Path(matrix_path).write_text("%%MatrixMarket matrix coordinate integer general\n%\n0 0 0\n", encoding="utf-8")
        Path(circ_index_path).write_text("", encoding="utf-8")
        Path(cell_index_path).write_text("", encoding="utf-8")

    monkeypatch.setattr("circyto.cli.smoke.prepare_alignment_cache", fake_prepare_alignment_cache)
    monkeypatch.setattr("circyto.cli.smoke.run_detector_alignment_manifest", fake_run_detector_alignment_manifest)
    monkeypatch.setattr("circyto.cli.smoke.collect_matrix", fake_collect_matrix)

    outdir = tmp_path / "out"
    result = runner.invoke(app, ["smoke", "--detector", "ciri3", "--aligner", "bwa-mem", "--outdir", str(outdir)])
    assert result.exit_code == 0, result.stdout
    summary = json.loads((outdir / "smoke_summary.json").read_text(encoding="utf-8"))
    assert summary["fixture"]["fixture_kind"] == "packaged-demo"
    assert summary["fixture"]["read_layout"] == "paired-end"
    assert summary["smoke_pass"] is True
    assert summary["generated_index_paths"]


def test_smoke_star_ciri3_uses_generated_genome_dir(tmp_path: Path, monkeypatch) -> None:
    demo_dir = _make_packaged_demo(tmp_path / "pkg")
    monkeypatch.setattr("circyto.cli.smoke.get_packaged_smoke_demo_dir", lambda: demo_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})
    monkeypatch.setattr("circyto.cli.smoke._ensure_bwa_index", lambda ref_fa: [])
    seen: dict[str, str] = {}

    def fake_ensure_star_genome_dir(ref_fa: Path, gtf: Path | None, genome_dir: Path) -> Path:
        seen["genome_dir"] = str(genome_dir)
        genome_dir.mkdir(parents=True, exist_ok=True)
        (genome_dir / "genomeParameters.txt").write_text("ok\n", encoding="utf-8")
        return genome_dir

    def fake_prepare_alignment_cache(**kwargs):
        seen["extra_flags"] = kwargs["extra_flags"]
        seen["star_tmpdir"] = os.environ.get("CIRCYTO_STAR_TMPDIR", "")
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / "alignment_manifest.tsv"
        manifest_path.write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\t"
            "chimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\tmapper_mode\tartifact_bucket\tsortedness\n"
            f"smoke_pe\t\t{outdir/'cache'/'alignment.sam'}\tsmoke\tpaired-end\tstar\t{kwargs['ref_fa']}\t"
            f"cache1\tsmoke_manifest.tsv\t{outdir/'cache'/'chim.txt'}\t{outdir/'cache'/'mate1.fq'}\t"
            f"{outdir/'cache'/'mate2.fq'}\t{outdir/'cache'/'bwa.sam'}\t1\tstar\tunsorted\n",
            encoding="utf-8",
        )
        (outdir / "alignment_prepare_summary.json").write_text(
            json.dumps({"aligner": "star", "cells": [{"cell_id": "smoke_pe", "status": "aligned"}]}),
            encoding="utf-8",
        )
        return manifest_path

    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        tsv = outdir / "smoke_pe.tsv"
        tsv.write_text("circ_id\tchr\tstart\tend\tstrand\tsupport\n", encoding="utf-8")
        (outdir / "detector_run_summary.json").write_text(
            json.dumps(
                {"cells": [{"cell_id": "smoke_pe", "status": "empty", "normalized_row_count": 0, "input_file_type": "STAR tuple", "mapper_mode": "1", "tsv_path": str(tsv)}]}
            ),
            encoding="utf-8",
        )
        return []

    def fake_collect_matrix(cirifull_dir, matrix_path, circ_index_path, cell_index_path, min_count_per_cell=1):
        Path(matrix_path).parent.mkdir(parents=True, exist_ok=True)
        Path(matrix_path).write_text("%%MatrixMarket matrix coordinate integer general\n%\n0 0 0\n", encoding="utf-8")
        Path(circ_index_path).write_text("", encoding="utf-8")
        Path(cell_index_path).write_text("", encoding="utf-8")

    monkeypatch.setattr("circyto.cli.smoke._ensure_star_genome_dir", fake_ensure_star_genome_dir)
    monkeypatch.setattr("circyto.cli.smoke.prepare_alignment_cache", fake_prepare_alignment_cache)
    monkeypatch.setattr("circyto.cli.smoke.run_detector_alignment_manifest", fake_run_detector_alignment_manifest)
    monkeypatch.setattr("circyto.cli.smoke.collect_matrix", fake_collect_matrix)

    star_tmp = tmp_path / "star_tmp"
    result = runner.invoke(
        app,
        [
            "smoke",
            "--detector",
            "ciri3",
            "--aligner",
            "star",
            "--outdir",
            str(tmp_path / "out"),
            "--tmpdir",
            str(star_tmp),
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert "--genomeDir" in seen["extra_flags"]
    assert seen["star_tmpdir"] == str(star_tmp)
    assert seen["genome_dir"].endswith("star_index")


def test_smoke_require_nonempty_fails_on_empty_matrix(tmp_path: Path, monkeypatch) -> None:
    demo_dir = _make_packaged_demo(tmp_path / "pkg")
    monkeypatch.setattr("circyto.cli.smoke.get_packaged_smoke_demo_dir", lambda: demo_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})
    monkeypatch.setattr("circyto.cli.smoke._ensure_bwa_index", lambda ref_fa: [])

    def fake_prepare_alignment_cache(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / "alignment_manifest.tsv"
        manifest_path.write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\t"
            "chimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\tmapper_mode\tartifact_bucket\tsortedness\n"
            f"smoke_pe\t\t{outdir/'cache'/'alignment.sam'}\tsmoke\tpaired-end\tbwa-mem\t{kwargs['ref_fa']}\t"
            "cache1\tsmoke_manifest.tsv\t\t\t\t\t0\tbwa_mem\tunsorted\n",
            encoding="utf-8",
        )
        (outdir / "alignment_prepare_summary.json").write_text(
            json.dumps({"aligner": "bwa-mem", "cells": [{"cell_id": "smoke_pe", "status": "aligned"}]}),
            encoding="utf-8",
        )
        return manifest_path

    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        tsv = outdir / "smoke_pe.tsv"
        tsv.write_text("circ_id\tchr\tstart\tend\tstrand\tsupport\n", encoding="utf-8")
        (outdir / "detector_run_summary.json").write_text(
            json.dumps(
                {"cells": [{"cell_id": "smoke_pe", "status": "empty", "normalized_row_count": 0, "input_file_type": "SAM", "mapper_mode": "0", "tsv_path": str(tsv)}]}
            ),
            encoding="utf-8",
        )
        return []

    def fake_collect_matrix(cirifull_dir, matrix_path, circ_index_path, cell_index_path, min_count_per_cell=1):
        Path(matrix_path).parent.mkdir(parents=True, exist_ok=True)
        Path(matrix_path).write_text("%%MatrixMarket matrix coordinate integer general\n%\n0 0 0\n", encoding="utf-8")
        Path(circ_index_path).write_text("", encoding="utf-8")
        Path(cell_index_path).write_text("", encoding="utf-8")

    monkeypatch.setattr("circyto.cli.smoke.prepare_alignment_cache", fake_prepare_alignment_cache)
    monkeypatch.setattr("circyto.cli.smoke.run_detector_alignment_manifest", fake_run_detector_alignment_manifest)
    monkeypatch.setattr("circyto.cli.smoke.collect_matrix", fake_collect_matrix)

    result = runner.invoke(
        app,
        ["smoke", "--detector", "ciri3", "--aligner", "bwa-mem", "--outdir", str(tmp_path / "out"), "--require-nonempty"],
    )
    assert result.exit_code == 1
    summary = json.loads((tmp_path / "out" / "smoke_summary.json").read_text(encoding="utf-8"))
    assert summary["smoke_pass"] is False


def test_smoke_reports_missing_ciri3_runtime(tmp_path: Path, monkeypatch) -> None:
    demo_dir = _make_packaged_demo(tmp_path / "pkg")
    monkeypatch.setattr("circyto.cli.smoke.get_packaged_smoke_demo_dir", lambda: demo_dir)

    class _BrokenCiri3(_FakeCiri3):
        def validate_runtime(self):
            return False, ["No CIRI3 jar or wrapper detected."], {}

    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _BrokenCiri3()})
    result = runner.invoke(app, ["smoke", "--detector", "ciri3", "--aligner", "bwa-mem", "--outdir", str(tmp_path / "out")])
    assert result.exit_code != 0
    assert "CIRI3 runtime not ready" in (result.stdout + (result.stderr or ""))


def test_smoke_fails_early_for_mixed_layout_manifest(tmp_path: Path, monkeypatch) -> None:
    demo_dir = _make_packaged_demo(tmp_path / "pkg")
    monkeypatch.setattr("circyto.cli.smoke.get_packaged_smoke_demo_dir", lambda: demo_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})

    se_r1 = tmp_path / "se.fastq"
    pe_r1 = tmp_path / "pe_R1.fastq"
    pe_r2 = tmp_path / "pe_R2.fastq"
    se_r1.write_text("@se\nACGT\n+\nIIII\n", encoding="utf-8")
    pe_r1.write_text("@pe1\nACGT\n+\nIIII\n", encoding="utf-8")
    pe_r2.write_text("@pe2\nTGCA\n+\nIIII\n", encoding="utf-8")
    manifest = tmp_path / "mixed.tsv"
    manifest.write_text(
        "\n".join(
            [
                "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id\tread_layout",
                f"se1\tsmartseq2\t{se_r1}\t\t\tlib1\t10\tgrp\tsingle-end",
                f"pe1\tsmartseq2\t{pe_r1}\t{pe_r2}\t\tlib1\t10\tgrp\tpaired-end",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "smoke",
            "--detector",
            "ciri3",
            "--aligner",
            "bwa-mem",
            "--outdir",
            str(tmp_path / "out"),
            "--use-local-manifest",
            "--manifest",
            str(manifest),
            "--subset-cells",
            "2",
        ],
    )
    assert result.exit_code != 0
    assert "mixes read layouts" in (result.stdout + (result.stderr or ""))


def test_smoke_fails_early_for_invalid_reference_path(tmp_path: Path, monkeypatch) -> None:
    demo_dir = _make_packaged_demo(tmp_path / "pkg")
    monkeypatch.setattr("circyto.cli.smoke.get_packaged_smoke_demo_dir", lambda: demo_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})
    missing_ref = tmp_path / "missing.fa"
    result = runner.invoke(
        app,
        [
            "smoke",
            "--detector",
            "ciri3",
            "--aligner",
            "bwa-mem",
            "--outdir",
            str(tmp_path / "out"),
            "--ref-fa",
            str(missing_ref),
        ],
    )
    assert result.exit_code != 0
    assert "Reference FASTA not found:" in (result.stdout + (result.stderr or ""))


def test_smoke_rejects_single_end_star_demo(tmp_path: Path, monkeypatch) -> None:
    demo_dir = _make_packaged_demo(tmp_path / "pkg")
    monkeypatch.setattr("circyto.cli.smoke.get_packaged_smoke_demo_dir", lambda: demo_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})
    result = runner.invoke(
        app,
        [
            "smoke",
            "--detector",
            "ciri3",
            "--aligner",
            "star",
            "--read-layout",
            "single-end",
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
    assert "STAR smoke currently requires --read-layout paired-end" in (result.stdout + (result.stderr or ""))
