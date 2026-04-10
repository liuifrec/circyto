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


def _make_repo_fixture(root: Path) -> tuple[Path, Path]:
    tools_dir = root / "tools"
    test_dir = tools_dir / "CIRI-full_v2.0" / "CIRI-full_test"
    test_dir.mkdir(parents=True, exist_ok=True)
    (test_dir / "test_1.fq.gz").write_text("r1", encoding="utf-8")
    (test_dir / "test_2.fq.gz").write_text("r2", encoding="utf-8")

    ref_dir = root / "ref"
    ref_dir.mkdir(parents=True, exist_ok=True)
    (ref_dir / "chr21.fa").write_text(">chr21\nACGT\n", encoding="utf-8")
    (ref_dir / "chr21.gtf").write_text(
        'chr21\ttest\texon\t1\t4\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
        encoding="utf-8",
    )
    (ref_dir / "star_index_chr21").mkdir(parents=True, exist_ok=True)
    return root, tools_dir


def test_smoke_root_help_mentions_alignment_first_options() -> None:
    result = runner.invoke(app, ["smoke", "--help"])
    assert result.exit_code == 0
    assert "--aligner" in result.stdout
    assert "--allow-empty" in result.stdout
    assert "--genome-dir" in result.stdout


def test_smoke_bwa_ciri3_writes_summary(tmp_path: Path, monkeypatch) -> None:
    repo_root, tools_dir = _make_repo_fixture(tmp_path / "repo")
    monkeypatch.setattr("circyto.cli.smoke.get_repo_root", lambda: repo_root)
    monkeypatch.setattr("circyto.cli.smoke.get_tools_dir", lambda: tools_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})

    def fake_prepare_alignment_cache(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / "alignment_manifest.tsv"
        manifest_path.write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\t"
            "chimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\tmapper_mode\tartifact_bucket\tsortedness\n"
            f"smoke1\t\t{outdir/'cache'/'alignment.sam'}\tsmoke\tpaired-end\t{kwargs['aligner']}\t{kwargs['ref_fa']}\t"
            "cache1\tsmoke_manifest.tsv\t\t\t\t\t0\tbwa_mem\tunsorted\n",
            encoding="utf-8",
        )
        (outdir / "alignment_prepare_summary.json").write_text(
            json.dumps(
                {
                    "aligner": kwargs["aligner"],
                    "status_counts": {"aligned": 1},
                    "cells": [{"cell_id": "smoke1", "status": "aligned"}],
                }
            ),
            encoding="utf-8",
        )
        return manifest_path

    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        tsv = outdir / "smoke1.tsv"
        tsv.write_text("circ_id\tchr\tstart\tend\tstrand\tsupport\n", encoding="utf-8")
        (outdir / "detector_run_summary.json").write_text(
            json.dumps(
                {
                    "cells": [
                        {
                            "cell_id": "smoke1",
                            "status": "empty",
                            "normalized_row_count": 0,
                            "input_file_type": "SAM",
                            "mapper_mode": "0",
                            "tsv_path": str(tsv),
                        }
                    ]
                }
            ),
            encoding="utf-8",
        )
        return []

    def fake_collect_matrix(cirifull_dir, matrix_path, circ_index_path, cell_index_path, min_count_per_cell=1):
        Path(matrix_path).write_text("%%MatrixMarket matrix coordinate real symmetric\n%\n0 0 0\n", encoding="utf-8")
        Path(circ_index_path).write_text("", encoding="utf-8")
        Path(cell_index_path).write_text("", encoding="utf-8")

    monkeypatch.setattr("circyto.cli.smoke.prepare_alignment_cache", fake_prepare_alignment_cache)
    monkeypatch.setattr("circyto.cli.smoke.run_detector_alignment_manifest", fake_run_detector_alignment_manifest)
    monkeypatch.setattr("circyto.cli.smoke.collect_matrix", fake_collect_matrix)

    outdir = tmp_path / "out"
    result = runner.invoke(app, ["smoke", "--detector", "ciri3", "--aligner", "bwa-mem", "--outdir", str(outdir)])
    assert result.exit_code == 0, result.stdout
    summary = json.loads((outdir / "smoke_summary.json").read_text(encoding="utf-8"))
    assert summary["smoke_pass"] is True
    assert summary["aligner"] == "bwa-mem"
    assert summary["allow_empty"] is True
    assert summary["matrix_nnz"] == 0


def test_smoke_star_ciri3_passes_tmpdir_to_prepare(tmp_path: Path, monkeypatch) -> None:
    repo_root, tools_dir = _make_repo_fixture(tmp_path / "repo")
    monkeypatch.setattr("circyto.cli.smoke.get_repo_root", lambda: repo_root)
    monkeypatch.setattr("circyto.cli.smoke.get_tools_dir", lambda: tools_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})
    seen: dict[str, str] = {}

    def fake_prepare_alignment_cache(**kwargs):
        seen["aligner"] = kwargs["aligner"]
        seen["extra_flags"] = kwargs["extra_flags"]
        seen["star_tmpdir"] = os.environ.get("CIRCYTO_STAR_TMPDIR", "")
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / "alignment_manifest.tsv"
        manifest_path.write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\t"
            "chimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\tmapper_mode\tartifact_bucket\tsortedness\n"
            f"smoke1\t\t{outdir/'cache'/'alignment.sam'}\tsmoke\tpaired-end\tstar\t{kwargs['ref_fa']}\t"
            f"cache1\tsmoke_manifest.tsv\t{outdir/'cache'/'chim.txt'}\t{outdir/'cache'/'mate1.fq'}\t"
            f"{outdir/'cache'/'mate2.fq'}\t{outdir/'cache'/'bwa.sam'}\t1\tstar\tunsorted\n",
            encoding="utf-8",
        )
        (outdir / "alignment_prepare_summary.json").write_text(
            json.dumps({"aligner": "star", "status_counts": {"aligned": 1}, "cells": [{"cell_id": "smoke1", "status": "aligned"}]}),
            encoding="utf-8",
        )
        return manifest_path

    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        tsv = outdir / "smoke1.tsv"
        tsv.write_text("circ_id\tchr\tstart\tend\tstrand\tsupport\n", encoding="utf-8")
        (outdir / "detector_run_summary.json").write_text(
            json.dumps(
                {
                    "cells": [
                        {
                            "cell_id": "smoke1",
                            "status": "empty",
                            "normalized_row_count": 0,
                            "input_file_type": "STAR tuple",
                            "mapper_mode": "1",
                            "tsv_path": str(tsv),
                        }
                    ]
                }
            ),
            encoding="utf-8",
        )
        return []

    def fake_collect_matrix(cirifull_dir, matrix_path, circ_index_path, cell_index_path, min_count_per_cell=1):
        Path(matrix_path).write_text("%%MatrixMarket matrix coordinate real symmetric\n%\n0 0 0\n", encoding="utf-8")
        Path(circ_index_path).write_text("", encoding="utf-8")
        Path(cell_index_path).write_text("", encoding="utf-8")

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
    assert seen["aligner"] == "star"
    assert "--genomeDir" in seen["extra_flags"]
    assert seen["star_tmpdir"] == str(star_tmp)


def test_smoke_require_nonempty_fails_on_empty_matrix(tmp_path: Path, monkeypatch) -> None:
    repo_root, tools_dir = _make_repo_fixture(tmp_path / "repo")
    monkeypatch.setattr("circyto.cli.smoke.get_repo_root", lambda: repo_root)
    monkeypatch.setattr("circyto.cli.smoke.get_tools_dir", lambda: tools_dir)
    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _FakeCiri3()})

    def fake_prepare_alignment_cache(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        manifest_path = outdir / "alignment_manifest.tsv"
        manifest_path.write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\t"
            "chimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\tmapper_mode\tartifact_bucket\tsortedness\n"
            f"smoke1\t\t{outdir/'cache'/'alignment.sam'}\tsmoke\tpaired-end\tbwa-mem\t{kwargs['ref_fa']}\t"
            "cache1\tsmoke_manifest.tsv\t\t\t\t\t0\tbwa_mem\tunsorted\n",
            encoding="utf-8",
        )
        (outdir / "alignment_prepare_summary.json").write_text(
            json.dumps({"aligner": "bwa-mem", "status_counts": {"aligned": 1}, "cells": [{"cell_id": "smoke1", "status": "aligned"}]}),
            encoding="utf-8",
        )
        return manifest_path

    def fake_run_detector_alignment_manifest(**kwargs):
        outdir = Path(kwargs["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        tsv = outdir / "smoke1.tsv"
        tsv.write_text("circ_id\tchr\tstart\tend\tstrand\tsupport\n", encoding="utf-8")
        (outdir / "detector_run_summary.json").write_text(
            json.dumps(
                {
                    "cells": [
                        {
                            "cell_id": "smoke1",
                            "status": "empty",
                            "normalized_row_count": 0,
                            "input_file_type": "SAM",
                            "mapper_mode": "0",
                            "tsv_path": str(tsv),
                        }
                    ]
                }
            ),
            encoding="utf-8",
        )
        return []

    def fake_collect_matrix(cirifull_dir, matrix_path, circ_index_path, cell_index_path, min_count_per_cell=1):
        Path(matrix_path).write_text("%%MatrixMarket matrix coordinate real symmetric\n%\n0 0 0\n", encoding="utf-8")
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
    assert summary["allow_empty"] is False


def test_smoke_reports_missing_ciri3_runtime(tmp_path: Path, monkeypatch) -> None:
    repo_root, tools_dir = _make_repo_fixture(tmp_path / "repo")
    monkeypatch.setattr("circyto.cli.smoke.get_repo_root", lambda: repo_root)
    monkeypatch.setattr("circyto.cli.smoke.get_tools_dir", lambda: tools_dir)

    class _BrokenCiri3(_FakeCiri3):
        def validate_runtime(self):
            return False, ["No CIRI3 jar or wrapper detected."], {}

    monkeypatch.setattr("circyto.cli.smoke.build_default_engines", lambda: {"ciri3": _BrokenCiri3()})

    result = runner.invoke(app, ["smoke", "--detector", "ciri3", "--aligner", "bwa-mem", "--outdir", str(tmp_path / "out")])
    assert result.exit_code != 0
    assert "CIRI3 runtime not ready" in (result.stdout + (result.stderr or ""))
