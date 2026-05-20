from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_summary(workdir: Path) -> None:
    paths = {
        "manifest": str((workdir / "manifests" / "single_end_manifest.tsv").resolve()),
        "alignment_manifest": str((workdir / "align" / "alignment_manifest.tsv").resolve()),
        "detector_summary": str((workdir / "ciri3" / "detector_run_summary.json").resolve()),
        "matrix": str((workdir / "matrix" / "circ_counts.mtx").resolve()),
        "cell_index": str((workdir / "matrix" / "cell_index.txt").resolve()),
        "circ_index": str((workdir / "matrix" / "circ_index.txt").resolve()),
        "h5ad": str((workdir / "anndata" / "circ_counts.h5ad").resolve()),
    }
    summary = {
        "workflow": "full-length-circrna",
        "workflow_type": "full-length-circrna",
        "command_name": "circyto workflow full-length-circrna",
        "circyto_version": "0.9.0",
        "workflow_uuid": "123e4567-e89b-12d3-a456-426614174000",
        "protocol": "ramda",
        "read_layout": "single-end",
        "genome_fasta": "/tmp/ref.fa",
        "gtf": "/tmp/genes.gtf",
        "detector_backend": "ciri3",
        "started_at": "2026-05-19T00:00:00+00:00",
        "completed_at": "2026-05-19T00:01:00+00:00",
        "elapsed_seconds": 60.0,
        "hostname": "host",
        "python_version": "3.12.0",
        "stage_graph": [
            {"stage": "alignment", "status": "completed"},
            {"stage": "detector", "status": "completed"},
            {"stage": "matrix", "status": "completed"},
            {"stage": "h5ad_export", "status": "completed"},
        ],
        "completed_stages": ["alignment", "detector", "matrix", "h5ad_export"],
        "skipped_stages": [],
        "failed_stages": [],
        "partial_outputs_detected": [],
        "cleanup_summary": {"enabled": False, "scope": None, "performed": False, "deleted_paths": [], "reclaimed_bytes": 0},
        "paths": paths,
    }
    (workdir / "workflow_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def _write_fixture_workdir(tmp_path: Path) -> Path:
    workdir = tmp_path / "workflow"
    for rel in ("manifests", "align", "ciri3", "matrix", "anndata", "qc", "logs"):
        (workdir / rel).mkdir(parents=True, exist_ok=True)
    (workdir / "manifests" / "single_end_manifest.tsv").write_text("cell_id\tread1\ncell1\t/tmp/a.fastq\n", encoding="utf-8")
    (workdir / "align" / "alignment_manifest.tsv").write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
        f"cell1\t\t{workdir / 'align' / 'cell1.sam'}\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted\n",
        encoding="utf-8",
    )
    (workdir / "align" / "cell1.sam").write_text("@SQ\tSN:chr21\tLN:1000\n", encoding="utf-8")
    (workdir / "align" / "cache.sam").write_text("sam\n", encoding="utf-8")
    (workdir / "demux").mkdir(parents=True, exist_ok=True)
    (workdir / "demux" / "cell1_R1.fastq.gz").write_text("demux\n", encoding="utf-8")
    (workdir / "ciri3" / "detector_run_summary.json").write_text(json.dumps({"status_counts": {"success": 1}}, indent=2) + "\n", encoding="utf-8")
    (workdir / "matrix" / "circ_counts.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n%\n1 1 0\n", encoding="utf-8")
    (workdir / "matrix" / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
    (workdir / "matrix" / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
    (workdir / "anndata" / "circ_counts.h5ad").write_text("placeholder", encoding="utf-8")
    (workdir / "ciri3" / "cell1.tsv").write_text("circ_id\tcount\n", encoding="utf-8")
    (workdir / "ciri3" / "cell1.tsv.provenance.json").write_text(json.dumps({"cell_id": "cell1"}, indent=2) + "\n", encoding="utf-8")
    _write_summary(workdir)
    return workdir


def test_cleanup_workflow_dry_run_reports_plan_without_deleting(tmp_path: Path) -> None:
    workdir = _write_fixture_workdir(tmp_path)
    cache_sam = workdir / "align" / "cache.sam"
    result = runner.invoke(app, ["cleanup-workflow", "--workdir", str(workdir), "--scope", "alignments", "--dry-run"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["dry_run"] is True
    assert payload["cleanup"]["mode"] == "planned-only"
    assert payload["estimated_reclaimed_bytes"] >= 4
    assert cache_sam.exists()
    assert any(item["path"].endswith("cache.sam") for item in payload["planned_deletions"])


def test_cleanup_workflow_real_run_deletes_alignment_only_and_updates_summary(tmp_path: Path) -> None:
    workdir = _write_fixture_workdir(tmp_path)
    cache_sam = workdir / "align" / "cache.sam"
    demux_fastq = workdir / "demux" / "cell1_R1.fastq.gz"
    matrix = workdir / "matrix" / "circ_counts.mtx"
    result = runner.invoke(app, ["cleanup-workflow", "--workdir", str(workdir), "--scope", "alignments"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["cleanup"]["cleanup_performed"] is True
    assert not cache_sam.exists()
    assert demux_fastq.exists()
    assert matrix.exists()
    summary = json.loads((workdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["cleanup_summary"]["performed"] is True
    assert summary["cleanup_summary"]["scope"] == "alignments"
    assert any(path.endswith("cache.sam") for path in summary["cleanup_summary"]["deleted_paths"])


def test_cleanup_workflow_refuses_when_integrity_fails_without_force(tmp_path: Path) -> None:
    workdir = _write_fixture_workdir(tmp_path)
    (workdir / "matrix" / "circ_counts.mtx").unlink()
    result = runner.invoke(app, ["cleanup-workflow", "--workdir", str(workdir), "--scope", "alignments"])
    assert result.exit_code != 0
    assert "Workflow integrity check failed" in result.output
    assert "--force" in result.output


def test_cleanup_workflow_force_allows_cleanup_despite_integrity_failure(tmp_path: Path) -> None:
    workdir = _write_fixture_workdir(tmp_path)
    cache_sam = workdir / "align" / "cache.sam"
    (workdir / "matrix" / "circ_counts.mtx").unlink()
    result = runner.invoke(app, ["cleanup-workflow", "--workdir", str(workdir), "--scope", "alignments", "--force"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["workflow_integrity_ok"] is False
    assert payload["cleanup"]["cleanup_performed"] is True
    assert not cache_sam.exists()


def test_cleanup_workflow_help_lists_expected_options() -> None:
    result = runner.invoke(app, ["cleanup-workflow", "--help"])
    assert result.exit_code == 0
    assert "--workdir" in result.output
    assert "--scope" in result.output
    assert "--dry-run" in result.output
    assert "--force" in result.output
