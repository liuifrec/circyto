from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_summary(workdir: Path, *, include_rna: bool = False) -> None:
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
    if include_rna:
        summary["rna_import"] = {
            "method": "simple-overlap",
            "n_cells": 1,
            "n_genes": 2,
        }
    (workdir / "workflow_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def _write_fixture_workdir(tmp_path: Path, *, include_rna: bool = False) -> Path:
    workdir = tmp_path / "workflow"
    for rel in (
        "manifests",
        "align",
        "ciri3",
        "matrix",
        "anndata",
        "qc",
        "logs",
    ):
        (workdir / rel).mkdir(parents=True, exist_ok=True)
    (workdir / "manifests" / "single_end_manifest.tsv").write_text("cell_id\tread1\ncell1\t/tmp/a.fastq\n", encoding="utf-8")
    (workdir / "align" / "alignment_manifest.tsv").write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
        f"cell1\t\t{workdir / 'align' / 'cell1.sam'}\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted\n",
        encoding="utf-8",
    )
    (workdir / "align" / "cell1.sam").write_text("@SQ\tSN:chr21\tLN:1000\n", encoding="utf-8")
    (workdir / "ciri3" / "detector_run_summary.json").write_text(json.dumps({"status_counts": {"success": 1}}, indent=2) + "\n", encoding="utf-8")
    (workdir / "matrix" / "circ_counts.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n%\n1 1 0\n", encoding="utf-8")
    (workdir / "matrix" / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
    (workdir / "matrix" / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
    (workdir / "anndata" / "circ_counts.h5ad").write_text("placeholder", encoding="utf-8")
    (workdir / "ciri3" / "cell1.tsv.provenance.json").write_text(json.dumps({"cell_id": "cell1"}, indent=2) + "\n", encoding="utf-8")
    if include_rna:
        (workdir / "rna").mkdir(parents=True, exist_ok=True)
        (workdir / "rna" / "gene_counts.tsv").write_text("gene_id\tgene_name\tcell1\nG1\tGENE1\t1\n", encoding="utf-8")
        (workdir / "rna" / "gene_feature_table.tsv").write_text("gene_id\tgene_name\nG1\tGENE1\n", encoding="utf-8")
        (workdir / "rna" / "rna_import_summary.json").write_text(json.dumps({"method": "simple-overlap"}, indent=2) + "\n", encoding="utf-8")
    _write_summary(workdir, include_rna=include_rna)
    return workdir


def test_check_workflow_passes_on_valid_fixture(tmp_path: Path) -> None:
    workdir = _write_fixture_workdir(tmp_path, include_rna=True)
    result = runner.invoke(app, ["check-workflow", "--workdir", str(workdir)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["ok"] is True
    assert payload["workflow_type"] == "full-length-circrna"


def test_check_workflow_fails_on_missing_required_output(tmp_path: Path) -> None:
    workdir = _write_fixture_workdir(tmp_path)
    (workdir / "matrix" / "circ_counts.mtx").unlink()
    result = runner.invoke(app, ["check-workflow", "--workdir", str(workdir)])
    assert result.exit_code != 0
    payload = json.loads(result.output)
    assert payload["ok"] is False
    assert any("Missing required output: matrix" in item for item in payload["errors"])


def test_check_workflow_fails_on_corrupted_summary(tmp_path: Path) -> None:
    workdir = tmp_path / "workflow"
    workdir.mkdir()
    (workdir / "workflow_summary.json").write_text("{bad json\n", encoding="utf-8")
    result = runner.invoke(app, ["check-workflow", "--workdir", str(workdir)])
    assert result.exit_code != 0
    payload = json.loads(result.output)
    assert payload["ok"] is False
    assert any("Corrupted JSON file" in item for item in payload["errors"])
