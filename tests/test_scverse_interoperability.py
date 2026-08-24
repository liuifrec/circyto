from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.multimodal.sync import read_h5mu
from circyto.pipeline.gene_expression_velocity import HAS_MUDATA


runner = CliRunner()


def _write_validate_fixture(tmp_path: Path, *, include_h5mu: bool = False) -> Path:
    workdir = tmp_path / "workflow"
    (workdir / "matrix").mkdir(parents=True)
    (workdir / "anndata").mkdir()
    (workdir / "qc").mkdir()
    (workdir / "rna").mkdir()
    (workdir / "workflow_summary.json").write_text(
        json.dumps(
            {
                "workflow": "full-length-circrna",
                "workflow_type": "full-length-circrna",
                "command_name": "circyto workflow full-length-circrna",
                "circyto_version": "0.9.0",
                "workflow_uuid": "uuid-1",
                "protocol": "ramda",
                "read_layout": "single-end",
                "genome_fasta": "/tmp/ref.fa",
                "gtf": "/tmp/genes.gtf",
                "detector_backend": "ciri3",
                "started_at": "2026-05-22T00:00:00+00:00",
                "completed_at": "2026-05-22T00:01:00+00:00",
                "elapsed_seconds": 60.0,
                "hostname": "host",
                "python_version": "3.12.0",
                "cleanup_summary": {"enabled": False, "scope": None, "performed": False, "deleted_paths": [], "reclaimed_bytes": 0},
                "completed_stages": ["workflow"],
                "skipped_stages": [],
                "failed_stages": [],
                "partial_outputs_detected": [],
                "paths": {
                    "manifest": str((workdir / "manifest.tsv").resolve()),
                    "alignment_manifest": str((workdir / "align.tsv").resolve()),
                    "detector_summary": str((workdir / "detector_run_summary.json").resolve()),
                    "matrix": str((workdir / "matrix" / "circ_counts.mtx").resolve()),
                    "cell_index": str((workdir / "matrix" / "cell_index.txt").resolve()),
                    "circ_index": str((workdir / "matrix" / "circ_index.txt").resolve()),
                    "h5ad": str((workdir / "anndata" / "circ_counts.h5ad").resolve()),
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (workdir / "manifest.tsv").write_text("sample_id\tfastq_1\nsample1\ta.fastq.gz\n", encoding="utf-8")
    (workdir / "align.tsv").write_text("cell_id\tsam\ncell1\ta.sam\n", encoding="utf-8")
    (workdir / "detector_run_summary.json").write_text("{}\n", encoding="utf-8")
    (workdir / "matrix" / "circ_counts.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n%\n1 1 0\n", encoding="utf-8")
    (workdir / "matrix" / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
    (workdir / "matrix" / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
    (workdir / "matrix" / "circ_feature_table.tsv").write_text(
        "circ_id\tchrom\tstart\tend\tstrand\thost_gene\ncirc1\tchr1\t1\t10\t+\tGENE1\n",
        encoding="utf-8",
    )
    (workdir / "anndata" / "circ_counts.h5ad").write_text("not-a-real-h5ad", encoding="utf-8")
    (workdir / "rna" / "gene_counts.tsv").write_text("gene_id\tgene_name\tcell1\nG1\tGENE1\t1\n", encoding="utf-8")
    (workdir / "rna" / "gene_feature_table.tsv").write_text("gene_id\tgene_name\nG1\tGENE1\n", encoding="utf-8")
    (workdir / "rna" / "rna_import_summary.json").write_text(json.dumps({"method": "simple-overlap"}) + "\n", encoding="utf-8")
    (workdir / "qc" / "rna_qc.tsv").write_text("cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\ncell1\t1\t1\t0\t0\t0\n", encoding="utf-8")
    (workdir / "qc" / "rna_gene_qc.tsv").write_text("gene_id\tgene_name\tn_cells_detected\ttotal_counts\tmean_counts_nonzero\tgene_biotype\nG1\tGENE1\t1\t1\t1\tprotein_coding\n", encoding="utf-8")
    if include_h5mu:
        (workdir / "mudata").mkdir()
        (workdir / "mudata" / "full_length.h5mu").write_text("not-a-real-h5mu", encoding="utf-8")
    return workdir


def test_print_environment_smoke() -> None:
    result = runner.invoke(app, ["print-environment"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert "circyto_version" in payload
    assert "python_version" in payload
    assert "timestamp" in payload


def test_validate_workdir_reports_success_or_clean_optional_warnings(tmp_path: Path) -> None:
    workdir = _write_validate_fixture(tmp_path)
    result = runner.invoke(app, ["validate-workdir", "--workdir", str(workdir), "--json"])
    assert result.exit_code != 0
    payload = json.loads(result.output)
    assert "h5ad unreadable" in " ".join(payload["errors"])


@pytest.mark.skipif(not HAS_MUDATA, reason="mudata not installed")
def test_mudata_environment_metadata_serialized(tmp_path: Path) -> None:
    import mudata as mu

    root = tmp_path / "workflow"
    (root / "rna").mkdir(parents=True)
    (root / "matrix").mkdir(parents=True)
    (root / "qc").mkdir(parents=True)
    (root / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\tcellA\tDIYHEK_192\n"
        "G1\tMT-CO1\t10\t1\n"
        "G2\tRPLP0\t0\t5\n"
        "G3\tGENE3\t3\t2\n",
        encoding="utf-8",
    )
    (root / "rna" / "gene_feature_table.tsv").write_text(
        "gene_id\tgene_name\tgene_biotype\n"
        "G1\tMT-CO1\tprotein_coding\n"
        "G2\tRPLP0\tprotein_coding\n"
        "G3\tGENE3\tlncRNA\n",
        encoding="utf-8",
    )
    (root / "rna" / "rna_import_summary.json").write_text(
        json.dumps({"method": "simple-overlap", "n_genes": 3, "n_cells": 2}, indent=2) + "\n",
        encoding="utf-8",
    )
    (root / "qc" / "rna_qc.tsv").write_text(
        "cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\n"
        "cellA\t13\t2\t0.7692307692\t0.0\t2\n"
        "DIYHEK_192\t8\t3\t0.125\t0.625\t0\n",
        encoding="utf-8",
    )
    (root / "qc" / "rna_circ_cell_summary.tsv").write_text(
        "cell_id\tmembership\ttotal_rna_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\tcircRNA_total_support\n"
        "cellA\tshared\t13\t2\t0.7692307692\t0.0\t2\t3\n"
        "DIYHEK_192\trna_only\t8\t3\t0.125\t0.625\t0\t0\n",
        encoding="utf-8",
    )
    (root / "qc" / "rna_circ_summary.json").write_text(
        json.dumps({"n_shared_cells": 1, "n_rna_only_cells": 1, "n_circ_only_cells": 0}, indent=2) + "\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n3 1 2\n1 1 2\n3 1 1\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_index.txt").write_text("circ1\ncirc2\ncirc3\n", encoding="utf-8")
    (root / "matrix" / "cell_index.txt").write_text("cellA\n", encoding="utf-8")
    (root / "matrix" / "circ_feature_table.tsv").write_text(
        "circ_id\tchrom\tstart\tend\tstrand\thost_gene\n"
        "circ1\tchr1\t1\t10\t+\tGENE1\n"
        "circ2\tchr1\t11\t20\t+\tGENE2\n"
        "circ3\tchr1\t21\t30\t-\tGENE3\n",
        encoding="utf-8",
    )
    (root / "workflow_summary.json").write_text(
        json.dumps(
            {
                "workflow": "full-length-circrna",
                "workflow_type": "full-length-circrna",
                "workflow_uuid": "123e4567-e89b-12d3-a456-426614174000",
                "protocol": "ramda",
                "read_layout": "single-end",
                "cleanup": {"delete_candidates": []},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    out_path = tmp_path / "full_length.h5mu"
    result = runner.invoke(app, ["export-mudata", "--workdir", str(root), "--output", str(out_path)])
    assert result.exit_code == 0, result.output
    mdata = read_h5mu(out_path)
    assert "environment" in mdata.uns["circyto"]
    assert mdata.uns["circyto"]["mudata_sync_policy"] == "pull_on_update=True"
    env = mdata.uns["circyto"]["environment"]
    assert "circyto_version" in env
    assert "export_timestamp" in env
