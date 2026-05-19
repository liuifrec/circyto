from __future__ import annotations

import json
from pathlib import Path

import pytest
from scipy import sparse as sp
from scipy.io import mmwrite
from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_fastq(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")


def _write_ref_and_gtf(tmp_path: Path) -> tuple[Path, Path]:
    ref = tmp_path / "chr21.fa"
    gtf = tmp_path / "chr21.gtf"
    ref.write_text(">chr21\nACGT\n", encoding="utf-8")
    gtf.write_text('chr21\ttest\texon\t1\t4\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")
    return ref, gtf


def _write_single_end_manifest(tmp_path: Path, *, protocol: str = "ramda") -> Path:
    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "sample_id\tfastq_1\tfastq_2\tprotocol\tstrandedness\tread_layout\tn_input_reads",
                f"cell1\t{fastq}\t\t{protocol}\tunstranded\tsingle\t12",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


def _write_paired_manifest(tmp_path: Path, *, protocol: str = "smartseq3") -> Path:
    r1 = tmp_path / "reads_1.fastq"
    r2 = tmp_path / "reads_2.fastq"
    _write_fastq(r1)
    _write_fastq(r2)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "sample_id\tfastq_1\tfastq_2\tprotocol\tstrandedness\tread_layout",
                f"cell1\t{r1}\t{r2}\t{protocol}\tforward\tpaired-end",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


def _assert_disk_usage_summary(summary: dict[str, object]) -> None:
    for key in (
        "workdir_size_bytes",
        "align_size_bytes",
        "ciri3_size_bytes",
        "matrix_size_bytes",
        "anndata_size_bytes",
    ):
        assert key in summary
        assert isinstance(summary[key], int)
        assert int(summary[key]) >= 0
    largest_files = summary.get("largest_files")
    assert isinstance(largest_files, list)
    assert largest_files
    assert {"path", "bytes"}.issubset(largest_files[0].keys())
    assert int(largest_files[0]["bytes"]) >= 0


def test_full_length_workflow_dry_run_ramda_lists_demux_matrix_and_h5ad(tmp_path: Path, monkeypatch) -> None:
    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_single_end_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)

    def fake_run_ciri3_workflow(params, *, progress):
        logs = params.outdir / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        (logs / "run_ciri3_summary.json").write_text(
            json.dumps({"workflow": "run-ciri3", "dry_run": True}, indent=2) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    outdir = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.stdout
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["workflow"] == "full-length-circrna"
    assert summary["dry_run"] is True
    for key in (
        "command_name",
        "circyto_version",
        "workflow_type",
        "workflow_uuid",
        "started_at",
        "completed_at",
        "elapsed_seconds",
        "hostname",
        "python_version",
        "cleanup_summary",
        "completed_stages",
        "skipped_stages",
        "failed_stages",
        "partial_outputs_detected",
    ):
        assert key in summary
    assert summary["skip_demux_effective"] is True
    stages = {item["stage"]: item["status"] for item in summary["stage_graph"]}
    assert stages["demux"] == "skipped"
    assert stages["alignment"] == "planned"
    assert stages["detector"] == "planned"
    assert stages["matrix"] == "planned"
    assert stages["h5ad_export"] == "planned"


def test_full_length_workflow_dry_run_paired_ramda_plans_star_ciri3_route(tmp_path: Path, monkeypatch) -> None:
    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_paired_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    star_index = tmp_path / "star_index"
    star_index.mkdir()
    seen: dict[str, object] = {}

    def fake_run_ciri3_workflow(params, *, progress):
        seen["star_index"] = params.star_index
        seen["dry_run"] = params.dry_run
        seen["protocol"] = params.protocol
        logs = params.outdir / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        (logs / "run_ciri3_summary.json").write_text(
            json.dumps({"workflow": "run-ciri3", "dry_run": True}, indent=2) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    outdir = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--star-index",
            str(star_index),
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.stdout
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["dry_run"] is True
    assert summary["read_layout_counts"]["paired-end"] == 1
    assert summary["allow_paired_ramda"] is False
    assert summary["experimental_paired_ramda"] is False
    assert any("validated STAR+CIRI3 paired-end route" in warning for warning in summary["warnings"])
    assert seen["dry_run"] is True
    assert seen["protocol"] == "ramda"
    assert seen["star_index"] == star_index


def test_full_length_workflow_real_run_paired_ramda_requires_allow_flag(tmp_path: Path) -> None:
    manifest = _write_paired_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    star_index = tmp_path / "star_index"
    star_index.mkdir()

    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(tmp_path / "run"),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--star-index",
            str(star_index),
        ],
    )
    assert result.exit_code != 0
    assert result.exception is not None
    assert "--allow-paired-ramda" in str(result.exception)


def test_full_length_workflow_future_gene_expression_flags_fail_clearly(tmp_path: Path) -> None:
    manifest = _write_single_end_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)

    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(tmp_path / "run"),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--gene-expression-method",
            "featurecounts",
        ],
    )
    assert result.exit_code != 0
    assert result.exception is not None
    assert "featureCounts/velocyto-based" in str(result.exception)


def test_full_length_workflow_gene_counts_import_writes_rna_snapshot(tmp_path: Path, monkeypatch) -> None:
    pytest.importorskip("anndata")
    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_single_end_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    gene_counts = tmp_path / "gene_counts.tsv"
    gene_counts.write_text(
        "gene_id\tgene_name\tcell1\n"
        "ENSG1\tGENE1\t5\n"
        "ENSG2\tGENE2\t0\n",
        encoding="utf-8",
    )

    def fake_run_ciri3_workflow(params, *, progress):
        align = params.outdir / "align"
        ciri3 = params.outdir / "ciri3"
        matrix = params.outdir / "matrix"
        logs = params.outdir / "logs"
        for path in (align, ciri3, matrix, logs):
            path.mkdir(parents=True, exist_ok=True)
        (align / "alignment_manifest.tsv").write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
            "cell1\t\t/tmp/cell1.sam\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted\n",
            encoding="utf-8",
        )
        (align / "alignment_prepare_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"aligned": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "cells": [{"cell_id": "cell1", "status": "aligned", "seconds": 1.2}],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (ciri3 / "cell1.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\ncirc1\tchr21\t10\t20\t+\t4\n",
            encoding="utf-8",
        )
        (ciri3 / "detector_run_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"success": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "elapsed_seconds": 2.4,
                    "cells": [
                        {
                            "cell_id": "cell1",
                            "status": "success",
                            "seconds": 2.4,
                            "raw_row_count": 1,
                            "normalized_row_count": 1,
                        }
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        X = sp.csr_matrix([[4]], dtype=int)
        mmwrite(matrix / "circ_counts.mtx", X)
        (matrix / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
        (matrix / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
        (matrix / "circ_feature_table.tsv").write_text(
            "circ_id\tchrom\tstart\tend\tstrand\thost_gene\ncirc1\tchr21\t10\t20\t+\tGENE1\n",
            encoding="utf-8",
        )
        (logs / "run_ciri3_summary.json").write_text(
            json.dumps({"workflow": "run-ciri3", "dry_run": False}, indent=2) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    outdir = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--gene-counts",
            str(gene_counts),
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert (outdir / "rna" / "gene_counts.tsv").exists()
    assert (outdir / "rna" / "gene_feature_table.tsv").exists()
    assert (outdir / "rna" / "rna_import_summary.json").exists()
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["rna_import"]["n_genes"] == 2
    assert summary["paths"]["rna_dir"] == str((outdir / "rna").resolve())


def test_full_length_workflow_cleanup_intermediates_dry_run_plans_without_deleting(tmp_path: Path, monkeypatch) -> None:
    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_single_end_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    align_dir = tmp_path / "run" / "align"
    align_dir.mkdir(parents=True, exist_ok=True)
    planned_sam = align_dir / "planned.sam"
    planned_sam.write_text("sam", encoding="utf-8")

    def fake_run_ciri3_workflow(params, *, progress):
        logs = params.outdir / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        (logs / "run_ciri3_summary.json").write_text(
            json.dumps({"workflow": "run-ciri3", "dry_run": True}, indent=2) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    outdir = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--cleanup-intermediates",
            "alignments",
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert planned_sam.exists()
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    cleanup = summary["cleanup"]
    assert cleanup["enabled"] is True
    assert cleanup["mode"] == "planned-only"
    assert cleanup["planned_scope"] == "alignments"
    assert cleanup["candidate_count"] >= 1
    assert any(item["path"].endswith("planned.sam") for item in cleanup["delete_candidates"])


def test_full_length_workflow_simple_overlap_writes_rna_snapshot(tmp_path: Path, monkeypatch) -> None:
    pytest.importorskip("anndata")
    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_single_end_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    gtf.write_text(
        'chr21\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE1";\n',
        encoding="utf-8",
    )
    sam = tmp_path / "cell1.sam"
    sam.write_text(
        "@SQ\tSN:chr21\tLN:1000\n"
        "read1\t0\tchr21\t1\t255\t4M\t*\t0\t0\tACGT\t!!!!\n",
        encoding="utf-8",
    )

    def fake_run_ciri3_workflow(params, *, progress):
        align = params.outdir / "align"
        ciri3 = params.outdir / "ciri3"
        matrix = params.outdir / "matrix"
        logs = params.outdir / "logs"
        for path in (align, ciri3, matrix, logs):
            path.mkdir(parents=True, exist_ok=True)
        (align / "alignment_manifest.tsv").write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
            f"cell1\t\t{sam}\tgrp\tsingle-end\tbwa-mem\t{ref}\tk1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted\n",
            encoding="utf-8",
        )
        (align / "alignment_prepare_summary.json").write_text(
            json.dumps(
                {"status_counts": {"aligned": 1}, "completed_cells": 1, "failed_cells": 0, "cells": [{"cell_id": "cell1", "status": "aligned"}]},
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (ciri3 / "cell1.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\ncirc1\tchr21\t1\t4\t+\t1\n",
            encoding="utf-8",
        )
        (ciri3 / "detector_run_summary.json").write_text(
            json.dumps(
                {"status_counts": {"success": 1}, "completed_cells": 1, "failed_cells": 0, "elapsed_seconds": 0.1, "cells": [{"cell_id": "cell1", "status": "success", "raw_row_count": 1, "normalized_row_count": 1}]},
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        X = sp.csr_matrix([[1]], dtype=int)
        mmwrite(matrix / "circ_counts.mtx", X)
        (matrix / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
        (matrix / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
        (matrix / "circ_feature_table.tsv").write_text(
            "circ_id\tchrom\tstart\tend\tstrand\thost_gene\ncirc1\tchr21\t1\t4\t+\tGENE1\n",
            encoding="utf-8",
        )
        (logs / "run_ciri3_summary.json").write_text(json.dumps({"workflow": "run-ciri3", "dry_run": False}) + "\n", encoding="utf-8")

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    outdir = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--gene-expression-method",
            "simple-overlap",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert (outdir / "rna" / "gene_counts.tsv").exists()
    assert (outdir / "rna" / "gene_feature_table.tsv").exists()
    assert (outdir / "rna" / "rna_import_summary.json").exists()
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["rna_import"]["method"] == "simple-overlap"
    assert summary["rna_import"]["assigned_templates"] == 1


def test_full_length_workflow_cleanup_execution_deletes_alignment_intermediates(tmp_path: Path, monkeypatch) -> None:
    pytest.importorskip("anndata")
    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_single_end_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)

    def fake_run_ciri3_workflow(params, *, progress):
        align = params.outdir / "align"
        ciri3 = params.outdir / "ciri3"
        matrix = params.outdir / "matrix"
        logs = params.outdir / "logs"
        for path in (align, ciri3, matrix, logs):
            path.mkdir(parents=True, exist_ok=True)
        doomed = align / "planned.sam"
        doomed.write_text("sam", encoding="utf-8")
        (align / "alignment_manifest.tsv").write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
            f"cell1\t\t{doomed}\tgrp\tsingle-end\tbwa-mem\t{ref}\tk1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted\n",
            encoding="utf-8",
        )
        (align / "alignment_prepare_summary.json").write_text(
            json.dumps({"status_counts": {"aligned": 1}, "completed_cells": 1, "failed_cells": 0, "cells": [{"cell_id": "cell1", "status": "aligned"}]}, indent=2) + "\n",
            encoding="utf-8",
        )
        (ciri3 / "cell1.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\ncirc1\tchr21\t1\t4\t+\t1\n",
            encoding="utf-8",
        )
        (ciri3 / "detector_run_summary.json").write_text(
            json.dumps({"status_counts": {"success": 1}, "completed_cells": 1, "failed_cells": 0, "elapsed_seconds": 0.1, "cells": [{"cell_id": "cell1", "status": "success", "raw_row_count": 1, "normalized_row_count": 1}]}, indent=2) + "\n",
            encoding="utf-8",
        )
        X = sp.csr_matrix([[1]], dtype=int)
        mmwrite(matrix / "circ_counts.mtx", X)
        (matrix / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
        (matrix / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
        (matrix / "circ_feature_table.tsv").write_text(
            "circ_id\tchrom\tstart\tend\tstrand\thost_gene\ncirc1\tchr21\t1\t4\t+\tGENE1\n",
            encoding="utf-8",
        )
        (logs / "run_ciri3_summary.json").write_text(json.dumps({"workflow": "run-ciri3", "dry_run": False}) + "\n", encoding="utf-8")

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    outdir = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--cleanup-intermediates",
            "alignments",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert not (outdir / "align" / "planned.sam").exists()
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["cleanup_performed"] is True
    assert summary["cleanup_scope"] == "alignments"
    assert summary["cleanup_reclaimed_bytes"] >= 3


def test_full_length_workflow_failed_run_does_not_delete_intermediates(tmp_path: Path, monkeypatch) -> None:
    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_single_end_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    doomed = tmp_path / "run" / "align" / "failed.sam"
    doomed.parent.mkdir(parents=True, exist_ok=True)
    doomed.write_text("sam", encoding="utf-8")

    def fake_run_ciri3_workflow(params, *, progress):
        raise RuntimeError("synthetic failure")

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(tmp_path / "run"),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--cleanup-intermediates",
            "alignments",
        ],
    )
    assert result.exit_code != 0
    assert doomed.exists()


def test_full_length_workflow_real_run_writes_h5ad_and_qc(tmp_path: Path, monkeypatch) -> None:
    pytest.importorskip("anndata")
    import anndata as ad

    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_single_end_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)

    def fake_run_ciri3_workflow(params, *, progress):
        align = params.outdir / "align"
        ciri3 = params.outdir / "ciri3"
        matrix = params.outdir / "matrix"
        logs = params.outdir / "logs"
        for path in (align, ciri3, matrix, logs):
            path.mkdir(parents=True, exist_ok=True)
        (align / "alignment_manifest.tsv").write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
            "cell1\t\t/tmp/cell1.sam\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted\n",
            encoding="utf-8",
        )
        (align / "alignment_prepare_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"aligned": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "cells": [{"cell_id": "cell1", "status": "aligned", "seconds": 1.2}],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (ciri3 / "cell1.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\ncirc1\tchr21\t10\t20\t+\t4\n",
            encoding="utf-8",
        )
        (ciri3 / "detector_run_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"success": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "elapsed_seconds": 2.4,
                    "cells": [
                        {
                            "cell_id": "cell1",
                            "status": "success",
                            "seconds": 2.4,
                            "raw_row_count": 1,
                            "normalized_row_count": 1,
                        }
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        X = sp.csr_matrix([[4]], dtype=int)
        mmwrite(matrix / "circ_counts.mtx", X)
        (matrix / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
        (matrix / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
        (matrix / "circ_feature_table.tsv").write_text(
            "circ_id\tchrom\tstart\tend\tstrand\thost_gene\ncirc1\tchr21\t10\t20\t+\tGENE1\n",
            encoding="utf-8",
        )
        (logs / "run_ciri3_summary.json").write_text(
            json.dumps({"workflow": "run-ciri3", "dry_run": False}, indent=2) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    outdir = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--export-h5ad",
        ],
    )
    assert result.exit_code == 0, result.stdout
    h5ad_path = outdir / "anndata" / "circ_counts.h5ad"
    assert h5ad_path.exists()
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["paths"]["h5ad"] == str(h5ad_path.resolve())
    assert (outdir / "qc" / "cell_qc.tsv").exists()
    assert (outdir / "qc" / "circ_qc.tsv").exists()
    _assert_disk_usage_summary(summary)
    assert summary["workflow_type"] == "full-length-circrna"
    assert summary["cleanup_summary"]["performed"] is False
    assert isinstance(summary["workflow_uuid"], str)
    assert "summary_qc" in summary["completed_stages"]
    assert summary["workdir_size_bytes"] >= summary["matrix_size_bytes"]
    assert summary["anndata_size_bytes"] > 0
    adata = ad.read_h5ad(h5ad_path)
    assert list(adata.obs_names) == ["cell1"]
    assert adata.obs["protocol"].tolist() == ["ramda"]
    assert adata.obs["read_layout"].tolist() == ["single-end"]
    assert list(adata.var_names) == ["circ1"]


def test_full_length_workflow_paired_ramda_with_allow_flag_writes_h5ad_and_qc(tmp_path: Path, monkeypatch) -> None:
    pytest.importorskip("anndata")
    import anndata as ad

    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_paired_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    star_index = tmp_path / "star_index"
    star_index.mkdir()
    seen: dict[str, object] = {}

    def fake_run_ciri3_workflow(params, *, progress):
        seen["star_index"] = params.star_index
        seen["dry_run"] = params.dry_run
        align = params.outdir / "align"
        ciri3 = params.outdir / "ciri3"
        matrix = params.outdir / "matrix"
        logs = params.outdir / "logs"
        for path in (align, ciri3, matrix, logs):
            path.mkdir(parents=True, exist_ok=True)
        (align / "alignment_manifest.tsv").write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\tchimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\n"
            "cell1\t\t/tmp/cell1.star.sam\tgrp\tpaired-end\tstar\t/tmp/ref.fa\tk1\t/tmp/source.tsv\t1\tstar\tunsorted\t/tmp/cell1.Chimeric.out.junction\t/tmp/cell1.Unmapped.out.mate1\t/tmp/cell1.Unmapped.out.mate2\t/tmp/cell1.bwa_rescue.sam\n",
            encoding="utf-8",
        )
        (align / "alignment_prepare_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"aligned": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "cells": [{"cell_id": "cell1", "status": "aligned", "seconds": 1.5}],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (ciri3 / "cell1.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\ncirc1\tchr21\t10\t20\t+\t4\n",
            encoding="utf-8",
        )
        (ciri3 / "detector_run_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"success": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "elapsed_seconds": 2.7,
                    "cells": [
                        {
                            "cell_id": "cell1",
                            "status": "success",
                            "seconds": 2.7,
                            "raw_row_count": 1,
                            "normalized_row_count": 1,
                        }
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        X = sp.csr_matrix([[4]], dtype=int)
        mmwrite(matrix / "circ_counts.mtx", X)
        (matrix / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
        (matrix / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
        (matrix / "circ_feature_table.tsv").write_text(
            "circ_id\tchrom\tstart\tend\tstrand\thost_gene\ncirc1\tchr21\t10\t20\t+\tGENE1\n",
            encoding="utf-8",
        )
        (logs / "run_ciri3_summary.json").write_text(
            json.dumps({"workflow": "run-ciri3", "dry_run": False}, indent=2) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    outdir = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--star-index",
            str(star_index),
            "--allow-paired-ramda",
            "--export-h5ad",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert seen["dry_run"] is False
    assert seen["star_index"] == star_index
    h5ad_path = outdir / "anndata" / "circ_counts.h5ad"
    assert h5ad_path.exists()
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["allow_paired_ramda"] is True
    assert summary["experimental_paired_ramda"] is True
    assert summary["paths"]["h5ad"] == str(h5ad_path.resolve())
    assert any("validated STAR+CIRI3 paired-end route" in warning for warning in summary["warnings"])
    _assert_disk_usage_summary(summary)
    assert summary["workdir_size_bytes"] >= summary["align_size_bytes"]
    assert summary["anndata_size_bytes"] > 0
    adata = ad.read_h5ad(h5ad_path)
    assert adata.obs["protocol"].tolist() == ["ramda"]
    assert adata.obs["read_layout"].tolist() == ["paired-end"]
    assert list(adata.var_names) == ["circ1"]


def test_full_length_workflow_paired_ramda_deprecated_flag_alias_still_works(tmp_path: Path, monkeypatch) -> None:
    from circyto.pipeline import workflow_full_length_circrna as workflow

    manifest = _write_paired_manifest(tmp_path, protocol="ramda")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    star_index = tmp_path / "star_index"
    star_index.mkdir()
    seen: dict[str, object] = {}

    def fake_run_ciri3_workflow(params, *, progress):
        seen["star_index"] = params.star_index
        align = params.outdir / "align"
        ciri3 = params.outdir / "ciri3"
        matrix = params.outdir / "matrix"
        logs = params.outdir / "logs"
        for path in (align, ciri3, matrix, logs):
            path.mkdir(parents=True, exist_ok=True)
        (align / "alignment_manifest.tsv").write_text(
            "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
            "cell1\t\t/tmp/cell1.star.sam\tgrp\tpaired-end\tstar\t/tmp/ref.fa\tk1\t/tmp/source.tsv\t1\tstar\tunsorted\n",
            encoding="utf-8",
        )
        (align / "alignment_prepare_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"aligned": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "cells": [{"cell_id": "cell1", "status": "aligned", "seconds": 1.5}],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        (ciri3 / "cell1.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\ncirc1\tchr21\t10\t20\t+\t4\n",
            encoding="utf-8",
        )
        (ciri3 / "detector_run_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"success": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "elapsed_seconds": 1.0,
                    "cells": [
                        {
                            "cell_id": "cell1",
                            "status": "success",
                            "seconds": 1.0,
                            "raw_row_count": 1,
                            "normalized_row_count": 1,
                        }
                    ],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        mmwrite(matrix / "circ_counts.mtx", sp.csr_matrix([[4]], dtype=int))
        (matrix / "circ_index.txt").write_text("circ1\n", encoding="utf-8")
        (matrix / "cell_index.txt").write_text("cell1\n", encoding="utf-8")
        (matrix / "circ_feature_table.tsv").write_text(
            "circ_id\tchrom\tstart\tend\tstrand\thost_gene\ncirc1\tchr21\t10\t20\t+\tGENE1\n",
            encoding="utf-8",
        )
        (logs / "run_ciri3_summary.json").write_text(
            json.dumps({"workflow": "run-ciri3", "dry_run": False}, indent=2) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(workflow, "run_ciri3_workflow", fake_run_ciri3_workflow)

    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(tmp_path / "run"),
            "--protocol",
            "ramda",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--star-index",
            str(star_index),
            "--experimental-paired-ramda",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert seen["star_index"] == star_index


def test_full_length_workflow_smartseq3_requires_skip_demux(tmp_path: Path) -> None:
    manifest = _write_paired_manifest(tmp_path, protocol="smartseq3")
    ref, gtf = _write_ref_and_gtf(tmp_path)
    star_index = tmp_path / "star_index"
    star_index.mkdir()

    result = runner.invoke(
        app,
        [
            "workflow",
            "full-length-circrna",
            "--manifest",
            str(manifest),
            "--outdir",
            str(tmp_path / "run"),
            "--protocol",
            "smartseq3",
            "--genome-fasta",
            str(ref),
            "--gtf",
            str(gtf),
            "--star-index",
            str(star_index),
        ],
    )
    assert result.exit_code != 0
    assert result.exception is not None
    assert "smartseq3-ciri3" in str(result.exception)


def test_existing_smartseq3_workflow_command_still_present() -> None:
    result = runner.invoke(app, ["workflow", "smartseq3-ciri3", "--help"])
    assert result.exit_code == 0
    assert "--read1" in result.stdout
    assert "--index1" in result.stdout
