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
