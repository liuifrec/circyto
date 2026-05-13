from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
import pytest
from scipy import io as scio
from scipy import sparse as sp
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.manifest.alignment import AlignmentManifestRow, write_alignment_manifest_tsv
from circyto.pipeline.collect import collect_matrix
from circyto.pipeline.workflow_reporting import (
    HAS_MUDATA,
    build_cell_qc_table,
    build_circ_qc_table,
    export_circ_h5ad,
    export_mudata_bundle,
    load_circ_feature_table,
    numeric_summary,
)
from circyto.pipeline.workflow_smartseq3_ciri3 import create_selected_manifest


runner = CliRunner()


def _open_maybe_gz(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with _open_maybe_gz(path, "wt") as handle:
        for read_id, seq in records:
            handle.write(f"@{read_id}\n{seq}\n+\n{'I' * len(seq)}\n")


def _build_cli_inputs(tmp_path: Path) -> dict[str, Path]:
    read1 = tmp_path / "R1.fastq.gz"
    read2 = tmp_path / "R2.fastq.gz"
    index1 = tmp_path / "I1.fastq.gz"
    index2 = tmp_path / "I2.fastq.gz"
    annotation = tmp_path / "annotation.tsv"
    ref_fa = tmp_path / "ref.fa"
    star_index = tmp_path / "star_index"

    _write_fastq(read1, [("read1/1", "ACGT"), ("read2/1", "TGCA"), ("read3/1", "CCCC")])
    _write_fastq(read2, [("read1/2", "TTTT"), ("read2/2", "GGGG"), ("read3/2", "AAAA")])
    _write_fastq(index1, [("read1", "AAAA"), ("read2", "CCCC"), ("read3", "AAAA")])
    _write_fastq(index2, [("read1", "TTTT"), ("read2", "GGGG"), ("read3", "TTTT")])
    annotation.write_text(
        "cell_id\tindex1\tindex2\n"
        "cellA\tAAAA\tTTTT\n"
        "cellB\tCCCC\tGGGG\n",
        encoding="utf-8",
    )
    ref_fa.write_text(">chr1\nACGT\n", encoding="utf-8")
    star_index.mkdir()
    return {
        "read1": read1,
        "read2": read2,
        "index1": index1,
        "index2": index2,
        "annotation": annotation,
        "ref_fa": ref_fa,
        "star_index": star_index,
    }


def _write_matrix_bundle(
    tmp_path: Path,
    *,
    circ_ids: list[str],
    cell_ids: list[str],
    rows: list[tuple[int, int, int]],
) -> tuple[Path, Path, Path, Path]:
    matrix_path = tmp_path / "circ_counts.mtx"
    circ_index_path = tmp_path / "circ_index.txt"
    cell_index_path = tmp_path / "cell_index.txt"
    feature_path = tmp_path / "circ_feature_table.tsv"
    shape = (len(circ_ids), len(cell_ids))
    if rows:
        data = [value for _, _, value in rows]
        row_idx = [row for row, _, _ in rows]
        col_idx = [col for _, col, _ in rows]
        X = sp.csr_matrix((data, (row_idx, col_idx)), shape=shape, dtype=int)
    else:
        X = sp.csr_matrix(shape, dtype=int)
    scio.mmwrite(str(matrix_path), X)
    circ_index_path.write_text("\n".join(circ_ids) + ("\n" if circ_ids else ""), encoding="utf-8")
    cell_index_path.write_text("\n".join(cell_ids) + ("\n" if cell_ids else ""), encoding="utf-8")
    pd.DataFrame(
        [
            {"circ_id": circ_id, "chrom": "chr1", "start": idx * 10 + 1, "end": idx * 10 + 5, "strand": "+", "host_gene": f"GENE{idx + 1}"}
            for idx, circ_id in enumerate(circ_ids)
        ]
    ).to_csv(feature_path, sep="\t", index=False)
    return matrix_path, circ_index_path, cell_index_path, feature_path


def _write_gene_counts_tsv(path: Path) -> None:
    path.write_text(
        "gene_id\tcellA\tcellB\tcellC\n"
        "GENE1\t10\t0\t5\n"
        "GENE2\t0\t3\t7\n",
        encoding="utf-8",
    )


def test_create_selected_manifest_top_n(tmp_path: Path) -> None:
    inputs = _build_cli_inputs(tmp_path)
    from circyto.demux.smartseq3 import SmartSeq3DemuxParams, demux_smartseq3

    demux_out = tmp_path / "demux"
    demux_smartseq3(
        SmartSeq3DemuxParams(
            read1=inputs["read1"],
            read2=inputs["read2"],
            index1=inputs["index1"],
            index2=inputs["index2"],
            annotation=inputs["annotation"],
            outdir=demux_out,
            cell_id_column="cell_id",
            index1_column="index1",
            index2_column="index2",
            write_sink=False,
        )
    )

    selection = create_selected_manifest(
        demux_manifest=demux_out / "manifest.tsv",
        demux_summary_path=demux_out / "demux_summary.json",
        output_path=tmp_path / "manifests" / "top1_manifest.tsv",
        top_n=1,
    )

    assert selection["selected_cell_count"] == 1
    assert selection["selected_cell_ids"] == ["cellA"]
    manifest_text = (tmp_path / "manifests" / "top1_manifest.tsv").read_text(encoding="utf-8")
    assert "cellA\tsmartseq3" in manifest_text
    assert "cellB\tsmartseq3" not in manifest_text


def test_workflow_resume_skips_completed_stages(tmp_path: Path, monkeypatch) -> None:
    from circyto.pipeline import workflow_smartseq3_ciri3 as workflow

    inputs = _build_cli_inputs(tmp_path)
    outdir = tmp_path / "workflow"
    demux_dir = outdir / "demux"
    manifests_dir = outdir / "manifests"
    align_dir = outdir / "align"
    ciri3_dir = outdir / "ciri3"
    matrix_dir = outdir / "matrix"
    for path in (demux_dir, manifests_dir, align_dir, ciri3_dir, matrix_dir, outdir / "logs"):
        path.mkdir(parents=True, exist_ok=True)

    cella_r1 = demux_dir / "fastq" / "cellA_R1.fastq.gz"
    cella_r2 = demux_dir / "fastq" / "cellA_R2.fastq.gz"
    _write_fastq(cella_r1, [("read1/1", "ACGT")])
    _write_fastq(cella_r2, [("read1/2", "TTTT")])
    (demux_dir / "manifest.tsv").write_text(
        "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tread_layout\n"
        f"cellA\tsmartseq3\t{cella_r1}\t{cella_r2}\t\tsmartseq3_demux\t2\tpaired-end\n",
        encoding="utf-8",
    )
    (demux_dir / "demux_summary.json").write_text(
        json.dumps(
            {
                "assigned_records": 2,
                "total_records": 2,
                "number_of_cells_detected": 1,
                "reads_per_cell": {"cellA": 2},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (manifests_dir / "top1_manifest.tsv").write_text((demux_dir / "manifest.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    bam = align_dir / "aligned" / "cellA.bam"
    bam.parent.mkdir(parents=True, exist_ok=True)
    bam.write_text("bam", encoding="utf-8")
    write_alignment_manifest_tsv(
        [
            AlignmentManifestRow(
                cell_id="cellA",
                bam=str(bam),
                group_id="smartseq3_demux",
                read_layout="paired-end",
                aligner="star",
                reference=str(inputs["ref_fa"]),
                cache_key="k1",
                    source_manifest=str((manifests_dir / "top1_manifest.tsv").resolve()),
            )
        ],
        align_dir / "alignment_manifest.tsv",
    )
    (align_dir / "alignment_prepare_summary.json").write_text(
        json.dumps(
            {
                "status_counts": {"aligned": 1},
                "completed_cells": 1,
                "failed_cells": 0,
                "sentinel_cells": 0,
                "elapsed_seconds": 1.2,
                "cells": [{"cell_id": "cellA", "status": "aligned", "seconds": 1.2}],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (ciri3_dir / "cellA.tsv").write_text(
        "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
        "circA\tchr1\t10\t20\t+\t3\n",
        encoding="utf-8",
    )
    (ciri3_dir / "detector_run_summary.json").write_text(
        json.dumps(
            {
                "status_counts": {"success": 1},
                "completed_cells": 1,
                "failed_cells": 0,
                "elapsed_seconds": 0.8,
                "cells": [{"cell_id": "cellA", "status": "success", "seconds": 0.8}],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (matrix_dir / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n"
        "1 1 1\n"
        "1 1 3\n",
        encoding="utf-8",
    )
    (matrix_dir / "circ_index.txt").write_text("circA\n", encoding="utf-8")
    (matrix_dir / "cell_index.txt").write_text("cellA\n", encoding="utf-8")

    def _fail(*args, **kwargs):
        raise AssertionError("resume should skip stage execution")

    monkeypatch.setattr(workflow, "demux_smartseq3", _fail)
    monkeypatch.setattr(workflow, "prepare_alignment_cache", _fail)
    monkeypatch.setattr(workflow, "run_detector_alignment_manifest", _fail)
    monkeypatch.setattr(workflow, "collect_matrix", _fail)

    result = runner.invoke(
        app,
        [
            "workflow",
            "smartseq3-ciri3",
            "--read1",
            str(inputs["read1"]),
            "--read2",
            str(inputs["read2"]),
            "--index1",
            str(inputs["index1"]),
            "--index2",
            str(inputs["index2"]),
            "--annotation",
            str(inputs["annotation"]),
            "--cell-id-column",
            "cell_id",
            "--index1-column",
            "index1",
            "--index2-column",
            "index2",
            "--ref-fa",
            str(inputs["ref_fa"]),
            "--star-index",
            str(inputs["star_index"]),
            "--outdir",
            str(outdir),
            "--top-n",
            "1",
            "--resume",
        ],
    )

    assert result.exit_code == 0
    assert "resume skip: demux" in result.stdout
    assert "resume skip: alignment prep" in result.stdout
    assert "resume skip: detector" in result.stdout
    assert "resume skip: matrix" in result.stdout


def test_workflow_writes_summary_json(tmp_path: Path, monkeypatch) -> None:
    from circyto.pipeline import workflow_smartseq3_ciri3 as workflow

    inputs = _build_cli_inputs(tmp_path)
    outdir = tmp_path / "workflow"

    def _fake_prepare_alignment_cache(**kwargs):
        align_outdir = kwargs["outdir"]
        manifest = kwargs["manifest"]
        bam = align_outdir / "aligned" / "cellA.bam"
        bam.parent.mkdir(parents=True, exist_ok=True)
        bam.write_text("bam", encoding="utf-8")
        write_alignment_manifest_tsv(
            [
                AlignmentManifestRow(
                    cell_id="cellA",
                    bam=str(bam),
                    group_id="smartseq3_demux",
                    read_layout="paired-end",
                    aligner="star",
                    reference=str(kwargs["ref_fa"]),
                    cache_key="k1",
                    source_manifest=str(manifest.resolve()),
                    chimeric_junction="",
                    unmapped_mate1="",
                    unmapped_mate2="",
                    bwa_sam="",
                    mapper_mode="1",
                    artifact_bucket="star_bam",
                    sortedness="sorted",
                )
            ],
            align_outdir / "alignment_manifest.tsv",
        )
        (align_outdir / "alignment_prepare_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"aligned": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "sentinel_cells": 0,
                    "elapsed_seconds": 1.5,
                    "cells": [{"cell_id": "cellA", "status": "aligned", "seconds": 1.5}],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return align_outdir / "alignment_manifest.tsv"

    def _fake_run_detector_alignment_manifest(**kwargs):
        detector_outdir = kwargs["outdir"]
        detector_outdir.mkdir(parents=True, exist_ok=True)
        (detector_outdir / "cellA.tsv").write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
            "circA\tchr1\t10\t20\t+\t5\n"
            "circB\tchr1\t30\t40\t-\t2\n",
            encoding="utf-8",
        )
        (detector_outdir / "detector_run_summary.json").write_text(
            json.dumps(
                {
                    "status_counts": {"success": 1},
                    "completed_cells": 1,
                    "failed_cells": 0,
                    "elapsed_seconds": 0.7,
                    "cells": [{"cell_id": "cellA", "status": "success", "seconds": 0.7}],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return []

    monkeypatch.setattr(workflow, "prepare_alignment_cache", _fake_prepare_alignment_cache)
    monkeypatch.setattr(workflow, "run_detector_alignment_manifest", _fake_run_detector_alignment_manifest)

    result = runner.invoke(
        app,
        [
            "workflow",
            "smartseq3-ciri3",
            "--read1",
            str(inputs["read1"]),
            "--read2",
            str(inputs["read2"]),
            "--index1",
            str(inputs["index1"]),
            "--index2",
            str(inputs["index2"]),
            "--annotation",
            str(inputs["annotation"]),
            "--cell-id-column",
            "cell_id",
            "--index1-column",
            "index1",
            "--index2-column",
            "index2",
            "--ref-fa",
            str(inputs["ref_fa"]),
            "--star-index",
            str(inputs["star_index"]),
            "--outdir",
            str(outdir),
            "--top-n",
            "1",
            "--no-write-sink",
        ],
    )

    assert result.exit_code == 0
    summary = json.loads((outdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["workflow"] == "smartseq3-ciri3"
    assert summary["experimental"] is True
    assert summary["demux"]["assigned_records"] == 3
    assert summary["demux"]["assignment_rate"] == 1.0
    assert summary["selected_cell_count"] == 1
    assert summary["alignment_status_counts"] == {"aligned": 1}
    assert summary["detector_status_counts"] == {"success": 1}
    assert summary["matrix"]["n_rows"] == 2
    assert summary["matrix"]["n_circRNAs"] == 2
    assert summary["matrix"]["n_cells"] == 1
    assert summary["matrix"]["nnz"] == 2
    assert summary["detector"]["success"] == 1
    assert summary["paths"]["cell_qc"].endswith("cell_qc.tsv")
    assert summary["paths"]["circ_qc"].endswith("circ_qc.tsv")
    assert (outdir / "qc" / "cell_qc.tsv").exists()
    assert (outdir / "qc" / "circ_qc.tsv").exists()
    assert (outdir / "anndata" / "circ_counts.h5ad").exists()
    assert summary["command_options"]["top_n"] == 1
    assert summary["paths"]["selected_manifest"].endswith("top1_manifest.tsv")


def test_qc_summary_generation_from_small_matrix(tmp_path: Path) -> None:
    circ_ids = ["circA", "circB"]
    cell_ids = ["cellA", "cellB", "cellC"]
    _, _, _, feature_path = _write_matrix_bundle(
        tmp_path,
        circ_ids=circ_ids,
        cell_ids=cell_ids,
        rows=[(0, 0, 5), (1, 0, 1), (0, 1, 2)],
    )
    feature_df = load_circ_feature_table(circ_ids, feature_path)
    X_cells_by_circ = sp.csr_matrix([[5, 1], [2, 0], [0, 0]], dtype=int)
    cell_qc = build_cell_qc_table(
        selected_cell_ids=cell_ids,
        assigned_reads={"cellA": 10, "cellB": 6, "cellC": 1},
        alignment_cells={"cellA": {"status": "aligned", "seconds": 1.0}, "cellB": {"status": "aligned", "seconds": 1.5}},
        detector_cells={
            "cellA": {"status": "success", "seconds": 0.5},
            "cellB": {"status": "success", "seconds": 0.6},
            "cellC": {"status": "empty", "seconds": 0.2},
        },
        X_cells_by_circ=X_cells_by_circ,
    )
    circ_qc = build_circ_qc_table(circ_ids=circ_ids, feature_df=feature_df, X_cells_by_circ=X_cells_by_circ)

    assert cell_qc.loc["cellA", "circRNA_count"] == 2
    assert cell_qc.loc["cellB", "circRNA_count"] == 1
    assert bool(cell_qc.loc["cellC", "empty_flag"]) is True
    assert circ_qc.loc["circA", "n_cells_detected"] == 2
    assert circ_qc.loc["circA", "total_support"] == 7
    assert circ_qc.loc["circB", "max_support"] == 1
    assert numeric_summary([2, 1, 0]) == {"min": 0.0, "median": 1.0, "mean": 1.0, "max": 2.0}


def test_circ_qc_blank_host_gene_preserves_columns_and_cell_bounds(tmp_path: Path) -> None:
    circ_ids = ["circA", "circB"]
    feature_path = tmp_path / "circ_feature_table.tsv"
    feature_path.write_text(
        "circ_id\tchrom\tstart\tend\tstrand\thost_gene\n"
        "circA\tchr1\t1\t5\t+\t\n"
        "circB\tchr1\t11\t15\t-\tGENE2\n",
        encoding="utf-8",
    )

    feature_df = load_circ_feature_table(circ_ids, feature_path)
    X_cells_by_circ = sp.csr_matrix([[5, 1], [2, 0], [0, 0]], dtype=int)
    circ_qc = build_circ_qc_table(circ_ids=circ_ids, feature_df=feature_df, X_cells_by_circ=X_cells_by_circ)

    assert circ_qc.columns.tolist() == [
        "chrom",
        "start",
        "end",
        "strand",
        "host_gene",
        "n_cells_detected",
        "total_support",
        "max_support",
        "mean_support_detected_cells",
    ]
    assert circ_qc.loc["circA", "host_gene"] == ""
    assert int(circ_qc["n_cells_detected"].max()) <= X_cells_by_circ.shape[0]

    out_path = tmp_path / "circ_qc.tsv"
    circ_qc.reset_index().to_csv(out_path, sep="\t", index=False)
    written = pd.read_csv(out_path, sep="\t", keep_default_na=False)

    assert written.columns.tolist() == [
        "circ_id",
        "chrom",
        "start",
        "end",
        "strand",
        "host_gene",
        "n_cells_detected",
        "total_support",
        "max_support",
        "mean_support_detected_cells",
    ]
    assert written.loc[0, "host_gene"] == ""
    assert written.loc[0, "n_cells_detected"] == 2
    assert written.loc[0, "total_support"] == 7
    assert written.loc[0, "max_support"] == 5


def test_collect_and_qc_preserve_string_host_gene_without_numeric_overwrite(tmp_path: Path) -> None:
    input_dir = tmp_path / "cirifull"
    input_dir.mkdir()
    (input_dir / "cellA.tsv").write_text(
        "circ_id\tchr\tstart\tend\tstrand\thost_gene\thost_gene_n\tsupport\n"
        "circA\tchr1\t1\t5\t-\tGENE1\t1\t4\n"
        "circB\tchr2\t10\t20\t+\t\t0\t1\n",
        encoding="utf-8",
    )
    (input_dir / "cellB.tsv").write_text(
        "circ_id\tchr\tstart\tend\tstrand\thost_gene\thost_gene_n\tsupport\n"
        "circA\tchr1\t1\t5\t-\tGENE1\t1\t2\n",
        encoding="utf-8",
    )

    matrix_dir = tmp_path / "matrix"
    matrix_dir.mkdir()
    collect_matrix(
        str(input_dir),
        str(matrix_dir / "circ_counts.mtx"),
        str(matrix_dir / "circ_index.txt"),
        str(matrix_dir / "cell_index.txt"),
    )

    feature_path = matrix_dir / "circ_feature_table.tsv"
    feature_tsv = pd.read_csv(feature_path, sep="\t", keep_default_na=False)
    assert feature_tsv["host_gene"].tolist() == ["GENE1", ""]
    assert not any(str(value).isdigit() for value in feature_tsv["host_gene"].tolist())

    feature_df = load_circ_feature_table(["circA", "circB"], feature_path)
    X_cells_by_circ = sp.csr_matrix([[4, 1], [2, 0]], dtype=int)
    circ_qc = build_circ_qc_table(circ_ids=["circA", "circB"], feature_df=feature_df, X_cells_by_circ=X_cells_by_circ)
    assert circ_qc["host_gene"].tolist() == ["GENE1", ""]
    assert all(isinstance(value, str) for value in circ_qc["host_gene"].tolist())
    assert not any(str(value).isdigit() for value in circ_qc["host_gene"].tolist())
    assert circ_qc.loc["circA", "n_cells_detected"] == 2
    assert circ_qc.loc["circA", "host_gene"] == "GENE1"

    out_path = tmp_path / "circ_qc.tsv"
    circ_qc.reset_index().to_csv(out_path, sep="\t", index=False)
    written = pd.read_csv(out_path, sep="\t", keep_default_na=False)
    assert written["host_gene"].tolist() == ["GENE1", ""]
    assert not any(str(value).isdigit() for value in written["host_gene"].tolist())


def test_h5ad_export_shape_orientation_and_indices(tmp_path: Path) -> None:
    import anndata as ad

    cell_qc = pd.DataFrame(
        {
            "assigned_reads": [10, 6, 1],
            "alignment_status": ["aligned", "aligned", "missing"],
            "detector_status": ["success", "success", "empty"],
            "circRNA_count": [2, 1, 0],
            "total_circRNA_support": [6, 2, 0],
            "empty_flag": [False, False, True],
            "detector_seconds": [0.5, 0.6, 0.2],
            "alignment_seconds": [1.0, 1.5, None],
        },
        index=["cellA", "cellB", "cellC"],
    )
    circ_qc = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [1, 11],
            "end": [5, 15],
            "strand": ["+", "+"],
            "host_gene": ["GENE1", "GENE2"],
            "n_cells_detected": [2, 1],
            "total_support": [7, 1],
            "max_support": [5, 1],
            "mean_support_detected_cells": [3.5, 1.0],
        },
        index=["circA", "circB"],
    )
    X_cells_by_circ = sp.csr_matrix([[5, 1], [2, 0], [0, 0]], dtype=int)
    out_path = tmp_path / "circ_counts.h5ad"
    export_circ_h5ad(
        out_path=out_path,
        X_cells_by_circ=X_cells_by_circ,
        cell_qc=cell_qc,
        circ_qc=circ_qc,
        uns_payload={"workflow_name": "smartseq3-ciri3"},
    )

    adata = ad.read_h5ad(out_path)
    assert adata.shape == (3, 2)
    assert list(adata.obs_names) == ["cellA", "cellB", "cellC"]
    assert list(adata.var_names) == ["circA", "circB"]
    assert adata.X.toarray().tolist() == [[5, 1], [2, 0], [0, 0]]
    assert adata.uns["circyto"]["workflow_name"] == "smartseq3-ciri3"


def test_empty_matrix_behavior_in_h5ad_export(tmp_path: Path) -> None:
    import anndata as ad

    cell_qc = pd.DataFrame(
        {
            "assigned_reads": [4, 2],
            "alignment_status": ["aligned", "aligned"],
            "detector_status": ["empty", "empty"],
            "circRNA_count": [0, 0],
            "total_circRNA_support": [0, 0],
            "empty_flag": [True, True],
            "detector_seconds": [0.2, 0.3],
            "alignment_seconds": [0.9, 1.1],
        },
        index=["cellA", "cellB"],
    )
    circ_qc = pd.DataFrame(columns=["chrom", "start", "end", "strand", "host_gene", "n_cells_detected", "total_support", "max_support", "mean_support_detected_cells"]).rename_axis("circ_id")
    out_path = tmp_path / "empty.h5ad"
    export_circ_h5ad(
        out_path=out_path,
        X_cells_by_circ=sp.csr_matrix((2, 0), dtype=int),
        cell_qc=cell_qc,
        circ_qc=circ_qc,
        uns_payload={"workflow_name": "smartseq3-ciri3"},
    )
    adata = ad.read_h5ad(out_path)
    assert adata.shape == (2, 0)
    assert list(adata.obs_names) == ["cellA", "cellB"]


def test_mudata_export_from_synthetic_gene_counts_tsv(tmp_path: Path) -> None:
    if not HAS_MUDATA:
        pytest.skip("mudata not installed")

    import mudata as mu

    gene_counts = tmp_path / "gene_counts.tsv"
    _write_gene_counts_tsv(gene_counts)
    cell_qc = pd.DataFrame(
        {
            "assigned_reads": [10, 6],
            "alignment_status": ["aligned", "aligned"],
            "detector_status": ["success", "success"],
            "circRNA_count": [2, 1],
            "total_circRNA_support": [6, 2],
            "empty_flag": [False, False],
            "detector_seconds": [0.5, 0.6],
            "alignment_seconds": [1.0, 1.5],
        },
        index=["cellA", "cellB"],
    )
    circ_qc = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [1, 11],
            "end": [5, 15],
            "strand": ["+", "+"],
            "host_gene": ["GENE1", "GENE2"],
            "n_cells_detected": [2, 1],
            "total_support": [7, 1],
            "max_support": [5, 1],
            "mean_support_detected_cells": [3.5, 1.0],
        },
        index=["circA", "circB"],
    )
    out_path = tmp_path / "circyto_multimodal.h5mu"
    export_mudata_bundle(
        out_path=out_path,
        circ_X=sp.csr_matrix([[5, 1], [2, 0]], dtype=int),
        circ_qc=circ_qc,
        cell_qc=cell_qc,
        uns_payload={"workflow_name": "smartseq3-ciri3"},
        gene_counts=gene_counts,
        gene_counts_format="tsv",
        cell_join="inner",
    )

    mdata = mu.read_h5mu(out_path)
    assert "rna" in mdata.mod
    assert "circ" in mdata.mod
    assert mdata.mod["rna"].shape == (2, 2)
    assert mdata.mod["circ"].shape == (2, 2)
    assert list(mdata.obs_names) == ["cellA", "cellB"]


def test_mudata_requested_without_package_raises(tmp_path: Path) -> None:
    if HAS_MUDATA:
        pytest.skip("mudata installed")

    gene_counts = tmp_path / "gene_counts.tsv"
    _write_gene_counts_tsv(gene_counts)
    with pytest.raises(RuntimeError, match="mudata is not installed"):
        export_mudata_bundle(
            out_path=tmp_path / "circyto_multimodal.h5mu",
            circ_X=sp.csr_matrix([[1]], dtype=int),
            circ_qc=pd.DataFrame(
                {
                    "chrom": ["chr1"],
                    "start": [1],
                    "end": [2],
                    "strand": ["+"],
                    "host_gene": ["GENE1"],
                    "n_cells_detected": [1],
                    "total_support": [1],
                    "max_support": [1],
                    "mean_support_detected_cells": [1.0],
                },
                index=["circA"],
            ),
            cell_qc=pd.DataFrame(
                {
                    "assigned_reads": [1],
                    "alignment_status": ["aligned"],
                    "detector_status": ["success"],
                    "circRNA_count": [1],
                    "total_circRNA_support": [1],
                    "empty_flag": [False],
                    "detector_seconds": [0.1],
                    "alignment_seconds": [0.2],
                },
                index=["cellA"],
            ),
            uns_payload={"workflow_name": "smartseq3-ciri3"},
            gene_counts=gene_counts,
            gene_counts_format="tsv",
            cell_join="inner",
        )


def test_workflow_help_mentions_experimental() -> None:
    result = runner.invoke(app, ["workflow", "smartseq3-ciri3", "--help"])
    assert result.exit_code == 0
    assert "Experimental end-to-end SMART-Seq3 to CIRI3 workflow" in result.stdout
    assert "--star-index" in result.stdout
    assert "--resume" in result.stdout
    assert "--no-resume" in result.stdout
    assert "--export-h5ad" in result.stdout
    assert "--export-mudata" in result.stdout
