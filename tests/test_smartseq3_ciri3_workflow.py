from __future__ import annotations

import gzip
import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.manifest.alignment import AlignmentManifestRow, write_alignment_manifest_tsv
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
        json.dumps({"status_counts": {"aligned": 1}}, indent=2) + "\n",
        encoding="utf-8",
    )
    (ciri3_dir / "cellA.tsv").write_text(
        "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
        "circA\tchr1\t10\t20\t+\t3\n",
        encoding="utf-8",
    )
    (ciri3_dir / "detector_run_summary.json").write_text(
        json.dumps({"status_counts": {"success": 1}}, indent=2) + "\n",
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
            json.dumps({"status_counts": {"aligned": 1}}, indent=2) + "\n",
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
            json.dumps({"status_counts": {"success": 1}}, indent=2) + "\n",
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
    assert summary["selected_cell_count"] == 1
    assert summary["alignment_status_counts"] == {"aligned": 1}
    assert summary["detector_status_counts"] == {"success": 1}
    assert summary["matrix"] == {"n_rows": 2, "n_cols": 1, "nnz": 2}
    assert summary["command_options"]["top_n"] == 1
    assert summary["paths"]["selected_manifest"].endswith("top1_manifest.tsv")


def test_workflow_help_mentions_experimental() -> None:
    result = runner.invoke(app, ["workflow", "smartseq3-ciri3", "--help"])
    assert result.exit_code == 0
    assert "Experimental end-to-end SMART-Seq3 to CIRI3 workflow" in result.stdout
    assert "--star-index" in result.stdout
    assert "--resume" in result.stdout
    assert "--no-resume" in result.stdout
