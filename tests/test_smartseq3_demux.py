from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.demux.smartseq3 import SmartSeq3DemuxParams, demux_smartseq3
from circyto.pipeline.align_manifest import read_source_manifest


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


def _read_fastq_headers(path: Path) -> list[str]:
    headers: list[str] = []
    with _open_maybe_gz(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            assert seq and plus and qual
            headers.append(header.strip())
    return headers


def _build_params(tmp_path: Path, *, max_records: int | None = None, write_sink: bool = True) -> SmartSeq3DemuxParams:
    read1 = tmp_path / "R1.fastq.gz"
    read2 = tmp_path / "R2.fastq.gz"
    index1 = tmp_path / "I1.fastq.gz"
    index2 = tmp_path / "I2.fastq.gz"
    annotation = tmp_path / "annotation.tsv"

    _write_fastq(read1, [("read1/1", "ACGT"), ("read2/1", "TGCA"), ("read3/1", "CCCC")])
    _write_fastq(read2, [("read1/2", "TTTT"), ("read2/2", "GGGG"), ("read3/2", "AAAA")])
    _write_fastq(index1, [("read1", "AAAA"), ("read2", "CCCC"), ("read3", "GGGG")])
    _write_fastq(index2, [("read1", "TTTT"), ("read2", "GGGG"), ("read3", "AAAA")])
    annotation.write_text(
        "cell_id\tindex1\tindex2\n"
        "cellA\tAAAA\tTTTT\n"
        "cellB\tCCCC\tGGGG\n",
        encoding="utf-8",
    )
    return SmartSeq3DemuxParams(
        read1=read1,
        read2=read2,
        index1=index1,
        index2=index2,
        annotation=annotation,
        outdir=tmp_path / "out",
        cell_id_column="cell_id",
        index1_column="index1",
        index2_column="index2",
        max_records=max_records,
        write_sink=write_sink,
    )


def test_smartseq3_demux_lockstep_writes_manifest_and_summary(tmp_path: Path) -> None:
    params = _build_params(tmp_path)

    summary = demux_smartseq3(params)

    assert summary["total_records"] == 3
    assert summary["assigned_records"] == 2
    assert summary["unmatched_records"] == 1
    assert summary["malformed_records"] == 0
    assert summary["number_of_cells_detected"] == 2
    assert summary["reads_per_cell"] == {"cellA": 1, "cellB": 1}

    r1_a = params.outdir / "fastq" / "cellA_R1.fastq.gz"
    r2_a = params.outdir / "fastq" / "cellA_R2.fastq.gz"
    r1_b = params.outdir / "fastq" / "cellB_R1.fastq.gz"
    r2_b = params.outdir / "fastq" / "cellB_R2.fastq.gz"
    assert r1_a.exists() and r2_a.exists() and r1_b.exists() and r2_b.exists()
    assert _read_fastq_headers(r1_a) == ["@read1/1"]
    assert _read_fastq_headers(r2_b) == ["@read2/2"]

    manifest = params.outdir / "manifest.tsv"
    manifest_text = manifest.read_text(encoding="utf-8")
    assert "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tread_layout" in manifest_text
    assert "cellA\tsmartseq3" in manifest_text
    assert "paired-end" in manifest_text

    source_rows = read_source_manifest(manifest, validate_files=True)
    assert [row.cell_id for row in source_rows] == ["cellA", "cellB"]
    assert all(row.read_layout == "paired-end" for row in source_rows)

    summary_json = json.loads((params.outdir / "demux_summary.json").read_text(encoding="utf-8"))
    assert summary_json["top_unmatched_index_pairs"] == [{"index_pair": "GGGG+AAAA", "count": 1}]


def test_smartseq3_demux_unmatched_sink_can_be_disabled(tmp_path: Path) -> None:
    params = _build_params(tmp_path, write_sink=False)

    demux_smartseq3(params)

    assert not (params.outdir / "sink").exists()


def test_smartseq3_demux_malformed_fastq_fails(tmp_path: Path) -> None:
    params = _build_params(tmp_path)
    with gzip.open(params.read1, "wt") as handle:
        handle.write("@read1\nACGT\n+\n")

    with pytest.raises(ValueError, match="Malformed FASTQ"):
        demux_smartseq3(params)


def test_smartseq3_demux_read_id_mismatch_fails_clearly(tmp_path: Path) -> None:
    params = _build_params(tmp_path)
    _write_fastq(params.index2, [("wrong1", "TTTT"), ("read2", "GGGG"), ("read3", "AAAA")])

    with pytest.raises(ValueError, match="FASTQ read ID mismatch"):
        demux_smartseq3(params)


def test_smartseq3_demux_max_records_limits_processing(tmp_path: Path) -> None:
    params = _build_params(tmp_path, max_records=2)

    summary = demux_smartseq3(params)

    assert summary["total_records"] == 2
    assert summary["assigned_records"] == 2
    assert summary["unmatched_records"] == 0
    assert not (params.outdir / "sink" / "unmatched_R1.fastq.gz").exists()


def test_smartseq3_demux_missing_annotation_columns_reports_available_columns(tmp_path: Path) -> None:
    params = _build_params(tmp_path)
    params.annotation.write_text("cell\ti7\ti5\ncellA\tAAAA\tTTTT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Available columns: cell, i7, i5"):
        demux_smartseq3(params)


def test_smartseq3_help_mentions_experimental() -> None:
    result = runner.invoke(app, ["demux", "smartseq3", "--help"])
    assert result.exit_code == 0
    assert "Experimental SMART-Seq3 pooled demultiplexing" in result.stdout
    assert "--cell-id-column" in result.stdout
