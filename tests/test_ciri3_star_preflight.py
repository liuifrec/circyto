from __future__ import annotations

from pathlib import Path

import pytest

from circyto.pipeline.align_manifest import plan_alignment_cache, prepare_alignment_cache


def _write_manifest(tmp_path: Path, *, read1: str, read2: str = "", read_layout: str = "") -> Path:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id\tread_layout",
                f"cell1\tsmartseq2\t{read1}\t{read2}\t\tlib1\t10\tgrp\t{read_layout}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


def test_plan_alignment_cache_rejects_single_end_star_ciri3(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read1 = tmp_path / "reads.fastq"
    read1.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    manifest = _write_manifest(tmp_path, read1=str(read1), read_layout="single-end")

    monkeypatch.setattr("circyto.pipeline.align_manifest._tool_exists", lambda name: True)

    plan = plan_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="star",
        detector_hint="ciri3",
        ref_fa=ref,
        extra_flags="--genomeDir /tmp/star-index",
    )

    message = " ".join(plan["errors"])
    assert "STAR + CIRI3 is currently supported only for paired-end data" in message
    assert "Single-end STAR + CIRI3 is not implemented in the current hybrid rescue path" in message
    assert "Use BWA + CIRI3" in message


def test_prepare_alignment_cache_allows_paired_end_star_ciri3(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read1 = tmp_path / "reads_1.fastq"
    read2 = tmp_path / "reads_2.fastq"
    read1.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    read2.write_text("@r2\nTGCA\n+\n!!!!\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    manifest = _write_manifest(tmp_path, read1=str(read1), read2=str(read2), read_layout="paired-end")

    monkeypatch.setattr("circyto.pipeline.align_manifest._tool_exists", lambda name: True)

    plan = plan_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="star",
        detector_hint="ciri3",
        ref_fa=ref,
        extra_flags="--genomeDir /tmp/star-index",
    )

    assert not plan["errors"]


def test_prepare_alignment_cache_rejects_gz_star_without_readfilescommand(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read1 = tmp_path / "reads.fastq.gz"
    read2 = tmp_path / "reads_2.fastq.gz"
    read1.write_text("gz-placeholder", encoding="utf-8")
    read2.write_text("gz-placeholder", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    manifest = _write_manifest(tmp_path, read1=str(read1), read2=str(read2), read_layout="paired-end")

    monkeypatch.setattr("circyto.pipeline.align_manifest._tool_exists", lambda name: True)

    with pytest.raises(RuntimeError, match="--readFilesCommand zcat"):
        prepare_alignment_cache(
            manifest=manifest,
            outdir=tmp_path / "align",
            aligner="star",
            detector_hint="ciri3",
            ref_fa=ref,
            extra_flags="--genomeDir /tmp/star-index",
        )
