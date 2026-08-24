from __future__ import annotations

from pathlib import Path

import pytest

from circyto.manifest.long_read import (
    BIOLOGICAL_INTERPRETATION_BOUNDARY,
    LONG_READ_SCHEMA_VERSION,
    LongReadManifestRow,
    read_long_read_manifest_tsv,
    validate_long_read_manifest_tsv,
    write_long_read_manifest_tsv,
)


def _write_fastq(path: Path, read_name: str = "read1") -> None:
    path.write_text(f"@{read_name}\nACGT\n+\nIIII\n", encoding="utf-8")


def _row(cell_id: str, fastq: Path, *, molecule_type: str = "cdna") -> LongReadManifestRow:
    return LongReadManifestRow(
        cell_id=cell_id,
        long_read_fastq=str(fastq),
        protocol="synthetic_nanopore_cdna",
        sequencing_platform="OXFORD_NANOPORE",
        archive_library_selection="RACE",
        library_preparation_summary="Synthetic cDNA fixture; not a biological circRNA validation.",
        molecule_type=molecule_type,
        barcode_status="not_applicable_physical_single_cell",
        source_accession=f"SYNTHETIC_{cell_id}",
        dataset_id="synthetic_dataset",
    )


def test_long_read_manifest_one_cell_roundtrip(tmp_path: Path) -> None:
    fastq = tmp_path / "cell_a.fastq"
    manifest = tmp_path / "long_read_manifest.tsv"
    _write_fastq(fastq)

    write_long_read_manifest_tsv([_row("cell_a", fastq)], manifest)
    rows = read_long_read_manifest_tsv(manifest, validate_files=True)

    assert len(rows) == 1
    assert rows[0].schema_version == LONG_READ_SCHEMA_VERSION
    assert rows[0].long_read_fastq == str(fastq.resolve())
    assert rows[0].biological_interpretation_boundary == BIOLOGICAL_INTERPRETATION_BOUNDARY
    ok, errors, summary = validate_long_read_manifest_tsv(manifest, strict=True)
    assert ok is True
    assert errors == []
    assert summary == {
        "cells": 1,
        "cdna_rows": 1,
        "direct_rna_rows": 0,
        "missing_files": 0,
    }


def test_long_read_manifest_two_cells_preserves_distinct_inputs(tmp_path: Path) -> None:
    first = tmp_path / "cell_a.fastq"
    second = tmp_path / "cell_b.fastq"
    manifest = tmp_path / "long_read_manifest.tsv"
    _write_fastq(first, "a")
    _write_fastq(second, "b")

    write_long_read_manifest_tsv(
        [_row("cell_b", second), _row("cell_a", first)],
        manifest,
    )
    rows = read_long_read_manifest_tsv(manifest, validate_files=True)

    assert [row.cell_id for row in rows] == ["cell_a", "cell_b"]
    assert len({row.long_read_fastq for row in rows}) == 2


def test_long_read_manifest_rejects_duplicate_and_unsafe_cell_ids(tmp_path: Path) -> None:
    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq)
    with pytest.raises(ValueError, match="duplicate cell_id"):
        write_long_read_manifest_tsv(
            [_row("same", fastq), _row("same", fastq)],
            tmp_path / "duplicate.tsv",
        )
    with pytest.raises(ValueError, match="unsafe cell_id"):
        write_long_read_manifest_tsv(
            [_row("../collision", fastq)],
            tmp_path / "unsafe.tsv",
        )


def test_long_read_manifest_rejects_boundary_reinterpretation(tmp_path: Path) -> None:
    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq)
    row = _row("cell", fastq)
    invalid = LongReadManifestRow(
        **{
            **row.to_dict(),
            "biological_interpretation_boundary": "circRNA_validation",
        }
    )
    with pytest.raises(ValueError, match="biological_interpretation_boundary"):
        write_long_read_manifest_tsv([invalid], tmp_path / "invalid.tsv")
