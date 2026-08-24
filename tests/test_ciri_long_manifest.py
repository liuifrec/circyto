from __future__ import annotations

import csv
from pathlib import Path

import pytest

from circyto.manifest.ciri_long import (
    CIRI_LONG_INTERPRETATION_BOUNDARY,
    CIRI_LONG_SCHEMA_VERSION,
    CiriLongManifestRow,
    read_ciri_long_manifest_tsv,
    validate_ciri_long_manifest_tsv,
    write_ciri_long_manifest_tsv,
)


def _write_fastq(path: Path) -> None:
    path.write_text("@r1\nACGTACGT\n+\nIIIIIIII\n", encoding="utf-8")


def _row(reads: Path, **overrides: object) -> CiriLongManifestRow:
    values: dict[str, object] = {
        "sample_id": "sampleA",
        "reads_path": str(reads),
        "source_accession": "SYNTHETIC_RCRT",
        "dataset_id": "synthetic_ciri_long",
        "reference_id": "synthetic_reference",
        "reference_build": "synthetic_v1",
    }
    values.update(overrides)
    return CiriLongManifestRow(**values)  # type: ignore[arg-type]


def test_ciri_long_manifest_roundtrip_and_gate(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    _write_fastq(reads)
    manifest = tmp_path / "ciri_long_manifest.tsv"
    write_ciri_long_manifest_tsv([_row(reads)], manifest)

    rows = read_ciri_long_manifest_tsv(manifest, validate_files=True)
    assert len(rows) == 1
    assert rows[0].schema_version == CIRI_LONG_SCHEMA_VERSION
    assert rows[0].reads_path == str(reads.resolve())
    assert rows[0].biological_interpretation_boundary == (
        CIRI_LONG_INTERPRETATION_BOUNDARY
    )
    ok, errors, summary = validate_ciri_long_manifest_tsv(manifest, strict=True)
    assert ok is True
    assert errors == []
    assert summary["chemistry_gate"] == "accepted_rcrt_circrna_enriched"
    assert summary["reference_identities"] == [
        "synthetic_reference:synthetic_v1"
    ]


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"sequencing_platform": "ILLUMINA"}, "sequencing_platform"),
        ({"molecule_type": "direct_rna"}, "direct RNA"),
        ({"molecule_type": "unknown"}, "unknown or unsupported"),
        ({"library_preparation": "ordinary_cdna"}, "Ordinary ONT cDNA"),
        ({"circRNA_enrichment": False}, "circRNA_enrichment must be true"),
        ({"reference_id": ""}, "empty reference_id"),
        ({"reference_build": ""}, "empty reference_build"),
        (
            {"biological_interpretation_boundary": "single_cell_validation"},
            "biological_interpretation_boundary",
        ),
    ],
)
def test_ciri_long_manifest_rejects_incompatible_inputs(
    tmp_path: Path,
    overrides: dict[str, object],
    message: str,
) -> None:
    reads = tmp_path / "reads.fastq"
    _write_fastq(reads)
    with pytest.raises(ValueError, match=message):
        write_ciri_long_manifest_tsv(
            [_row(reads, **overrides)],
            tmp_path / "invalid.tsv",
        )


def test_ciri_long_manifest_rejects_alignment_and_polya_selection(
    tmp_path: Path,
) -> None:
    bam = tmp_path / "generic_minimap2.bam"
    bam.write_bytes(b"BAM")
    with pytest.raises(ValueError, match="not a generic alignment"):
        write_ciri_long_manifest_tsv(
            [_row(bam)],
            tmp_path / "bam.tsv",
        )

    reads = tmp_path / "reads.fastq"
    _write_fastq(reads)
    with pytest.raises(ValueError, match=r"poly\(A\)/oligo"):
        write_ciri_long_manifest_tsv(
            [
                _row(
                    reads,
                    extra={"library_selection": "poly(A) selected"},
                )
            ],
            tmp_path / "polya.tsv",
        )


def test_ciri_long_manifest_requires_source_and_unique_safe_samples(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "reads.fastq"
    _write_fastq(reads)
    with pytest.raises(ValueError, match="at least one of source_accession"):
        write_ciri_long_manifest_tsv(
            [_row(reads, source_accession="", dataset_id="")],
            tmp_path / "source.tsv",
        )
    with pytest.raises(ValueError, match="duplicate sample_id"):
        write_ciri_long_manifest_tsv(
            [_row(reads), _row(reads)],
            tmp_path / "duplicate.tsv",
        )
    with pytest.raises(ValueError, match="unsafe sample_id"):
        write_ciri_long_manifest_tsv(
            [_row(reads, sample_id="../sample")],
            tmp_path / "unsafe.tsv",
        )


def test_ciri_long_manifest_reports_false_enrichment_from_tsv(
    tmp_path: Path,
) -> None:
    reads = tmp_path / "reads.fastq"
    _write_fastq(reads)
    manifest = tmp_path / "manual.tsv"
    row = _row(reads).to_dict()
    row["circRNA_enrichment"] = "false"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    ok, errors, summary = validate_ciri_long_manifest_tsv(manifest, strict=True)
    assert ok is False
    assert "circRNA_enrichment must be true" in errors[0]
    assert summary["chemistry_gate"] == "rejected"
