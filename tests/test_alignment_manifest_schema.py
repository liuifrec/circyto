from __future__ import annotations

from pathlib import Path

import pytest

from circyto.manifest.alignment import read_alignment_manifest_tsv, write_alignment_manifest_tsv, AlignmentManifestRow


def test_alignment_manifest_read_normalizes_paths_and_requires_layout(tmp_path: Path) -> None:
    sam = tmp_path / "alignments" / "sample.sam"
    sam.parent.mkdir(parents=True, exist_ok=True)
    sam.write_text("@SQ\tSN:chr1\tLN:4\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    source_manifest = tmp_path / "source.tsv"
    source_manifest.write_text("cell_id\tread1\n", encoding="utf-8")
    manifest = tmp_path / "alignment_manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tchimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\tmapper_mode\tartifact_bucket\tsortedness",
                f"cell1\t\talignments/sample.sam\tgrp\tsingle-end\tbwa-mem\t{ref}\tcache1\t{source_manifest}\t\t\t\t\t0\tbwa_mem\tunsorted",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = read_alignment_manifest_tsv(manifest, validate_files=True)
    assert len(rows) == 1
    row = rows[0]
    assert row.read_layout == "single-end"
    assert row.sam == str(sam.resolve())
    assert row.reference == str(ref.resolve())
    assert row.source_manifest == str(source_manifest.resolve())


def test_alignment_manifest_read_rejects_missing_layout(tmp_path: Path) -> None:
    sam = tmp_path / "sample.sam"
    sam.write_text("@SQ\tSN:chr1\tLN:4\n", encoding="utf-8")
    manifest = tmp_path / "alignment_manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tchimeric_junction\tunmapped_mate1\tunmapped_mate2\tbwa_sam\tmapper_mode\tartifact_bucket\tsortedness",
                f"cell1\t\t{sam}\tgrp\t\tbwa-mem\t\tcache1\t\t\t\t\t\t0\tbwa_mem\tunsorted",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing required read_layout"):
        read_alignment_manifest_tsv(manifest, validate_files=False)


def test_alignment_manifest_write_normalizes_paths(tmp_path: Path) -> None:
    sam = tmp_path / "sample.sam"
    sam.write_text("@SQ\tSN:chr1\tLN:4\n", encoding="utf-8")
    manifest = tmp_path / "alignment_manifest.tsv"
    write_alignment_manifest_tsv(
        [
            AlignmentManifestRow(
                cell_id="cell1",
                sam=str(sam),
                group_id="grp",
                read_layout="paired-end",
                aligner="bwa-mem",
                reference=str(tmp_path / "ref.fa"),
                cache_key="cache1",
                source_manifest=str(tmp_path / "source.tsv"),
                mapper_mode="0",
                artifact_bucket="bwa_mem",
                sortedness="unsorted",
            )
        ],
        manifest,
    )
    text = manifest.read_text(encoding="utf-8")
    assert str(sam.resolve()) in text
