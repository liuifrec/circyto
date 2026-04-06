from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from circyto.demux.smartseq2 import SmartSeq2DemuxParams, demux_smartseq2_pooled


def _write_fastq_gz(path: Path, seqs: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for i, seq in enumerate(seqs):
            handle.write(f"@read{i}\n")
            handle.write(f"{seq}\n")
            handle.write("+\n")
            handle.write(f"{'I' * len(seq)}\n")


def test_demux_smartseq2_rejects_mismatched_fastq_lengths(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    barcodes = tmp_path / "barcodes.tsv"

    _write_fastq_gz(r1, ["ACGTACGTAAAA"])
    _write_fastq_gz(r2, ["AACGTGATTTTT", "AACGTGATCCCC"])
    barcodes.write_text("sc01\tAACGTGAT\n", encoding="utf-8")

    params = SmartSeq2DemuxParams(
        r1=r1,
        r2=r2,
        barcodes_tsv=barcodes,
        outdir=tmp_path / "out",
        manifest_path=tmp_path / "out" / "manifest.tsv",
        library_id="LIB",
        barcode_read="R2",
        barcode_start=0,
        barcode_length=8,
    )

    with pytest.raises(ValueError) as excinfo:
        demux_smartseq2_pooled(params)

    assert "ended at different lengths" in str(excinfo.value)
