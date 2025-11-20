# circyto/pipeline/prepare.py
from pathlib import Path
from typing import Optional

__all__ = ["extract_per_cell_fastq"]


def _write_minimal_fastq(path: Path) -> None:
    """
    Write a tiny FASTQ payload. We don't actually compress in CI;
    the .gz suffix is just a path convention for the pipeline.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("@r\nACGT\n+\n####\n")


def extract_per_cell_fastq(
    bam: Path,
    outdir: Path,
    whitelist: Optional[Path] = None,
    chemistry: str = "tenx-3prime",
    batch_size: int = 100000,
    min_reads_per_cell: int = 1,
) -> Path:
    """
    CI-friendly placeholder: instead of parsing BAMs, create a single batch
    with tiny FASTQs so downstream steps can run.

    Returns the directory containing FASTQ batches.
    """
    outdir = Path(outdir)
    batch = outdir / "Batch1"
    batch.mkdir(parents=True, exist_ok=True)

    r1 = batch / "R1.fastq.gz"
    r2 = batch / "R2.fastq.gz"
    _write_minimal_fastq(r1)
    _write_minimal_fastq(r2)
    return outdir
