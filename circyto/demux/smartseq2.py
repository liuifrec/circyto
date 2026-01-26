from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple
import gzip
import json


def _open_text_maybe_gz(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def load_barcodes_tsv(path: Path) -> Dict[str, str]:
    """
    Accepts either:
      1) one barcode per line -> cell_id = barcode
      2) cell_id<TAB>barcode
    Returns dict: barcode -> cell_id
    """
    m: Dict[str, str] = {}
    with _open_text_maybe_gz(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) == 1:
                bc = parts[0]
                cell_id = bc
            else:
                cell_id, bc = parts[0], parts[1]
            m[bc] = cell_id
    return m


@dataclass
class SmartSeq2DemuxParams:
    r1: Path
    r2: Path
    barcodes_tsv: Path
    outdir: Path
    manifest_path: Path
    library_id: str
    barcode_read: str  # "R1" or "R2"
    barcode_start: int
    barcode_length: int
    umi_start: Optional[int] = None
    umi_length: Optional[int] = None
    trim_barcode: bool = True
    trim_umi: bool = False
    max_mismatch: int = 0  # v1: only exact match
    limit_reads: Optional[int] = None
    overwrite: bool = False


def _fastq_iter(handle):
    while True:
        h = handle.readline()
        if not h:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            return
        yield h.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def _trim(seq: str, qual: str, cut_start: int, cut_len: int) -> Tuple[str, str]:
    if cut_len <= 0:
        return seq, qual
    left = cut_start
    right = cut_start + cut_len
    return seq[:left] + seq[right:], qual[:left] + qual[right:]


def demux_smartseq2_pooled(p: SmartSeq2DemuxParams) -> Dict[str, Dict[str, int]]:
    out_fastq = p.outdir / "fastq"
    sink_dir = p.outdir / "sink"
    out_fastq.mkdir(parents=True, exist_ok=True)
    sink_dir.mkdir(parents=True, exist_ok=True)

    if p.max_mismatch != 0:
        raise ValueError("v1 only supports exact barcode match (max_mismatch=0)")

    if not p.overwrite:
        if p.manifest_path.exists():
            raise FileExistsError(f"Refusing to overwrite existing manifest: {p.manifest_path}")
        if any(out_fastq.glob("*_R1.fastq*")) or any(out_fastq.glob("*_R2.fastq*")):
            raise FileExistsError(f"Refusing to overwrite existing FASTQs in: {out_fastq}")

    barcode_to_cell = load_barcodes_tsv(p.barcodes_tsv)

    writers: Dict[str, Tuple[gzip.GzipFile, gzip.GzipFile]] = {}
    stats: Dict[str, Dict[str, int]] = {}

    unknown_r1 = gzip.open(sink_dir / "unknown_barcode_R1.fastq.gz", "wt")
    unknown_r2 = gzip.open(sink_dir / "unknown_barcode_R2.fastq.gz", "wt")

    def get_writer(cell_id: str):
        if cell_id in writers:
            return writers[cell_id]
        r1_path = out_fastq / f"{cell_id}_R1.fastq.gz"
        r2_path = out_fastq / f"{cell_id}_R2.fastq.gz"
        w1 = gzip.open(r1_path, "wt")
        w2 = gzip.open(r2_path, "wt")
        writers[cell_id] = (w1, w2)
        stats.setdefault(cell_id, {"n_reads": 0})
        return w1, w2

    total_seen = 0
    assigned = 0
    unknown = 0

    with _open_text_maybe_gz(p.r1, "rt") as f1, _open_text_maybe_gz(p.r2, "rt") as f2:
        it1 = _fastq_iter(f1)
        it2 = _fastq_iter(f2)

        for rec1, rec2 in zip(it1, it2):
            total_seen += 1
            if p.limit_reads and total_seen > p.limit_reads:
                break

            h1, s1, pl1, q1 = rec1
            h2, s2, pl2, q2 = rec2

            seq_bc = s1 if p.barcode_read.upper() == "R1" else s2
            bc = seq_bc[p.barcode_start : p.barcode_start + p.barcode_length]
            cell_id = barcode_to_cell.get(bc)

            if p.trim_barcode:
                if p.barcode_read.upper() == "R1":
                    s1, q1 = _trim(s1, q1, p.barcode_start, p.barcode_length)
                else:
                    s2, q2 = _trim(s2, q2, p.barcode_start, p.barcode_length)

            if p.umi_start is not None and p.umi_length is not None and p.trim_umi:
                if p.barcode_read.upper() == "R1":
                    s1, q1 = _trim(s1, q1, p.umi_start, p.umi_length)
                else:
                    s2, q2 = _trim(s2, q2, p.umi_start, p.umi_length)

            if cell_id is None:
                unknown += 1
                unknown_r1.write(f"{h1}\n{s1}\n{pl1}\n{q1}\n")
                unknown_r2.write(f"{h2}\n{s2}\n{pl2}\n{q2}\n")
                continue

            w1, w2 = get_writer(cell_id)
            w1.write(f"{h1}\n{s1}\n{pl1}\n{q1}\n")
            w2.write(f"{h2}\n{s2}\n{pl2}\n{q2}\n")
            stats[cell_id]["n_reads"] += 1
            assigned += 1

    for w1, w2 in writers.values():
        w1.close()
        w2.close()
    unknown_r1.close()
    unknown_r2.close()

    report = {
        "platform": "smartseq2",
        "inputs": {"r1": str(p.r1), "r2": str(p.r2), "barcodes": str(p.barcodes_tsv)},
        "params": {
            "barcode_read": p.barcode_read,
            "barcode_start": p.barcode_start,
            "barcode_length": p.barcode_length,
            "umi_start": p.umi_start,
            "umi_length": p.umi_length,
            "trim_barcode": p.trim_barcode,
            "trim_umi": p.trim_umi,
            "max_mismatch": p.max_mismatch,
            "limit_reads": p.limit_reads,
        },
        "summary": {
            "total_reads_seen": total_seen,
            "assigned_reads": assigned,
            "unknown_barcode_reads": unknown,
            "cell_count": len(stats),
        },
        "outputs": {
            "fastq_dir": str((p.outdir / "fastq").resolve()),
            "sink_dir": str((p.outdir / "sink").resolve()),
            "manifest": str(p.manifest_path.resolve()),
        },
    }
    with (p.outdir / "demux_report.json").open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    return stats
