from pathlib import Path
from typing import Optional
from rich.console import Console
import gzip, pysam, json
from collections import defaultdict

console = Console()

def _read_whitelist(path: Optional[Path]):
if path is None:
return None
op = gzip.open if str(path).endswith(".gz") else open
with op(path, "rt") as f:
return {line.strip().split()[0] for line in f if line.strip()}

def extract_per_cell_fastq(
bam: Path,
outdir: Path,
whitelist: Optional[Path],
chemistry: str,
batch_size: int,
min_reads_per_cell: int,
):
"""
BAM → batched FASTQs. Heuristic reconstruction:
- tenx-3p: write aligned read seq to R1, dummy R2
- tenx-5p: write aligned read seq to R1 (assuming alignment from cDNA mate), dummy R2
Real datasets should prefer manifest mode with original R1/R2 FASTQs.
"""
chemistry = chemistry.lower()
if chemistry not in ("tenx-3p","tenx-5p","tenx-3","tenx-5"):
raise ValueError("chemistry must be tenx-3p or tenx-5p")

outdir.mkdir(parents=True, exist_ok=True)
wl = _read_whitelist(whitelist)
console.print(f"[bold]Reading BAM[/bold]: {bam}  [italic](chemistry={chemistry})[/italic]")
bamf = pysam.AlignmentFile(bam, "rb")
cell_counts = defaultdict(int)

def new_batch_writers(bid: int):
    d = outdir / f"batch_{bid:05d}"
    d.mkdir(parents=True, exist_ok=True)
    return (
        gzip.open(d / "R1.fastq.gz", "wt"),
        gzip.open(d / "R2.fastq.gz", "wt"),
        open(d / "cells.txt", "w"),
        d,
    )

batch_id = 0
r1, r2, cf, batch_dir = new_batch_writers(batch_id)
n_written = 0
seen_cells = set()

for r in bamf.fetch(until_eof=True):
    if r.is_unmapped:
        continue
    tags = dict(r.tags)
    cb = tags.get("CB", None)
    if cb is None:
        continue
    if wl is not None and cb not in wl:
        continue
    cell_counts[cb] += 1
    seen_cells.add(cb)

    rname = r.query_name
    seq = r.query_sequence or "N"
    qual = "".join(chr((q if q is not None else 30)+33) for q in (r.query_qualities or [30]*len(seq)))

    # For both 3' and 5' heuristics, place aligned sequence in R1.
    r1.write(f"@{rname}\n{seq}\n+\n{qual}\n")
    r2.write(f"@{rname}\n{'N'*len(seq)}\n+\n{'I'*len(seq)}\n")
    n_written += 1

    if n_written >= batch_size * 1000:
        r1.close(); r2.close()
        with cf:
            cf.write("\n".join(sorted(seen_cells)))
        batch_id += 1
        r1, r2, cf, batch_dir = new_batch_writers(batch_id)
        n_written = 0
        seen_cells = set()

r1.close(); r2.close()
with cf:
    cf.write("\n".join(sorted(seen_cells)))

qc = outdir / "prepare_qc.json"
with open(qc, "w") as f:
    json.dump({"cells": len(cell_counts), "reads_by_cell": dict(cell_counts)}, f, indent=2)
console.print(f"[green]Prepared[/green] batches in {outdir}. Cells observed: {len(cell_counts)}")


