from pathlib import Path
from typing import Optional, List
from joblib import Parallel, delayed
from rich.console import Console
import pandas as pd
import subprocess
import shlex

from ..detectors import get_engine

console = Console()

def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    return p

def _norm_path(x: str) -> str:
    return str(Path(x).resolve())

def _run_one(engine_name: str, r1: str, r2: Optional[str], ref_fa: str, gtf: str, out_tsv: str, threads: int):
    detector = get_engine(engine_name)
    detector(r1=_norm_path(r1),
             r2=_norm_path(r2) if r2 else "",
             ref_fa=_norm_path(ref_fa),
             gtf=_norm_path(gtf),
             out_tsv=_norm_path(out_tsv),
             threads=threads)

def run_detector(
    engine: str,
    manifest: str,
    ref_fa: str,
    gtf: str,
    outdir: str,
    threads: int = 8,
):
    """
    Run a detector for each row in a manifest TSV OR for a single BAM/FASTQ entry.

    Manifest format (tab-separated):
        cell_id    r1    [r2]

    Notes:
    - For BAM-based engines (e.g., CIRI-long), put the BAM path in 'r1' and leave r2 blank.
    - Output files are named as <cell_id>.tsv in outdir.
    """
    out_dir = _ensure_dir(Path(outdir))
    man = Path(manifest)
    if man.exists() and man.is_file():
        df = pd.read_csv(man, sep="\t")
        if "cell_id" not in df.columns or "r1" not in df.columns:
            raise ValueError("Manifest must have columns: cell_id, r1, [r2]")
        console.print(f"[bold]{engine}[/bold] jobs: {len(df)} → {out_dir}")
        Parallel(n_jobs=threads)(
            delayed(_run_one)(
                engine, row["r1"], row.get("r2", None),
                ref_fa, gtf,
                str(out_dir / f"{row['cell_id']}.tsv"),
                1  # each process runs its own internal threads if supported
            )
            for _, row in df.iterrows()
        )
    else:
        # Treat 'manifest' arg as a single file (BAM or FASTQ), write a single output
        cid = Path(manifest).stem
        console.print(f"[bold]{engine}[/bold] single run for: {manifest}")
        _run_one(engine, manifest, None, ref_fa, gtf, str(out_dir / f"{cid}.tsv"), threads)
    console.print(f"[green]Done[/green]: outputs in {out_dir}")
