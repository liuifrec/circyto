from __future__ import annotations
from pathlib import Path
from typing import Optional, Iterable, Tuple
import subprocess
import pandas as pd
from joblib import Parallel, delayed

def _ensure_dir(p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)

def _run_one(cmd_template: str, ref_fa: Path, gtf: Path,
             r1: Path, r2: Optional[Path], out_tsv: Path) -> Tuple[str, int]:
    _ensure_dir(out_tsv)
    cmd = cmd_template.format(
        ref_fa=str(ref_fa), gtf=str(gtf), r1=str(r1), r2=str(r2 or ""), out_tsv=str(out_tsv)
    )
    cp = subprocess.run(cmd, shell=True, check=True)
    return (out_tsv.name, cp.returncode)

def _iter_fastq_batches(fastq_dir: Path) -> Iterable[Tuple[str, Path, Optional[Path]]]:
    for sub in sorted(p for p in Path(fastq_dir).iterdir() if p.is_dir()):
        r1 = sub / "R1.fastq.gz"
        r2 = sub / "R2.fastq.gz"
        if r1.exists():
            yield (sub.name, r1, r2 if r2.exists() else None)

def run_cirifull_over_fastqs(
    fastq_dir: Path, outdir: Path, cmd_template: str, ref_fa: Path, gtf: Path, threads: int = 8
) -> None:
    jobs = list(_iter_fastq_batches(Path(fastq_dir)))
    if not jobs:
        print(f"[circyto] No FASTQ batches under {fastq_dir}")
        return
    outdir = Path(outdir)
    print(f"[circyto] ciri-full jobs: {len(jobs)} → {outdir}")
    Parallel(n_jobs=threads)(
        delayed(_run_one)(
            cmd_template, Path(ref_fa), Path(gtf), r1, r2, outdir / f"{bid}.tsv"
        ) for (bid, r1, r2) in jobs
    )
    print("[circyto] run_cirifull_over_fastqs done.")

def _iter_manifest(manifest_tsv: Path) -> Iterable[Tuple[str, Path, Optional[Path]]]:
    df = pd.read_csv(manifest_tsv, sep="\t")
    if "cell_id" not in df.columns or "r1" not in df.columns:
        raise ValueError("Manifest must have columns: cell_id, r1, [r2]")
    for _, row in df.iterrows():
        cid = str(row["cell_id"])
        r1 = Path(str(row["r1"]))
        r2 = Path(str(row["r2"])) if "r2" in df.columns and pd.notna(row["r2"]) else None
        yield (cid, r1, r2)

def run_cirifull_with_manifest(
    manifest_tsv: Path, outdir: Path, cmd_template: str, ref_fa: Path, gtf: Path, threads: int = 8
) -> None:
    jobs = list(_iter_manifest(Path(manifest_tsv)))
    if not jobs:
        print(f"[circyto] No rows in manifest: {manifest_tsv}")
        return
    outdir = Path(outdir)
    print(f"[circyto] ciri-full(manifest) jobs: {len(jobs)} → {outdir}")
    Parallel(n_jobs=threads)(
        delayed(_run_one)(
            cmd_template, Path(ref_fa), Path(gtf), r1, r2, outdir / f"{cid}.tsv"
        ) for (cid, r1, r2) in jobs
    )
    print("[circyto] run_cirifull_with_manifest done.")
