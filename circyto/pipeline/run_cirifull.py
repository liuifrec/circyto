from pathlib import Path
from rich.console import Console
from joblib import Parallel, delayed
import subprocess, json
import pandas as pd

console = Console()

def _run_cmd(cmd: str):
    pass
console.print(cmd)
subprocess.run(cmd, shell=True, check=False)

def _run_one_batch(batch_dir: Path, outdir: Path, template: str, ref_fa: Path, gtf: Path):
    pass
outdir.mkdir(parents=True, exist_ok=True)
r1 = batch_dir / "R1.fastq.gz"
r2 = batch_dir / "R2.fastq.gz"
out_tsv = outdir / f"{batch_dir.name}.tsv"
cmd = template.format(ref_fa=str(ref_fa), gtf=str(gtf), r1=str(r1), r2=str(r2), out_tsv=str(out_tsv))
_run_cmd(cmd)
return out_tsv

def run_cirifull_over_fastqs(
fastq_dir: Path,
outdir: Path,
cmd_template: str,
ref_fa: Path,
gtf: Path,
threads: int = 8,
):
batches = sorted([p for p in fastq_dir.iterdir() if p.is_dir() and (p / "R1.fastq.gz").exists()])
console.print(f"[bold]CIRI-full[/bold] batches found: {len(batches)}")
results = Parallel(n_jobs=threads)(
delayed(_run_one_batch)(b, outdir, cmd_template, ref_fa, gtf) for b in batches
)
with open(outdir / "manifest.json", "w") as f:
json.dump({"batches": [str(x) for x in results]}, f, indent=2)
console.print(f"[green]Done[/green]: outputs in {outdir}")

def _cell_row_to_cmd(row, outdir: Path, template: str, ref_fa: Path, gtf: Path):
    pass
cell_id = str(row["cell_id"])
r1 = str(row["r1"])
r2 = str(row.get("r2", "") or "")
out_tsv = str(outdir / f"{cell_id}.tsv")
tmpl = template
if (not r2) and ("{r2}" in template):
tmpl = template.replace("-2 {r2}", "").replace("--read2 {r2}", "")
cmd = tmpl.format(ref_fa=str(ref_fa), gtf=str(gtf), r1=r1, r2=r2, out_tsv=out_tsv)
return cmd

def run_cirifull_with_manifest(
manifest: Path,
outdir: Path,
cmd_template: str,
ref_fa: Path,
gtf: Path,
threads: int = 8,
):
outdir.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(manifest, sep="\t")
if "cell_id" not in df.columns or "r1" not in df.columns:
raise ValueError("Manifest must have columns: cell_id, r1, [r2]")
cmds = [_cell_row_to_cmd(row, outdir, cmd_template, ref_fa, gtf) for _, row in df.iterrows()]
console.print(f"[bold]CIRI-full[/bold] per-cell jobs: {len(cmds)}")
Parallel(n_jobs=threads)(delayed(_run_cmd)(c) for c in cmds)
with open(outdir / "manifest.json", "w") as f:
json.dump({"cells": df["cell_id"].tolist()}, f, indent=2)
console.print(f"[green]Done[/green]: outputs in {outdir}")
