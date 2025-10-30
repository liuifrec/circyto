from pathlib import Path
from rich.console import Console
from joblib import Parallel, delayed
import subprocess, json
console = Console()
def _run_one(batch_dir: Path, outdir: Path, template: str, ref_fa: Path, gtf: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    r1 = batch_dir / 'R1.fastq.gz'; r2 = batch_dir / 'R2.fastq.gz'; out_tsv = outdir / f'{batch_dir.name}.tsv'
    cmd = template.format(ref_fa=str(ref_fa), gtf=str(gtf), r1=str(r1), r2=str(r2), out_tsv=str(out_tsv))
    console.print(f'[bold]CIRI-full[/bold] on {batch_dir.name}: {cmd}')
    subprocess.run(cmd, shell=True, check=False)
    return out_tsv
def run_cirifull_over_fastqs(fastq_dir: Path, outdir: Path, cmd_template: str, ref_fa: Path, gtf: Path, threads: int = 8):
    batches = sorted([p for p in fastq_dir.iterdir() if p.is_dir() and (p / 'R1.fastq.gz').exists()])
    console.print(f'Found {len(batches)} batches.')
    results = Parallel(n_jobs=threads)(delayed(_run_one)(b, outdir, cmd_template, ref_fa, gtf) for b in batches)
    with open(outdir / 'manifest.json', 'w') as f: json.dump({'batches': [str(x) for x in results]}, f, indent=2)
    console.print(f'[green]Done[/green]: outputs in {outdir}')
