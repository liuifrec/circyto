from pathlib import Path
from typing import Optional
import typer
from rich.console import Console
from ..pipeline.prepare import extract_per_cell_fastq
from ..pipeline.run_cirifull import run_cirifull_over_fastqs
from ..pipeline.collect import collect_matrix
from ..writers.convert import convert_matrix_files
app = typer.Typer(add_completion=False, help="circyto.py — CIRI-full wrapper for scRNA-seq circRNA quantification")
console = Console()
@app.command()
def prepare(bam: Path = typer.Option(..., exists=True), outdir: Path = typer.Option(...), whitelist: Optional[Path] = typer.Option(None), batch_size: int = typer.Option(100), min_reads_per_cell: int = typer.Option(200)):
    extract_per_cell_fastq(bam, outdir, whitelist, batch_size, min_reads_per_cell)
@app.command()
def run(fastq_dir: Path = typer.Option(..., exists=True, file_okay=False), outdir: Path = typer.Option(...), cmd_template: str = typer.Option(...), ref_fa: Path = typer.Option(..., exists=True), gtf: Path = typer.Option(..., exists=True), threads: int = typer.Option(8)):
    run_cirifull_over_fastqs(fastq_dir, outdir, cmd_template, ref_fa, gtf, threads)
@app.command()
def collect(cirifull_dir: Path = typer.Option(..., exists=True), matrix: Path = typer.Option(...), circ_index: Path = typer.Option(...), cell_index: Path = typer.Option(...), min_count_per_cell: int = typer.Option(1)):
    collect_matrix(cirifull_dir, matrix, circ_index, cell_index, min_count_per_cell)
@app.command()
def convert(matrix: Path = typer.Option(..., exists=True), circ_index: Path = typer.Option(..., exists=True), cell_index: Path = typer.Option(..., exists=True), loom: Optional[Path] = typer.Option(None), h5ad: Optional[Path] = typer.Option(None)):
    convert_matrix_files(matrix, circ_index, cell_index, loom, h5ad)
@app.command()
def make(bam: Path = typer.Option(..., exists=True), whitelist: Optional[Path] = typer.Option(None), ref_fa: Path = typer.Option(..., exists=True), gtf: Path = typer.Option(..., exists=True), outdir: Path = typer.Option(...), threads: int = typer.Option(8), cmd_template: str = typer.Option(...)):
    fq_dir = outdir / "fastq_by_cell"; ciri_dir = outdir / "cirifull_out"; mat = outdir / "circ_counts.mtx"; circ_idx = outdir / "circ_index.txt"; cell_idx = outdir / "cell_index.txt"
    extract_per_cell_fastq(bam, fq_dir, whitelist, 100, 200)
    run_cirifull_over_fastqs(fq_dir, ciri_dir, cmd_template, ref_fa, gtf, threads)
    collect_matrix(ciri_dir, mat, circ_idx, cell_idx, 1)
    convert_matrix_files(mat, circ_idx, cell_idx, outdir / "circ.loom", outdir / "circ.h5ad")
