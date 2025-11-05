from pathlib import Path
from typing import Optional
import typer
from rich.console import Console

from ..pipeline.prepare import extract_per_cell_fastq
from ..pipeline.run_cirifull import run_cirifull_over_fastqs, run_cirifull_with_manifest
from ..pipeline.collect import collect_matrix
from ..writers.convert import convert_matrix_files

app = typer.Typer(add_completion=False, help="circyto.py — CIRI-full wrapper for plate (manifest) and 10x (3'/5') scRNA-seq")
console = Console()

@app.command()
def prepare(
bam: Path = typer.Option(..., exists=True, help="10x-style BAM with CB/UB tags"),
outdir: Path = typer.Option(..., help="Output directory for batched FASTQs"),
whitelist: Optional[Path] = typer.Option(None, help="Optional whitelist of barcodes (gz/tsv)"),
chemistry: str = typer.Option("tenx-3p", help="tenx-3p or tenx-5p"),
batch_size: int = typer.Option(100, help="Cells per batch FASTQ pair"),
min_reads_per_cell: int = typer.Option(200, help="Discard low-read cells"),
):
    extract_per_cell_fastq(bam, outdir, whitelist, chemistry, batch_size, min_reads_per_cell)

@app.command()
def run(
fastq_dir: Path = typer.Option(..., exists=True, file_okay=False, help="Directory from prepare with FASTQ batches"),
outdir: Path = typer.Option(..., help="Output dir for CIRI-full results"),
cmd_template: str = typer.Option(..., help="CIRI-full command template, uses {ref_fa},{gtf},{r1},{r2},{out_tsv}"),
ref_fa: Path = typer.Option(..., exists=True, help="Reference FASTA"),
gtf: Path = typer.Option(..., exists=True, help="GTF/GFF"),
threads: int = typer.Option(8, help="Parallel workers"),
):
    run_cirifull_over_fastqs(fastq_dir, outdir, cmd_template, ref_fa, gtf, threads)

@app.command("run-manifest")
def run_manifest(
manifest: Path = typer.Option(..., exists=True, help="TSV with columns: cell_id, r1, [r2]"),
outdir: Path = typer.Option(..., help="Output dir for CIRI-full results"),
cmd_template: str = typer.Option(..., help="CIRI-full command template"),
ref_fa: Path = typer.Option(..., exists=True, help="Reference FASTA"),
gtf: Path = typer.Option(..., exists=True, help="GTF/GFF"),
threads: int = typer.Option(8, help="Parallel workers"),
):
    run_cirifull_with_manifest(manifest, outdir, cmd_template, ref_fa, gtf, threads)

@app.command()
def collect(
cirifull_dir: Path = typer.Option(..., exists=True, help="Directory with CIRI-full outputs"),
matrix: Path = typer.Option(..., help="Output sparse matrix .mtx"),
circ_index: Path = typer.Option(..., help="Output circ index (rows)"),
cell_index: Path = typer.Option(..., help="Output cell index (cols)"),
min_count_per_cell: int = typer.Option(1, help="Filter cells with total counts < threshold"),
):
    collect_matrix(cirifull_dir, matrix, circ_index, cell_index, min_count_per_cell)

@app.command()
def convert(
matrix: Path = typer.Option(..., exists=True, help="Sparse matrix .mtx (rows=circ, cols=cells)"),
circ_index: Path = typer.Option(..., exists=True, help="Text file of circ IDs (one per line)"),
cell_index: Path = typer.Option(..., exists=True, help="Text file of cell barcodes"),
loom: Optional[Path] = typer.Option(None, help="Path to write .loom"),
h5ad: Optional[Path] = typer.Option(None, help="Path to write .h5ad"),
):
    convert_matrix_files(matrix, circ_index, cell_index, loom, h5ad)

@app.command()
def make(
outdir: Path = typer.Option(..., help="Work output directory"),
cmd_template: str = typer.Option(..., help="CIRI-full command template"),
ref_fa: Path = typer.Option(..., exists=True, help="Reference FASTA"),
gtf: Path = typer.Option(..., exists=True, help="GTF/GFF"),
# plate/manifest path:
manifest: Optional[Path] = typer.Option(None, help="TSV listing cells and FASTQs (plate/full-length)"),
# 10x path:
bam: Optional[Path] = typer.Option(None, exists=True, help="10x-style BAM with CB/UB tags"),
whitelist: Optional[Path] = typer.Option(None, help="Optional barcode whitelist (gz/tsv)"),
chemistry: str = typer.Option("tenx-3p", help="tenx-3p or tenx-5p (only if using --bam)"),
threads: int = typer.Option(8, help="Parallel workers"),
):
    if manifest:
        run_cirifull_with_manifest(manifest, outdir / "cirifull_out", cmd_template, ref_fa, gtf, threads)
    elif bam:
        fq_dir = outdir / "fastq_by_cell"
        ciri_dir = outdir / "cirifull_out"
        extract_per_cell_fastq(bam, fq_dir, whitelist, chemistry, 100, 200)
        run_cirifull_over_fastqs(fq_dir, ciri_dir, cmd_template, ref_fa, gtf, threads)
    else:
        raise typer.BadParameter("Provide either --manifest (plate) or --bam (10x).")
# Collect + convert (h5ad by default)
    mat = outdir / "circ_counts.mtx"; circ_idx = outdir / "circ_index.txt"; cell_idx = outdir / "cell_index.txt"
    collect_matrix(outdir / "cirifull_out", mat, circ_idx, cell_idx, 1)
    convert_matrix_files(mat, circ_idx, cell_idx, h5ad=outdir / "circ.h5ad")
