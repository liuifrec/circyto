from pathlib import Path
from typing import Optional
import typer
from rich.console import Console
from typing import List

from ..pipeline.prepare import extract_per_cell_fastq
from ..pipeline.run_cirifull import run_cirifull_over_fastqs, run_cirifull_with_manifest
from ..pipeline.collect import collect_matrix
from ..writers.convert import convert_matrix_files
from ..pipeline.export_multimodal import export_multimodal as _export_multimodal
from ..pipeline.annotate_host_gene import annotate_host_genes
from circyto.detectors import available_detectors
from circyto.pipeline.run_detector import run_detector_manifest, run_multidetector
from circyto.compare import compare_detectors_from_root  # already added earlier



app = typer.Typer(
    add_completion=False,
    help="circyto.py — CIRI-full wrapper for plate (manifest) and 10x (3'/5') scRNA-seq",
)
console = Console()


@app.command()
def prepare(
    bam: Path = typer.Option(..., exists=True, help="10x-style BAM with CB/UB tags"),
    outdir: Path = typer.Option(..., help="Output directory for batched FASTQs"),
    whitelist: Optional[Path] = typer.Option(
        None, help="Optional whitelist of barcodes (gz/tsv)"
    ),
    chemistry: str = typer.Option("tenx-3p", help="tenx-3p or tenx-5p"),
    batch_size: int = typer.Option(100, help="Cells per batch FASTQ pair"),
    min_reads_per_cell: int = typer.Option(200, help="Discard low-read cells"),
):
    extract_per_cell_fastq(
        bam, outdir, whitelist, chemistry, batch_size, min_reads_per_cell
    )


@app.command()
def run(
    fastq_dir: Path = typer.Option(
        ...,
        exists=True,
        file_okay=False,
        help="Directory from prepare with FASTQ batches",
    ),
    outdir: Path = typer.Option(..., help="Output dir for CIRI-full results"),
    cmd_template: str = typer.Option(
        ..., help="CIRI-full command template, uses {ref_fa},{gtf},{r1},{r2},{out_tsv}"
    ),
    ref_fa: Path = typer.Option(..., exists=True, help="Reference FASTA"),
    gtf: Path = typer.Option(..., exists=True, help="GTF/GFF"),
    threads: int = typer.Option(8, help="Parallel workers"),
):
    run_cirifull_over_fastqs(fastq_dir, outdir, cmd_template, ref_fa, gtf, threads)


@app.command("run-manifest")
def run_manifest(
    manifest: Path = typer.Option(
        ..., exists=True, help="TSV with columns: cell_id, r1, [r2]"
    ),
    outdir: Path = typer.Option(..., help="Output dir for CIRI-full results"),
    cmd_template: str = typer.Option(..., help="CIRI-full command template"),
    ref_fa: Path = typer.Option(..., exists=True, help="Reference FASTA"),
    gtf: Path = typer.Option(..., exists=True, help="GTF/GFF"),
    threads: int = typer.Option(8, help="Parallel workers"),
):
    run_cirifull_with_manifest(manifest, outdir, cmd_template, ref_fa, gtf, threads)


@app.command()
def collect(
    cirifull_dir: Path = typer.Option(
        ..., exists=True, help="Directory with CIRI-full outputs"
    ),
    matrix: Path = typer.Option(..., help="Output sparse matrix .mtx"),
    circ_index: Path = typer.Option(..., help="Output circ index (rows)"),
    cell_index: Path = typer.Option(..., help="Output cell index (cols)"),
    min_count_per_cell: int = typer.Option(
        1, help="Filter cells with total counts < threshold"
    ),
):
    collect_matrix(cirifull_dir, matrix, circ_index, cell_index, min_count_per_cell)


@app.command()
def convert(
    matrix: Path = typer.Option(
        ..., exists=True, help="Sparse matrix .mtx (rows=circ, cols=cells)"
    ),
    circ_index: Path = typer.Option(
        ..., exists=True, help="Text file of circ IDs (one per line)"
    ),
    cell_index: Path = typer.Option(
        ..., exists=True, help="Text file of cell barcodes"
    ),
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
    manifest: Optional[Path] = typer.Option(
        None, help="TSV listing cells and FASTQs (plate/full-length)"
    ),
    # 10x path:
    bam: Optional[Path] = typer.Option(
        None, exists=True, help="10x-style BAM with CB/UB tags"
    ),
    whitelist: Optional[Path] = typer.Option(
        None, help="Optional barcode whitelist (gz/tsv)"
    ),
    chemistry: str = typer.Option(
        "tenx-3p", help="tenx-3p or tenx-5p (only if using --bam)"
    ),
    threads: int = typer.Option(8, help="Parallel workers"),
):
    if manifest:
        run_cirifull_with_manifest(
            manifest, outdir / "cirifull_out", cmd_template, ref_fa, gtf, threads
        )
    elif bam:
        fq_dir = outdir / "fastq_by_cell"
        ciri_dir = outdir / "cirifull_out"
        extract_per_cell_fastq(bam, fq_dir, whitelist, chemistry, 100, 200)
        run_cirifull_over_fastqs(fq_dir, ciri_dir, cmd_template, ref_fa, gtf, threads)
    else:
        raise typer.BadParameter("Provide either --manifest (plate) or --bam (10x).")
    # Collect + convert (h5ad by default)
    mat = outdir / "circ_counts.mtx"
    circ_idx = outdir / "circ_index.txt"
    cell_idx = outdir / "cell_index.txt"
    collect_matrix(outdir / "cirifull_out", mat, circ_idx, cell_idx, 1)
    convert_matrix_files(mat, circ_idx, cell_idx, h5ad=outdir / "circ.h5ad")
@app.command("export-multimodal")
def export_multimodal_cmd(
    genes_h5ad: Path = typer.Option(..., exists=True, help="Base gene expression .h5ad"),
    circ_matrix: Path = typer.Option(..., exists=True, help="circRNA MatrixMarket .mtx"),
    circ_index: Path = typer.Option(..., exists=True, help="circRNA index (rows)"),
    cell_index: Path = typer.Option(..., exists=True, help="Cell index (columns)"),
    out: Path = typer.Option(..., help="Output multimodal .h5ad"),
    circ_feature_table: Optional[Path] = typer.Option(
        None,
        help="Optional circ_feature_table.tsv with host gene annotations",
    ),
):
    """
    Attach circRNA counts as a separate modality (obsm['X_circ']) to an existing
    gene-expression AnnData. If a feature table is provided, include circRNA
    annotations and host-gene mapping.
    """
    _export_multimodal(
        genes_h5ad=genes_h5ad,
        circ_matrix=circ_matrix,
        circ_index=circ_index,
        cell_index=cell_index,
        out=out,
        circ_feature_table=circ_feature_table,
    )

@app.command("annotate-host-genes")
def annotate_host_genes_cmd(
    circ_feature_table: Path = typer.Option(
        ..., exists=True, help="circ_feature_table.tsv from circyto collect"
    ),
    gtf: Path = typer.Option(..., exists=True, help="Reference GTF used for circ calling"),
    out: Optional[Path] = typer.Option(
        None,
        help="Output feature table (TSV). If omitted, overwrites circ_feature_table.tsv",
    ),
    max_genes_per_circ: int = typer.Option(
        5, help="Maximum number of host genes to record per circRNA"
    ),
):
    """
    Annotate circRNAs with host gene(s) using a reference GTF.

    Adds columns:
      - host_gene
      - host_gene_id
      - host_genes_multi
      - host_gene_ids_multi
      - host_gene_n
    """
    annotate_host_genes(
        circ_feature_table=circ_feature_table,
        gtf_path=gtf,
        out=out,
        max_genes_per_circ=max_genes_per_circ,
    )
@ app.command()
def run_detector(
    detector: str = typer.Argument(..., help="Detector name (e.g., ciri-full)"),
    manifest: Path = typer.Option(..., help="Manifest TSV"),
    outdir: Path = typer.Option(..., help="Output directory"),
    ref_fa: Optional[Path] = typer.Option(None, help="Reference (FA)"),
    gtf: Optional[Path] = typer.Option(None, help="Annotation (GTF/GFF)"),
    threads: int = typer.Option(8),
    parallel: int = typer.Option(4),
):


    engines = available_detectors()
    if detector not in engines:
        raise typer.Exit(f"Detector {detector} not available. Available: {list(engines)}")

    det = engines[detector]
    print(f"[circyto] Running detector '{det.name}' (version={det.version()})")

    results = run_detector_manifest(
        detector=det,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
        parallel=parallel,
    )

    print(f"[circyto] Completed {len(results)} jobs.")





# ... existing @app.command() functions ...


@app.command()
def compare_detectors(
    root_dir: Path = typer.Option(..., help="Root directory containing per-detector subdirectories"),
    detectors: str = typer.Option(..., help="Comma-separated detector names (subdir names under root_dir)"),
    outdir: Path = typer.Option(..., help="Output directory for comparison results"),
):
    """
    Compare multiple detectors by their circRNA calls.

    Expects a layout like:

      root_dir/
        ciri-full/
          cell1.tsv
          cell2.tsv
        ciri-long/
          cell1.tsv
          cell2.tsv

    and will write circ × detector presence matrix, union/intersection lists,
    and a summary JSON into `outdir`.
    """
    det_list = [d.strip() for d in detectors.split(",") if d.strip()]
    if not det_list:
        raise typer.Exit("No detectors specified.")

    print(f"[circyto] Comparing detectors: {det_list}")
    summary = compare_detectors_from_root(root_dir, det_list, outdir)
    print("[circyto] Comparison complete.")
    print(f"[circyto] Detectors: {summary.get('n_detectors', 0)}")
    print(f"[circyto] Union circRNAs: {summary.get('union_size', 0)}")
    print(f"[circyto] Intersection circRNAs: {summary.get('intersection_size', 0)}")
@app.command()
def run_multidetector_cli(
    detectors: str = typer.Option(..., help="Comma-separated detector names (e.g. 'ciri-full,ciri-long')"),
    manifest: Path = typer.Option(..., help="Manifest TSV with cell_id, r1, r2"),
    root_outdir: Path = typer.Option(..., help="Root output directory (per-detector subdirs will be created)"),
    ref_fa: Path = typer.Option(..., help="Reference FASTA"),
    gtf: Path = typer.Option(..., help="Annotation GTF/GFF"),
    threads: int = typer.Option(8, help="Threads per detector"),
    parallel: int = typer.Option(4, help="Parallel jobs per detector"),
):
    """
    Run multiple detectors over the same manifest.

    Output layout:

      root_outdir/
        <detector_name>/
          <cell>.tsv
    """
    det_list = [d.strip() for d in detectors.split(",") if d.strip()]
    if not det_list:
        raise typer.Exit("No detectors specified.")

    engines_all = available_detectors()
    missing = [d for d in det_list if d not in engines_all]
    if missing:
        raise typer.Exit(f"Detectors not available: {missing}. Available: {list(engines_all)}")

    engines = {d: engines_all[d] for d in det_list}

    print(f"[circyto] Running multi-detector pipeline: {det_list}")
    run_multidetector(
        detectors=engines,
        manifest=manifest,
        root_outdir=root_outdir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
        parallel=parallel,
    )
    print("[circyto] Multi-detector run complete.")
@ app.command()
def run_multidetector(
    detectors: List[str] = typer.Argument(..., help="List of detectors, e.g. ciri-full ciri2"),
    manifest: Path = typer.Option(..., exists=True),
    outdir: Path = typer.Option(...),
    ref_fa: Path = typer.Option(None),
    gtf: Path = typer.Option(None),
    threads: int = typer.Option(8),
    parallel: int = typer.Option(1),
):
    """
    Run multiple detectors on the same manifest.

    Example:
        circyto run-multidetector ciri-full ciri2 \
            --manifest manifest.tsv \
            --outdir work/multi \
            --ref-fa ref.fa --gtf genes.gtf
    """
    from circyto.detectors import build_default_engines
    from circyto.pipeline.run_detector import run_detector_manifest
    from circyto.utils import ensure_dir

    engines = build_default_engines()
    for det in detectors:
        if det not in engines:
            raise ValueError(f"Detector '{det}' is not available. Available: {list(engines.keys())}")

    ensure_dir(outdir)

    metadata = {}

    for det in detectors:
        det_engine = engines[det]
        det_out = outdir / det
        ensure_dir(det_out)

        print(f"[multidetector] Running detector: {det}")
        results = run_detector_manifest(
            detector=det_engine,
            manifest=manifest,
            outdir=det_out,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
            parallel=parallel,
        )

        metadata[det] = {
            "n_cells": len(results),
            "detector": det,
            "outdir": str(det_out),
        }

    # write summary json
    import json
    summary_path = outdir / "summary.json"
    summary_path.write_text(json.dumps(metadata, indent=2))

    print(f"[multidetector] Completed. Summary at {summary_path}")
@ app.command()
def collect_multidetector(
    multi_out: Path = typer.Argument(..., help="Output dir from run-multidetector"),
):
    """
    Build per-detector circRNA sparse matrices.
    """
    from circyto.pipeline.multidetector_collect import (
        build_detector_matrix,
        write_matrix,
    )
    from circyto.detectors import build_default_engines
    from circyto.utils import ensure_dir

    ensure_dir(multi_out)
    det_dirs = [d for d in multi_out.iterdir() if d.is_dir()]

    for det_dir in det_dirs:
        tsv_dir = det_dir
        det_name = det_dir.name

        print(f"[collect-multidetector] Building matrix for {det_name}")

        X, circ_ids, cell_ids = build_detector_matrix(tsv_dir)

        matrices_dir = multi_out / "matrices"
        ensure_dir(matrices_dir)

        prefix = matrices_dir / det_name
        write_matrix(X, circ_ids, cell_ids, prefix)

        print(f"[collect-multidetector] WROTE: {prefix}.mtx")
