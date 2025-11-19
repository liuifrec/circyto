# circyto/cli/circyto.py
from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import typer
from rich.console import Console

from ..pipeline.prepare import extract_per_cell_fastq
from ..pipeline.run_cirifull import (
    run_cirifull_over_fastqs,
    run_cirifull_with_manifest,
)
from ..pipeline.collect import collect_matrix
from ..writers.convert import convert_matrix_files
from ..pipeline.export_multimodal import export_multimodal as _export_multimodal
from ..pipeline.annotate_host_gene import annotate_host_genes
from ..pipeline.run_detector import run_detector_manifest
from ..pipeline.run_multidetector import run_multidetector_pipeline
from ..pipeline.merge_detectors import merge_detectors as _merge_detectors
from ..pipeline.compare_detectors import compare_detectors as _compare_detectors
from ..pipeline.multidetector_collect import (
    build_detector_matrix,
    write_matrix,
)
from ..utils import ensure_dir
from ..detectors import available_detectors, build_default_engines


app = typer.Typer(
    add_completion=False,
    help=(
        "circyto — CLI toolkit for single-cell circRNA detection and integration.\n\n"
        "Includes:\n"
        "  - legacy CIRI-full wrappers (prepare/run/run-manifest/collect/make)\n"
        "  - detector-API based runners (run-detector / run-multidetector)\n"
        "  - multi-detector merge/compare utilities\n"
        "  - multimodal export and host-gene annotation\n"
    ),
)
console = Console()


# --------------------------------------------------------------------------------------
# Core 10x / CIRI-full pipeline commands (legacy, still supported)
# --------------------------------------------------------------------------------------


@app.command()
def prepare(
    bam: Path = typer.Option(
        ...,
        exists=True,
        help="10x-style BAM with CB/UB tags",
    ),
    outdir: Path = typer.Option(
        ...,
        help="Output directory for per-cell FASTQs",
    ),
    whitelist: Optional[Path] = typer.Option(
        None,
        help="Optional whitelist of barcodes (gz/tsv)",
    ),
    chemistry: str = typer.Option(
        "tenx-3p",
        help="Library type: tenx-3p or tenx-5p",
    ),
    batch_size: int = typer.Option(
        100,
        help="Cells per FASTQ batch",
    ),
    min_reads_per_cell: int = typer.Option(
        200,
        help="Discard cells with < N reads",
    ),
) -> None:
    """
    Extract per-cell FASTQs from a 10x-style BAM.
    """
    extract_per_cell_fastq(
        bam=bam,
        outdir=outdir,
        whitelist=whitelist,
        chemistry=chemistry,
        batch_size=batch_size,
        min_reads_per_cell=min_reads_per_cell,
    )


@app.command()
def run(
    fastq_dir: Path = typer.Option(
        ...,
        exists=True,
        help="Directory of FASTQ batches (from `prepare`)",
    ),
    outdir: Path = typer.Option(
        ...,
        help="Output directory for CIRI-full results",
    ),
    cmd_template: str = typer.Option(
        ...,
        help=(
            "CIRI-full command template; uses "
            "{ref_fa}, {gtf}, {r1}, {r2}, {out_tsv}"
        ),
    ),
    ref_fa: Path = typer.Option(
        ...,
        exists=True,
        help="Reference FASTA",
    ),
    gtf: Path = typer.Option(
        ...,
        exists=True,
        help="GTF/GFF annotation",
    ),
    threads: int = typer.Option(
        8,
        help="Parallel workers",
    ),
) -> None:
    """
    Run CIRI-full over batched FASTQs.
    """
    run_cirifull_over_fastqs(
        fastq_dir=fastq_dir,
        outdir=outdir,
        cmd_template=cmd_template,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
    )


@app.command("run-manifest")
def run_manifest(
    manifest: Path = typer.Option(
        ...,
        exists=True,
        help="TSV with columns: cell_id, r1, [r2]",
    ),
    outdir: Path = typer.Option(
        ...,
        help="Output directory for CIRI-full results",
    ),
    cmd_template: str = typer.Option(
        ...,
        help="CIRI-full command template",
    ),
    ref_fa: Path = typer.Option(
        ...,
        exists=True,
        help="Reference FASTA",
    ),
    gtf: Path = typer.Option(
        ...,
        exists=True,
        help="GTF/GFF annotation",
    ),
    threads: int = typer.Option(
        8,
        help="Parallel workers",
    ),
) -> None:
    """
    Run CIRI-full over a plate/full-length manifest.
    """
    run_cirifull_with_manifest(
        manifest=manifest,
        outdir=outdir,
        cmd_template=cmd_template,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
    )


@app.command()
def collect(
    cirifull_dir: Path = typer.Option(
        ...,
        exists=True,
        help="Directory with CIRI-full per-cell outputs",
    ),
    matrix: Path = typer.Option(
        ...,
        help="Output sparse matrix (.mtx)",
    ),
    circ_index: Path = typer.Option(
        ...,
        help="Output circ index (rows)",
    ),
    cell_index: Path = typer.Option(
        ...,
        help="Output cell index (columns)",
    ),
    min_count_per_cell: int = typer.Option(
        1,
        help="Drop cells with total circ counts < threshold",
    ),
) -> None:
    """
    Collect CIRI-full per-cell TSVs into a circ × cell MatrixMarket matrix.
    """
    collect_matrix(
        cirifull_dir=cirifull_dir,
        matrix=matrix,
        circ_index=circ_index,
        cell_index=cell_index,
        min_count_per_cell=min_count_per_cell,
    )


@app.command()
def convert(
    matrix: Path = typer.Option(
        ...,
        exists=True,
        help="Sparse matrix .mtx (rows=circ, cols=cells)",
    ),
    circ_index: Path = typer.Option(
        ...,
        exists=True,
        help="Text file of circ IDs (one per line)",
    ),
    cell_index: Path = typer.Option(
        ...,
        exists=True,
        help="Text file of cell barcodes",
    ),
    loom: Optional[Path] = typer.Option(
        None,
        help="Optional path to write .loom",
    ),
    h5ad: Optional[Path] = typer.Option(
        None,
        help="Optional path to write .h5ad",
    ),
) -> None:
    """
    Convert circ × cell matrix and index files to loom/h5ad.
    """
    convert_matrix_files(
        matrix=matrix,
        circ_index=circ_index,
        cell_index=cell_index,
        loom=loom,
        h5ad=h5ad,
    )


@app.command()
def make(
    outdir: Path = typer.Option(
        ...,
        help="Work output directory",
    ),
    cmd_template: str = typer.Option(
        ...,
        help="CIRI-full command template",
    ),
    ref_fa: Path = typer.Option(
        ...,
        exists=True,
        help="Reference FASTA",
    ),
    gtf: Path = typer.Option(
        ...,
        exists=True,
        help="GTF/GFF annotation",
    ),
    manifest: Optional[Path] = typer.Option(
        None,
        help="Plate-style TSV listing cells and FASTQs",
    ),
    bam: Optional[Path] = typer.Option(
        None,
        exists=True,
        help="10x-style BAM with CB/UB tags",
    ),
    whitelist: Optional[Path] = typer.Option(
        None,
        help="Optional barcode whitelist (gz/tsv) for 10x path",
    ),
    chemistry: str = typer.Option(
        "tenx-3p",
        help="tenx-3p or tenx-5p (for 10x path)",
    ),
    threads: int = typer.Option(
        8,
        help="Parallel workers",
    ),
) -> None:
    """
    Convenience wrapper:

      - Plate path: --manifest
      - 10x path:   --bam (+ whitelist / chemistry)

    Then runs collect + convert to produce circ_counts.mtx and circ.h5ad.
    """
    if manifest:
        ciri_dir = outdir / "cirifull_out"
        run_cirifull_with_manifest(
            manifest=manifest,
            outdir=ciri_dir,
            cmd_template=cmd_template,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
        )
    elif bam:
        fq_dir = outdir / "fastq_by_cell"
        ciri_dir = outdir / "cirifull_out"
        extract_per_cell_fastq(
            bam=bam,
            outdir=fq_dir,
            whitelist=whitelist,
            chemistry=chemistry,
            batch_size=100,
            min_reads_per_cell=200,
        )
        run_cirifull_over_fastqs(
            fastq_dir=fq_dir,
            outdir=ciri_dir,
            cmd_template=cmd_template,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
        )
    else:
        raise typer.BadParameter("Provide either --manifest (plate) or --bam (10x).")

    mat = outdir / "circ_counts.mtx"
    circ_idx = outdir / "circ_index.txt"
    cell_idx = outdir / "cell_index.txt"

    collect_matrix(
        cirifull_dir=outdir / "cirifull_out",
        matrix=mat,
        circ_index=circ_idx,
        cell_index=cell_idx,
        min_count_per_cell=1,
    )
    convert_matrix_files(
        matrix=mat,
        circ_index=circ_idx,
        cell_index=cell_idx,
        h5ad=outdir / "circ.h5ad",
    )


# --------------------------------------------------------------------------------------
# Multimodal export and host-gene annotation
# --------------------------------------------------------------------------------------


@app.command("export-multimodal")
def export_multimodal_cmd(
    genes_h5ad: Path = typer.Option(
        ...,
        exists=True,
        help="Base gene-expression .h5ad",
    ),
    circ_matrix: Path = typer.Option(
        ...,
        exists=True,
        help="circRNA MatrixMarket .mtx",
    ),
    circ_index: Path = typer.Option(
        ...,
        exists=True,
        help="circRNA index (rows)",
    ),
    cell_index: Path = typer.Option(
        ...,
        exists=True,
        help="Cell index (columns)",
    ),
    out: Path = typer.Option(
        ...,
        help="Output multimodal .h5ad",
    ),
    circ_feature_table: Optional[Path] = typer.Option(
        None,
        help="Optional circ_feature_table.tsv with host-gene annotations",
    ),
) -> None:
    """
    Attach circRNA counts as obsm['X_circ'] to an existing gene-expression AnnData.
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
        ...,
        exists=True,
        help="circ_feature_table.tsv from `circyto collect`",
    ),
    gtf: Path = typer.Option(
        ...,
        exists=True,
        help="Reference GTF used for circ calling",
    ),
    out: Optional[Path] = typer.Option(
        None,
        help="Output TSV (default: overwrite circ_feature_table)",
    ),
    max_genes_per_circ: int = typer.Option(
        5,
        help="Maximum number of host genes to record per circRNA",
    ),
) -> None:
    """
    Annotate circRNAs with host gene(s) using a reference GTF.
    """
    annotate_host_genes(
        circ_feature_table=circ_feature_table,
        gtf_path=gtf,
        out=out,
        max_genes_per_circ=max_genes_per_circ,
    )


# --------------------------------------------------------------------------------------
# Detector API: single-detector and multi-detector runners
# --------------------------------------------------------------------------------------


@app.command("run-detector")
def run_detector_cmd(
    detector: str = typer.Argument(
        ...,
        help="Detector name (e.g. 'ciri-full', 'ciri2')",
    ),
    manifest: Path = typer.Option(
        ...,
        exists=True,
        help="Manifest TSV (columns: cell_id, r1, [r2])",
    ),
    outdir: Path = typer.Option(
        ...,
        help="Output directory for per-cell TSVs",
    ),
    ref_fa: Optional[Path] = typer.Option(
        None,
        help="Reference FASTA (required for most detectors)",
    ),
    gtf: Optional[Path] = typer.Option(
        None,
        help="Annotation GTF/GFF (required for most detectors)",
    ),
    threads: int = typer.Option(
        8,
        help="Threads per detector process",
    ),
    parallel: int = typer.Option(
        4,
        help="Number of cells to run in parallel",
    ),
) -> None:
    """
    Run a single detector over a manifest using the detector API.
    """
    engines = build_default_engines()
    if detector not in engines:
        available = ", ".join(sorted(engines.keys()))
        raise typer.Exit(f"Detector '{detector}' not available. Available: {available}")

    det_engine = engines[detector]
    console.print(
        f"[bold cyan][circyto][/bold cyan] Running detector '{det_engine.name}' "
        f"(version={det_engine.version()})"
    )

    outdir.mkdir(parents=True, exist_ok=True)

    results = run_detector_manifest(
        detector=det_engine,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
        parallel=parallel,
    )

    console.print(
        f"[bold cyan][circyto][/bold cyan] Completed {len(results)} jobs into {outdir}"
    )

@app.command("run-multidetector")
def run_multidetector_cmd(
    detectors: List[str] = typer.Argument(..., help="List of detectors to run"),
    manifest: Path = typer.Option(..., exists=True, dir_okay=False),
    outdir: Path = typer.Argument(...),
    ref_fa: Path = typer.Option(None),
    gtf: Path = typer.Option(None),
    threads: int = 8,
    parallel: int = 1,
):
    """
    Run multiple detectors on the same manifest.
    """
    from circyto.pipeline.run_multidetector import run_multidetector_pipeline

    outdir.mkdir(parents=True, exist_ok=True)

    result = run_multidetector_pipeline(
        detectors=detectors,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
        parallel=parallel,
    )

    typer.echo(f"[run-multidetector] Completed. Summary at {outdir/'summary.json'}")


# --------------------------------------------------------------------------------------
# Multi-detector: merge, collect matrices, compare
# --------------------------------------------------------------------------------------


@app.command("merge-detectors")
def merge_detectors_cmd(
    indir: Path = typer.Argument(
        ...,
        exists=True,
        help="Output dir from `run-multidetector` (must contain summary.json)",
    ),
    outdir: Path = typer.Argument(
        ...,
        help="Output directory for merged TSVs (union/long) and metadata",
    ),
) -> None:
    """
    Merge per-detector per-cell TSVs into union and long-format tables.

    Produces:

      - circ_union.tsv        (one row per circRNA, per-detector support columns)
      - circ_by_detector.tsv  (long-form circ × detector × cell)
      - metadata.json         (basic bookkeeping)
    """
    outdir.mkdir(parents=True, exist_ok=True)
    res = _merge_detectors(indir=indir, outdir=outdir)

    msg = "[merge-detectors] finished."
    if isinstance(res, dict):
        union = res.get("union_path")
        long = res.get("long_path")
        meta = res.get("metadata_path")
        if union:
            msg += f" union={union}"
        if long:
            msg += f" long={long}"
        if meta:
            msg += f" meta={meta}"
    console.print(msg)


@app.command("collect-multidetector")
def collect_multidetector_cmd(
    multi_out: Path = typer.Argument(
        ...,
        exists=True,
        help="Output dir from `run-multidetector` (per-detector subdirs with TSVs)",
    ),
) -> None:
    """
    Build per-detector circRNA sparse matrices from multi-detector TSVs.

    For each detector subdir under multi_out, writes:

      multi_out/matrices/<detector>.mtx
      multi_out/matrices/<detector>.circ.txt
      multi_out/matrices/<detector>.cell.txt
    """
    ensure_dir(multi_out)
    matrices_dir = multi_out / "matrices"
    ensure_dir(matrices_dir)

    det_dirs = [d for d in multi_out.iterdir() if d.is_dir()]

    for det_dir in sorted(det_dirs):
        det_name = det_dir.name
        console.print(
            f"[bold cyan][collect-multidetector][/bold cyan] Building matrix for {det_name}"
        )

        X, circ_ids, cell_ids = build_detector_matrix(det_dir)
        prefix = matrices_dir / det_name
        write_matrix(X, circ_ids, cell_ids, prefix)

        console.print(
            f"[bold cyan][collect-multidetector][/bold cyan] Wrote {prefix}.mtx"
        )


@app.command("compare-detectors")
def compare_detectors_cmd(
    indir: Path = typer.Argument(
        ...,
        exists=True,
        help="Directory containing circ_union.tsv (output of `merge-detectors`)",
    ),
    outdir: Path = typer.Argument(
        ...,
        help="Output directory for Jaccard matrix, summary, metadata",
    ),
) -> None:
    """
    Compare detectors using a merged circ_union.tsv.

    Produces:

      - jaccard.tsv           (detector × detector Jaccard similarity)
      - detector_summary.tsv  (per-detector stats)
      - compare_metadata.json
    """
    outdir.mkdir(parents=True, exist_ok=True)
    res = _compare_detectors(indir=indir, outdir=outdir)

    msg = "[compare-detectors] finished."
    if isinstance(res, dict):
        jac = res.get("jaccard_path")
        summ = res.get("summary_path")
        meta = res.get("metadata_path")
        if jac:
            msg += f" jaccard={jac}"
        if summ:
            msg += f" summary={summ}"
        if meta:
            msg += f" meta={meta}"
    console.print(msg)
