from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import typer
from rich.console import Console

from circyto.pipeline.prepare import extract_per_cell_fastq
from circyto.pipeline.run_cirifull import (
    run_cirifull_over_fastqs,
    run_cirifull_with_manifest,
)
from circyto.pipeline.collect import collect_matrix
from circyto.writers.convert import convert_matrix_files
from circyto.pipeline.export_multimodal import export_multimodal as _export_multimodal
from circyto.pipeline.annotate_host_gene import annotate_host_genes
from circyto.pipeline.run_detector import run_detector_manifest
from circyto.pipeline.run_multidetector import run_multidetector_pipeline
from circyto.pipeline.merge_detectors import merge_detectors as _merge_detectors
from circyto.pipeline.compare_detectors import compare_detectors as _compare_detectors
from circyto.pipeline.multidetector_collect import (
    build_detector_matrix,
    write_matrix,
)
from circyto.pipeline.collect_find_circ3 import collect_find_circ3_matrix
from circyto.pipeline.collect_circexplorer2_matrix import (
    collect_circexplorer2_matrix as collect_circexplorer2_matrix_from_dir,
)
from circyto.utils import ensure_dir
from circyto.detectors import build_default_engines

app = typer.Typer(
    add_completion=False,
    help=(
        "circyto — CLI toolkit for single-cell circRNA detection and integration.\n\n"
        "Includes:\n"
        "  [LEGACY] CIRI-full wrappers (prepare/run/run-manifest/collect/make)\n"
        "  [NEW]    Detector-API runners (run-detector / run-batch / run-multidetector)\n"
        "  [MATRIX] collect-matrix, collect-find-circ3, collect-circexplorer2,\n"
        "           collect-multidetector\n"
        "  [IO]     export-multimodal, annotate-host-genes, convert\n"
        "  [COMPARE] merge-detectors, compare-detectors\n"
    ),
)
console = Console()

# --------------------------------------------------------------------------------------
# LEGACY 10x / CIRI-full pipeline commands (kept for backwards compatibility)
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
        help="[LEGACY] Output directory for per-cell FASTQs",
    ),
    whitelist: Optional[Path] = typer.Option(
        None,
        help="[LEGACY] Optional whitelist of barcodes (gz/tsv)",
    ),
    chemistry: str = typer.Option(
        "tenx-3p",
        help="[LEGACY] Library type: tenx-3p or tenx-5p",
    ),
    batch_size: int = typer.Option(
        100,
        help="[LEGACY] Cells per FASTQ batch",
    ),
    min_reads_per_cell: int = typer.Option(
        200,
        help="[LEGACY] Discard cells with < N reads",
    ),
) -> None:
    """
    [LEGACY] Extract per-cell FASTQs from a 10x-style BAM.

    This helper is kept for backwards compatibility. New workflows should usually
    go through detector-API commands (run-detector / run-batch) instead.
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
        help="[LEGACY] Directory of FASTQ batches (from `prepare`)",
    ),
    outdir: Path = typer.Option(
        ...,
        help="[LEGACY] Output directory for CIRI-full results",
    ),
    cmd_template: str = typer.Option(
        ...,
        help=(
            "[LEGACY] CIRI-full command template; uses placeholders:\n"
            "  {ref_fa}, {gtf}, {r1}, {r2}, {out_tsv}\n"
            "Example:\n"
            "  'bash -lc \"export R1={r1} R2={r2} REF_FA={ref_fa} GTF={gtf} "
            "OUT_TSV={out_tsv} ; tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh\"'"
        ),
    ),
    ref_fa: Path = typer.Option(
        ...,
        exists=True,
        help="[LEGACY] Reference FASTA",
    ),
    gtf: Path = typer.Option(
        ...,
        exists=True,
        help="[LEGACY] GTF/GFF annotation",
    ),
    threads: int = typer.Option(
        8,
        help="[LEGACY] Number of batches to run in parallel",
    ),
) -> None:
    """
    [LEGACY] Run CIRI-full over batched FASTQs using a shell command template.

    New workflows should prefer `run-detector ciri-full` / `run-batch --detector ciri-full`.
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
        help="[LEGACY] Output directory for CIRI-full results",
    ),
    cmd_template: str = typer.Option(
        ...,
        help=(
            "[LEGACY] CIRI-full command template. Valid placeholders:\n"
            "  {ref_fa}  reference FASTA path\n"
            "  {gtf}     GTF/GFF annotation path\n"
            "  {r1}      read 1 FASTQ\n"
            "  {r2}      read 2 FASTQ (empty if missing)\n"
            "  {out_tsv} per-cell output TSV path\n"
        ),
    ),
    ref_fa: Path = typer.Option(
        ...,
        exists=True,
        help="[LEGACY] Reference FASTA",
    ),
    gtf: Path = typer.Option(
        ...,
        exists=True,
        help="[LEGACY] GTF/GFF annotation",
    ),
    threads: int = typer.Option(
        8,
        help="[LEGACY] Number of cells to run in parallel",
    ),
) -> None:
    """
    [LEGACY] Run CIRI-full over a plate/full-length manifest using a shell template.

    This command exists for backwards compatibility and for complex shell-based
    CIRI-full setups. For new work, prefer:

        circyto run-batch --detector ciri-full --manifest ...

    or:

        circyto run-detector ciri-full --manifest ...

    which use the detector API instead of a raw shell template.
    """
    console.print(
        "[yellow][LEGACY][/yellow] `run-manifest` uses a CIRI-full shell template.\n"
        "New projects should consider `run-batch --detector ciri-full` instead."
    )
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
        help="[LEGACY] Directory with CIRI-full per-cell outputs (*.tsv)",
    ),
    matrix: Path = typer.Option(
        ...,
        help="Output sparse matrix (.mtx, rows=circ, cols=cells)",
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
    [LEGACY] Collect CIRI-full per-cell TSVs into a circ × cell MatrixMarket matrix.

    New code should prefer `collect-matrix --detector ciri-full` as the unified collector.
    """
    collect_matrix(
        cirifull_dir=str(cirifull_dir),
        matrix_path=str(matrix),
        circ_index_path=str(circ_index),
        cell_index_path=str(cell_index),
        min_count_per_cell=min_count_per_cell,
    )


@app.command("collect-find-circ3")
def collect_find_circ3_cmd(
    findcirc3_dir: Path = typer.Option(
        ...,
        exists=True,
        help=(
            "Directory with find_circ3 per-cell outputs "
            "(<cell_id>/<cell_id>_splice_sites.bed)."
        ),
    ),
    matrix: Path = typer.Option(
        ...,
        help="Output sparse matrix (.mtx, rows=circ, cols=cells)",
    ),
    circ_index: Path = typer.Option(
        ...,
        help="Output circ index (rows, one circ_id per line)",
    ),
    cell_index: Path = typer.Option(
        ...,
        help="Output cell index (columns, one cell_id per line)",
    ),
    min_count_per_cell: int = typer.Option(
        1,
        help="Minimum total n_reads per cell to keep that column.",
    ),
) -> None:
    """
    Collect find_circ3 per-cell splice_sites.bed files into a circ × cell
    MatrixMarket matrix + circ/cell index files.
    """
    collect_find_circ3_matrix(
        findcirc3_dir=str(findcirc3_dir),
        matrix_path=str(matrix),
        circ_index_path=str(circ_index),
        cell_index_path=str(cell_index),
        min_count_per_cell=min_count_per_cell,
    )


@app.command("collect-circexplorer2")
def collect_circexplorer2_cmd(
    circexplorer2_dir: Path = typer.Option(
        ...,
        exists=True,
        help=(
            "Directory with CIRCexplorer2 per-cell outputs "
            "(<cell_id>/circularRNA_known.txt)."
        ),
    ),
    matrix: Path = typer.Option(
        ...,
        help="Output sparse matrix (.mtx, rows=circ, cols=cells)",
    ),
    circ_index: Path = typer.Option(
        ...,
        help="Output circ index (rows, one circ_id per line)",
    ),
    cell_index: Path = typer.Option(
        ...,
        help="Output cell index (columns, one cell_id per line)",
    ),
    min_support: int = typer.Option(
        1,
        help="Minimum readNumber (support) per circRNA per cell.",
    ),
) -> None:
    """
    Collect CIRCexplorer2 per-cell circularRNA_known.txt files into a circ × cell
    MatrixMarket matrix + circ/cell index files.
    """
    outdir = matrix.parent
    outdir.mkdir(parents=True, exist_ok=True)

    collect_circexplorer2_matrix_from_dir(
        indir=circexplorer2_dir,
        outdir=outdir,
        min_support=min_support,
    )

    # Move/rename default outputs to requested paths if needed
    default_matrix = outdir / "circ_counts.mtx"
    default_circ = outdir / "circ_index.txt"
    default_cell = outdir / "cell_index.txt"

    if default_matrix.exists() and matrix != default_matrix:
        default_matrix.replace(matrix)
    if default_circ.exists() and circ_index != default_circ:
        default_circ.replace(circ_index)
    if default_cell.exists() and cell_index != default_cell:
        default_cell.replace(cell_index)


# --------------------------------------------------------------------------------------
# Conversion, one-shot make
# --------------------------------------------------------------------------------------


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
        matrix_path=matrix,
        circ_index_path=circ_index,
        cell_index_path=cell_index,
        loom=loom,
        h5ad=h5ad,
    )


@app.command()
def make(
    outdir: Path = typer.Option(
        ...,
        help="[LEGACY] Work output directory",
    ),
    cmd_template: str = typer.Option(
        ...,
        help="[LEGACY] CIRI-full command template (see run-manifest help for placeholders)",
    ),
    ref_fa: Path = typer.Option(
        ...,
        exists=True,
        help="[LEGACY] Reference FASTA",
    ),
    gtf: Path = typer.Option(
        ...,
        exists=True,
        help="[LEGACY] GTF/GFF annotation",
    ),
    manifest: Optional[Path] = typer.Option(
        None,
        help="[LEGACY] Plate-style TSV listing cells and FASTQs",
    ),
    bam: Optional[Path] = typer.Option(
        None,
        exists=True,
        help="[LEGACY] 10x-style BAM with CB/UB tags",
    ),
    whitelist: Optional[Path] = typer.Option(
        None,
        help="[LEGACY] Optional barcode whitelist (gz/tsv) for 10x path",
    ),
    chemistry: str = typer.Option(
        "tenx-3p",
        help="[LEGACY] tenx-3p or tenx-5p (for 10x path)",
    ),
    threads: int = typer.Option(
        8,
        help="[LEGACY] Number of cells/batches to run in parallel",
    ),
) -> None:
    """
    [LEGACY] Convenience wrapper combining CIRI-full calling + collect + convert.

    - Plate path: --manifest
    - 10x path:   --bam (+ whitelist / chemistry)

    Writes circ_counts.mtx, circ_index.txt, cell_index.txt, and circ.h5ad into outdir.
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
        cirifull_dir=str(outdir / "cirifull_out"),
        matrix_path=str(mat),
        circ_index_path=str(circ_idx),
        cell_index_path=str(cell_idx),
        min_count_per_cell=1,
    )
    convert_matrix_files(
        matrix_path=mat,
        circ_index_path=circ_idx,
        cell_index_path=cell_idx,
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
        help="circ_feature_table.tsv from `circyto collect` or collect-matrix",
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
# Detector API: single-detector and NEW unified batch runner
# --------------------------------------------------------------------------------------


@app.command("run-detector")
def run_detector_cmd(
    detector: str = typer.Argument(
        ...,
        help="Detector name (e.g. 'ciri-full', 'ciri2', 'find-circ3', 'circexplorer2')",
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

    This is the preferred interface for new workflows (single-detector).
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


@app.command("run-batch")
def run_batch_cmd(
    detector: str = typer.Option(
        ...,
        "--detector",
        "-d",
        help="Detector name (e.g. 'ciri-full', 'find-circ3', 'circexplorer2')",
    ),
    manifest: Path = typer.Option(
        ...,
        exists=True,
        help="Manifest TSV (columns: cell_id, r1, [r2])",
    ),
    outdir: Path = typer.Option(
        ...,
        help="Output directory for per-cell detector outputs",
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
    Unified batch runner for any detector.

    This is the RECOMMENDED entry point for multi-cell runs.

        circyto run-batch --detector ciri-full --manifest manifest.tsv --outdir out/ ...

    is equivalent to:

        circyto run-detector ciri-full --manifest manifest.tsv --outdir out/ ...
    """
    engines = build_default_engines()
    if detector not in engines:
        available = ", ".join(sorted(engines.keys()))
        raise typer.Exit(f"Detector '{detector}' not available. Available: {available}")

    det_engine = engines[detector]
    console.print(
        f"[bold cyan][circyto][/bold cyan] [run-batch] detector='{det_engine.name}' "
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
        f"[bold cyan][circyto][/bold cyan] [run-batch] Completed {len(results)} jobs into {outdir}"
    )


@app.command("run-multidetector")
def run_multidetector_cmd(
    detectors: List[str] = typer.Argument(..., help="List of detectors to run"),
    manifest: Path = typer.Option(..., exists=True, dir_okay=False),
    outdir: Path = typer.Argument(...),
    ref_fa: Optional[Path] = typer.Option(None),
    gtf: Optional[Path] = typer.Option(None),
    threads: int = typer.Option(8),
    parallel: int = typer.Option(1),
) -> None:
    """
    Run multiple detectors on the same manifest.

    Example:

        circyto run-multidetector ciri-full find-circ3 \\
          --manifest manifest.tsv \\
          --ref-fa ref.fa --gtf ref.gtf \\
          out/multidetector_run
    """
    outdir.mkdir(parents=True, exist_ok=True)

    run_multidetector_pipeline(
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


@app.command("collect-matrix")
def collect_matrix_cmd(
    detector: str = typer.Option(
        ...,
        "--detector",
        "-d",
        help="Detector name for matrix collection (e.g. 'ciri-full', 'find-circ3', 'circexplorer2')",
    ),
    indir: Path = typer.Option(
        ...,
        exists=True,
        help="Directory with per-cell outputs for the detector",
    ),
    matrix: Path = typer.Option(
        ...,
        help="Output sparse matrix (.mtx, rows=circ, cols=cells)",
    ),
    circ_index: Path = typer.Option(
        ...,
        help="Output circ index (rows, one circ_id per line)",
    ),
    cell_index: Path = typer.Option(
        ...,
        help="Output cell index (columns, one cell_id per line)",
    ),
    min_count_per_cell: int = typer.Option(
        1,
        help="Minimum total counts per cell to keep that cell.",
    ),
) -> None:
    """
    Unified matrix collector.

    Dispatches to the appropriate collector based on --detector:

      - ciri-full      → legacy `collect` (CIRI-full TSVs)
      - find-circ3     → `collect-find-circ3`
      - circexplorer2  → `collect-circexplorer2`
    """
    detector = detector.lower()
    if detector == "ciri-full":
        collect_matrix(
            cirifull_dir=str(indir),
            matrix_path=str(matrix),
            circ_index_path=str(circ_index),
            cell_index_path=str(cell_index),
            min_count_per_cell=min_count_per_cell,
        )
    elif detector == "find-circ3":
        collect_find_circ3_matrix(
            findcirc3_dir=str(indir),
            matrix_path=str(matrix),
            circ_index_path=str(circ_index),
            cell_index_path=str(cell_index),
            min_count_per_cell=min_count_per_cell,
        )
    elif detector == "circexplorer2":
        outdir = matrix.parent
        outdir.mkdir(parents=True, exist_ok=True)
        collect_circexplorer2_matrix_from_dir(
            indir=indir,
            outdir=outdir,
            min_support=min_count_per_cell,
        )
        default_matrix = outdir / "circ_counts.mtx"
        default_circ = outdir / "circ_index.txt"
        default_cell = outdir / "cell_index.txt"
        if default_matrix.exists() and matrix != default_matrix:
            default_matrix.replace(matrix)
        if default_circ.exists() and circ_index != default_circ:
            default_circ.replace(circ_index)
        if default_cell.exists() and cell_index != default_cell:
            default_cell.replace(cell_index)
    else:
        raise typer.Exit(
            f"Unknown detector '{detector}' for collect-matrix. "
            "Supported: ciri-full, find-circ3, circexplorer2."
        )


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
    Build per-detector circRNA sparse matrices from multi-detector outputs.

    For each detector listed in multi_out/summary.json, writes:

      multi_out/matrices/<detector>.mtx
      multi_out/matrices/<detector>.circ.txt
      multi_out/matrices/<detector>.cell.txt

    Uses detector-specific collectors:
      - ciri-full      → legacy CIRI-full collector
      - find-circ3     → find_circ3 splice_sites.bed collector
      - circexplorer2  → CIRCexplorer2 circularRNA_known.txt collector
    """
    import json  # at top of file is also fine; here is okay too.

    console = Console()
    multi_out = Path(multi_out)
    summary_path = multi_out / "summary.json"

    if not summary_path.exists():
        raise typer.Exit(
            f"No summary.json found in {multi_out}. "
            "Run `circyto run-multidetector` first."
        )

    with summary_path.open() as f:
        payload = json.load(f)

    results = payload.get("results", {})
    if not results:
        raise typer.Exit(f"No 'results' key in {summary_path}.")

    matrices_dir = multi_out / "matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)

    # Fixed default for now (same as other collectors)
    min_count_per_cell = 1

    for det_name, det_results in results.items():
        console.print(
            f"[bold cyan][collect-multidetector][/bold cyan] "
            f"Building matrix for detector {det_name}"
        )

        det_outdir = multi_out / det_name
        if not det_outdir.exists():
            console.print(
                f"[yellow][collect-multidetector][/yellow] "
                f"Detector dir missing for {det_name}: {det_outdir} (skipping)"
            )
            continue

        prefix = matrices_dir / det_name
        matrix_path = prefix.with_suffix(".mtx")
        circ_index_path = matrices_dir / f"{det_name}.circ.txt"
        cell_index_path = matrices_dir / f"{det_name}.cell.txt"

        if det_name == "ciri-full":
            # legacy CIRI-full TSV collector
            collect_matrix(
                cirifull_dir=str(det_outdir),
                matrix_path=str(matrix_path),
                circ_index_path=str(circ_index_path),
                cell_index_path=str(cell_index_path),
                min_count_per_cell=min_count_per_cell,
            )

        elif det_name == "find-circ3":
            # per-cell <cell_id>/<cell_id>_splice_sites.bed
            collect_find_circ3_matrix(
                findcirc3_dir=str(det_outdir),
                matrix_path=str(matrix_path),
                circ_index_path=str(circ_index_path),
                cell_index_path=str(cell_index_path),
                min_count_per_cell=min_count_per_cell,
            )

        elif det_name == "circexplorer2":
            collect_circexplorer2_matrix(
                circexplorer2_dir=str(det_outdir),
                matrix_path=str(matrix_path),
                circ_index_path=str(circ_index_path),
                cell_index_path=str(cell_index_path),
                min_count_per_cell=min_count_per_cell,
            )

        else:
            console.print(
                f"[yellow][collect-multidetector][/yellow] "
                f"No collector wired for detector {det_name}; skipping."
            )
            continue

        console.print(
            f"[bold cyan][collect-multidetector][/bold cyan] "
            f"Wrote {matrix_path.name} "
            f"(circ_index={circ_index_path.name}, cell_index={cell_index_path.name})"
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


if __name__ == "__main__":
    app()
