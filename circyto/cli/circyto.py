from __future__ import annotations

import shutil
from pathlib import Path
from typing import List, Optional, Tuple

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
from circyto.pipeline.collect_find_circ3 import collect_find_circ3_matrix
from circyto.pipeline.collect_circexplorer2_matrix import (
    collect_circexplorer2_matrix as collect_circexplorer2_matrix_from_dir,
)
from circyto.detectors import build_default_engines

app = typer.Typer(
    add_completion=False,
    help=(
        "circyto — CLI toolkit for single-cell circRNA detection and integration.\n\n"
        "Conventions (locked):\n"
        "  - Output directories use:  --outdir / -o\n"
        "  - Input directories use:   --indir\n"
        "  - Most commands accept either:\n"
        "        COMMAND INDIR OUTDIR\n"
        "    or  COMMAND --indir INDIR --outdir OUTDIR\n"
        "  - Collectors default to writing:\n"
        "        OUTDIR/circ_counts.mtx\n"
        "        OUTDIR/circ_index.txt\n"
        "        OUTDIR/cell_index.txt\n\n"
        "Includes:\n"
        "  [LEGACY] CIRI-full wrappers (prepare/run/run-manifest/make)\n"
        "  [RUN]    run-detector / run-batch / run-multidetector\n"
        "  [MATRIX] collect-matrix (+ per-detector collectors)\n"
        "  [MERGE]  merge-detectors\n"
        "  [COMPARE] compare-ids (fuzzy/exact), compare-detectors (merged outputs)\n"
    ),
)
console = Console()


# --------------------------------------------------------------------------------------
# Helpers: consistent INDIR/OUTDIR + default output naming
# --------------------------------------------------------------------------------------


def _pick_one(
    pos: Optional[Path],
    opt: Optional[Path],
    *,
    name: str,
    allow_both_equal: bool = True,
) -> Optional[Path]:
    if pos is None and opt is None:
        return None
    if pos is not None and opt is not None:
        if allow_both_equal and pos == opt:
            return pos
        raise typer.BadParameter(f"Provide {name} positionally OR via option, not both.")
    return pos if pos is not None else opt


def _auto_outdir(*parts: str) -> Path:
    safe = "_".join([p for p in parts if p]).replace("/", "_").replace(" ", "_")
    return Path("work") / safe


def _default_matrix_paths(outdir: Path) -> Tuple[Path, Path, Path]:
    return (
        outdir / "circ_counts.mtx",
        outdir / "circ_index.txt",
        outdir / "cell_index.txt",
    )


def _resolve_collect_paths(
    outdir: Optional[Path],
    matrix: Optional[Path],
    circ_index: Optional[Path],
    cell_index: Optional[Path],
) -> Tuple[Optional[Path], Optional[Path], Optional[Path], Optional[Path]]:
    """
    If outdir is provided, fill in default matrix/index paths unless explicitly set.
    Returns (outdir, matrix, circ_index, cell_index).
    """
    if outdir is None:
        return (None, matrix, circ_index, cell_index)

    outdir.mkdir(parents=True, exist_ok=True)
    dmat, dcirc, dcell = _default_matrix_paths(outdir)
    if matrix is None:
        matrix = dmat
    if circ_index is None:
        circ_index = dcirc
    if cell_index is None:
        cell_index = dcell
    return (outdir, matrix, circ_index, cell_index)


def _require_paths(*paths: Optional[Path], names: List[str]) -> None:
    missing = [n for p, n in zip(paths, names) if p is None]
    if missing:
        raise typer.BadParameter(
            "Missing required outputs: "
            + ", ".join(missing)
            + ". Provide --outdir/-o or explicit --matrix/--circ-index/--cell-index."
        )


# --------------------------------------------------------------------------------------
# LEGACY 10x / CIRI-full pipeline commands (kept for backwards compatibility)
# --------------------------------------------------------------------------------------


@app.command()
def prepare(
    bam: Path = typer.Option(..., exists=True, help="10x-style BAM with CB/UB tags"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="[LEGACY] Output directory"),
    whitelist: Optional[Path] = typer.Option(None, help="[LEGACY] Optional whitelist"),
    chemistry: str = typer.Option("tenx-3p", help="[LEGACY] tenx-3p or tenx-5p"),
    batch_size: int = typer.Option(100, help="[LEGACY] Cells per FASTQ batch"),
    min_reads_per_cell: int = typer.Option(200, help="[LEGACY] Discard cells with < N reads"),
) -> None:
    """
    [LEGACY] Extract per-cell FASTQs from a 10x-style BAM.
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
    fastq_dir: Path = typer.Option(..., exists=True, help="[LEGACY] FASTQ batches dir"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="[LEGACY] Output directory"),
    cmd_template: str = typer.Option(..., help="[LEGACY] CIRI-full command template"),
    ref_fa: Path = typer.Option(..., exists=True, help="[LEGACY] Reference FASTA"),
    gtf: Path = typer.Option(..., exists=True, help="[LEGACY] GTF/GFF annotation"),
    threads: int = typer.Option(8, help="[LEGACY] Number of batches to run in parallel"),
) -> None:
    """
    [LEGACY] Run CIRI-full over batched FASTQs using a shell command template.
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
    manifest: Path = typer.Option(..., exists=True, help="TSV: cell_id, r1, [r2]"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="[LEGACY] Output directory"),
    cmd_template: str = typer.Option(..., help="[LEGACY] CIRI-full command template"),
    ref_fa: Path = typer.Option(..., exists=True, help="[LEGACY] Reference FASTA"),
    gtf: Path = typer.Option(..., exists=True, help="[LEGACY] GTF/GFF annotation"),
    threads: int = typer.Option(8, help="[LEGACY] Number of cells to run in parallel"),
) -> None:
    """
    [LEGACY] Run CIRI-full over a plate/full-length manifest using a shell template.
    """
    console.print(
        "[yellow][LEGACY][/yellow] `run-manifest` uses a CIRI-full shell template.\n"
        "Prefer `run-detector` / `run-batch` for new workflows."
    )
    run_cirifull_with_manifest(
        manifest=manifest,
        outdir=outdir,
        cmd_template=cmd_template,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
    )


# --------------------------------------------------------------------------------------
# Collectors (now consistent): accept INDIR OUTDIR OR --indir/--outdir OR explicit file paths
# --------------------------------------------------------------------------------------


@app.command()
def collect(
    indir_pos: Optional[Path] = typer.Argument(
        None, metavar="INDIR", help="[LEGACY] CIRI-full per-cell outputs dir (*.tsv)"
    ),
    outdir_pos: Optional[Path] = typer.Argument(
        None, metavar="OUTDIR", help="Output directory (defaults to work/...)"
    ),
    cirifull_dir: Optional[Path] = typer.Option(
        None,
        "--cirifull-dir",
        "--indir",
        exists=True,
        help="[LEGACY] Directory with CIRI-full per-cell outputs (*.tsv)",
    ),
    outdir_opt: Optional[Path] = typer.Option(
        None, "--outdir", "-o", help="Output directory (writes circ_counts.mtx + indexes)"
    ),
    matrix: Optional[Path] = typer.Option(
        None, "--matrix", help="Output sparse matrix (.mtx). Optional if --outdir is used."
    ),
    circ_index: Optional[Path] = typer.Option(
        None, "--circ-index", help="Output circ index. Optional if --outdir is used."
    ),
    cell_index: Optional[Path] = typer.Option(
        None, "--cell-index", help="Output cell index. Optional if --outdir is used."
    ),
    min_count_per_cell: int = typer.Option(
        1, "--min-count-per-cell", help="Drop cells with total circ counts < threshold"
    ),
) -> None:
    """
    [LEGACY] Collect CIRI-full per-cell TSVs into a circ × cell MatrixMarket matrix.

    Consistent usage:
      circyto collect INDIR OUTDIR
      circyto collect --indir INDIR --outdir OUTDIR
    """
    indir = _pick_one(indir_pos, cirifull_dir, name="INDIR/--cirifull-dir/--indir")
    if indir is None:
        raise typer.BadParameter("Provide INDIR positionally or via --cirifull-dir/--indir.")
    if not indir.exists():
        raise typer.BadParameter(f"INDIR does not exist: {indir}")

    outdir = _pick_one(outdir_pos, outdir_opt, name="OUTDIR/--outdir")
    if outdir is None and (matrix is None or circ_index is None or cell_index is None):
        outdir = _auto_outdir("collect_ciri-full", indir.name)

    outdir, matrix, circ_index, cell_index = _resolve_collect_paths(
        outdir, matrix, circ_index, cell_index
    )
    _require_paths(matrix, circ_index, cell_index, names=["--matrix", "--circ-index", "--cell-index"])

    assert matrix and circ_index and cell_index
    matrix.parent.mkdir(parents=True, exist_ok=True)
    circ_index.parent.mkdir(parents=True, exist_ok=True)
    cell_index.parent.mkdir(parents=True, exist_ok=True)

    collect_matrix(
        cirifull_dir=str(indir),
        matrix_path=str(matrix),
        circ_index_path=str(circ_index),
        cell_index_path=str(cell_index),
        min_count_per_cell=min_count_per_cell,
    )
    console.print(f"[bold cyan][collect][/bold cyan] Wrote: {matrix} {circ_index} {cell_index}")


@app.command("collect-find-circ3")
def collect_find_circ3_cmd(
    indir_pos: Optional[Path] = typer.Argument(
        None, metavar="INDIR", help="find_circ3 outputs dir (<cell>/<cell>_splice_sites.bed)"
    ),
    outdir_pos: Optional[Path] = typer.Argument(
        None, metavar="OUTDIR", help="Output directory (defaults to work/...)"
    ),
    findcirc3_dir: Optional[Path] = typer.Option(
        None,
        "--findcirc3-dir",
        "--indir",
        exists=True,
        help="Directory with find_circ3 per-cell outputs (<cell>/<cell>_splice_sites.bed).",
    ),
    outdir_opt: Optional[Path] = typer.Option(
        None, "--outdir", "-o", help="Output directory (writes circ_counts.mtx + indexes)"
    ),
    matrix: Optional[Path] = typer.Option(None, "--matrix"),
    circ_index: Optional[Path] = typer.Option(None, "--circ-index"),
    cell_index: Optional[Path] = typer.Option(None, "--cell-index"),
    min_count_per_cell: int = typer.Option(
        1, "--min-count-per-cell", help="Minimum total n_reads per cell to keep that column."
    ),
) -> None:
    """
    Collect find_circ3 per-cell splice_sites.bed into circ × cell matrix.

    Consistent usage:
      circyto collect-find-circ3 INDIR OUTDIR
      circyto collect-find-circ3 --indir INDIR --outdir OUTDIR
    """
    indir = _pick_one(indir_pos, findcirc3_dir, name="INDIR/--findcirc3-dir/--indir")
    if indir is None:
        raise typer.BadParameter("Provide INDIR positionally or via --findcirc3-dir/--indir.")
    if not indir.exists():
        raise typer.BadParameter(f"INDIR does not exist: {indir}")

    outdir = _pick_one(outdir_pos, outdir_opt, name="OUTDIR/--outdir")
    if outdir is None and (matrix is None or circ_index is None or cell_index is None):
        outdir = _auto_outdir("collect_find-circ3", indir.name)

    outdir, matrix, circ_index, cell_index = _resolve_collect_paths(
        outdir, matrix, circ_index, cell_index
    )
    _require_paths(matrix, circ_index, cell_index, names=["--matrix", "--circ-index", "--cell-index"])

    assert matrix and circ_index and cell_index
    matrix.parent.mkdir(parents=True, exist_ok=True)
    circ_index.parent.mkdir(parents=True, exist_ok=True)
    cell_index.parent.mkdir(parents=True, exist_ok=True)

    collect_find_circ3_matrix(
        findcirc3_dir=str(indir),
        matrix_path=str(matrix),
        circ_index_path=str(circ_index),
        cell_index_path=str(cell_index),
        min_count_per_cell=min_count_per_cell,
    )
    console.print(
        f"[bold cyan][collect-find-circ3][/bold cyan] Wrote: {matrix} {circ_index} {cell_index}"
    )


@app.command("collect-circexplorer2")
def collect_circexplorer2_cmd(
    indir_pos: Optional[Path] = typer.Argument(
        None, metavar="INDIR", help="CIRCexplorer2 outputs dir (<cell>/circularRNA_known.txt)"
    ),
    outdir_pos: Optional[Path] = typer.Argument(
        None, metavar="OUTDIR", help="Output directory (defaults to work/...)"
    ),
    circexplorer2_dir: Optional[Path] = typer.Option(
        None,
        "--circexplorer2-dir",
        "--indir",
        exists=True,
        help="Directory with CIRCexplorer2 per-cell outputs (<cell>/circularRNA_known.txt).",
    ),
    outdir_opt: Optional[Path] = typer.Option(
        None, "--outdir", "-o", help="Output directory (writes circ_counts.mtx + indexes)"
    ),
    matrix: Optional[Path] = typer.Option(None, "--matrix"),
    circ_index: Optional[Path] = typer.Option(None, "--circ-index"),
    cell_index: Optional[Path] = typer.Option(None, "--cell-index"),
    min_support: int = typer.Option(
        1, "--min-support", help="Minimum readNumber (support) per circRNA per cell."
    ),
) -> None:
    """
    Collect CIRCexplorer2 outputs into circ × cell matrix.

    Consistent usage:
      circyto collect-circexplorer2 INDIR OUTDIR
      circyto collect-circexplorer2 --indir INDIR --outdir OUTDIR
    """
    indir = _pick_one(indir_pos, circexplorer2_dir, name="INDIR/--circexplorer2-dir/--indir")
    if indir is None:
        raise typer.BadParameter("Provide INDIR positionally or via --circexplorer2-dir/--indir.")
    if not indir.exists():
        raise typer.BadParameter(f"INDIR does not exist: {indir}")

    outdir = _pick_one(outdir_pos, outdir_opt, name="OUTDIR/--outdir")
    if outdir is None and (matrix is None or circ_index is None or cell_index is None):
        outdir = _auto_outdir("collect_circexplorer2", indir.name)

    outdir, matrix, circ_index, cell_index = _resolve_collect_paths(
        outdir, matrix, circ_index, cell_index
    )
    _require_paths(matrix, circ_index, cell_index, names=["--matrix", "--circ-index", "--cell-index"])

    assert outdir and matrix and circ_index and cell_index
    outdir.mkdir(parents=True, exist_ok=True)

    # Write defaults into a temp folder to avoid collisions, then rename/move to requested paths.
    tmp = outdir / "_tmp_circexplorer2_collect"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True, exist_ok=True)

    collect_circexplorer2_matrix_from_dir(indir=indir, outdir=tmp, min_support=min_support)

    default_matrix = tmp / "circ_counts.mtx"
    default_circ = tmp / "circ_index.txt"
    default_cell = tmp / "cell_index.txt"

    matrix.parent.mkdir(parents=True, exist_ok=True)
    circ_index.parent.mkdir(parents=True, exist_ok=True)
    cell_index.parent.mkdir(parents=True, exist_ok=True)

    if default_matrix.exists():
        default_matrix.replace(matrix)
    if default_circ.exists():
        default_circ.replace(circ_index)
    if default_cell.exists():
        default_cell.replace(cell_index)

    shutil.rmtree(tmp, ignore_errors=True)
    console.print(
        f"[bold cyan][collect-circexplorer2][/bold cyan] Wrote: {matrix} {circ_index} {cell_index}"
    )


# --------------------------------------------------------------------------------------
# Conversion, one-shot make
# --------------------------------------------------------------------------------------


@app.command()
def convert(
    matrix: Path = typer.Option(..., exists=True, help="Sparse matrix .mtx (rows=circ, cols=cells)"),
    circ_index: Path = typer.Option(..., exists=True, help="Text file of circ IDs (one per line)"),
    cell_index: Path = typer.Option(..., exists=True, help="Text file of cell IDs (one per line)"),
    loom: Optional[Path] = typer.Option(None, help="Optional path to write .loom"),
    h5ad: Optional[Path] = typer.Option(None, help="Optional path to write .h5ad"),
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
    outdir: Path = typer.Option(..., "--outdir", "-o", help="[LEGACY] Work output directory"),
    cmd_template: str = typer.Option(..., help="[LEGACY] CIRI-full command template"),
    ref_fa: Path = typer.Option(..., exists=True, help="[LEGACY] Reference FASTA"),
    gtf: Path = typer.Option(..., exists=True, help="[LEGACY] GTF/GFF annotation"),
    manifest: Optional[Path] = typer.Option(None, help="[LEGACY] Plate-style TSV listing cells and FASTQs"),
    bam: Optional[Path] = typer.Option(None, exists=True, help="[LEGACY] 10x-style BAM with CB/UB tags"),
    whitelist: Optional[Path] = typer.Option(None, help="[LEGACY] Optional barcode whitelist"),
    chemistry: str = typer.Option("tenx-3p", help="[LEGACY] tenx-3p or tenx-5p"),
    threads: int = typer.Option(8, help="[LEGACY] Number of cells/batches to run in parallel"),
) -> None:
    """
    [LEGACY] Convenience wrapper combining CIRI-full calling + collect + convert.
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
    genes_h5ad: Path = typer.Option(..., exists=True, help="Base gene-expression .h5ad"),
    circ_matrix: Path = typer.Option(..., exists=True, help="circRNA MatrixMarket .mtx"),
    circ_index: Path = typer.Option(..., exists=True, help="circRNA index (rows)"),
    cell_index: Path = typer.Option(..., exists=True, help="Cell index (columns)"),
    out: Path = typer.Option(..., help="Output multimodal .h5ad"),
    circ_feature_table: Optional[Path] = typer.Option(None, help="Optional circ_feature_table.tsv"),
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
    circ_feature_table: Path = typer.Option(..., exists=True, help="circ_feature_table.tsv"),
    gtf: Path = typer.Option(..., exists=True, help="Reference GTF used for circ calling"),
    out: Optional[Path] = typer.Option(None, help="Output TSV (default: overwrite circ_feature_table)"),
    max_genes_per_circ: int = typer.Option(5, help="Maximum number of host genes to record per circRNA"),
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
# Detector API runners (consistent): --outdir everywhere; detector can be arg OR --detector
# --------------------------------------------------------------------------------------


def _run_detector_impl(
    *,
    detector: str,
    manifest: Path,
    outdir: Path,
    ref_fa: Optional[Path],
    gtf: Optional[Path],
    threads: int,
    parallel: int,
) -> None:
    engines = build_default_engines()
    if detector not in engines:
        available = ", ".join(sorted(engines.keys()))
        raise typer.Exit(f"Detector '{detector}' not available. Available: {available}")

    det_engine = engines[detector]
    console.print(
        f"[bold cyan][circyto][/bold cyan] Running detector='{det_engine.name}' "
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


@app.command("run-detector")
def run_detector_cmd(
    detector_pos: Optional[str] = typer.Argument(
        None, metavar="DETECTOR", help="Detector name (e.g. ciri-full, find-circ3, circexplorer2)"
    ),
    detector_opt: Optional[str] = typer.Option(
        None, "--detector", "-d", help="Detector name (alias of positional DETECTOR; e.g. circexplorer2)"
    ),
    manifest: Path = typer.Option(..., exists=True, help="Manifest TSV (cell_id, r1, [r2])"),
    outdir: Optional[Path] = typer.Option(
        None, "--outdir", "-o", help="Output directory for per-cell detector outputs"
    ),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa", help="Reference FASTA"),
    gtf: Optional[Path] = typer.Option(None, "--gtf", help="Annotation GTF/GFF"),
    threads: int = typer.Option(8, "--threads", help="Threads per detector process"),
    parallel: int = typer.Option(4, "--parallel", help="Number of cells to run in parallel"),
) -> None:
    """
    Run a single detector over a manifest using the detector API.

    Examples:
      circyto run-detector find-circ3 --manifest m.tsv --outdir out/ --ref-fa ref.fa
      circyto run-detector --detector find-circ3 --manifest m.tsv --outdir out/ --ref-fa ref.fa
    """
    detector = detector_pos or detector_opt
    if detector_pos and detector_opt and detector_pos != detector_opt:
        raise typer.BadParameter("Provide DETECTOR positionally OR via --detector, not both.")
    if detector is None:
        raise typer.BadParameter("Missing detector. Provide DETECTOR or --detector/-d.")

    if outdir is None:
        outdir = _auto_outdir("run-detector", detector, manifest.stem)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

    _run_detector_impl(
        detector=detector,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
        parallel=parallel,
    )


@app.command("run-batch")
def run_batch_cmd(
    detector: str = typer.Option(..., "--detector", "-d", help="Detector name"),
    manifest: Path = typer.Option(..., exists=True, help="Manifest TSV (cell_id, r1, [r2])"),
    outdir: Optional[Path] = typer.Option(None, "--outdir", "-o", help="Output directory"),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa", help="Reference FASTA"),
    gtf: Optional[Path] = typer.Option(None, "--gtf", help="Annotation GTF/GFF"),
    threads: int = typer.Option(8, "--threads", help="Threads per detector process"),
    parallel: int = typer.Option(4, "--parallel", help="Number of cells to run in parallel"),
) -> None:
    """
    Unified batch runner for any detector.

    NOTE: This is effectively an alias of run-detector with a strict `--detector`.
    """
    if outdir is None:
        outdir = _auto_outdir("run-batch", detector, manifest.stem)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

    _run_detector_impl(
        detector=detector,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
        parallel=parallel,
    )


@app.command("run-multidetector")
def run_multidetector_cmd(
    detectors: List[str] = typer.Argument(..., help="List of detectors to run"),
    manifest: Path = typer.Option(..., exists=True, dir_okay=False, help="Manifest TSV"),
    outdir: Optional[Path] = typer.Option(
        None, "--outdir", "-o", help="Output directory for multi-detector run"
    ),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa", help="Reference FASTA"),
    gtf: Optional[Path] = typer.Option(None, "--gtf", help="Annotation GTF/GFF"),
    threads: int = typer.Option(8, "--threads"),
    parallel: int = typer.Option(1, "--parallel"),
) -> None:
    """
    Run multiple detectors on the same manifest.

    Consistent usage:
      circyto run-multidetector ciri-full find-circ3 --manifest m.tsv --outdir out/
    """
    if outdir is None:
        outdir = _auto_outdir("run-multidetector", manifest.stem)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

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
# Unified matrix collector (consistent): --detector --indir --outdir
# --------------------------------------------------------------------------------------


@app.command("collect-matrix")
def collect_matrix_cmd(
    detector: str = typer.Option(
        ...,
        "--detector",
        "-d",
        help="Detector name (ciri-full, find-circ3, circexplorer2)",
    ),
    indir: Path = typer.Option(..., "--indir", exists=True, help="Directory with per-cell outputs"),
    outdir: Optional[Path] = typer.Option(
        None, "--outdir", "-o", help="Output directory (writes circ_counts.mtx + indexes)"
    ),
    matrix: Optional[Path] = typer.Option(None, "--matrix"),
    circ_index: Optional[Path] = typer.Option(None, "--circ-index"),
    cell_index: Optional[Path] = typer.Option(None, "--cell-index"),
    min_count_per_cell: int = typer.Option(1, "--min-count-per-cell"),
) -> None:
    """
    Unified matrix collector.

    Preferred usage:
      circyto collect-matrix -d find-circ3 --indir work/.../find-circ3 --outdir out/
    """
    detector = detector.lower()
    if outdir is None and (matrix is None or circ_index is None or cell_index is None):
        outdir = _auto_outdir("collect-matrix", detector, indir.name)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

    outdir, matrix, circ_index, cell_index = _resolve_collect_paths(
        outdir, matrix, circ_index, cell_index
    )
    _require_paths(matrix, circ_index, cell_index, names=["--matrix", "--circ-index", "--cell-index"])
    assert matrix and circ_index and cell_index

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
        # Use temp dir to avoid default-name collisions.
        tmp = matrix.parent / "_tmp_collect_matrix_circexplorer2"
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir(parents=True, exist_ok=True)
        collect_circexplorer2_matrix_from_dir(indir=indir, outdir=tmp, min_support=min_count_per_cell)

        dm = tmp / "circ_counts.mtx"
        dc = tmp / "circ_index.txt"
        dl = tmp / "cell_index.txt"
        if dm.exists():
            dm.replace(matrix)
        if dc.exists():
            dc.replace(circ_index)
        if dl.exists():
            dl.replace(cell_index)
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        raise typer.Exit(
            f"Unknown detector '{detector}' for collect-matrix. "
            "Supported: ciri-full, find-circ3, circexplorer2."
        )

    console.print(f"[bold cyan][collect-matrix][/bold cyan] Wrote: {matrix} {circ_index} {cell_index}")


# --------------------------------------------------------------------------------------
# Multi-detector: merge, compare (consistent): accept positional OR --indir/--outdir
# --------------------------------------------------------------------------------------


@app.command("merge-detectors")
def merge_detectors_cmd(
    indir_pos: Optional[Path] = typer.Argument(None, metavar="INDIR", help="Output dir from run-multidetector"),
    outdir_pos: Optional[Path] = typer.Argument(None, metavar="OUTDIR", help="Output directory for merged tables"),
    indir_opt: Optional[Path] = typer.Option(None, "--indir", exists=True, help="Input directory (alias)"),
    outdir_opt: Optional[Path] = typer.Option(None, "--outdir", "-o", help="Output directory (alias)"),
) -> None:
    """
    Merge per-detector per-cell TSVs into union + long-format tables.
    """
    indir = _pick_one(indir_pos, indir_opt, name="INDIR/--indir")
    if indir is None:
        raise typer.BadParameter("Provide INDIR positionally or via --indir.")
    if not indir.exists():
        raise typer.BadParameter(f"INDIR does not exist: {indir}")

    outdir = _pick_one(outdir_pos, outdir_opt, name="OUTDIR/--outdir")
    if outdir is None:
        outdir = _auto_outdir("merge-detectors", indir.name)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

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


@app.command("compare-ids")
def compare_ids_cmd(
    a: Path = typer.Option(..., "--a", help="circ_index.txt or TSV with circ_id column"),
    b: Path = typer.Option(..., "--b", help="circ_index.txt or TSV with circ_id column"),
    window: int = typer.Option(5, "--window", help="Fuzzy window in bp (±window)"),
    col: str = typer.Option("circ_id", "--col", help="TSV column name for circ_id"),
) -> None:
    """
    Compare two circ ID lists (exact + fuzzy ±window bp).

    This is the “sanity compare” tool you used for CIRI vs find-circ3.
    """
    from circyto.analysis.compare_detectors import load_keys, fuzzy_hits, fuzzy_jaccard

    A = load_keys(a, col=col)
    B = load_keys(b, col=col)

    exact_i = len(A & B)
    exact_u = len(A) + len(B) - exact_i

    ha = fuzzy_hits(A, B, window=window)
    hb = fuzzy_hits(B, A, window=window)
    fj = fuzzy_jaccard(A, B, window=window)

    typer.echo(f"A: {len(A)}")
    typer.echo(f"B: {len(B)}")
    typer.echo(f"Exact intersect: {exact_i}")
    typer.echo(f"Exact Jaccard:   {exact_i/exact_u if exact_u else 0.0:.6f}")
    typer.echo(f"Fuzzy hits A→B (±{window}): {ha}")
    typer.echo(f"Fuzzy hits B→A (±{window}): {hb}")
    typer.echo(f"Fuzzy Jaccard (sym): {fj:.6f}")


@app.command("collect-multidetector")
def collect_multidetector_cmd(
    indir_pos: Optional[Path] = typer.Argument(
        None, metavar="INDIR", help="Output dir from run-multidetector (contains summary.json)"
    ),
    indir_opt: Optional[Path] = typer.Option(
        None, "--indir", exists=True, help="Input directory (alias of positional INDIR)"
    ),
) -> None:
    """
    Build per-detector circRNA matrices from a multi-detector run into:
      INDIR/matrices/<detector>.mtx
      INDIR/matrices/<detector>.circ.txt
      INDIR/matrices/<detector>.cell.txt
    """
    import json

    multi_out = _pick_one(indir_pos, indir_opt, name="INDIR/--indir")
    if multi_out is None:
        raise typer.BadParameter("Provide INDIR positionally or via --indir.")
    if not multi_out.exists():
        raise typer.BadParameter(f"INDIR does not exist: {multi_out}")

    summary_path = multi_out / "summary.json"
    if not summary_path.exists():
        raise typer.Exit(
            f"No summary.json found in {multi_out}. Run `circyto run-multidetector` first."
        )

    with summary_path.open() as f:
        payload = json.load(f)

    results = payload.get("results", {})
    if not results:
        raise typer.Exit(f"No 'results' key in {summary_path}.")

    matrices_dir = multi_out / "matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)

    min_count_per_cell = 1

    for det_name in results.keys():
        console.print(
            f"[bold cyan][collect-multidetector][/bold cyan] Building matrix for detector {det_name}"
        )

        det_outdir = multi_out / det_name
        if not det_outdir.exists():
            console.print(
                f"[yellow][collect-multidetector][/yellow] Missing detector dir: {det_outdir} (skipping)"
            )
            continue

        prefix = matrices_dir / det_name
        matrix_path = prefix.with_suffix(".mtx")
        circ_index_path = matrices_dir / f"{det_name}.circ.txt"
        cell_index_path = matrices_dir / f"{det_name}.cell.txt"

        if det_name == "ciri-full":
            collect_matrix(
                cirifull_dir=str(det_outdir),
                matrix_path=str(matrix_path),
                circ_index_path=str(circ_index_path),
                cell_index_path=str(cell_index_path),
                min_count_per_cell=min_count_per_cell,
            )

        elif det_name == "find-circ3":
            collect_find_circ3_matrix(
                findcirc3_dir=str(det_outdir),
                matrix_path=str(matrix_path),
                circ_index_path=str(circ_index_path),
                cell_index_path=str(cell_index_path),
                min_count_per_cell=min_count_per_cell,
            )

        elif det_name == "circexplorer2":
            # Important: circexplorer2 collector writes fixed filenames;
            # use a per-detector temp directory to avoid collisions.
            tmp = matrices_dir / f"_tmp_{det_name}"
            if tmp.exists():
                shutil.rmtree(tmp)
            tmp.mkdir(parents=True, exist_ok=True)

            collect_circexplorer2_matrix_from_dir(
                indir=det_outdir,
                outdir=tmp,
                min_support=min_count_per_cell,
            )

            dm = tmp / "circ_counts.mtx"
            dc = tmp / "circ_index.txt"
            dl = tmp / "cell_index.txt"
            if dm.exists():
                dm.replace(matrix_path)
            if dc.exists():
                dc.replace(circ_index_path)
            if dl.exists():
                dl.replace(cell_index_path)

            shutil.rmtree(tmp, ignore_errors=True)

        else:
            console.print(
                f"[yellow][collect-multidetector][/yellow] No collector wired for {det_name}; skipping."
            )
            continue

        console.print(
            f"[bold cyan][collect-multidetector][/bold cyan] Wrote {matrix_path.name} "
            f"(circ_index={circ_index_path.name}, cell_index={cell_index_path.name})"
        )


@app.command("compare-detectors")
def compare_detectors_merged_cmd(
    indir_pos: Optional[Path] = typer.Argument(None, metavar="INDIR", help="Directory containing circ_union.tsv"),
    outdir_pos: Optional[Path] = typer.Argument(None, metavar="OUTDIR", help="Output directory for compare outputs"),
    indir_opt: Optional[Path] = typer.Option(None, "--indir", exists=True, help="Input directory (alias)"),
    outdir_opt: Optional[Path] = typer.Option(None, "--outdir", "-o", help="Output directory (alias)"),
) -> None:
    """
    Compare detectors using a merged circ_union.tsv (output of merge-detectors).

    Produces:
      - jaccard.tsv
      - detector_summary.tsv
      - compare_metadata.json
    """
    indir = _pick_one(indir_pos, indir_opt, name="INDIR/--indir")
    if indir is None:
        raise typer.BadParameter("Provide INDIR positionally or via --indir.")
    if not indir.exists():
        raise typer.BadParameter(f"INDIR does not exist: {indir}")

    outdir = _pick_one(outdir_pos, outdir_opt, name="OUTDIR/--outdir")
    if outdir is None:
        outdir = _auto_outdir("compare-detectors", indir.name)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

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
