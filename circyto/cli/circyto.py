from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import List, Optional, Tuple

import typer
from rich.console import Console
from circyto import __version__
from circyto.cli.analyze import analyze_app
from circyto.cli.demux import demux_app
from circyto.cli.manifest import manifest_app
from circyto.cli.workflow import workflow_app

from circyto.pipeline.prepare import extract_per_cell_fastq
from circyto.pipeline.run_cirifull import (
    run_cirifull_over_fastqs,
    run_cirifull_with_manifest,
)
from circyto.pipeline.collect import collect_matrix
from circyto.pipeline.annotate_circs import (
    annotate_circ_table,
    parse_annotation_db_spec,
)
from circyto.pipeline.annotate_host_gene import annotate_host_genes, repair_host_gene_file
from circyto.pipeline.align_manifest import (
    export_manifest_subset,
    plan_alignment_cache,
    prepare_alignment_cache,
    run_detector_alignment_manifest,
    summarize_alignment_chunks,
    summarize_run_state,
)
from circyto.pipeline.run_detector import run_detector_manifest
from circyto.pipeline.run_multidetector import run_multidetector_pipeline
from circyto.pipeline.run_ciri3 import RunCiri3Params, run_ciri3_workflow
from circyto.pipeline.scomatic_interop import (
    export_scomatic_inputs,
    join_circ_snv_summary,
)
from circyto.pipeline.scomatic_full_length_adapter import (
    import_scomatic,
    merge_scomatic,
    normalize_full_length_protocol,
    prepare_scomatic_input,
    run_scomatic,
)
from circyto.pipeline.scomatic_normalization import normalize_scomatic_results
from circyto.pipeline.scrr_cell_mapping import (
    build_scrr_cell_map,
    merge_scrr_cnv,
    remap_scrr_mudata_obs,
)
from circyto.pipeline.scrr_cnv_import import import_scrr_cnv
from circyto.pipeline.scrr_rt_import import import_scrr_rt, merge_scrr_rt
from circyto.pipeline.gene_expression_velocity import (
    SUPPORTED_CLEANUP_SCOPES,
    add_posthoc_rna_profile,
    cleanup_completed_workflow,
    export_completed_workflow_mudata,
    export_circ_bed,
    get_environment_summary,
    import_dna_snv_summary,
    inspect_completed_workdir,
    inspect_mudata_file,
    refresh_rna_qc_from_existing_outputs,
    summarize_benchmark_workdirs,
    summarize_circ_host_genes,
    summarize_dna_rna_circ,
    summarize_mudata_qc,
    summarize_rna_circ_integration,
    validate_completed_workdir,
)
from circyto.pipeline.workflow_integrity import check_workflow_integrity
from circyto.pipeline.merge_detectors import merge_detectors as _merge_detectors
from circyto.pipeline.compare_detectors import compare_detectors as _compare_detectors
from circyto.pipeline.collect_find_circ3 import collect_find_circ3_matrix
from circyto.pipeline.collect_circexplorer2_matrix import (
    collect_circexplorer2_matrix as collect_circexplorer2_matrix_from_dir,
)
from circyto.pipeline.public_dataset_prepare import prepare_public_dataset
from circyto.pipeline.scanpy_downstream import scanpy_pca_workflow, scanpy_qc_report
from circyto.detectors import build_default_engines, get_detector_capabilities
from circyto.detectors.ciri3 import Ciri3Detector
from circyto.paths import get_repo_root, get_tools_dir
from circyto.cli.doctor import doctor_app
from circyto.cli.detectors import detectors_app
from circyto.cli.smoke import smoke_app

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
        "  [ALIGN]  prepare-alignment-cache / plan-alignment-cache / align-manifest / run-detector-from-alignments\n"
        "  [MATRIX] collect-matrix (+ per-detector collectors)\n"
        "  [WORKFLOW] workflow smartseq3-ciri3 / full-length-circrna (experimental)\n"
        "  [ANNOTATE] annotate-circs\n"
        "  [INTEROP] prepare-scomatic-input, run-scomatic, import-scomatic, merge-scomatic\n"
        "  [INTEROP] export-scomatic-inputs, join-circ-snv-summary\n"
        "  [ANALYZE] analyze summarize-h5ad\n"
        "  [MERGE]  merge-detectors\n"
        "  [COMPARE] compare-ids (fuzzy/exact), compare-detectors (merged outputs)\n"
    ),
)
console = Console()


def _version_callback(value: bool) -> None:
    if not value:
        return
    typer.echo(__version__)
    raise typer.Exit()


@app.callback()
def app_callback(
    version: bool = typer.Option(
        False,
        "--version",
        help="Print the installed circyto version and exit.",
        callback=_version_callback,
        is_eager=True,
    ),
) -> None:
    """
    Top-level CLI callback for global options.
    """

app.add_typer(doctor_app, name="doctor")
app.add_typer(detectors_app, name="detectors")

app.add_typer(demux_app, name="demux")
app.add_typer(manifest_app, name="manifest")
app.add_typer(workflow_app, name="workflow")
app.add_typer(analyze_app, name="analyze")
app.add_typer(smoke_app, name="smoke")
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


@app.command("prepare-public-dataset")
def prepare_public_dataset_command(
    dataset_id: str = typer.Option(..., "--dataset-id", help="Public dataset identifier to plan."),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory for planning artifacts."),
    max_runs: Optional[int] = typer.Option(None, "--max-runs", help="Limit the number of selected runs."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Write planning artifacts only; do not attempt downloads."),
    protocol: str = typer.Option(
        ...,
        "--protocol",
        help="Expected protocol for the requested dataset: smartseq3, ramda, shin-ramda, or scrr.",
    ),
    download_method: str = typer.Option(
        "sra",
        "--download-method",
        help="How to render the download plan: sra, ena, or none.",
    ),
) -> None:
    """
    Prepare a lightweight public-dataset download plan without fetching large files.
    """
    try:
        summary = prepare_public_dataset(
            dataset_id=dataset_id,
            outdir=outdir,
            max_runs=max_runs,
            dry_run=dry_run,
            protocol=protocol,
            download_method=download_method,
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    warnings = list(summary.get("warnings", []))
    if dry_run and warnings:
        typer.echo("WARNING:")
        for warning in warnings:
            typer.echo(warning)
        typer.echo("")

    console.print(
        json.dumps(
            {
                "dataset_id": summary["dataset_id"],
                "protocol": summary["protocol"],
                "organism": summary["organism"],
                "expected_read_layout": summary["expected_read_layout"],
                "expected_reference": summary["expected_reference"],
                "recommended_route": summary["recommended_route"],
                "download_method": summary["download_method"],
                "dry_run": summary["dry_run"],
                "row_mode": summary["row_mode"],
                "warnings": summary["warnings"],
                "selected_run_count": summary["selected_run_count"],
                "selected_runs_tsv": str(summary["selected_runs_path"]),
                "download_plan_sh": str(summary["download_plan_path"]),
                "readme_next_steps_md": str(summary["readme_path"]),
            },
            indent=2,
        )
    )


@app.command("add-rna-profile")
def add_rna_profile_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory."),
    gtf: Path = typer.Option(..., "--gtf", exists=True, help="Annotation GTF for simple-overlap profiling."),
    method: str = typer.Option("simple-overlap", "--method", help="RNA profiling method. Currently only simple-overlap is supported."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Discover the existing alignment manifest and planned RNA outputs without writing files."),
) -> None:
    """
    Add a lightweight post-hoc RNA profile to an already completed workflow directory.
    """
    try:
        summary = add_posthoc_rna_profile(
            workdir=workdir,
            gtf_path=gtf,
            method=method,
            dry_run=dry_run,
        )
    except (FileNotFoundError, ValueError, NotImplementedError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("check-workflow")
def check_workflow_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Workflow directory to inspect."),
) -> None:
    """
    Read-only integrity check for a workflow directory.
    """
    summary = check_workflow_integrity(workdir)
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))
    if not summary.get("ok", False):
        raise typer.Exit(code=1)


@app.command("validate-workdir")
def validate_workdir_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory to validate for interoperability."),
    json_output: bool = typer.Option(False, "--json", help="Emit the full validation payload as JSON."),
) -> None:
    """
    Validate a completed workdir for scverse-facing interoperability and expected artifacts.
    """
    summary = validate_completed_workdir(workdir)
    if json_output:
        typer.echo(json.dumps(summary, indent=2, sort_keys=True))
    else:
        typer.echo(json.dumps(summary, indent=2, sort_keys=True))
    if not summary.get("ok", False):
        raise typer.Exit(code=1)


@app.command("print-environment")
def print_environment_command() -> None:
    """
    Print the current circyto and scverse package environment summary.
    """
    typer.echo(json.dumps(get_environment_summary(), indent=2, sort_keys=True))


@app.command("inspect-workdir")
def inspect_workdir_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory to inspect."),
    json_output: bool = typer.Option(False, "--json", help="Emit the full inspection payload as JSON."),
) -> None:
    """
    Read-only inspection of a completed workflow directory and its available modalities.
    """
    try:
        summary = inspect_completed_workdir(workdir)
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    if json_output:
        typer.echo(json.dumps(summary, indent=2, sort_keys=True))
        return
    typer.echo(
        "\n".join(
            [
                f"modalities={','.join(summary['available_modalities']) if summary['available_modalities'] else '-'}",
                f"workflow_type={summary['source_workflow_type']} protocol={summary['source_protocol']} read_layout={summary['source_read_layout']}",
                f"mudata_present={summary['mudata_present']} qc_present={summary['qc_present']}",
                f"matrices_present={json.dumps(summary['matrices_present'], sort_keys=True)}",
            ]
        )
    )


@app.command("cleanup-workflow")
def cleanup_workflow_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory to clean."),
    scope: str = typer.Option("alignments", "--scope", help="Cleanup scope: alignments, demux, or all."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Report planned deletions and estimated reclaimed bytes without deleting anything."),
    force: bool = typer.Option(False, "--force", help="Allow cleanup even if `circyto check-workflow` reports integrity problems."),
) -> None:
    """
    Clean regenerable workflow-owned intermediates from a completed workflow directory.
    """
    if scope not in SUPPORTED_CLEANUP_SCOPES:
        raise typer.BadParameter(
            f"Unsupported cleanup scope: {scope}. Choose from: {', '.join(SUPPORTED_CLEANUP_SCOPES)}"
        )
    try:
        summary = cleanup_completed_workflow(
            workdir=workdir,
            scope=scope,
            dry_run=dry_run,
            force=force,
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("refresh-rna-qc")
def refresh_rna_qc_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory with existing WORKDIR/rna outputs."),
) -> None:
    """
    Regenerate RNA QC summaries from existing RNA profile outputs without requiring alignments.
    """
    try:
        summary = refresh_rna_qc_from_existing_outputs(workdir=workdir)
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("summarize-rna-circ")
def summarize_rna_circ_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory with RNA and circ outputs."),
    write_summary: bool = typer.Option(False, "--write-summary", help="Write qc/rna_circ_cell_summary.tsv and qc/rna_circ_summary.json."),
    json_output: bool = typer.Option(False, "--json", help="Emit the full summary payload as JSON."),
) -> None:
    """
    Summarize RNA and circRNA overlap across cells in a completed workflow directory.
    """
    try:
        summary = summarize_rna_circ_integration(
            workdir=workdir,
            write_summary=write_summary,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    if json_output:
        typer.echo(json.dumps(summary, indent=2, sort_keys=True))
        return

    typer.echo(
        "\n".join(
            [
                f"RNA cells={summary['n_rna_cells']} circ cells={summary['n_circ_cells']} shared={summary['n_shared_cells']}",
                f"RNA-only={summary['n_rna_only_cells']} circ-only={summary['n_circ_only_cells']}",
                f"RNA-only cell IDs: {', '.join(summary['rna_only_cells']) if summary['rna_only_cells'] else '-'}",
                f"circ-only cell IDs: {', '.join(summary['circ_only_cells']) if summary['circ_only_cells'] else '-'}",
                f"Relationship: shared_cells={summary['rna_total_count_vs_circRNA_count_relationship']['shared_cells_considered']} "
                f"pearson={summary['rna_total_count_vs_circRNA_count_relationship']['pearson_correlation_total_rna_vs_circ_count']}",
            ]
        )
    )


@app.command("summarize-benchmark")
def summarize_benchmark_command(
    workdirs: list[Path] = typer.Option(..., "--workdirs", exists=True, file_okay=False, dir_okay=True, help="One or more completed workflow directories."),
    output: Optional[Path] = typer.Option(None, "--output", help="Output TSV path for the benchmark summary table."),
    json_output_path: Optional[Path] = typer.Option(None, "--json", help="Optional JSON summary path."),
) -> None:
    """
    Summarize completed workflow directories into a manuscript-facing benchmark table.
    """
    try:
        summary = summarize_benchmark_workdirs(workdirs=workdirs, output_tsv=output, output_json=json_output_path)
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("export-mudata")
def export_mudata_completed_workflow_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory with RNA and circ outputs."),
    output: Optional[Path] = typer.Option(None, "--output", help="Output h5mu path. Defaults to WORKDIR/mudata/full_length.h5mu."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite an existing h5mu output path."),
) -> None:
    """
    Export a completed RNA+circ workflow directory as a MuData multimodal bundle.
    """
    try:
        summary = export_completed_workflow_mudata(
            workdir=workdir,
            out_path=output,
            overwrite=overwrite,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("inspect-mudata")
def inspect_mudata_command(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="MuData .h5mu file to inspect."),
    json_output: bool = typer.Option(False, "--json", help="Emit the full inspection payload as JSON."),
) -> None:
    """
    Read-only structural inspection of a circyto MuData output.
    """
    try:
        summary = inspect_mudata_file(input)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    if json_output:
        typer.echo(json.dumps(summary, indent=2, sort_keys=True))
        return
    typer.echo(
        "\n".join(
            [
                f"modalities={','.join(summary['modalities']) if summary['modalities'] else '-'} n_obs={summary['n_obs']}",
                f"rna_shape={summary['rna_shape']} circ_shape={summary['circ_shape']}",
                f"obs_columns={', '.join(summary['obs_columns']) if summary['obs_columns'] else '-'}",
                f"rna_var_columns={', '.join(summary['rna_var_columns']) if summary['rna_var_columns'] else '-'}",
                f"circ_var_columns={', '.join(summary['circ_var_columns']) if summary['circ_var_columns'] else '-'}",
                f"circyto_uns_keys={', '.join(summary['circyto_uns_keys']) if summary['circyto_uns_keys'] else '-'}",
                f"membership: shared={summary['n_shared_cells']} rna_only={summary['n_rna_only_cells']} circ_only={summary['n_circ_only_cells']}",
            ]
        )
    )


@app.command("summarize-mudata-qc")
def summarize_mudata_qc_command(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="MuData .h5mu file to summarize."),
    json_output: bool = typer.Option(False, "--json", help="Emit the full QC summary payload as JSON."),
) -> None:
    """
    Read-only QC summary over a circyto MuData output.
    """
    try:
        summary = summarize_mudata_qc(input)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    if json_output:
        typer.echo(json.dumps(summary, indent=2, sort_keys=True))
        return
    lines = [f"n_obs={summary['n_obs']}"]
    for column in (
        "total_rna_counts",
        "detected_genes",
        "mitochondrial_fraction",
        "ribosomal_fraction",
        "circRNA_count",
        "circRNA_total_support",
    ):
        value = summary.get(column)
        if value is None:
            lines.append(f"{column}=NA")
        else:
            lines.append(
                f"{column}: min={value['min']} median={value['median']} mean={value['mean']} max={value['max']}"
            )
    lines.append(f"pearson_total_rna_vs_circRNA_count={summary['pearson_total_rna_vs_circRNA_count']}")
    typer.echo("\n".join(lines))


@app.command("import-dna-snv-summary")
def import_dna_snv_summary_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory."),
    dna_cell_summary: Path = typer.Option(..., "--dna-cell-summary", exists=True, dir_okay=False, help="dna_cell_summary.tsv"),
    dna_variant_summary: Optional[Path] = typer.Option(None, "--dna-variant-summary", exists=True, dir_okay=False, help="Optional dna_variant_summary.tsv"),
    scomatic_candidate_summary: Optional[Path] = typer.Option(None, "--scomatic-candidate-summary", exists=True, dir_okay=False, help="Optional scomatic_candidate_summary.tsv"),
) -> None:
    """
    Import lightweight DNA/CNV and RNA-derived candidate variant summaries into a completed workflow directory.
    """
    try:
        summary = import_dna_snv_summary(
            workdir=workdir,
            dna_cell_summary_path=dna_cell_summary,
            dna_variant_summary_path=dna_variant_summary,
            scomatic_candidate_summary_path=scomatic_candidate_summary,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("import-scrr-cnv")
def import_scrr_cnv_command(
    cnv_states: Path = typer.Option(..., "--cnv-states", exists=True, dir_okay=False, help="Processed scRR CNV states table, e.g. summary_CNV_states_bin50kb.txt.gz."),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False, help="Output directory for cnv.h5ad, cnv_cells.tsv, cnv_bins.tsv, and summary JSON."),
    cnv_mappabilitynorm: Optional[Path] = typer.Option(None, "--cnv-mappabilitynorm", exists=True, dir_okay=False, help="Optional processed scRR mappability-normalized CNV signal table matching --cnv-states bins."),
    cell_mapping: Optional[Path] = typer.Option(None, "--cell-mapping", exists=True, dir_okay=False, help="Optional TSV mapping DNA IDs to RNA/canonical cell IDs."),
    obs_id_strategy: str = typer.Option("canonical", "--obs-id-strategy", help="AnnData obs_names strategy: canonical, dna, or rna."),
    no_h5ad: bool = typer.Option(False, "--no-h5ad", help="Write only cnv_cells.tsv, cnv_bins.tsv, and summary JSON."),
) -> None:
    """
    Import processed scRR-seq DNA CNV summaries as a bin-level CNV modality.
    """
    try:
        summary = import_scrr_cnv(
            cnv_states_path=cnv_states,
            cnv_mappabilitynorm_path=cnv_mappabilitynorm,
            cell_mapping_path=cell_mapping,
            outdir=outdir,
            obs_id_strategy=obs_id_strategy,
            write_h5ad=not no_h5ad,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("build-scrr-cell-map")
def build_scrr_cell_map_command(
    soft: Path = typer.Option(..., "--soft", exists=True, dir_okay=False, help="GEO family.soft or family.soft.gz file."),
    out: Path = typer.Option(..., "--out", "-o", dir_okay=False, help="Output scRR GSM-to-cell mapping TSV."),
) -> None:
    """
    Build a scRR GSM -> RNA/DNA title -> canonical biological-cell mapping from GEO SOFT metadata.
    """
    try:
        summary = build_scrr_cell_map(soft_path=soft, output_path=out)
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("remap-scrr-mudata-obs")
def remap_scrr_mudata_obs_command(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="Existing RNA+circ MuData h5mu with GEO GSM obs IDs."),
    cell_map: Path = typer.Option(..., "--cell-map", exists=True, dir_okay=False, help="Mapping TSV from build-scrr-cell-map."),
    output: Path = typer.Option(..., "--output", "-o", dir_okay=False, help="Output remapped RNA+circ h5mu."),
    target_id: str = typer.Option("canonical_cell_id", "--target-id", help="Remapped obs ID target: canonical_cell_id or rna_cell_id."),
    allow_partial: bool = typer.Option(False, "--allow-partial", help="Allow obs IDs missing from the mapping and leave them unchanged."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite an existing output h5mu."),
) -> None:
    """
    Remap scRR RNA/circ MuData obs IDs from GSM accessions to biological cell IDs.
    """
    try:
        summary = remap_scrr_mudata_obs(
            input_h5mu=input,
            cell_map_path=cell_map,
            output_h5mu=output,
            target_id=target_id,
            allow_partial=allow_partial,
            overwrite=overwrite,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("merge-scrr-cnv")
def merge_scrr_cnv_command(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="Remapped RNA+circ MuData h5mu."),
    cnv: Path = typer.Option(..., "--cnv", exists=True, dir_okay=False, help="CNV AnnData h5ad from import-scrr-cnv."),
    output: Path = typer.Option(..., "--output", "-o", dir_okay=False, help="Output RNA+circ+CNV h5mu."),
    summary_json: Optional[Path] = typer.Option(None, "--summary-json", dir_okay=False, help="Output summary JSON. Defaults to OUTPUT with .summary.json suffix."),
    allow_partial: bool = typer.Option(False, "--allow-partial", help="Allow non-identical modality obs sets."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing output files."),
) -> None:
    """
    Merge remapped scRR RNA+circ MuData with a CNV AnnData modality.
    """
    try:
        summary = merge_scrr_cnv(
            input_h5mu=input,
            cnv_h5ad=cnv,
            output_h5mu=output,
            summary_json=summary_json,
            allow_partial=allow_partial,
            overwrite=overwrite,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("import-scrr-rt")
def import_scrr_rt_command(
    rt_table: Path = typer.Option(..., "--rt-table", exists=True, dir_okay=False, help="Processed scRR HAP1 replication timing/state table."),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False, help="Output directory for rt.h5ad, rt_cells.tsv, rt_features.tsv, and summary JSON."),
    avg_rt_bedgraph: Optional[Path] = typer.Option(None, "--avg-rt-bedgraph", exists=True, dir_okay=False, help="Optional processed average RT bedGraph. Stored as var['avg_rt'] only when coordinates match --rt-table features."),
    cell_mapping: Optional[Path] = typer.Option(None, "--cell-mapping", exists=True, dir_okay=False, help="Optional TSV mapping DNA/RNA IDs to canonical cell IDs."),
    obs_id_strategy: str = typer.Option("canonical", "--obs-id-strategy", help="AnnData obs_names strategy: canonical, dna, or rna."),
    no_h5ad: bool = typer.Option(False, "--no-h5ad", help="Write only rt_cells.tsv, rt_features.tsv, and summary JSON."),
) -> None:
    """
    Import processed scRR-seq DNA replication timing/state summaries as an RT modality.
    """
    try:
        summary = import_scrr_rt(
            rt_table_path=rt_table,
            avg_rt_bedgraph_path=avg_rt_bedgraph,
            cell_mapping_path=cell_mapping,
            outdir=outdir,
            obs_id_strategy=obs_id_strategy,
            write_h5ad=not no_h5ad,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("merge-scrr-rt")
def merge_scrr_rt_command(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="Remapped RNA+circ MuData h5mu."),
    rt: Path = typer.Option(..., "--rt", exists=True, dir_okay=False, help="Replication timing/state AnnData h5ad from import-scrr-rt."),
    output: Path = typer.Option(..., "--output", "-o", dir_okay=False, help="Output RNA+circ+RT h5mu."),
    summary_json: Optional[Path] = typer.Option(None, "--summary-json", dir_okay=False, help="Output summary JSON. Defaults to OUTPUT with .summary.json suffix."),
    allow_partial: bool = typer.Option(False, "--allow-partial", help="Allow non-identical modality obs sets."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing output files."),
) -> None:
    """
    Merge remapped scRR RNA+circ MuData with a replication timing/state AnnData modality.
    """
    try:
        summary = merge_scrr_rt(
            input_h5mu=input,
            rt_h5ad=rt,
            output_h5mu=output,
            summary_json=summary_json,
            allow_partial=allow_partial,
            overwrite=overwrite,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("prepare-scomatic-input")
def prepare_scomatic_input_command(
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False, help="Output directory for SComatic-ready inputs and provenance."),
    protocol: str = typer.Option(..., "--protocol", help="Full-length RNA protocol: smartseq2, smartseq3, ramda, shinramda, or scrr-rna."),
    alignment_manifest: Optional[Path] = typer.Option(None, "--alignment-manifest", exists=True, dir_okay=False, help="One-cell-per-BAM/SAM alignment manifest with cell_id and bam or sam columns."),
    merged_bam: Optional[Path] = typer.Option(None, "--merged-bam", exists=True, dir_okay=False, help="Already merged multi-cell BAM with CB tags."),
    cell_metadata: Optional[Path] = typer.Option(None, "--cell-metadata", exists=True, dir_okay=False, help="Optional cell metadata TSV. Required for --merged-bam mode."),
    reference_fasta: Optional[Path] = typer.Option(None, "--reference-fasta", exists=True, dir_okay=False, help="Reference FASTA used for provenance and later SComatic runs."),
    sample_id: str = typer.Option("scomatic", "--sample-id", help="Sample prefix for merged SComatic BAM outputs."),
    threads: int = typer.Option(1, "--threads", help="Threads for samtools sort/merge."),
    cell_type_column: Optional[str] = typer.Option(None, "--cell-type-column", help="Manifest or metadata column to use as SComatic Cell_type."),
    default_cell_type: Optional[str] = typer.Option(None, "--default-cell-type", help="Fallback SComatic Cell_type when no column value is available."),
    samtools: str = typer.Option("samtools", "--samtools", help="samtools executable for one-cell-per-BAM/SAM preparation."),
) -> None:
    """
    Prepare full-length RNA alignments and metadata for external SComatic interoperability.
    """
    try:
        summary = prepare_scomatic_input(
            outdir=outdir,
            protocol=protocol,
            alignment_manifest_path=alignment_manifest,
            merged_bam_path=merged_bam,
            cell_metadata_path=cell_metadata,
            reference_fasta_path=reference_fasta,
            sample_id=sample_id,
            threads=threads,
            cell_type_column=cell_type_column,
            default_cell_type=default_cell_type,
            samtools=samtools,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("run-scomatic")
def run_scomatic_command(
    prepared_dir: Path = typer.Option(..., "--prepared-dir", exists=True, file_okay=False, help="Output directory from prepare-scomatic-input."),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False, help="Output directory for SComatic command plan, logs, and summary."),
    scomatic_dir: Path = typer.Option(..., "--scomatic-dir", exists=True, file_okay=False, help="SComatic checkout directory."),
    reference_fasta: Path = typer.Option(..., "--reference-fasta", exists=True, dir_okay=False, help="Reference FASTA for SComatic."),
    threads: int = typer.Option(1, "--threads", help="Thread count recorded in the BaseCellCounter command plan."),
    python_executable: str = typer.Option("python", "--python-executable", help="Python executable for SComatic scripts."),
    basecellcounter_script: Optional[Path] = typer.Option(None, "--basecellcounter-script", dir_okay=False, help="Override BaseCellCounter.py path."),
    step1_script: Optional[Path] = typer.Option(None, "--step1-script", dir_okay=False, help="Override BaseCellCalling.step1.py path."),
    step2_script: Optional[Path] = typer.Option(None, "--step2-script", dir_okay=False, help="Override BaseCellCalling.step2.py path."),
    basecellcounter_args: Optional[str] = typer.Option(None, "--basecellcounter-args", help="Shell-style BaseCellCounter args. Supports placeholders like {merged_bam}, {reference_fasta}, {basecellcounter_out}, {threads}."),
    step1_args: Optional[str] = typer.Option(None, "--step1-args", help="Shell-style Step1 args. Supports placeholders like {basecellcounter_out} and {step1_out}."),
    step2_args: Optional[str] = typer.Option(None, "--step2-args", help="Shell-style Step2 args. Supports placeholders like {step1_out} and {step2_out}."),
    run_step1: bool = typer.Option(False, "--run-step1", help="Include SComatic Step1 in the command plan."),
    run_step2: bool = typer.Option(False, "--run-step2", help="Include SComatic Step2 in the command plan."),
    execute: bool = typer.Option(False, "--execute", help="Actually execute the generated external SComatic commands."),
) -> None:
    """
    Generate, and optionally execute, an external SComatic command plan.
    """
    try:
        summary = run_scomatic(
            prepared_dir=prepared_dir,
            outdir=outdir,
            scomatic_dir=scomatic_dir,
            reference_fasta_path=reference_fasta,
            threads=threads,
            python_executable=python_executable,
            basecellcounter_script=basecellcounter_script,
            step1_script=step1_script,
            step2_script=step2_script,
            basecellcounter_args=basecellcounter_args,
            step1_args=step1_args,
            step2_args=step2_args,
            run_step1=run_step1,
            run_step2=run_step2,
            execute=execute,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("import-scomatic")
def import_scomatic_command(
    scomatic_output: List[Path] = typer.Option(..., "--scomatic-output", exists=True, dir_okay=False, help="SComatic Step1/Step2 or candidate TSV/CSV output. Repeat for multiple files."),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False, help="Directory for scomatic_candidate_summary.tsv and import summaries."),
    cell_annotations: Optional[Path] = typer.Option(None, "--cell-annotations", exists=True, dir_okay=False, help="Optional SComatic Index/Cell_type annotation table."),
    provenance_metadata: Optional[Path] = typer.Option(None, "--provenance-metadata", exists=True, dir_okay=False, help="Optional JSON or text SComatic run metadata."),
) -> None:
    """
    Import external SComatic Step1/Step2 outputs into circyto's candidate-signal schema.
    """
    try:
        summary = import_scomatic(
            scomatic_output_paths=scomatic_output,
            outdir=outdir,
            cell_annotation_path=cell_annotations,
            provenance_metadata_path=provenance_metadata,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("merge-scomatic")
def merge_scomatic_command(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="Input RNA/circ or multimodal MuData h5mu."),
    scomatic_candidates: Path = typer.Option(..., "--scomatic-candidates", exists=True, dir_okay=False, help="scomatic_candidate_summary.tsv from import-scomatic or normalize-scomatic-results."),
    output: Path = typer.Option(..., "--output", "-o", dir_okay=False, help="Output MuData h5mu with candidate_snv modality."),
    summary_json: Optional[Path] = typer.Option(None, "--summary-json", dir_okay=False, help="Output summary JSON. Defaults to OUTPUT with .summary.json suffix."),
    allow_partial: bool = typer.Option(False, "--allow-partial", help="Allow candidate cells absent from input MuData obs, useful for cell-type-level Step1/Step2 tables."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing output files."),
) -> None:
    """
    Merge RNA-derived SComatic candidate signals into MuData as candidate_snv.
    """
    try:
        summary = merge_scomatic(
            input_h5mu=input,
            scomatic_candidates_path=scomatic_candidates,
            output_h5mu=output,
            summary_json=summary_json,
            allow_partial=allow_partial,
            overwrite=overwrite,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("normalize-scomatic-results")
def normalize_scomatic_results_command(
    scomatic_output: List[Path] = typer.Option(
        ...,
        "--scomatic-output",
        exists=True,
        dir_okay=False,
        help="External SComatic TSV/CSV output file. Repeat this option for multiple files.",
    ),
    outdir: Path = typer.Option(
        ...,
        "--outdir",
        "-o",
        file_okay=False,
        help="Directory for scomatic_candidate_summary.tsv and normalize_scomatic_results_summary.json.",
    ),
    cell_annotations: Optional[Path] = typer.Option(
        None,
        "--cell-annotations",
        exists=True,
        dir_okay=False,
        help="Optional SComatic cell annotation table with Index and Cell_type columns.",
    ),
    provenance_metadata: Optional[Path] = typer.Option(
        None,
        "--provenance-metadata",
        exists=True,
        dir_okay=False,
        help="Optional JSON or text provenance metadata to preserve in the normalization summary.",
    ),
) -> None:
    """
    Normalize external SComatic outputs into circyto's candidate-signal schema without running SComatic.
    """
    try:
        summary = normalize_scomatic_results(
            scomatic_output_paths=scomatic_output,
            outdir=outdir,
            cell_annotation_path=cell_annotations,
            provenance_metadata_path=provenance_metadata,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("summarize-dna-rna-circ")
def summarize_dna_rna_circ_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory with RNA/circ and imported DNA summaries."),
    write_summary: bool = typer.Option(False, "--write-summary", help="Write qc/dna_rna_circ_cell_summary.tsv and qc/dna_rna_circ_summary.json."),
) -> None:
    """
    Summarize joined DNA, RNA, and circRNA per-cell summaries for scRR-style workflows.
    """
    try:
        summary = summarize_dna_rna_circ(workdir=workdir, write_summary=write_summary)
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("summarize-circ-host-genes")
def summarize_circ_host_genes_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory with circ matrix outputs."),
    output: Optional[Path] = typer.Option(None, "--output", help="Output TSV path. Defaults to WORKDIR/qc/circ_host_gene_summary.tsv."),
    json_output: bool = typer.Option(False, "--json", help="Emit the full summary payload as JSON."),
) -> None:
    """
    Summarize circRNA host-gene recurrence and support across a completed workflow directory.
    """
    try:
        summary = summarize_circ_host_genes(workdir=workdir, output_path=output)
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    if json_output:
        typer.echo(json.dumps(summary, indent=2, sort_keys=True))
        return
    typer.echo(
        "\n".join(
            [
                f"n_host_genes={summary['n_host_genes']} n_circ_with_host_gene={summary['n_circ_with_host_gene']}",
                f"output={summary['output_path']}",
            ]
        )
    )


@app.command("export-circ-bed")
def export_circ_bed_command(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False, dir_okay=True, help="Completed workflow directory with circ matrix outputs."),
    output: Optional[Path] = typer.Option(None, "--output", help="Output BED-like path. Defaults to WORKDIR/qc/circs.bed."),
) -> None:
    """
    Export circRNA coordinates and total support in a BED-like format for genome-browser interoperability.
    """
    try:
        summary = export_circ_bed(workdir=workdir, output_path=output)
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("scanpy-qc-report")
def scanpy_qc_report_command(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="MuData .h5mu file to analyze."),
    output_dir: Path = typer.Option(..., "--output-dir", file_okay=False, help="Directory for exploratory Scanpy QC outputs."),
) -> None:
    """
    Exploratory Scanpy QC reporting over the RNA modality of a circyto MuData file.
    """
    try:
        summary = scanpy_qc_report(input_path=input, output_dir=output_dir)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("scanpy-pca")
def scanpy_pca_command(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="MuData .h5mu file to analyze."),
    output_dir: Path = typer.Option(..., "--output-dir", file_okay=False, help="Directory for exploratory Scanpy PCA/UMAP outputs."),
) -> None:
    """
    Exploratory Scanpy PCA/UMAP/Leiden workflow over the RNA modality of a circyto MuData file.
    """
    try:
        summary = scanpy_pca_workflow(input_path=input, output_dir=output_dir)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


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
    from circyto.writers.convert import convert_matrix_files

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
    from circyto.writers.convert import convert_matrix_files

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
    from circyto.pipeline.export_multimodal import export_multimodal as _export_multimodal

    _export_multimodal(
        genes_h5ad=genes_h5ad,
        circ_matrix=circ_matrix,
        circ_index=circ_index,
        cell_index=cell_index,
        out=out,
        circ_feature_table=circ_feature_table,
    )


@app.command("export-scomatic-inputs")
def export_scomatic_inputs_cmd(
    bam_manifest: Path = typer.Option(..., "--bam-manifest", exists=True, help="TSV with per-cell BAM paths"),
    cell_metadata: Path = typer.Option(..., "--cell-metadata", exists=True, help="TSV with per-cell metadata"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="Output directory"),
    reference_fasta: Path = typer.Option(..., "--reference-fasta", exists=True, help="Reference FASTA for external SComatic runs"),
    protocol: str = typer.Option(..., "--protocol", help="One of: smartseq3, ramda, shin-ramda"),
) -> None:
    """
    Export lightweight interoperability tables for an external SComatic run.
    """
    try:
        normalized_protocol = normalize_full_length_protocol(protocol)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    export_scomatic_inputs(
        bam_manifest=bam_manifest,
        cell_metadata=cell_metadata,
        outdir=outdir,
        reference_fasta=reference_fasta,
        protocol=normalized_protocol,
    )
    typer.echo(f"Wrote SComatic interoperability scaffold to {outdir}")


@app.command("join-circ-snv-summary")
def join_circ_snv_summary_cmd(
    circ_matrix: Path = typer.Option(..., "--circ-matrix", exists=True, help="Wide circRNA-by-cell TSV"),
    circ_feature_table: Path = typer.Option(..., "--circ-feature-table", exists=True, help="circ_feature_table.tsv or similar"),
    scomatic_candidates: Path = typer.Option(..., "--scomatic-candidates", exists=True, help="Exploratory SComatic candidate SNV TSV"),
    cell_metadata: Path = typer.Option(..., "--cell-metadata", exists=True, help="TSV with per-cell metadata"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="Output directory"),
) -> None:
    """
    Join circRNA summaries with exploratory SComatic candidate SNV tables.
    """
    join_circ_snv_summary(
        circ_matrix=circ_matrix,
        circ_feature_table=circ_feature_table,
        scomatic_candidates=scomatic_candidates,
        cell_metadata=cell_metadata,
        outdir=outdir,
    )
    typer.echo(f"Wrote circ/SNV summary tables to {outdir}")


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


@app.command("repair-host-genes")
def repair_host_genes_cmd(
    input: Path = typer.Option(..., "--input", exists=True, dir_okay=False, help="Input circRNA .h5ad or MuData .h5mu file."),
    output: Path = typer.Option(..., "--output", "-o", dir_okay=False, help="Output repaired .h5ad or .h5mu file."),
    circ_mod: str = typer.Option("circ", "--circ-mod", help="CircRNA modality name for .h5mu inputs."),
    gtf: Optional[Path] = typer.Option(None, "--gtf", exists=True, dir_okay=False, help="Optional GTF/GFF for primary coordinate-based host-gene annotation."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite an existing output path."),
) -> None:
    """
    Repair host-gene provenance columns in existing circRNA h5ad or h5mu outputs.
    """
    try:
        summary = repair_host_gene_file(
            input_path=input,
            output_path=output,
            circ_mod=circ_mod,
            gtf_path=gtf,
            overwrite=overwrite,
        )
    except (FileNotFoundError, ValueError, RuntimeError, ImportError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


@app.command("annotate-circs")
def annotate_circs_cmd(
    circ_table: Path = typer.Option(..., "--circ-table", exists=True, help="Input circ_qc.tsv-like table"),
    annotation_db: List[str] = typer.Option(
        [],
        "--annotation-db",
        help=(
            "Repeated DB spec. Use 'name:path' for canonical columns "
            "(chrom,start,end,strand,id,host_gene) or "
            "'name=...;path=...;chrom=...;start=...;end=...;strand=...;id=...;host_gene=...;extra=col1,col2'"
        ),
    ),
    out: Path = typer.Option(..., "--out", help="Output annotated circ table TSV"),
    summary_out: Optional[Path] = typer.Option(None, "--summary-out", help="Optional output path for annotation_summary.json"),
    max_bsj_distance: int = typer.Option(0, "--max-bsj-distance", help="Maximum absolute start/end distance for near-BSJ matches"),
    enable_host_gene_match: bool = typer.Option(
        False,
        "--enable-host-gene-match/--no-enable-host-gene-match",
        help="Allow lower-confidence host-gene plus overlapping-locus matches",
    ),
    update_h5ad: Optional[Path] = typer.Option(None, "--update-h5ad", exists=True, help="Optional circ_counts.h5ad to update in-place"),
) -> None:
    """
    Annotate circRNAs against one or more known circRNA databases.
    """
    if not annotation_db:
        raise typer.BadParameter("Provide at least one --annotation-db spec.")
    try:
        db_specs = [parse_annotation_db_spec(item) for item in annotation_db]
    except Exception as exc:
        raise typer.BadParameter(str(exc)) from exc

    summary = annotate_circ_table(
        circ_table=circ_table,
        db_specs=db_specs,
        out_path=out,
        summary_path=summary_out,
        max_bsj_distance=max_bsj_distance,
        enable_host_gene_match=enable_host_gene_match,
        update_h5ad=update_h5ad,
    )
    typer.echo(f"Annotated circ table: {out}")
    if summary_out is None:
        summary_out = out.with_name("annotation_summary.json")
    typer.echo(f"Annotation summary: {summary_out}")
    typer.echo(
        f"Status counts: novel={summary['novel_count']} total={summary['total_circRNAs']}"
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


def _get_detector_engine(detector: str):
    engines = build_default_engines()
    if detector not in engines:
        available = ", ".join(sorted(engines.keys()))
        raise typer.Exit(f"Detector '{detector}' not available. Available: {available}")
    return engines[detector]


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


@app.command("run-ciri3")
def run_ciri3_cmd(
    manifest: Path = typer.Option(..., "--manifest", exists=True, help="Manifest TSV with FASTQ rows"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="Workflow output directory"),
    genome_fasta: Path = typer.Option(..., "--genome-fasta", exists=True, help="Reference FASTA"),
    gtf: Path = typer.Option(..., "--gtf", exists=True, help="Annotation GTF"),
    star_index: Path | None = typer.Option(None, "--star-index", exists=True, help="STAR genomeDir for paired-end rows"),
    protocol: str | None = typer.Option(
        None,
        "--protocol",
        help="Protocol preset: smartseq3, ramda, shin-ramda. If omitted, use manifest protocol/platform.",
    ),
    strandedness: str | None = typer.Option(
        None,
        "--strandedness",
        help="Explicit strandedness override/preset metadata: forward, reverse, or unstranded.",
    ),
    threads: int = typer.Option(8, "--threads", help="Threads per alignment/detector task"),
    parallel: int = typer.Option(1, "--parallel", help="Concurrent alignment/detector tasks"),
    chunk_size: int = typer.Option(1, "--chunk-size", help="Chunk size for resumable execution"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print planned STAR/BWA and CIRI3 commands without executing"),
    fail_fast: bool = typer.Option(False, "--fail-fast", help="Stop after the first failed chunk"),
    command_template: str | None = typer.Option(None, "--command-template", help="Optional CIRI3 command template override"),
) -> None:
    """
    Run a protocol-aware CIRI3 alignment-first workflow.

    Paired-end rows use STAR+CIRI3 hybrid preparation.
    Single-end rows use BWA+CIRI3 while preserving the same output layout.
    """
    run_ciri3_workflow(
        RunCiri3Params(
            manifest=manifest,
            outdir=outdir,
            genome_fasta=genome_fasta,
            gtf=gtf,
            star_index=star_index,
            protocol=protocol,
            strandedness=strandedness,
            threads=threads,
            parallel=parallel,
            chunk_size=chunk_size,
            dry_run=dry_run,
            fail_fast=fail_fast,
            command_template=command_template,
        ),
        progress=typer.echo,
    )


@app.command("prepare-alignment-cache")
def prepare_alignment_cache_cmd(
    manifest: Path = typer.Option(..., exists=True, help="Input manifest TSV with FASTQ and/or BAM columns"),
    outdir: Optional[Path] = typer.Option(None, "--outdir", "-o", help="Output directory for cached alignments"),
    aligner: str = typer.Option(
        "reuse-existing",
        "--aligner",
        help="Alignment strategy. For --aligner star with --detector ciri3, circyto runs the official STAR+BWA hybrid prep and records aligned.sam + chimeric junctions + unmapped mates + bwa rescue SAM.",
    ),
    detector: Optional[str] = typer.Option(None, "--detector", "-d", help="Optional detector hint for cache keying"),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa", help="Reference FASTA for alignment"),
    threads: int = typer.Option(8, "--threads", help="Threads per alignment task"),
    parallel: int = typer.Option(4, "--parallel", help="Concurrent alignment tasks"),
    sentinel_cells: int = typer.Option(0, "--sentinel-cells", help="Run the first N cells before the rest"),
    chunk_size: int = typer.Option(25, "--chunk-size", help="Chunk size for resumable alignment execution"),
    command_template: Optional[str] = typer.Option(
        None,
        "--command-template",
        help="Shell template for custom aligners; placeholders: {cell_id} {read1} {read2} {ref_fa} {out_path} {threads} {extra_flags}",
    ),
    extra_flags: str = typer.Option("", "--extra-flags", help="Extra aligner flags included in cache key"),
    link_mode: str = typer.Option("symlink", "--link-mode", help="symlink or copy staged alignments into the working directory"),
    index_bam: bool = typer.Option(False, "--index-bam", help="Index BAM outputs when samtools is available"),
    output_format: str = typer.Option("bam", "--output-format", help="bam or sam"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Preflight only; write a plan JSON without running alignments"),
    fail_fast: bool = typer.Option(False, "--fail-fast", help="Stop after the first failed chunk instead of continuing"),
) -> None:
    """
    Prepare reusable alignment artifacts from a source manifest.

    This is the first stage of the alignment-first execution track.

    Detector-specific note:
    - `--aligner star --detector ciri3` uses the official CIRI3 hybrid contract:
      STAR writes `Aligned.out.sam`, `Chimeric.out.junction`, and unmapped mates,
      then circyto runs `bwa mem -T 19` on the unmapped mates and records `bwa_sam`.
    """
    if outdir is None:
        outdir = _auto_outdir("prepare-alignment-cache", detector or aligner, manifest.stem)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

    outdir.mkdir(parents=True, exist_ok=True)
    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=outdir,
        aligner=aligner,
        ref_fa=ref_fa,
        detector_hint=detector,
        threads=threads,
        parallel=parallel,
        sentinel_cells=sentinel_cells,
        chunk_size=chunk_size,
        command_template=command_template,
        extra_flags=extra_flags,
        link_mode=link_mode,
        index_bam=index_bam,
        output_format=output_format,
        dry_run=dry_run,
        fail_fast=fail_fast,
    )
    if dry_run:
        typer.echo(f"[prepare-alignment-cache] Wrote plan: {outdir/'alignment_prepare_plan.json'}")
    else:
        typer.echo(f"[prepare-alignment-cache] Wrote alignment manifest: {alignment_manifest}")


@app.command("align-manifest")
def align_manifest_cmd(
    manifest: Path = typer.Option(..., exists=True, help="Input manifest TSV with FASTQ and/or BAM columns"),
    outdir: Optional[Path] = typer.Option(None, "--outdir", "-o", help="Output directory"),
    aligner: str = typer.Option("reuse-existing", "--aligner"),
    detector: Optional[str] = typer.Option(None, "--detector", "-d"),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa"),
    threads: int = typer.Option(8, "--threads"),
    parallel: int = typer.Option(4, "--parallel"),
    sentinel_cells: int = typer.Option(0, "--sentinel-cells"),
    chunk_size: int = typer.Option(25, "--chunk-size"),
    command_template: Optional[str] = typer.Option(None, "--command-template"),
    extra_flags: str = typer.Option("", "--extra-flags"),
    output_format: str = typer.Option("bam", "--output-format"),
    dry_run: bool = typer.Option(False, "--dry-run"),
    fail_fast: bool = typer.Option(False, "--fail-fast"),
) -> None:
    """
    Alias for prepare-alignment-cache focused on manifest generation.
    """
    if outdir is None:
        outdir = _auto_outdir("align-manifest", detector or aligner, manifest.stem)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=outdir,
        aligner=aligner,
        ref_fa=ref_fa,
        detector_hint=detector,
        threads=threads,
        parallel=parallel,
        sentinel_cells=sentinel_cells,
        chunk_size=chunk_size,
        command_template=command_template,
        extra_flags=extra_flags,
        output_format=output_format,
        dry_run=dry_run,
        fail_fast=fail_fast,
    )
    if dry_run:
        typer.echo(outdir / "alignment_prepare_plan.json")
    else:
        typer.echo(alignment_manifest)


@app.command("plan-alignment-cache")
def plan_alignment_cache_cmd(
    manifest: Path = typer.Option(..., exists=True, help="Input manifest TSV with FASTQ and/or BAM columns"),
    outdir: Optional[Path] = typer.Option(None, "--outdir", "-o", help="Output directory"),
    aligner: str = typer.Option("reuse-existing", "--aligner"),
    detector: Optional[str] = typer.Option(None, "--detector", "-d"),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa"),
    threads: int = typer.Option(8, "--threads"),
    parallel: int = typer.Option(4, "--parallel"),
    sentinel_cells: int = typer.Option(0, "--sentinel-cells"),
    chunk_size: int = typer.Option(25, "--chunk-size"),
    command_template: Optional[str] = typer.Option(None, "--command-template"),
    extra_flags: str = typer.Option("", "--extra-flags"),
    output_format: str = typer.Option("bam", "--output-format"),
    preview_rows: int = typer.Option(3, "--preview-rows", help="Show exact first N commands in the plan"),
) -> None:
    """
    Print and persist a preflight plan for alignment cache preparation.
    """
    if outdir is None:
        outdir = _auto_outdir("plan-alignment-cache", detector or aligner, manifest.stem)
    payload = plan_alignment_cache(
        manifest=manifest,
        outdir=outdir,
        aligner=aligner,
        ref_fa=ref_fa,
        detector_hint=detector,
        threads=threads,
        parallel=parallel,
        sentinel_cells=sentinel_cells,
        chunk_size=chunk_size,
        command_template=command_template,
        extra_flags=extra_flags,
        output_format=output_format,
        preview_rows=preview_rows,
    )
    outdir.mkdir(parents=True, exist_ok=True)
    plan_path = outdir / "alignment_prepare_plan.json"
    plan_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    typer.echo(f"rows={payload['n_rows']} chunks={payload['chunk_count']} errors={len(payload['errors'])}")
    typer.echo(plan_path)


@app.command("run-detector-from-alignments")
def run_detector_from_alignments_cmd(
    detector: str = typer.Option(..., "--detector", "-d", help="Alignment-capable detector name"),
    manifest: Path = typer.Option(..., exists=True, help="Alignment manifest TSV"),
    outdir: Optional[Path] = typer.Option(None, "--outdir", "-o", help="Detector output directory"),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa", help="Reference FASTA"),
    gtf: Optional[Path] = typer.Option(None, "--gtf", help="Annotation GTF/GFF"),
    command_template: Optional[str] = typer.Option(
        None,
        "--command-template",
        help="Optional detector command template override. For ciri3, an explicit template always overrides direct java -jar execution.",
    ),
    threads: int = typer.Option(8, "--threads", help="Threads per detector process"),
    parallel: int = typer.Option(4, "--parallel", help="Number of alignment rows to run in parallel"),
    chunk_size: int = typer.Option(50, "--chunk-size", help="Chunk size for detector execution"),
    sentinel_cells: int = typer.Option(0, "--sentinel-cells", help="Run the first N cells before the rest"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Preflight only; write a detector plan JSON"),
    fail_fast: bool = typer.Option(False, "--fail-fast", help="Stop after the first failed chunk instead of continuing"),
) -> None:
    """
    Run an alignment-capable detector over a reusable alignment manifest.

    For CIRI3 STAR rows (`mapper_mode=1`), the manifest must include:
    - STAR aligned SAM
    - `chimeric_junction`
    - `bwa_sam` rescue alignment
    """
    det_engine = _get_detector_engine(detector)
    if detector == "ciri3" and command_template is not None and isinstance(det_engine, Ciri3Detector):
        det_engine.command_template = command_template
    caps = get_detector_capabilities(det_engine)
    if not caps.accepts_alignment:
        raise typer.BadParameter(
            f"Detector '{detector}' does not accept alignment inputs. Use run-detector with FASTQ manifests instead."
        )
    if outdir is None:
        outdir = _auto_outdir("run-detector-from-alignments", detector, manifest.stem)
        console.print(f"[yellow]--outdir not provided; using[/yellow] {outdir}")

    results = run_detector_alignment_manifest(
        detector=det_engine,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
        parallel=parallel,
        chunk_size=chunk_size,
        sentinel_cells=sentinel_cells,
        dry_run=dry_run,
        fail_fast=fail_fast,
    )
    if dry_run:
        typer.echo(f"[run-detector-from-alignments] Wrote plan: {outdir/'detector_alignment_plan.json'}")
    else:
        typer.echo(f"[run-detector-from-alignments] Completed {len(results)} jobs into {outdir}")


@app.command("validate-ciri3-template")
def validate_ciri3_template_cmd(
    template: Optional[str] = typer.Option(None, "--template", help="CIRI3 command template to validate"),
) -> None:
    """
    Validate the configured CIRI3 command template before a large cluster run.
    """
    det = Ciri3Detector(command_template=template) if template else Ciri3Detector()
    ok, errors, details = det.validate_runtime(template=template)
    if ok:
        typer.echo("OK")
        typer.echo(str(details))
        return
    for err in errors:
        typer.echo(f"ERROR: {err}", err=True)
    typer.echo(str(details), err=True)
    raise typer.Exit(code=1)


@app.command("alignment-first-smoke")
def alignment_first_smoke_cmd(
    outdir: Path = typer.Option(Path("work/alignment_first_smoke"), "--outdir", "-o"),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa", help="Reference FASTA; defaults to repo ref/chr21.fa when present"),
    r1: Optional[Path] = typer.Option(None, "--r1", help="Override paired-end smoke read1 FASTQ"),
    r2: Optional[Path] = typer.Option(None, "--r2", help="Override paired-end smoke read2 FASTQ"),
    threads: int = typer.Option(2, "--threads"),
    parallel: int = typer.Option(1, "--parallel"),
    chunk_size: int = typer.Option(1, "--chunk-size"),
    dry_run: bool = typer.Option(False, "--dry-run"),
) -> None:
    """
    Run a small local alignment-first smoke path using repo-shipped FASTQs plus a bundled CIRI3-compatible template.

    This validates the workflow plumbing, not biological CIRI3 correctness.
    """
    repo_root = get_repo_root()
    tools_dir = get_tools_dir()
    if repo_root is None or tools_dir is None:
        raise typer.BadParameter("Could not resolve repo-local smoke assets.")
    default_r1 = tools_dir / "CIRI-full_v2.0" / "CIRI-full_test" / "test_1.fq.gz"
    default_r2 = tools_dir / "CIRI-full_v2.0" / "CIRI-full_test" / "test_2.fq.gz"
    default_ref = repo_root / "ref" / "chr21.fa"
    smoke_template = tools_dir / "ciri3_smoke_template.sh"
    r1 = r1 or default_r1
    r2 = r2 or default_r2
    ref_fa = ref_fa or default_ref
    if not r1.exists() or not r2.exists():
        raise typer.BadParameter(f"Smoke FASTQs not found: {r1} {r2}")
    if not ref_fa.exists():
        raise typer.BadParameter(f"Reference FASTA not found: {ref_fa}")
    if not smoke_template.exists():
        raise typer.BadParameter(f"Smoke template not found: {smoke_template}")

    outdir.mkdir(parents=True, exist_ok=True)
    manifest = outdir / "smoke_manifest.tsv"
    manifest.write_text(
        "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id\n"
        f"smoke1\tsmartseq2\t{r1}\t{r2}\t\talignment-first-smoke\t0\tsmoke\n",
        encoding="utf-8",
    )
    prepare_dir = outdir / "align"
    detector_dir = outdir / "ciri3"
    matrix_dir = outdir / "matrix"
    template = (
        f"bash {smoke_template} {{alignment}} {{raw_output}} {{cell_id}}"
        " # {alignment_format} {threads} {outdir}"
    )
    if dry_run:
        prepare_alignment_cache(
            manifest=manifest,
            outdir=prepare_dir,
            aligner="bwa-mem",
            ref_fa=ref_fa,
            detector_hint="ciri3",
            threads=threads,
            parallel=parallel,
            chunk_size=chunk_size,
            dry_run=True,
        )
        det = Ciri3Detector(command_template=template)
        run_detector_alignment_manifest(
            detector=det,
            manifest=prepare_dir / "alignment_manifest.tsv",
            outdir=detector_dir,
            ref_fa=ref_fa,
            threads=threads,
            parallel=parallel,
            chunk_size=chunk_size,
            dry_run=True,
        )
        typer.echo(f"[alignment-first-smoke] Wrote plans under {outdir}")
        return

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=prepare_dir,
        aligner="bwa-mem",
        ref_fa=ref_fa,
        detector_hint="ciri3",
        threads=threads,
        parallel=parallel,
        chunk_size=chunk_size,
        index_bam=True,
    )
    det = Ciri3Detector(command_template=template)
    run_detector_alignment_manifest(
        detector=det,
        manifest=alignment_manifest,
        outdir=detector_dir,
        ref_fa=ref_fa,
        threads=threads,
        parallel=parallel,
        chunk_size=chunk_size,
    )
    collect_matrix_cmd(
        detector="ciri3",
        indir=detector_dir,
        outdir=matrix_dir,
        matrix=None,
        circ_index=None,
        cell_index=None,
        min_count_per_cell=1,
    )
    typer.echo(f"[alignment-first-smoke] Completed under {outdir}")


@app.command("summarize-chunks")
def summarize_chunks_cmd(
    indir: Path = typer.Option(..., "--indir", exists=True, help="Alignment or detector output directory with chunks/"),
    json_out: bool = typer.Option(False, "--json", help="Emit JSON"),
) -> None:
    """
    Summarize per-chunk checkpoint files for cluster logs and resume inspection.
    """
    payload = summarize_alignment_chunks(indir)
    if json_out:
        typer.echo(json.dumps(payload, indent=2, sort_keys=True))
        return
    typer.echo(f"chunks={payload['n_chunks']} status_counts={payload['status_counts']}")
    if payload["failed_chunks"]:
        typer.echo(f"failed_chunks={payload['failed_chunks']}")
    if payload["failed_cells"]:
        typer.echo(f"failed_cells={payload['failed_cells']}")
    for chunk in payload["chunks"]:
        typer.echo(
            f"chunk={chunk['chunk_index']} status={chunk['status']} size={chunk['chunk_size']} "
            f"elapsed={chunk['elapsed_seconds']}"
        )


@app.command("summarize-run-state")
def summarize_run_state_cmd(
    manifest: Path = typer.Option(..., "--manifest", exists=True, help="Source or alignment manifest"),
    run_dir: Path = typer.Option(..., "--run-dir", exists=True, help="Prepare or detector run directory"),
    mode: str = typer.Option(..., "--mode", help="prepare or detector"),
    json_out: bool = typer.Option(False, "--json", help="Emit JSON"),
) -> None:
    """
    Summarize a run against its manifest, including missing and stale cells.
    """
    payload = summarize_run_state(manifest=manifest, run_dir=run_dir, mode=mode)
    if json_out:
        typer.echo(json.dumps(payload, indent=2, sort_keys=True))
        return
    typer.echo(
        f"planned={payload['planned_cells']} completed={payload['completed_cells']} "
        f"failed={payload['failed_cells']} missing={payload['missing_cells']} stale={len(payload['stale_cells'])}"
    )
    if payload["failed_cell_ids"]:
        typer.echo(f"failed_cells={payload['failed_cell_ids']}")
    if payload["missing_cell_ids"]:
        typer.echo(f"missing_cells={payload['missing_cell_ids']}")
    if payload["stale_cells"]:
        typer.echo(f"stale_cells={payload['stale_cells']}")


@app.command("export-run-subset")
def export_run_subset_cmd(
    manifest: Path = typer.Option(..., "--manifest", exists=True, help="Source or alignment manifest to filter"),
    run_dir: Path = typer.Option(..., "--run-dir", exists=True, help="Prepare or detector run directory"),
    out: Path = typer.Option(..., "--out", help="Output subset manifest"),
    subset: str = typer.Option("failed", "--subset", help="failed, missing, stale, incomplete, all-failed-chunks"),
    chunk_index: Optional[int] = typer.Option(None, "--chunk-index", help="Export one specific chunk as a manifest"),
) -> None:
    """
    Export failed, missing, stale, incomplete, or chunk-specific rows as a new manifest for reruns.
    """
    out_path = export_manifest_subset(
        manifest=manifest,
        run_dir=run_dir,
        out_path=out,
        subset=subset,
        chunk_index=chunk_index,
    )
    typer.echo(out_path)


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
        help="Detector name (ciri-full, ciri2, ciri3, find-circ3, circexplorer2)",
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

    if detector in {"ciri-full", "ciri2", "ciri3"}:
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
            "Supported: ciri-full, ciri2, ciri3, find-circ3, circexplorer2."
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

        if det_name in {"ciri-full", "ciri2", "ciri3"}:
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
