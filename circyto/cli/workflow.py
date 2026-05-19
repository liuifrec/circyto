from __future__ import annotations

from pathlib import Path

import typer

from circyto.pipeline.workflow_full_length_circrna import (
    FullLengthCircRNAWorkflowParams,
    run_full_length_circrna_workflow,
)
from circyto.pipeline.workflow_smartseq3_ciri3 import (
    SmartSeq3Ciri3WorkflowParams,
    run_smartseq3_ciri3_workflow,
)


workflow_app = typer.Typer(help="Experimental high-level workflows.")


@workflow_app.command("full-length-circrna")
def full_length_circrna(
    manifest: Path = typer.Option(..., "--manifest", exists=True, help="[EXPERIMENTAL] Per-cell FASTQ manifest TSV"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="[EXPERIMENTAL] Workflow output directory"),
    protocol: str = typer.Option(..., "--protocol", help="[EXPERIMENTAL] One of: smartseq3, ramda, shin-ramda"),
    genome_fasta: Path = typer.Option(..., "--genome-fasta", exists=True, help="[EXPERIMENTAL] Reference FASTA"),
    gtf: Path = typer.Option(..., "--gtf", exists=True, help="[EXPERIMENTAL] Annotation GTF"),
    detector: str = typer.Option("ciri3", "--detector", help="[EXPERIMENTAL] Detector backend; currently only ciri3"),
    star_index: Path | None = typer.Option(None, "--star-index", exists=True, help="[EXPERIMENTAL] STAR genomeDir for paired-end rows"),
    threads: int = typer.Option(8, "--threads", help="[EXPERIMENTAL] Threads per alignment/detector task"),
    parallel: int = typer.Option(1, "--parallel", help="[EXPERIMENTAL] Concurrent alignment/detector tasks"),
    chunk_size: int = typer.Option(1, "--chunk-size", help="[EXPERIMENTAL] Resumable chunk size"),
    skip_demux: bool = typer.Option(False, "--skip-demux", help="[EXPERIMENTAL] Skip demux stage. Implied for ramda/shin-ramda."),
    export_h5ad: bool = typer.Option(True, "--export-h5ad/--no-export-h5ad", help="[EXPERIMENTAL] Write OUTDIR/anndata/circ_counts.h5ad"),
    gene_counts: Path | None = typer.Option(None, "--gene-counts", exists=True, help="[EXPERIMENTAL] Optional normalized gene_counts.tsv import for RNA-side snapshotting under OUTDIR/rna/"),
    export_mudata: bool = typer.Option(False, "--export-mudata/--no-export-mudata", help="[PLANNED] Future multimodal output; not implemented yet for full-length workflows."),
    gene_expression_method: str = typer.Option("none", "--gene-expression-method", help="[PLANNED] Future gene-expression method: featurecounts|velocyto|none"),
    velocity_layers: str = typer.Option("none", "--velocity-layers", help="[PLANNED] Future velocity-layer mode: none|velocyto"),
    cleanup_intermediates: bool = typer.Option(False, "--cleanup-intermediates/--no-cleanup-intermediates", help="[PLANNED] After successful completion, optionally clean regenerable alignment intermediates. Dry-run planning only for now."),
    dry_run: bool = typer.Option(False, "--dry-run", help="[EXPERIMENTAL] Plan stages and underlying commands without executing"),
    fail_fast: bool = typer.Option(False, "--fail-fast", help="[EXPERIMENTAL] Stop after the first failed chunk"),
    command_template: str | None = typer.Option(None, "--command-template", help="[EXPERIMENTAL] Optional CIRI3 command template override"),
    allow_paired_ramda: bool = typer.Option(
        False,
        "--allow-paired-ramda",
        "--experimental-paired-ramda",
        help="Allow real paired-end ramda/shin-ramda execution through the validated STAR+CIRI3 path. --experimental-paired-ramda is accepted as a deprecated alias.",
    ),
) -> None:
    """
    Experimental manifest-driven full-length circRNA workflow.
    """
    run_full_length_circrna_workflow(
        FullLengthCircRNAWorkflowParams(
            manifest=manifest,
            outdir=outdir,
            protocol=protocol,
            genome_fasta=genome_fasta,
            gtf=gtf,
            detector=detector,
            star_index=star_index,
            threads=threads,
            parallel=parallel,
            chunk_size=chunk_size,
            skip_demux=skip_demux,
            export_h5ad=export_h5ad,
            gene_counts=gene_counts,
            export_mudata=export_mudata,
            gene_expression_method=gene_expression_method,
            velocity_layers=velocity_layers,
            cleanup_intermediates=cleanup_intermediates,
            dry_run=dry_run,
            fail_fast=fail_fast,
            command_template=command_template,
            experimental_paired_ramda=allow_paired_ramda,
        ),
        progress=typer.echo,
    )


@workflow_app.command("smartseq3-ciri3")
def smartseq3_ciri3(
    read1: Path = typer.Option(..., "--read1", exists=True, help="[EXPERIMENTAL] SMART-Seq3 transcript read1 FASTQ(.gz)"),
    read2: Path = typer.Option(..., "--read2", exists=True, help="[EXPERIMENTAL] SMART-Seq3 transcript read2 FASTQ(.gz)"),
    index1: Path = typer.Option(..., "--index1", exists=True, help="[EXPERIMENTAL] SMART-Seq3 index1 FASTQ(.gz)"),
    index2: Path = typer.Option(..., "--index2", exists=True, help="[EXPERIMENTAL] SMART-Seq3 index2 FASTQ(.gz)"),
    annotation: Path = typer.Option(..., "--annotation", exists=True, help="[EXPERIMENTAL] Plate annotation TSV/CSV"),
    cell_id_column: str = typer.Option(..., "--cell-id-column", help="[EXPERIMENTAL] Annotation cell ID column"),
    index1_column: str = typer.Option(..., "--index1-column", help="[EXPERIMENTAL] Annotation index1 column"),
    index2_column: str = typer.Option(..., "--index2-column", help="[EXPERIMENTAL] Annotation index2 column"),
    ref_fa: Path = typer.Option(..., "--ref-fa", exists=True, help="[EXPERIMENTAL] Reference FASTA"),
    star_index: Path = typer.Option(..., "--star-index", exists=True, help="[EXPERIMENTAL] STAR genomeDir"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="[EXPERIMENTAL] Workflow output directory"),
    top_n: int | None = typer.Option(None, "--top-n", help="[EXPERIMENTAL] Keep only the top N cells by reads_per_cell"),
    max_records: int | None = typer.Option(None, "--max-records", help="[EXPERIMENTAL] Stop demux after N records"),
    threads: int = typer.Option(8, "--threads", help="[EXPERIMENTAL] Threads per alignment/detector task"),
    parallel: int = typer.Option(1, "--parallel", help="[EXPERIMENTAL] Concurrent alignment/detector tasks"),
    chunk_size: int = typer.Option(1, "--chunk-size", help="[EXPERIMENTAL] Resumable chunk size"),
    max_mismatch: int = typer.Option(0, "--max-mismatch", help="[EXPERIMENTAL] Combined I1+I2 mismatches allowed"),
    write_sink: bool = typer.Option(True, "--write-sink/--no-write-sink", help="[EXPERIMENTAL] Write unmatched transcript reads to OUTDIR/demux/sink/"),
    resume: bool = typer.Option(False, "--resume/--no-resume", help="[EXPERIMENTAL] Skip completed workflow stages after validating expected outputs"),
    export_h5ad: bool = typer.Option(True, "--export-h5ad/--no-export-h5ad", help="[EXPERIMENTAL] Write OUTDIR/anndata/circ_counts.h5ad"),
    gene_counts: Path | None = typer.Option(None, "--gene-counts", exists=True, help="[EXPERIMENTAL] Optional gene-count input for multimodal export"),
    gene_counts_format: str = typer.Option("tsv", "--gene-counts-format", help="[EXPERIMENTAL] Gene-count format: tsv or mtx-dir"),
    export_mudata: bool = typer.Option(False, "--export-mudata/--no-export-mudata", help="[EXPERIMENTAL] Write OUTDIR/mudata/circyto_multimodal.h5mu"),
    cell_join: str = typer.Option("inner", "--cell-join", help="[EXPERIMENTAL] How to align RNA and circ cells: inner or outer"),
) -> None:
    """
    Experimental end-to-end SMART-Seq3 to CIRI3 workflow.
    """
    run_smartseq3_ciri3_workflow(
        SmartSeq3Ciri3WorkflowParams(
            read1=read1,
            read2=read2,
            index1=index1,
            index2=index2,
            annotation=annotation,
            cell_id_column=cell_id_column,
            index1_column=index1_column,
            index2_column=index2_column,
            ref_fa=ref_fa,
            star_index=star_index,
            outdir=outdir,
            top_n=top_n,
            max_records=max_records,
            threads=threads,
            parallel=parallel,
            chunk_size=chunk_size,
            max_mismatch=max_mismatch,
            write_sink=write_sink,
            resume=resume,
            export_h5ad=export_h5ad,
            gene_counts=gene_counts,
            gene_counts_format=gene_counts_format,
            export_mudata=export_mudata,
            cell_join=cell_join,
        ),
        progress=typer.echo,
    )
