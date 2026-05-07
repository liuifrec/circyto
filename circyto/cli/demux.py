from __future__ import annotations

from pathlib import Path
import typer

from circyto.demux.smartseq2 import SmartSeq2DemuxParams, demux_smartseq2_pooled
from circyto.demux.smartseq3 import SmartSeq3DemuxParams, demux_smartseq3

demux_app = typer.Typer(help="Resolve pooled reads into cell-resolved inputs and emit a manifest.")


@demux_app.command("smartseq2")
def smartseq2(
    r1: Path = typer.Option(..., "--r1", help="Pooled R1 FASTQ(.gz)"),
    r2: Path = typer.Option(..., "--r2", help="Pooled R2 FASTQ(.gz)"),
    barcodes: Path = typer.Option(..., "--barcodes", help="Barcode whitelist TSV (barcode per line OR cell_id<TAB>barcode)"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="Output directory"),
    manifest: Path = typer.Option(None, "--manifest", help="Manifest output path (default: OUTDIR/manifest.tsv)"),
    library_id: str = typer.Option("libA", "--library-id", help="Library/pool/run identifier"),
    barcode_read: str = typer.Option("R1", "--barcode-read", help="Which read contains barcode: R1 or R2"),
    barcode_start: int = typer.Option(..., "--barcode-start", help="0-based start position of barcode"),
    barcode_length: int = typer.Option(..., "--barcode-length", help="Length of barcode"),
    trim_barcode: bool = typer.Option(True, "--trim-barcode/--no-trim-barcode", help="Trim barcode bases from sequences"),
    max_mismatch: int = typer.Option(0, "--max-mismatch", help="Barcode mismatches allowed (v1: only 0 supported)"),
    limit_reads: int = typer.Option(None, "--limit-reads", help="Stop after N read-pairs (smoke test)"),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing outputs"),
):
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest if manifest is not None else (outdir / "manifest.tsv")

    params = SmartSeq2DemuxParams(
        r1=r1,
        r2=r2,
        barcodes_tsv=barcodes,
        outdir=outdir,
        manifest_path=manifest_path,
        library_id=library_id,
        barcode_read=barcode_read,
        barcode_start=barcode_start,
        barcode_length=barcode_length,
        trim_barcode=trim_barcode,
        max_mismatch=max_mismatch,
        limit_reads=limit_reads,
        overwrite=overwrite,
    )

    stats = demux_smartseq2_pooled(params)
    typer.echo(f"Demux complete: cells={len(stats)}")
    typer.echo(f"Manifest: {manifest_path}")
    typer.echo(f"Report:   {outdir / 'demux_report.json'}")


@demux_app.command("smartseq3")
def smartseq3(
    read1: Path = typer.Option(..., "--read1", help="[EXPERIMENTAL] Transcript R1 FASTQ(.gz)"),
    read2: Path = typer.Option(..., "--read2", help="[EXPERIMENTAL] Transcript R2 FASTQ(.gz)"),
    index1: Path = typer.Option(..., "--index1", help="[EXPERIMENTAL] Index I1 FASTQ(.gz)"),
    index2: Path = typer.Option(..., "--index2", help="[EXPERIMENTAL] Index I2 FASTQ(.gz)"),
    annotation: Path = typer.Option(
        ...,
        "--annotation",
        help="[EXPERIMENTAL] TSV/CSV annotation with explicit cell and index columns",
    ),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="[EXPERIMENTAL] Output directory"),
    cell_id_column: str = typer.Option(..., "--cell-id-column", help="[EXPERIMENTAL] Annotation cell ID column"),
    index1_column: str = typer.Option(..., "--index1-column", help="[EXPERIMENTAL] Annotation I1 column"),
    index2_column: str = typer.Option(..., "--index2-column", help="[EXPERIMENTAL] Annotation I2 column"),
    max_mismatch: int = typer.Option(
        0,
        "--max-mismatch",
        help="[EXPERIMENTAL] Maximum combined I1+I2 mismatches for assignment",
    ),
    max_records: int = typer.Option(
        None,
        "--max-records",
        help="[EXPERIMENTAL] Stop after N records for smoke/preflight runs",
    ),
    write_sink: bool = typer.Option(
        True,
        "--write-sink/--no-write-sink",
        help="[EXPERIMENTAL] Write unmatched transcript reads to OUTDIR/sink/",
    ),
    compress_output: bool = typer.Option(
        True,
        "--compress-output/--no-compress-output",
        help="[EXPERIMENTAL] Write per-cell FASTQs as .fastq.gz or plain .fastq",
    ),
    emit_manifest: bool = typer.Option(
        True,
        "--emit-manifest/--no-emit-manifest",
        help="[EXPERIMENTAL] Emit OUTDIR/manifest.tsv for downstream circyto alignment workflows",
    ),
):
    """
    Experimental SMART-Seq3 pooled demultiplexing using transcript FASTQs plus dual index FASTQs.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    params = SmartSeq3DemuxParams(
        read1=read1,
        read2=read2,
        index1=index1,
        index2=index2,
        annotation=annotation,
        outdir=outdir,
        cell_id_column=cell_id_column,
        index1_column=index1_column,
        index2_column=index2_column,
        max_mismatch=max_mismatch,
        max_records=max_records,
        write_sink=write_sink,
        compress_output=compress_output,
        emit_manifest=emit_manifest,
    )
    summary = demux_smartseq3(params)
    typer.echo(
        "Experimental SMART-Seq3 demux complete: "
        f"assigned={summary['assigned_records']} total={summary['total_records']} cells={summary['number_of_cells_detected']}"
    )
    if emit_manifest:
        typer.echo(f"Manifest: {outdir / 'manifest.tsv'}")
    typer.echo(f"Summary:  {outdir / 'demux_summary.json'}")
