from __future__ import annotations

from pathlib import Path
import typer

from circyto.demux.smartseq2 import SmartSeq2DemuxParams, demux_smartseq2_pooled

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
