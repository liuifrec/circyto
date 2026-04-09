from __future__ import annotations

from pathlib import Path
import typer

from circyto.manifest.alignment import validate_alignment_manifest_tsv
from circyto.manifest.v1 import validate_manifest_tsv

manifest_app = typer.Typer(help="Manifest utilities (validate, summarize).")


@manifest_app.command("validate")
def validate(
    manifest: Path = typer.Argument(..., help="Manifest v1 TSV"),
    strict: bool = typer.Option(False, "--strict", help="Require referenced files to exist"),
):
    ok, errs, summary = validate_manifest_tsv(manifest, strict=strict)
    typer.echo(f"cells={summary['cells']} fastq_rows={summary['fastq_rows']} bam_rows={summary['bam_rows']}")
    if strict:
        typer.echo(f"missing_files={summary['missing_files']}")
    if not ok:
        for e in errs:
            typer.echo(f"ERROR: {e}", err=True)
        raise typer.Exit(code=1)
    typer.echo("OK")


@manifest_app.command("validate-alignment")
def validate_alignment(
    manifest: Path = typer.Argument(..., help="Alignment manifest TSV"),
    strict: bool = typer.Option(False, "--strict", help="Require referenced files to exist"),
):
    ok, errs, summary = validate_alignment_manifest_tsv(manifest, strict=strict)
    typer.echo(f"cells={summary['cells']} bam_rows={summary['bam_rows']} sam_rows={summary['sam_rows']}")
    if strict:
        typer.echo(f"missing_files={summary['missing_files']}")
    if not ok:
        for e in errs:
            typer.echo(f"ERROR: {e}", err=True)
        raise typer.Exit(code=1)
    typer.echo("OK")
