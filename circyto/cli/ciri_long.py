from __future__ import annotations

import json
from pathlib import Path

import typer

from circyto.manifest.ciri_long import validate_ciri_long_manifest_tsv
from circyto.pipeline.ciri_long_adapter import (
    check_ciri_long_readiness,
    normalize_ciri_long_outputs,
    run_ciri_long_call_stage,
    run_ciri_long_collapse_stage,
)


ciri_long_app = typer.Typer(
    add_completion=False,
    help=(
        "Optional CIRI-long adapter for circRNA-enriched rolling-circle "
        "reverse-transcription Nanopore libraries."
    ),
)


def _echo(payload: dict[str, object]) -> None:
    typer.echo(json.dumps(payload, indent=2, sort_keys=True))


@ciri_long_app.command("doctor")
def doctor_command(
    ciri_long: str = typer.Option("CIRI-long", "--ciri-long"),
    bwa: str = typer.Option("bwa", "--bwa"),
    reference: Path | None = typer.Option(
        None, "--reference", dir_okay=False, help="Reference FASTA to check and index."
    ),
    gtf: Path | None = typer.Option(
        None, "--gtf", dir_okay=False, help="Optional GTF annotation."
    ),
    circ_annotation: Path | None = typer.Option(
        None,
        "--circ-annotation",
        dir_okay=False,
        help="Optional known-circRNA BED/GTF annotation.",
    ),
) -> None:
    """Check CIRI-long, BWA, reference indexes, and optional annotations."""
    payload = check_ciri_long_readiness(
        ciri_long=ciri_long,
        bwa=bwa,
        reference_fasta=reference,
        gtf=gtf,
        circ_annotation=circ_annotation,
    )
    _echo(payload)
    if not payload["ok"]:
        raise typer.Exit(code=2)


@ciri_long_app.command("validate-manifest")
def validate_manifest_command(
    manifest: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    strict: bool = typer.Option(
        False, "--strict", help="Require every raw read path to exist."
    ),
) -> None:
    """Validate the versioned RCRT/circRNA-enriched input manifest."""
    ok, errors, summary = validate_ciri_long_manifest_tsv(
        manifest, strict=strict
    )
    payload: dict[str, object] = {"ok": ok, "errors": errors, **summary}
    _echo(payload)
    if not ok:
        raise typer.Exit(code=2)


def _call_options(
    *,
    manifest: Path,
    reference: Path,
    outdir: Path,
    gtf: Path | None,
    circ_annotation: Path | None,
    threads: int,
    ciri_long: str,
    bwa: str,
    execute: bool,
) -> None:
    try:
        payload = run_ciri_long_call_stage(
            manifest_path=manifest,
            reference_fasta=reference,
            outdir=outdir,
            gtf=gtf,
            circ_annotation=circ_annotation,
            threads=threads,
            ciri_long=ciri_long,
            bwa=bwa,
            execute=execute,
        )
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    _echo(payload)


@ciri_long_app.command("plan")
def plan_command(
    manifest: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    reference: Path = typer.Option(..., "--reference", exists=True, dir_okay=False),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False),
    gtf: Path | None = typer.Option(None, "--gtf", exists=True, dir_okay=False),
    circ_annotation: Path | None = typer.Option(
        None, "--circ-annotation", exists=True, dir_okay=False
    ),
    threads: int = typer.Option(1, "--threads", min=1),
    ciri_long: str = typer.Option("CIRI-long", "--ciri-long"),
    bwa: str = typer.Option("bwa", "--bwa"),
) -> None:
    """Validate readiness and write shell-free CIRI-long call plans."""
    _call_options(
        manifest=manifest,
        reference=reference,
        outdir=outdir,
        gtf=gtf,
        circ_annotation=circ_annotation,
        threads=threads,
        ciri_long=ciri_long,
        bwa=bwa,
        execute=False,
    )


@ciri_long_app.command("call")
def call_command(
    manifest: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    reference: Path = typer.Option(..., "--reference", exists=True, dir_okay=False),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False),
    gtf: Path | None = typer.Option(None, "--gtf", exists=True, dir_okay=False),
    circ_annotation: Path | None = typer.Option(
        None, "--circ-annotation", exists=True, dir_okay=False
    ),
    threads: int = typer.Option(1, "--threads", min=1),
    ciri_long: str = typer.Option("CIRI-long", "--ciri-long"),
    bwa: str = typer.Option("bwa", "--bwa"),
    execute: bool = typer.Option(
        False,
        "--execute",
        help="Execute the planned call stage. Without this flag, write a dry-run plan.",
    ),
) -> None:
    """Plan or execute `CIRI-long call` independently for each manifest sample."""
    _call_options(
        manifest=manifest,
        reference=reference,
        outdir=outdir,
        gtf=gtf,
        circ_annotation=circ_annotation,
        threads=threads,
        ciri_long=ciri_long,
        bwa=bwa,
        execute=execute,
    )


@ciri_long_app.command("collapse")
def collapse_command(
    call_manifest: Path = typer.Option(
        ..., "--call-manifest", exists=True, dir_okay=False
    ),
    reference: Path = typer.Option(..., "--reference", exists=True, dir_okay=False),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False),
    prefix: str = typer.Option("cohort", "--prefix"),
    gtf: Path | None = typer.Option(None, "--gtf", exists=True, dir_okay=False),
    circ_annotation: Path | None = typer.Option(
        None, "--circ-annotation", exists=True, dir_okay=False
    ),
    threads: int = typer.Option(1, "--threads", min=1),
    ciri_long: str = typer.Option("CIRI-long", "--ciri-long"),
    bwa: str = typer.Option("bwa", "--bwa"),
    execute: bool = typer.Option(
        False,
        "--execute",
        help="Execute collapse. Without this flag, write a shell-free dry-run plan.",
    ),
) -> None:
    """Plan or execute `CIRI-long collapse` from a call-stage manifest."""
    try:
        payload = run_ciri_long_collapse_stage(
            call_manifest_path=call_manifest,
            reference_fasta=reference,
            outdir=outdir,
            prefix=prefix,
            gtf=gtf,
            circ_annotation=circ_annotation,
            threads=threads,
            ciri_long=ciri_long,
            bwa=bwa,
            execute=execute,
        )
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    _echo(payload)


@ciri_long_app.command("import")
def import_command(
    collapse_dir: Path = typer.Option(
        ..., "--collapse-dir", exists=True, file_okay=False
    ),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False),
    prefix: str = typer.Option("cohort", "--prefix"),
) -> None:
    """Normalize official collapse outputs into separate BSJ/isoform/evidence tables."""
    try:
        payload = normalize_ciri_long_outputs(
            collapse_dir=collapse_dir,
            outdir=outdir,
            prefix=prefix,
        )
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    _echo(payload)
