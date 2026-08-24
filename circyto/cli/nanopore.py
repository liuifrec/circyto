from __future__ import annotations

import json
from pathlib import Path

import typer

from circyto.manifest.long_read import validate_long_read_manifest_tsv
from circyto.paths import get_package_root
from circyto.pipeline.nanopore_archive import (
    load_expected_run,
    query_ena_run,
    validate_ena_run_identity,
)
from circyto.pipeline.nanopore_interop import (
    inspect_alignment_for_bsj,
    prepare_nanopore_alignments,
)
from circyto.pipeline.workflow_reporting import utc_now_iso, write_json


nanopore_app = typer.Typer(
    add_completion=False,
    help="Experimental single-cell Oxford Nanopore ingestion and alignment interoperability.",
)


def expected_run_path(name: str) -> Path:
    normalized = name.strip().lower().replace("_", "")
    if normalized not in {"srr4048177", "mandalorion-srr4048177"}:
        raise ValueError(
            f"Unknown expected run profile {name!r}; available profile: srr4048177"
        )
    return get_package_root() / "resources" / "nanopore" / "srr4048177_expected.json"


@nanopore_app.command("query-run")
def query_run_command(
    accession: str = typer.Option(..., "--accession", help="ENA/SRA run accession."),
    expected: str = typer.Option(
        "srr4048177",
        "--expected",
        help="Checked-in hard identity profile.",
    ),
    output: Path = typer.Option(..., "--output", "-o", dir_okay=False),
    timeout_seconds: float = typer.Option(30.0, "--timeout-seconds", min=0.1),
) -> None:
    """Query ENA metadata and enforce the documented run identity without downloading reads."""
    try:
        expectation = load_expected_run(expected_run_path(expected))
        metadata = query_ena_run(accession, timeout_seconds=timeout_seconds)
        warnings = validate_ena_run_identity(metadata, expectation)
    except (ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    payload = {
        "schema_version": "circyto.ena_metadata_snapshot.v1",
        "metadata_retrieved_at": utc_now_iso(),
        "hard_identity_check": "passed",
        "warning_level_comparisons": warnings,
        "run": metadata.to_dict(),
    }
    write_json(output, payload)
    typer.echo(json.dumps(payload, indent=2, sort_keys=True))


@nanopore_app.command("validate-manifest")
def validate_manifest_command(
    manifest: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    strict: bool = typer.Option(False, "--strict", help="Require every FASTQ path to exist."),
) -> None:
    """Validate the versioned single-cell long-read manifest."""
    ok, errors, summary = validate_long_read_manifest_tsv(manifest, strict=strict)
    payload = {"ok": ok, "errors": errors, **summary}
    typer.echo(json.dumps(payload, indent=2, sort_keys=True))
    if not ok:
        raise typer.Exit(code=2)


@nanopore_app.command("align")
def align_command(
    manifest: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    reference: Path = typer.Option(..., "--reference", exists=True, dir_okay=False),
    reference_id: str = typer.Option(..., "--reference-id"),
    reference_build: str = typer.Option(..., "--reference-build"),
    reference_sha256: str = typer.Option(..., "--reference-sha256"),
    outdir: Path = typer.Option(..., "--outdir", "-o", file_okay=False),
    threads: int = typer.Option(8, "--threads", min=1),
    minimap2: str = typer.Option("minimap2", "--minimap2"),
    samtools: str = typer.Option("samtools", "--samtools"),
    keep_sam: bool = typer.Option(
        False,
        "--keep-sam",
        help="Retain a full SAM only for debugging; disabled by default.",
    ),
    archive_metadata: Path | None = typer.Option(
        None,
        "--archive-metadata",
        exists=True,
        dir_okay=False,
    ),
) -> None:
    """Stream minimap2 into coordinate-sorting samtools and emit cell-scoped outputs."""
    try:
        alignment_manifest = prepare_nanopore_alignments(
            manifest_path=manifest,
            reference_fasta=reference,
            reference_id=reference_id,
            reference_build=reference_build,
            reference_sha256=reference_sha256,
            outdir=outdir,
            threads=threads,
            minimap2=minimap2,
            samtools=samtools,
            keep_sam=keep_sam,
            archive_metadata_path=archive_metadata,
        )
    except (FileNotFoundError, FileExistsError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(alignment_manifest)


@nanopore_app.command("inspect-bsj")
def inspect_bsj_command(
    alignment: Path = typer.Option(..., "--alignment", exists=True, dir_okay=False),
    cell_id: str = typer.Option(..., "--cell-id"),
    input_query_count: int = typer.Option(..., "--input-query-count", min=0),
    output: Path = typer.Option(..., "--output", "-o", dir_okay=False),
    qc_output: Path | None = typer.Option(None, "--qc-output", dir_okay=False),
    samtools: str = typer.Option("samtools", "--samtools"),
) -> None:
    """Report conservative alignment-order patterns; never emit circRNA calls."""
    try:
        qc = inspect_alignment_for_bsj(
            alignment_path=alignment,
            cell_id=cell_id,
            input_query_count=input_query_count,
            output_path=output,
            qc_output_path=qc_output,
            samtools=samtools,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(qc, indent=2, sort_keys=True))
