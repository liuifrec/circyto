from pathlib import Path

import typer

from circyto.pipeline.qc_report import generate_report


def report(
    workdir: Path = typer.Option(..., "--workdir", exists=True, file_okay=False,
                                help="Existing workflow results directory; source artifacts are read-only."),
) -> None:
    """Create an offline HTML QC report and metrics.json from existing outputs."""
    try:
        payload = generate_report(workdir)
    except (OSError, ValueError) as exc:
        typer.echo(f"Report generation failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"HTML report: {workdir / 'qc' / 'report.html'}")
    typer.echo(f"Metrics: {workdir / 'qc' / 'metrics.json'}")
    typer.echo(f"Recorded workflow status: {payload['workflow']['status']}")
    if payload["workflow"]["status"] != "completed":
        typer.echo("WARNING: Workflow completion is not established. Inspect the report before interpreting counts.")
    if payload["warnings"]:
        typer.echo(f"Report warnings: {len(payload['warnings'])}; see Warnings and availability in the HTML report.")
