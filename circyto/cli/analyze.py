from __future__ import annotations

from pathlib import Path

import typer

from circyto.pipeline.workflow_reporting import summarize_h5ad


analyze_app = typer.Typer(help="Lightweight downstream summaries for analysis-ready outputs.")


@analyze_app.command("summarize-h5ad")
def summarize_h5ad_cmd(
    input: Path = typer.Option(..., "--input", exists=True, help="Input circRNA .h5ad"),
    outdir: Path = typer.Option(..., "--outdir", "-o", help="Output summary directory"),
) -> None:
    """
    Summarize a circRNA h5ad and export small downstream-ready tables.
    """
    summarize_h5ad(input, outdir)
