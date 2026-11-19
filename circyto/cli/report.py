from pathlib import Path
from typing import List, Optional

import typer

from circyto.pipeline.qc_report import generate_report


def report(
    workdir: Optional[Path] = typer.Option(None, "--workdir", exists=True, file_okay=False,
                                help="Existing workflow results directory; source artifacts are read-only."),
    explorer: Optional[Path] = typer.Option(None, "--explorer", exists=True, dir_okay=False,
                                help="Opt in: circ AnnData .h5ad or MuData .h5mu containing circ/rna modalities."),
    rna: Optional[Path] = typer.Option(None, "--rna", exists=True, dir_okay=False,
                                help="Separate RNA .h5ad supplies expression, annotations and existing embedding."),
    output_dir: Optional[Path] = typer.Option(None, "--output-dir", file_okay=False,
                                help="Report destination; required without --workdir. Defaults to WORKDIR/qc."),
    title: str = typer.Option("Single-cell circRNA explorer", "--title", help="Dataset label for the optional explorer."),
    cells: Optional[Path] = typer.Option(None, "--cells", exists=True, dir_okay=False,
                                help="Explicit report cohort: one exact cell ID per line, no header. Source objects are unchanged."),
    umap_key: str = typer.Option("X_umap", "--umap-key", help="Existing RNA obsm key; never recomputed."),
    rna_layer: str = typer.Option("X", "--rna-layer", help="Stored RNA X or layer name; no transformation."),
    expression_unit: str = typer.Option("stored RNA values (units not declared)", "--expression-unit",
                                help="Explicit label for the selected RNA matrix (e.g. log1p CP10K)."),
    support_unit: Optional[str] = typer.Option(None, "--support-unit", help="Native support label; defaults to recorded support_definition."),
    cell_id_key: Optional[str] = typer.Option(None, "--cell-id-key", help="Exact biological cell ID column in every obs; defaults to obs_names."),
    gene: Optional[List[str]] = typer.Option(None, "--gene", help="Exact stable ID or symbol; repeatable. Duplicate symbols retain all matched IDs."),
    circ: Optional[List[str]] = typer.Option(None, "--circ", help="Exact BSJ, catalogue alias or host gene; repeatable."),
    obs_column: Optional[List[str]] = typer.Option(None, "--obs-column", help="Observation field to embed, e.g. rna:condition; repeatable."),
    candidate_metadata: Optional[Path] = typer.Option(None, "--candidate-metadata", exists=True, dir_okay=False,
                                help="Optional TSV keyed by circ_id: aliases, annotation sources, evidence URLs/DOIs. Conflicts rejected."),
    max_genes: int = typer.Option(32, "--max-genes", min=0, help="Maximum embedded gene features."),
    max_circs: int = typer.Option(32, "--max-circs", min=0, help="Maximum embedded circRNA features."),
    max_cells: int = typer.Option(50000, "--max-cells", min=1, help="Cell bound; oversized objects fail without sampling."),
    max_values: int = typer.Option(2000000, "--max-values", min=1, help="Bound on embedded cell × feature measurements."),
    max_bytes: int = typer.Option(32000000, "--max-bytes", min=1, help="Bound on serialized explorer JSON bytes."),
) -> None:
    """Create offline QC reports, optionally with an exact-cell RNA/circRNA explorer."""
    try:
        if explorer is None:
            if workdir is None:
                raise ValueError("Supply --workdir for QC or --explorer and --output-dir for a standalone explorer")
            if (rna is not None or gene or circ or candidate_metadata is not None or obs_column or cells is not None
                    or umap_key != "X_umap" or rna_layer != "X" or cell_id_key is not None or support_unit is not None
                    or expression_unit != "stored RNA values (units not declared)" or title != "Single-cell circRNA explorer"
                    or (max_genes, max_circs, max_cells, max_values, max_bytes) != (32, 32, 50000, 2000000, 32000000)):
                raise ValueError("RNA/feature options require --explorer")
            data = None
        else:
            from circyto.pipeline.single_cell_explorer import ExplorerOptions, load_explorer

            if workdir is None and output_dir is None:
                raise ValueError("Standalone exploration requires --output-dir")
            options = ExplorerOptions(title=title, umap_key=umap_key, rna_layer=rna_layer,
                expression_unit=expression_unit, support_unit=support_unit, cell_id_key=cell_id_key,
                genes=tuple(gene or ()), circs=tuple(circ or ()), obs_columns=tuple(obs_column or ()),
                max_genes=max_genes, max_circs=max_circs, max_cells=max_cells, max_values=max_values, max_bytes=max_bytes)
            data = load_explorer(explorer, rna, options=options, candidate_metadata_path=candidate_metadata, cells_path=cells)
        payload = generate_report(workdir, explorer=data, output_dir=output_dir,
                                  protected_sources=tuple(p for p in (explorer, rna, candidate_metadata, cells) if p is not None))
    except (OSError, ValueError) as exc:
        typer.echo(f"Report generation failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    destination = output_dir if output_dir is not None else workdir / "qc"
    typer.echo(f"HTML report: {destination / 'report.html'}")
    if data is not None:
        typer.echo(f"Explorer data: {destination / 'explorer.json'}")
        typer.echo(f"Embedded: {len(data['cells'])} cells, {len(data['genes'])} genes, {len(data['circs'])} circRNAs; {data['embedding']['plotted_cells']} RNA coordinates")
        for warning in data["warnings"]:
            typer.echo(f"Explorer note: {warning}")
    if workdir is None:
        return
    typer.echo(f"Metrics: {destination / 'metrics.json'}")
    typer.echo(f"Recorded workflow status: {payload['workflow']['status']}")
    if payload["workflow"]["status"] != "completed":
        typer.echo("WARNING: Workflow completion is not established. Inspect the report before interpreting counts.")
    if payload["warnings"]:
        typer.echo(f"Report warnings: {len(payload['warnings'])}; see Warnings and availability in the HTML report.")
