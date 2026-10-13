"""Public detector-free interoperability entry point."""

from pathlib import Path

import typer


def import_circ(
    cells: Path = typer.Option(
        ...,
        "--cells",
        exists=True,
        dir_okay=False,
        help="TSV: cell_id, optional detector_status and metadata; row order is retained.",
    ),
    outdir: Path = typer.Option(
        ...,
        "--outdir",
        help="New output directory; never overwrites an existing directory.",
    ),
    genome_build: str = typer.Option(
        ...,
        "--genome-build",
        help="Exact assembly identifier; no liftover or inferred assembly.",
    ),
    coordinate_convention: str = typer.Option(
        ...,
        "--coordinate-convention",
        help="0-based-half-open or 1-based-closed; coordinates are retained unchanged.",
    ),
    support_definition: str = typer.Option(
        ...,
        "--support-definition",
        help="One explicit count unit, e.g. BSJ_reads or unique_bridges; never pooled across units.",
    ),
    detector: str = typer.Option(
        ...,
        "--detector",
        help="Source detector/assay identity, recorded without running it.",
    ),
    calls: Path | None = typer.Option(
        None,
        "--calls",
        exists=True,
        dir_okay=False,
        help="Long TSV: cell_id,circ_id,chrom,start,end,strand,support. Mutually exclusive with --matrix-dir.",
    ),
    matrix_dir: Path | None = typer.Option(
        None,
        "--matrix-dir",
        exists=True,
        file_okay=False,
        help="Existing circ_counts.mtx (circ x cell), indexes and circ_feature_table.tsv.",
    ),
    ciri3_manifest: Path | None = typer.Option(
        None,
        "--ciri3-manifest",
        exists=True,
        dir_okay=False,
        help="Native 13-column CIRI3 single-sample results: TSV mapping cell_id,sample_id,path,detector_version,result_scope. Exclusive with --calls/--matrix-dir.",
    ),
    rna: Path | None = typer.Option(
        None,
        "--rna",
        exists=True,
        dir_okay=False,
        help="Optional RNA AnnData, matched by exact cell IDs; requires circyto[mudata].",
    ),
) -> None:
    """Import existing cell-resolved support into sparse AnnData, optional MuData and QC HTML.

    Missing detector_status is unknown. Accepted states: completed, failed,
    skipped, unprocessed, unknown. A completed zero-support cell is not evidence
    of biological absence. No detector, alignment or coordinate conversion runs.
    """
    from circyto.pipeline.import_circ import import_circ as run

    try:
        result = run(
            cells=cells,
            outdir=outdir,
            genome_build=genome_build,
            coordinate_convention=coordinate_convention,
            support_definition=support_definition,
            detector=detector,
            calls=calls,
            matrix_dir=matrix_dir,
            rna=rna,
            ciri3_manifest=ciri3_manifest,
        )
    except (ValueError, OSError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(
        f"Imported {result['n_cells']} cells x {result['n_candidates']} candidates. Report: {outdir / 'qc/report.html'}"
    )
