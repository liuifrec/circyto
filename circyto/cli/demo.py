from __future__ import annotations

import json
import shutil
from pathlib import Path

import typer

from circyto.pipeline.collect import collect_matrix

demo_app = typer.Typer(help="Tiny self-contained demos for reviewer smoke tests.")


def _matrix_dims(matrix_path: Path) -> tuple[int, int, int]:
    with matrix_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("%"):
                continue
            rows, cols, nnz = text.split()[:3]
            return int(rows), int(cols), int(nnz)
    return 0, 0, 0


@demo_app.command("mini")
def mini_demo(
    outdir: Path = typer.Option(Path("demo_out"), "--out", "--outdir", "-o", help="Output directory."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Replace an existing demo output directory."),
) -> None:
    """
    Generate a tiny normalized circRNA fixture and collect it into a matrix.

    This command intentionally avoids external detector binaries and references.
    It validates package import, CLI wiring, normalized parser input, matrix
    collection, and host-gene provenance plumbing.
    """
    outdir = outdir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        if not overwrite:
            raise typer.BadParameter(f"Output directory already exists and is not empty: {outdir}")
        shutil.rmtree(outdir)

    detector_dir = outdir / "ciri3"
    matrix_dir = outdir / "matrix"
    detector_dir.mkdir(parents=True, exist_ok=True)
    matrix_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = outdir / "manifest.tsv"
    manifest_path.write_text(
        "cell_id\tdetector_output\n"
        f"demo_cell_1\t{detector_dir / 'demo_cell_1.tsv'}\n"
        f"demo_cell_2\t{detector_dir / 'demo_cell_2.tsv'}\n",
        encoding="utf-8",
    )
    (detector_dir / "demo_cell_1.tsv").write_text(
        "circ_id\tchr\tstart\tend\tstrand\tsupport\thost_gene\n"
        "chr21:100|200|+\tchr21\t100\t200\t+\t3\tSMOKE1\n"
        "chr21:300|420|-\tchr21\t300\t420\t-\t1\tSMOKE2\n",
        encoding="utf-8",
    )
    (detector_dir / "demo_cell_2.tsv").write_text(
        "circ_id\tchr\tstart\tend\tstrand\tsupport\thost_gene\n"
        "chr21:100|200|+\tchr21\t100\t200\t+\t2\tSMOKE1\n",
        encoding="utf-8",
    )

    matrix_path = matrix_dir / "circ_counts.mtx"
    circ_index_path = matrix_dir / "circ_index.txt"
    cell_index_path = matrix_dir / "cell_index.txt"
    collect_matrix(
        cirifull_dir=str(detector_dir),
        matrix_path=str(matrix_path),
        circ_index_path=str(circ_index_path),
        cell_index_path=str(cell_index_path),
        min_count_per_cell=1,
    )
    n_circs, n_cells, nnz = _matrix_dims(matrix_path)
    summary = {
        "demo": "mini",
        "purpose": "self-contained parser and matrix smoke test",
        "detector_fixture": "normalized CIRI3-style TSV",
        "external_tools_required": False,
        "manifest": str(manifest_path),
        "detector_dir": str(detector_dir),
        "matrix_path": str(matrix_path),
        "circ_index_path": str(circ_index_path),
        "cell_index_path": str(cell_index_path),
        "circ_feature_table": str(matrix_dir / "circ_feature_table.tsv"),
        "n_circs": n_circs,
        "n_cells": n_cells,
        "nnz": nnz,
    }
    summary_path = outdir / "demo_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))
