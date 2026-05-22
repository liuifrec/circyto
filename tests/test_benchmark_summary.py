from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_benchmark_workdir(root: Path, *, name: str, protocol: str, read_layout: str, circ_cells: list[str], include_h5mu: bool) -> Path:
    workdir = root / name
    (workdir / "rna").mkdir(parents=True)
    (workdir / "matrix").mkdir()
    (workdir / "qc").mkdir()
    if include_h5mu:
        (workdir / "mudata").mkdir()
        (workdir / "mudata" / "full_length.h5mu").write_text("dummy", encoding="utf-8")
    (workdir / "workflow_summary.json").write_text(
        json.dumps(
            {
                "workflow_type": "full-length-circrna",
                "protocol": protocol,
                "read_layout": read_layout,
                "cleanup_summary": {"performed": False},
                "stage_graph": [{"stage": "workflow", "status": "completed"}],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (workdir / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\tcellA\tcellB\nG1\tGENE1\t10\t20\nG2\tGENE2\t1\t2\n",
        encoding="utf-8",
    )
    (workdir / "rna" / "gene_feature_table.tsv").write_text(
        "gene_id\tgene_name\nG1\tGENE1\nG2\tGENE2\n",
        encoding="utf-8",
    )
    (workdir / "qc" / "rna_qc.tsv").write_text(
        "cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\n"
        "cellA\t11\t2\t0.1\t0.0\t2\n"
        "cellB\t22\t2\t0.0\t0.0\t0\n",
        encoding="utf-8",
    )
    (workdir / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 2\n1 1 3\n2 1 1\n",
        encoding="utf-8",
    )
    (workdir / "matrix" / "circ_index.txt").write_text("circ1\ncirc2\n", encoding="utf-8")
    (workdir / "matrix" / "cell_index.txt").write_text("\n".join(circ_cells) + "\n", encoding="utf-8")
    return workdir


def test_summarize_benchmark_writes_outputs(tmp_path: Path) -> None:
    w1 = _write_benchmark_workdir(tmp_path, name="scrr_imr90", protocol="ramda", read_layout="single", circ_cells=["cellA"], include_h5mu=True)
    w2 = _write_benchmark_workdir(tmp_path, name="scrr_hap1", protocol="ramda", read_layout="paired", circ_cells=["cellA"], include_h5mu=False)
    out_tsv = tmp_path / "benchmark_summary.tsv"
    out_json = tmp_path / "benchmark_summary.json"
    result = runner.invoke(
        app,
        [
            "summarize-benchmark",
            "--workdirs",
            str(w1),
            "--workdirs",
            str(w2),
            "--output",
            str(out_tsv),
            "--json",
            str(out_json),
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_workdirs"] == 2
    df = pd.read_csv(out_tsv, sep="\t")
    assert set(df.columns) >= {"dataset_name", "n_cells", "n_rna_features", "n_circ_features", "workflow_succeeded"}
    assert out_json.exists()
    assert bool(df.loc[df["dataset_name"] == "GSE278958 IMR90 scRR", "h5mu_exists"].iloc[0]) is True

