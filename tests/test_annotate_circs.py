from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from scipy import sparse as sp
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.annotate_circs import (
    annotate_circ_table,
    parse_annotation_db_spec,
)
from circyto.pipeline.workflow_reporting import export_circ_h5ad


runner = CliRunner()


def _write_circ_table(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "circ_id": "circA",
                "chrom": "chr1",
                "start": 100,
                "end": 200,
                "strand": "+",
                "host_gene": "GENE1",
                "n_cells_detected": 4,
                "total_support": 9,
            },
            {
                "circ_id": "circB",
                "chrom": "chr1",
                "start": 300,
                "end": 400,
                "strand": "+",
                "host_gene": "GENE2",
                "n_cells_detected": 3,
                "total_support": 5,
            },
            {
                "circ_id": "circC",
                "chrom": "chr2",
                "start": 500,
                "end": 600,
                "strand": "-",
                "host_gene": "GENE3",
                "n_cells_detected": 2,
                "total_support": 7,
            },
            {
                "circ_id": "circD",
                "chrom": "chr3",
                "start": 700,
                "end": 800,
                "strand": "+",
                "host_gene": "",
                "n_cells_detected": 1,
                "total_support": 1,
            },
        ]
    ).to_csv(path, sep="\t", index=False)


def _write_circatlas_style(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "seqnames": "chr1",
                "donor": 100,
                "acceptor": 200,
                "bsj_strand": "+",
                "atlas_id": "atlas_exact",
                "gene_symbol": "GENE1",
                "tissue": "brain",
            },
            {
                "seqnames": "chr1",
                "donor": 302,
                "acceptor": 398,
                "bsj_strand": "+",
                "atlas_id": "atlas_near",
                "gene_symbol": "GENE2",
                "tissue": "heart",
            },
            {
                "seqnames": "chr2",
                "donor": 500,
                "acceptor": 600,
                "bsj_strand": "+",
                "atlas_id": "atlas_wrong_strand",
                "gene_symbol": "GENE3",
                "tissue": "liver",
            },
        ]
    ).to_csv(path, sep="\t", index=False)


def _write_circsc_style(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "chromosome": "chr2",
                "start_1based": 500,
                "end_1based": 600,
                "circsc_id": "circsc_exact",
                "gene_label": "GENE3",
                "celltype": "neuron",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def test_annotate_circs_exact_near_and_novel(tmp_path: Path) -> None:
    circ_table = tmp_path / "circ_qc.tsv"
    db1 = tmp_path / "circatlas.tsv"
    db2 = tmp_path / "circsc.tsv"
    out = tmp_path / "circ_qc.annotated.tsv"

    _write_circ_table(circ_table)
    _write_circatlas_style(db1)
    _write_circsc_style(db2)

    summary = annotate_circ_table(
        circ_table=circ_table,
        db_specs=[
            parse_annotation_db_spec(
                f"name=circatlas;path={db1};chrom=seqnames;start=donor;end=acceptor;strand=bsj_strand;id=atlas_id;host_gene=gene_symbol;extra=tissue"
            ),
            parse_annotation_db_spec(
                f"name=circsc;path={db2};chrom=chromosome;start=start_1based;end=end_1based;id=circsc_id;host_gene=gene_label;extra=celltype"
            ),
        ],
        out_path=out,
        max_bsj_distance=3,
        enable_host_gene_match=False,
    )

    annotated = pd.read_csv(out, sep="\t", keep_default_na=False)
    by_id = annotated.set_index("circ_id")

    assert by_id.loc["circA", "circatlas_known"] in {True, "True", 1}
    assert by_id.loc["circA", "circatlas_match_type"] == "known_exact_bsj"
    assert by_id.loc["circA", "circatlas_ids"] == "atlas_exact"
    assert by_id.loc["circA", "best_annotation_status"] == "known_exact_bsj"

    assert by_id.loc["circB", "circatlas_match_type"] == "known_near_bsj"
    assert by_id.loc["circB", "circatlas_extra_summary"] == "tissue=heart"
    assert by_id.loc["circB", "best_annotation_status"] == "known_near_bsj"

    assert by_id.loc["circC", "circatlas_match_type"] == "novel"
    assert by_id.loc["circC", "circsc_match_type"] == "known_exact_bsj"
    assert by_id.loc["circC", "known_database_count"] == 1
    assert by_id.loc["circC", "known_databases"] == "circsc"

    assert by_id.loc["circD", "best_annotation_status"] == "novel"
    assert by_id.loc["circD", "known_database_count"] == 0

    assert summary["exact_matches_per_database"]["circatlas"] == 1
    assert summary["near_matches_per_database"]["circatlas"] == 1
    assert summary["exact_matches_per_database"]["circsc"] == 1
    assert summary["novel_count"] == 1
    assert summary["recurrent_known_count"] == 3
    assert summary["recurrent_novel_count"] == 0


def test_annotate_circs_strand_mismatch_and_host_gene_only(tmp_path: Path) -> None:
    circ_table = tmp_path / "circ_qc.tsv"
    db = tmp_path / "hostgene.tsv"
    out = tmp_path / "annotated.tsv"

    pd.DataFrame(
        [
            {
                "circ_id": "circX",
                "chrom": "chr4",
                "start": 1000,
                "end": 1100,
                "strand": "-",
                "host_gene": "GENEX",
            }
        ]
    ).to_csv(circ_table, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "chrom": "chr4",
                "start": 1000,
                "end": 1100,
                "strand": "+",
                "db_id": "wrong_strand",
                "host_gene": "GENEX",
            },
            {
                "chrom": "chr4",
                "start": 1030,
                "end": 1080,
                "strand": "-",
                "db_id": "host_only",
                "host_gene": "GENEX",
            },
        ]
    ).to_csv(db, sep="\t", index=False)

    annotate_circ_table(
        circ_table=circ_table,
        db_specs=[
            parse_annotation_db_spec(
                f"name=testdb;path={db};chrom=chrom;start=start;end=end;strand=strand;id=db_id;host_gene=host_gene"
            )
        ],
        out_path=out,
        enable_host_gene_match=True,
    )
    annotated = pd.read_csv(out, sep="\t", keep_default_na=False).set_index("circ_id")
    assert annotated.loc["circX", "testdb_match_type"] == "known_host_gene_only"
    assert annotated.loc["circX", "best_annotation_status"] == "known_host_gene_only"
    assert annotated.loc["circX", "testdb_ids"] == "host_only"

    annotate_circ_table(
        circ_table=circ_table,
        db_specs=[
            parse_annotation_db_spec(
                f"name=testdb;path={db};chrom=chrom;start=start;end=end;strand=strand;id=db_id;host_gene=host_gene"
            )
        ],
        out_path=out,
        enable_host_gene_match=False,
    )
    annotated_disabled = pd.read_csv(out, sep="\t", keep_default_na=False).set_index("circ_id")
    assert annotated_disabled.loc["circX", "testdb_match_type"] == "novel"
    assert annotated_disabled.loc["circX", "best_annotation_status"] == "novel"


def test_annotate_circs_updates_h5ad_var_preserving_shape_and_indices(tmp_path: Path) -> None:
    ad = pytest.importorskip("anndata")
    circ_table = tmp_path / "circ_qc.tsv"
    db1 = tmp_path / "circatlas.tsv"
    out = tmp_path / "annotated.tsv"
    h5ad_path = tmp_path / "circ_counts.h5ad"

    _write_circ_table(circ_table)
    _write_circatlas_style(db1)

    cell_qc = pd.DataFrame({"assigned_reads": [10, 8]}, index=["cellA", "cellB"])
    circ_qc = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr3"],
            "start": [100, 300, 500, 700],
            "end": [200, 400, 600, 800],
            "strand": ["+", "+", "-", "+"],
            "host_gene": ["GENE1", "GENE2", "GENE3", ""],
            "n_cells_detected": [2, 1, 1, 0],
            "total_support": [6, 1, 2, 0],
            "max_support": [4, 1, 2, 0],
            "mean_support_detected_cells": [3.0, 1.0, 2.0, 0.0],
        },
        index=["circA", "circB", "circC", "circD"],
    )
    export_circ_h5ad(
        out_path=h5ad_path,
        X_cells_by_circ=sp.csr_matrix([[4, 1, 0, 0], [2, 0, 2, 0]], dtype=int),
        cell_qc=cell_qc,
        circ_qc=circ_qc,
        uns_payload={"workflow_name": "smartseq3-ciri3"},
    )

    annotate_circ_table(
        circ_table=circ_table,
        db_specs=[
            parse_annotation_db_spec(
                f"name=circatlas;path={db1};chrom=seqnames;start=donor;end=acceptor;strand=bsj_strand;id=atlas_id;host_gene=gene_symbol;extra=tissue"
            )
        ],
        out_path=out,
        update_h5ad=h5ad_path,
    )

    adata = ad.read_h5ad(h5ad_path)
    assert adata.shape == (2, 4)
    assert list(adata.var_names) == ["circA", "circB", "circC", "circD"]
    assert adata.var.loc["circA", "circatlas_match_type"] == "known_exact_bsj"
    assert adata.var.loc["circD", "best_annotation_status"] == "novel"
    assert adata.uns["circyto"]["circ_annotation"]["databases"] == ["circatlas"]


def test_annotate_circs_cli_writes_outputs(tmp_path: Path) -> None:
    circ_table = tmp_path / "circ_qc.tsv"
    db = tmp_path / "db.tsv"
    out = tmp_path / "annotated.tsv"
    summary_path = tmp_path / "annotation_summary.json"

    pd.DataFrame(
        [
            {
                "circ_id": "circA",
                "chrom": "chr1",
                "start": 10,
                "end": 20,
                "strand": "+",
                "host_gene": "GENE1",
            }
        ]
    ).to_csv(circ_table, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "start": 10,
                "end": 20,
                "strand": "+",
                "id": "db1",
                "host_gene": "GENE1",
            }
        ]
    ).to_csv(db, sep="\t", index=False)

    result = runner.invoke(
        app,
        [
            "annotate-circs",
            "--circ-table",
            str(circ_table),
            "--annotation-db",
            f"testdb:{db}",
            "--out",
            str(out),
            "--summary-out",
            str(summary_path),
        ],
    )
    assert result.exit_code == 0, result.stdout
    annotated = pd.read_csv(out, sep="\t", keep_default_na=False)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert annotated.loc[0, "testdb_match_type"] == "known_exact_bsj"
    assert summary["exact_matches_per_database"]["testdb"] == 1
