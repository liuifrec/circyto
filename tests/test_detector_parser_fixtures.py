from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
from scipy.io import mmread

from circyto.detectors.ciri3 import _normalize_ciri3_output
from circyto.pipeline.collect import collect_matrix
from circyto.pipeline.collect_circexplorer2_matrix import (
    _parse_circexplorer2_file,
    collect_circexplorer2_matrix,
)
from circyto.pipeline.collect_find_circ3 import _parse_splice_sites, collect_find_circ3_matrix

FIXTURES = Path(__file__).parent / "fixtures" / "detectors"


def test_ciri3_raw_fixture_normalizes_to_stable_schema(tmp_path: Path) -> None:
    out = tmp_path / "ciri3.tsv"
    _normalize_ciri3_output(FIXTURES / "ciri3_raw.tsv", out)

    df = pd.read_csv(out, sep="\t", keep_default_na=False)
    assert df.columns.tolist() == ["circ_id", "chr", "start", "end", "strand", "support"]
    assert df["circ_id"].tolist() == ["chr21:100|200|+", "chr21:300|420|-"]
    assert df.loc[0, "support"] == 3


def test_ciri_style_fixtures_collect_matrix_and_host_gene_provenance(tmp_path: Path) -> None:
    indir = tmp_path / "detector_tsvs"
    indir.mkdir()
    shutil.copy2(FIXTURES / "ciri2_normalized.tsv", indir / "cell_a.tsv")
    shutil.copy2(FIXTURES / "cirifull_normalized.tsv", indir / "cell_b.tsv")

    matrix = tmp_path / "matrix" / "circ_counts.mtx"
    circ_index = tmp_path / "matrix" / "circ_index.txt"
    cell_index = tmp_path / "matrix" / "cell_index.txt"
    collect_matrix(
        cirifull_dir=str(indir),
        matrix_path=str(matrix),
        circ_index_path=str(circ_index),
        cell_index_path=str(cell_index),
        min_count_per_cell=1,
    )

    mat = mmread(matrix).tocsr()
    assert mat.shape == (2, 2)
    assert mat.nnz == 3
    assert circ_index.read_text(encoding="utf-8").splitlines() == [
        "chr21:100|200|+",
        "chr21:300|420|-",
    ]
    features = pd.read_csv(matrix.with_name("circ_feature_table.tsv"), sep="\t", keep_default_na=False)
    by_id = features.set_index("circ_id")
    assert by_id.loc["chr21:100|200|+", "host_gene"] == "SMOKE1"
    assert by_id.loc["chr21:100|200|+", "host_gene_source"] == "gtf"


def test_find_circ3_bed_fixture_parser_and_matrix(tmp_path: Path) -> None:
    parsed = _parse_splice_sites(FIXTURES / "find_circ3_splice_sites.bed", "cell_a")
    assert parsed.to_dict("records") == [
        {"circ_id": "chr21:100|200|+", "cell_id": "cell_a", "support": 2}
    ]

    root = tmp_path / "find_circ3"
    for cell_id in ("cell_a", "cell_b"):
        cell_dir = root / cell_id
        cell_dir.mkdir(parents=True)
        shutil.copy2(FIXTURES / "find_circ3_splice_sites.bed", cell_dir / f"{cell_id}_splice_sites.bed")

    matrix = tmp_path / "find_circ3.mtx"
    circ_index = tmp_path / "find_circ3.circ.txt"
    cell_index = tmp_path / "find_circ3.cell.txt"
    collect_find_circ3_matrix(
        findcirc3_dir=str(root),
        matrix_path=str(matrix),
        circ_index_path=str(circ_index),
        cell_index_path=str(cell_index),
        min_count_per_cell=1,
    )

    mat = mmread(matrix).tocsr()
    assert mat.shape == (1, 2)
    assert mat.nnz == 2
    assert circ_index.read_text(encoding="utf-8").strip() == "chr21:100|200|+"
    assert cell_index.read_text(encoding="utf-8").splitlines() == ["cell_a", "cell_b"]


def test_circexplorer2_fixture_parser_and_matrix(tmp_path: Path) -> None:
    events = _parse_circexplorer2_file(FIXTURES / "circexplorer2_circ.txt")
    assert events == [("chr21:100|200|+", 1), ("chr21:300|420|-", 1)]

    root = tmp_path / "circexplorer2"
    for cell_id in ("cell_a", "cell_b"):
        cell_dir = root / cell_id
        cell_dir.mkdir(parents=True)
        shutil.copy2(FIXTURES / "circexplorer2_circ.txt", cell_dir / f"{cell_id}_CIRCexplorer2_circ.txt")

    matrix, circ_index, cell_index = collect_circexplorer2_matrix(root, tmp_path / "matrix")
    mat = mmread(matrix).tocsr()
    assert mat.shape == (2, 2)
    assert mat.nnz == 4
    assert circ_index.read_text(encoding="utf-8").splitlines() == [
        "chr21:100|200|+",
        "chr21:300|420|-",
    ]
    assert cell_index.read_text(encoding="utf-8").splitlines() == ["cell_a", "cell_b"]
