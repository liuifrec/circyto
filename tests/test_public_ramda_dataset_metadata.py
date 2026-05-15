from __future__ import annotations

import csv
from pathlib import Path


REQUIRED_COLUMNS = [
    "dataset_id",
    "protocol",
    "organism",
    "cell_type_or_system",
    "source",
    "suggested_use",
    "notes",
]


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = [{key: "" if value is None else str(value) for key, value in row.items()} for row in reader]
    return header, rows


def test_public_ramda_dataset_metadata_tsv_exists_and_has_required_columns() -> None:
    path = Path("testdata/public_datasets/ramda_datasets.tsv")
    assert path.exists(), f"Missing metadata TSV: {path}"

    header, rows = _read_rows(path)
    assert rows, "Metadata TSV must contain at least one data row"
    for column in REQUIRED_COLUMNS:
        assert column in header, f"Missing required column: {column}"


def test_public_ramda_dataset_metadata_covers_ramda_shin_ramda_and_scrr() -> None:
    path = Path("testdata/public_datasets/ramda_datasets.tsv")
    _, rows = _read_rows(path)

    protocols = {row.get("protocol", "").strip().lower() for row in rows}
    dataset_ids = {row.get("dataset_id", "").strip().lower() for row in rows}

    assert "ramda" in protocols
    assert "shin-ramda" in protocols
    assert "screpli-ramda-seq" in protocols
    assert any(dataset_id.startswith("gse278") for dataset_id in dataset_ids)


def test_human_ramda_candidate_runs_fixture_exists() -> None:
    path = Path("testdata/public_datasets/human_ramda_candidate_runs.tsv")
    assert path.exists(), f"Missing candidate fixture TSV: {path}"

    header, rows = _read_rows(path)
    assert rows, "Human candidate fixture must contain at least one data row"
    for column in [
        "dataset_id",
        "accession_kind",
        "accession",
        "sample_or_system",
        "protocol",
        "organism",
        "likely_read_layout",
        "public_source",
        "notes",
    ]:
        assert column in header, f"Missing required column: {column}"
