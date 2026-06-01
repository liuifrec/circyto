from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def test_normalize_synthetic_scomatic_like_fixture(tmp_path: Path) -> None:
    raw = tmp_path / "external_scomatic.tsv"
    raw.write_text(
        "CB\tCHROM\tPOS\tREF\tALT\tGene\tSComatic_filter\tRead_count_ALT\tVAF\n"
        "cellA\tchr21\t101\tA\tG\tGENE1\tPASS\t4\t0.25\n",
        encoding="utf-8",
    )
    annotation = tmp_path / "cell_annotations.tsv"
    annotation.write_text("Index\tCell_type\ncellA\tRamDA\n", encoding="utf-8")
    provenance = tmp_path / "scomatic_run_metadata.json"
    provenance.write_text(
        json.dumps({"scomatic_version": "synthetic-test", "reference": "chr21-mini"}),
        encoding="utf-8",
    )
    outdir = tmp_path / "normalized"

    result = runner.invoke(
        app,
        [
            "normalize-scomatic-results",
            "--scomatic-output",
            str(raw),
            "--cell-annotations",
            str(annotation),
            "--provenance-metadata",
            str(provenance),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["description"] == "Normalized external SComatic output files; SComatic was not run."
    assert "RNA-derived candidate variant signals" in payload["terminology_note"]
    assert payload["provenance_metadata"]["content"]["scomatic_version"] == "synthetic-test"
    assert payload["cell_annotation"]["rows"] == 1

    normalized = pd.read_csv(outdir / "scomatic_candidate_summary.tsv", sep="\t", keep_default_na=False)
    assert list(normalized.columns) == [
        "variant_id",
        "cell_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene",
        "filter_status",
        "candidate_variant_class",
        "read_support",
        "vaf",
        "caller",
    ]
    row = normalized.iloc[0].to_dict()
    assert row["variant_id"] == "chr21:101:A>G:cellA"
    assert row["cell_id"] == "cellA"
    assert row["filter_status"] == "PASS"
    assert row["candidate_variant_class"] == "RNA-derived candidate variant signal"
    assert row["read_support"] == 4
    assert row["vaf"] == 0.25
    assert row["caller"] == "SComatic"

    summary = json.loads((outdir / "normalize_scomatic_results_summary.json").read_text(encoding="utf-8"))
    assert summary["source_files"][0]["resolved_required_columns"]["cell_id"] == "CB"
    assert summary["output_scomatic_candidate_summary"].endswith("scomatic_candidate_summary.tsv")


def test_normalize_multiple_external_files(tmp_path: Path) -> None:
    raw_a = tmp_path / "scomatic_a.tsv"
    raw_b = tmp_path / "scomatic_b.tsv"
    header = "cell_id\tchrom\tpos\tref\talt\tgene\tfilter_status\tread_support\tvaf\n"
    raw_a.write_text(header + "cellA\tchr21\t101\tA\tG\tGENE1\tPASS\t4\t0.25\n", encoding="utf-8")
    raw_b.write_text(header + "cellB\tchr21\t102\tC\tT\tGENE2\tPASS\t2\t0.10\n", encoding="utf-8")
    outdir = tmp_path / "normalized"

    result = runner.invoke(
        app,
        [
            "normalize-scomatic-results",
            "--scomatic-output",
            str(raw_a),
            "--scomatic-output",
            str(raw_b),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_input_files"] == 2
    assert payload["n_candidates"] == 2
    normalized = pd.read_csv(outdir / "scomatic_candidate_summary.tsv", sep="\t")
    assert normalized["cell_id"].tolist() == ["cellA", "cellB"]


def test_normalize_fails_on_missing_required_columns(tmp_path: Path) -> None:
    raw = tmp_path / "external_scomatic.tsv"
    raw.write_text(
        "CB\tCHROM\tPOS\tREF\tALT\tGene\tSComatic_filter\tRead_count_ALT\n"
        "cellA\tchr21\t101\tA\tG\tGENE1\tPASS\t4\n",
        encoding="utf-8",
    )
    result = runner.invoke(
        app,
        [
            "normalize-scomatic-results",
            "--scomatic-output",
            str(raw),
            "--outdir",
            str(tmp_path / "normalized"),
        ],
    )

    assert result.exit_code != 0
    assert "missing required SComatic columns" in result.output
    assert "vaf" in result.output
    assert "Available columns" in result.output


def test_normalize_fails_when_cell_annotation_does_not_match(tmp_path: Path) -> None:
    raw = tmp_path / "external_scomatic.tsv"
    raw.write_text(
        "cell_id\tchrom\tpos\tref\talt\tgene\tfilter_status\tread_support\tvaf\n"
        "cellA\tchr21\t101\tA\tG\tGENE1\tPASS\t4\t0.25\n",
        encoding="utf-8",
    )
    annotation = tmp_path / "cell_annotations.tsv"
    annotation.write_text("Index\tCell_type\ncellB\tRamDA\n", encoding="utf-8")
    result = runner.invoke(
        app,
        [
            "normalize-scomatic-results",
            "--scomatic-output",
            str(raw),
            "--cell-annotations",
            str(annotation),
            "--outdir",
            str(tmp_path / "normalized"),
        ],
    )

    assert result.exit_code != 0
    assert "missing SComatic candidate cell IDs: cellA" in result.output
