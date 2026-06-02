from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.scomatic_normalization import normalize_scomatic_calling_table


runner = CliRunner()


SCOMATIC_CALLING_HEADER = (
    "#CHROM\tStart\tEnd\tREF\tALT\tFILTER\tCell_types\tUp_context\tDown_context\t"
    "N_ALT\tDp\tNc\tBc\tCc\tVAF\tCCF\tBCp\tCCp\tCell_types_min_BC\t"
    "Cell_types_min_CC\tRest_BC\tRest_CC\tFisher_p\tCell_type_Filter\tINFO\tHAP1\n"
)


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


def test_normalize_real_scomatic_step2_calling_output(tmp_path: Path) -> None:
    step2 = tmp_path / "hap1_batch10.HAP1.step2.calling.step2.tsv"
    step2.write_text(
        "##INFO=ALT,Description=Alternative alleles found\n"
        "##INFO=FILTER,Description=Filter status of the variant site\n"
        + SCOMATIC_CALLING_HEADER
        + "chr1\t100\t100\tA\tG\tPASS\tHAP1\tNNNNN\tNNNNN\t1\t20\t5\t6\t2\t0.3\t0.4\t0.0\t0.0\t1\t1\t0;14;1\t0;3;1\t.\tPASS\tDP|NC|CC|BC|BQ|BCf|BCr\t20|5|2:0:0:0:0:0|6:0:0:14:0:0|200:0:0:400:0:0|3:0:0:3:0:0|3:0:0:3:0:0\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "normalized"

    result = runner.invoke(
        app,
        [
            "normalize-scomatic-results",
            "--scomatic-output",
            str(step2),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_candidates"] == 1
    assert payload["source_files"][0]["parser"] == "scomatic-calling-step"
    assert payload["source_files"][0]["scomatic_stage"] == "step2"
    assert payload["source_files"][0]["metadata_line_count"] == 2
    assert "Cell_types" in payload["source_files"][0]["cell_id_semantics"]

    normalized = pd.read_csv(outdir / "scomatic_candidate_summary.tsv", sep="\t", keep_default_na=False)
    row = normalized.iloc[0].to_dict()
    assert row["variant_id"] == "chr1:100:A>G:HAP1"
    assert row["cell_id"] == "HAP1"
    assert row["filter_status"] == "PASS"
    assert row["gene"] == "."
    assert row["read_support"] == 6
    assert row["vaf"] == 0.3
    assert row["caller"] == "SComatic-step2"

    provenance = json.loads(
        (outdir / "scomatic_candidate_summary.tsv.provenance.json").read_text(encoding="utf-8")
    )
    assert provenance["output_scomatic_candidate_summary"].endswith("scomatic_candidate_summary.tsv")
    assert provenance["source_files"][0]["path"].endswith("hap1_batch10.HAP1.step2.calling.step2.tsv")


def test_normalize_real_scomatic_step2_validates_annotations_by_cell_type(tmp_path: Path) -> None:
    step2 = tmp_path / "hap1_batch10.HAP1.step2.calling.step2.tsv"
    step2.write_text(
        "##INFO=ALT,Description=Alternative alleles found\n"
        + SCOMATIC_CALLING_HEADER
        + "chr1\t100\t100\tA\tG\tPASS\tHAP1\tNNNNN\tNNNNN\t1\t20\t5\t6\t2\t0.3\t0.4\t0.0\t0.0\t1\t1\t0;14;1\t0;3;1\t.\tPASS\tDP|NC|CC|BC|BQ|BCf|BCr\t20|5|2:0:0:0:0:0|6:0:0:14:0:0|200:0:0:400:0:0|3:0:0:3:0:0|3:0:0:3:0:0\n",
        encoding="utf-8",
    )
    annotation = tmp_path / "cell_annotations.tsv"
    annotation.write_text("Index\tCell_type\ncellA\tHAP1\ncellB\tHAP1\n", encoding="utf-8")
    outdir = tmp_path / "normalized"

    result = runner.invoke(
        app,
        [
            "normalize-scomatic-results",
            "--scomatic-output",
            str(step2),
            "--cell-annotations",
            str(annotation),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["candidate_cell_semantics"] == "cell_type"
    assert payload["cell_annotation"]["validated_against"] == "Cell_type"


def test_normalize_real_scomatic_step1_skips_non_candidate_rows(tmp_path: Path) -> None:
    step1 = tmp_path / "hap1_batch10.HAP1.step1.calling.step1.tsv"
    step1.write_text(
        "##INFO=ALT,Description=Alternative alleles found\n"
        + SCOMATIC_CALLING_HEADER
        + "chr1\t100\t100\tA\tG\tPASS\tHAP1\tNNNNN\tNNNNN\t1\t20\t5\t6\t2\t0.3\t0.4\t0.0\t0.0\t1\t1\t0;14;1\t0;3;1\t.\tPASS\tDP|NC|CC|BC|BQ|BCf|BCr\t20|5|2:0:0:0:0:0|6:0:0:14:0:0|200:0:0:400:0:0|3:0:0:3:0:0|3:0:0:3:0:0\n"
        + "chr1\t101\t101\tC\t.\t.\t.\tNNNNN\tNNNNN\t.\t.\t.\t.\t.\t.\t.\t.\t.\t1\t1\t0;0;1\t0;0;1\t.\t.\tDP|NC|CC|BC|BQ|BCf|BCr\t10|4|0:4:0:0:0:0|0:10:0:0:0:0|0:400:0:0:0:0|0:5:0:0:0:0|0:5:0:0:0:0\n",
        encoding="utf-8",
    )

    normalized, summary = normalize_scomatic_calling_table(step1)

    assert summary["rows"] == 2
    assert summary["candidate_rows"] == 1
    assert summary["normalized_candidate_rows"] == 1
    assert summary["scomatic_stage"] == "step1"
    assert normalized["variant_id"].tolist() == ["chr1:100:A>G:HAP1"]
    assert normalized["caller"].tolist() == ["SComatic-step1"]


def test_normalize_real_scomatic_calling_output_flattens_multiallelic_cell_type(tmp_path: Path) -> None:
    step2 = tmp_path / "hap1_batch10.HAP1.step2.calling.step2.tsv"
    step2.write_text(
        "##INFO=ALT,Description=Alternative alleles found\n"
        + SCOMATIC_CALLING_HEADER
        + "chr2\t200\t200\tC\tA|T\tMulti-allelic\tHAP1\tNNNNN\tNNNNN\t2\t30\t8\t5|4\t2|2\t0.1667|0.1333\t0.25|0.25\t0.0|0.0\t0.0|0.0\t1\t1\t0;21;1\t0;6;1\t.\tMulti-allelic\tDP|NC|CC|BC|BQ|BCf|BCr\t30|8|0:2:2:0:0:0|5:21:4:0:0:0|150:630:120:0:0:0|2:10:2:0:0:0|3:11:2:0:0:0\n",
        encoding="utf-8",
    )

    normalized, summary = normalize_scomatic_calling_table(step2)

    assert summary["normalized_candidate_rows"] == 2
    assert normalized["alt"].tolist() == ["A", "T"]
    assert normalized["read_support"].tolist() == [5, 4]
    assert normalized["vaf"].tolist() == [0.1667, 0.1333]
