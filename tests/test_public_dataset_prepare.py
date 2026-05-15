from __future__ import annotations

import csv
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline import public_dataset_prepare as prep
from circyto.pipeline.public_dataset_prepare import get_public_dataset_spec, prepare_public_dataset


runner = CliRunner()


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_prepare_public_dataset_dry_run_creates_planning_files(tmp_path: Path) -> None:
    outdir = tmp_path / "gse98664_plan"

    result = runner.invoke(
        app,
        [
            "prepare-public-dataset",
            "--dataset-id",
            "GSE98664",
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--download-method",
            "sra",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert (outdir / "selected_runs.tsv").exists()
    assert (outdir / "download_plan.sh").exists()
    assert (outdir / "README_next_steps.md").exists()


def test_prepare_public_dataset_max_runs_limits_selected_runs(tmp_path: Path) -> None:
    outdir = tmp_path / "gse98664_plan"

    result = runner.invoke(
        app,
        [
            "prepare-public-dataset",
            "--dataset-id",
            "GSE98664",
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--download-method",
            "sra",
            "--dry-run",
            "--max-runs",
            "2",
        ],
    )

    assert result.exit_code == 0, result.stdout
    rows = _read_tsv_rows(outdir / "selected_runs.tsv")
    assert len(rows) == 2
    assert [row["run_id"] for row in rows] == [
        "SRR5516584",
        "SRR5516900",
    ]


def test_prepare_public_dataset_uses_curated_table_when_present_for_e_mtab_8735(tmp_path: Path) -> None:
    outdir = tmp_path / "emtab8735_plan"

    result = runner.invoke(
        app,
        [
            "prepare-public-dataset",
            "--dataset-id",
            "E-MTAB-8735",
            "--outdir",
            str(outdir),
            "--protocol",
            "smartseq3",
            "--download-method",
            "ena",
            "--dry-run",
            "--max-runs",
            "1",
        ],
    )

    assert result.exit_code == 0, result.stdout
    rows = _read_tsv_rows(outdir / "selected_runs.tsv")
    assert rows[0]["dataset_id"] == "E-MTAB-8735"
    assert rows[0]["protocol"] == "smartseq3"
    assert rows[0]["organism"] == "Homo sapiens"
    assert rows[0]["expected_reference"] == "hg38"
    assert rows[0]["expected_read_layout"] == "paired-end"
    assert rows[0]["recommended_route"] == "Smart-seq3 paired-end workflow"
    assert rows[0]["run_id"] == "TODO_REAL_ACCESSION_EMTAB8735_001"
    assert rows[0]["source"] == "Local curated planning table for E-MTAB-8735"
    assert "TODO_REAL_ACCESSION" in rows[0]["notes"]
    plan_text = (outdir / "download_plan.sh").read_text(encoding="utf-8")
    assert "Curated ENA/ArrayExpress-style URLs" in plan_text


def test_prepare_public_dataset_uses_curated_table_when_present_for_gse98664(tmp_path: Path) -> None:
    outdir = tmp_path / "gse98664_plan"

    result = runner.invoke(
        app,
        [
            "prepare-public-dataset",
            "--dataset-id",
            "GSE98664",
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--download-method",
            "sra",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert "WARNING:" in result.stdout
    assert "Dataset GSE98664 is annotated as Mus musculus." in result.stdout
    assert "Recommended references: mm10/mm39." in result.stdout
    assert "Do not use hg38 for biological validation." in result.stdout
    rows = _read_tsv_rows(outdir / "selected_runs.tsv")
    assert rows[0]["dataset_id"] == "GSE98664"
    assert rows[0]["organism"] == "Mus musculus"
    assert rows[0]["expected_reference"] == "mm10/mm39"
    assert rows[0]["expected_read_layout"] == "single-end"
    assert rows[0]["recommended_route"] == "BWA+CIRI3 single-end"
    assert rows[0]["run_id"] == "SRR5516584"
    assert "Real public SRA run" in rows[0]["notes"]
    plan_text = (outdir / "download_plan.sh").read_text(encoding="utf-8")
    assert "# Organism: Mus musculus" in plan_text
    assert "# Expected reference: mm10/mm39" in plan_text
    assert "# WARNING:" in plan_text
    assert "prefetch SRR5516584" in plan_text
    assert "fasterq-dump SRR5516584 --split-files --threads 8" in plan_text
    readme_text = (outdir / "README_next_steps.md").read_text(encoding="utf-8")
    assert "- Organism: Mus musculus" in readme_text
    assert "- Expected reference: mm10/mm39" in readme_text
    assert "## Warning" in readme_text


def test_prepare_public_dataset_placeholders_are_clearly_marked_when_fallback_is_used(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    outdir = tmp_path / "gse278944_plan"

    def fake_curated_path(spec):
        return tmp_path / "missing_curated.tsv"

    monkeypatch.setattr(prep, "_curated_run_table_path", fake_curated_path)

    summary = prepare_public_dataset(
        dataset_id="GSE278944",
        outdir=outdir,
        max_runs=1,
        dry_run=True,
        protocol="scrr",
        download_method="sra",
    )

    rows = _read_tsv_rows(outdir / "selected_runs.tsv")
    assert summary["row_mode"] == "placeholder"
    assert rows[0]["run_id"] == "TODO_REAL_ACCESSION_GSE278944_001"
    assert "TODO_REAL_ACCESSION" in rows[0]["notes"]
    readme_text = (outdir / "README_next_steps.md").read_text(encoding="utf-8")
    assert "Embedded placeholder rows were used instead" in readme_text


def test_prepare_public_dataset_unknown_dataset_id_gives_clear_error(tmp_path: Path) -> None:
    outdir = tmp_path / "unknown_plan"

    result = runner.invoke(
        app,
        [
            "prepare-public-dataset",
            "--dataset-id",
            "GSE000000",
            "--outdir",
            str(outdir),
            "--protocol",
            "ramda",
            "--download-method",
            "sra",
            "--dry-run",
        ],
    )

    assert result.exit_code != 0
    rendered = result.stdout + (result.stderr or "")
    assert "Unknown dataset id 'GSE000000'" in rendered
    assert "E-MTAB-8735" in rendered
    assert "GSE98664" in rendered


def test_prepare_public_dataset_does_not_attempt_real_downloads(tmp_path: Path) -> None:
    outdir = tmp_path / "shin_ramda_plan"

    summary = prepare_public_dataset(
        dataset_id="shin-ramda-riken",
        outdir=outdir,
        max_runs=1,
        dry_run=True,
        protocol="shin-ramda",
        download_method="none",
    )

    assert summary["selected_run_count"] == 1
    assert (outdir / "selected_runs.tsv").exists()
    assert (outdir / "download_plan.sh").exists()
    assert (outdir / "README_next_steps.md").exists()
    assert not (outdir / "fastq").exists()
    assert not (outdir / "sra").exists()
    plan_text = (outdir / "download_plan.sh").read_text(encoding="utf-8")
    assert "No download commands are generated" in plan_text


def test_public_dataset_specs_surface_scientific_metadata() -> None:
    mouse = get_public_dataset_spec("GSE98664")
    human = get_public_dataset_spec("E-MTAB-8735")

    assert mouse.organism == "Mus musculus"
    assert mouse.expected_reference == "mm10/mm39"
    assert mouse.expected_read_layout == "single-end"
    assert mouse.recommended_route == "BWA+CIRI3 single-end"

    assert human.organism == "Homo sapiens"
    assert human.expected_reference == "hg38"
    assert human.expected_read_layout == "paired-end"
    assert human.recommended_route == "Smart-seq3 paired-end workflow"
