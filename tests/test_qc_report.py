from __future__ import annotations

import csv
import hashlib
import json
from html.parser import HTMLParser
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.qc_report import collect_report_metrics, generate_report


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    root = tmp_path / "results with spaces"
    (root / "qc").mkdir(parents=True)
    (root / "matrix").mkdir()
    (root / "qc/cell_qc.tsv").write_text(
        "cell_id\tcircRNA_count\ttotal_circRNA_support\tdetector_status\talignment_status\n"
        "cellA\t2\t7\tsuccess\taligned\ncellB\t1\t3\tsuccess\taligned\ncellC\t0\t0\tempty\taligned\n",
        encoding="utf-8")
    (root / "qc/circ_qc.tsv").write_text(
        "circ_id\tn_cells_detected\ttotal_support\thost_gene\n"
        "chr1:10|20\t2\t8\tGENE1\nchr2:30|40\t1\t2\t\n", encoding="utf-8")
    (root / "matrix/circ_counts.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 5\n1 2 3\n2 1 2\n")
    (root / "matrix/circ_index.txt").write_text("chr1:10|20\nchr2:30|40\n")
    (root / "matrix/cell_index.txt").write_text("cellA\ncellB\n")
    write_json(root / "workflow_summary.json", {
        "workflow_type": "full-length-circrna", "protocol": "ramda", "detector_backend": "ciri3",
        "circyto_version": "0.10.0", "completed_at": "2026-01-01T00:01:00Z", "dry_run": False,
        "matrix": {"n_cells": 3, "n_circRNAs": 2}, "planned_cells": 3,
        "stage_graph": [{"stage": "alignment", "status": "completed"}, {"stage": "summary_qc", "status": "completed"}],
        "alignment_status_counts": {"aligned": 3}, "detector_status_counts": {"success": 2, "empty": 1},
        "hostname": "private-host.example", "command_options": {"token": "DO-NOT-EXPORT"},
        "genome_fasta": "/private/user/reference/genome.fa", "paths": {"matrix": "/old/machine/results/matrix/circ_counts.mtx"},
    })
    return root


def snapshot(root: Path) -> dict[str, str]:
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file() and path.name not in {"report.html", "metrics.json"}}


def test_cli_exact_metrics_preserves_inputs_and_is_reproducible(workdir: Path, monkeypatch):
    # Reporting must not invoke any subprocess, aligner, or detector.
    import subprocess
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: pytest.fail("report launched a subprocess"))
    before = snapshot(workdir)
    result = CliRunner().invoke(app, ["report", "--workdir", str(workdir)])
    assert result.exit_code == 0, result.output
    assert "Recorded workflow status: completed" in result.output
    payload = json.loads((workdir / "qc/metrics.json").read_text())
    metrics = payload["metrics"]
    assert metrics["cells"]["value"] == 3
    assert metrics["candidates"]["value"] == 2
    assert metrics["total_support"]["value"] == metrics["candidate_total_support"]["value"] == 10
    assert metrics["median_candidates_per_cell"]["value"] == 1
    assert metrics["zero_candidate_fraction"]["numerator"] == 1
    assert metrics["zero_candidate_fraction"]["denominator"] == 3
    assert metrics["zero_candidate_fraction"]["value"] == 1 / 3
    assert metrics["host_gene_coverage"]["value"] == 0.5
    assert payload["per_candidate"][0]["prevalence"] == 2 / 3
    with (workdir / "qc/cell_qc.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["candidate_count"] for row in payload["per_cell"]] == [int(row["circRNA_count"]) for row in rows]
    assert [row["total_support"] for row in payload["per_cell"]] == [int(row["total_circRNA_support"]) for row in rows]
    for source in payload["sources"]:
        if source["available"]:
            assert source["sha256"] == before[source["path"]]
    assert snapshot(workdir) == before
    first = [(workdir / name).read_bytes() for name in ("qc/report.html", "qc/metrics.json")]
    generate_report(workdir)
    assert first == [(workdir / name).read_bytes() for name in ("qc/report.html", "qc/metrics.json")]
    exported = "".join(part.decode() for part in first)
    assert "private-host" not in exported and "DO-NOT-EXPORT" not in exported
    assert "/private/" not in exported and "/old/machine" not in exported
    assert "genome.fa" in exported


@pytest.mark.parametrize("name", ["workflow_summary.json", "qc/cell_qc.tsv", "qc/circ_qc.tsv"])
def test_missing_sources_are_not_zero(workdir, name):
    (workdir / name).unlink()
    payload = generate_report(workdir)
    assert payload["workflow"]["status"] != "completed"
    if name.endswith("cell_qc.tsv"):
        assert payload["metrics"]["cells"]["value"] is None
        assert payload["metrics"]["total_support"]["value"] is None
        assert all(row["prevalence"] is None for row in payload["per_candidate"])
    if name.endswith("circ_qc.tsv"):
        assert payload["metrics"]["candidates"]["value"] is None
        assert payload["metrics"]["host_gene_coverage"]["value"] is None
    assert "Not available" in (workdir / "qc/report.html").read_text()


@pytest.mark.parametrize("value", ["", "NA", "NaN", "inf", "-1", "1.25", "nonsense"])
def test_missing_or_invalid_numeric_cell_is_never_imputed(workdir, value):
    path = workdir / "qc/cell_qc.tsv"
    path.write_text(path.read_text().replace("cellA\t2\t7", f"cellA\t{value}\t7"))
    payload = generate_report(workdir)
    assert payload["per_cell"][0]["candidate_count"] is None
    assert payload["metrics"]["zero_candidate_fraction"]["value"] is None
    assert payload["metrics"]["median_candidates_per_cell"]["value"] is None
    assert payload["metrics"]["total_support"]["value"] == 10


def test_missing_host_column_not_invented(workdir):
    (workdir / "qc/circ_qc.tsv").write_text("circ_id\tn_cells_detected\ttotal_support\nchr1:10|20\t2\t8\nchr2:30|40\t1\t2\n")
    payload = generate_report(workdir)
    assert payload["metrics"]["host_gene_annotated_candidates"]["value"] is None
    assert payload["metrics"]["host_gene_coverage"]["value"] is None


@pytest.mark.parametrize("source", ["cell", "stage", "summary", "optional_detector", "optional_alignment"])
def test_failures_override_completed_stages(workdir, source):
    path = workdir / "workflow_summary.json"
    summary = json.loads(path.read_text())
    if source == "cell":
        cell = workdir / "qc/cell_qc.tsv"
        cell.write_text(cell.read_text().replace("empty", "failed"))
    elif source == "stage":
        summary["failed_stages"] = ["detector"]
    elif source == "summary":
        summary["detector"] = {"failed": 1}
    else:
        name = "ciri3/detector_run_summary.json" if source == "optional_detector" else "align/alignment_prepare_summary.json"
        write_json(workdir / name, {"status_counts": {"success": 2, "failed": 1}})
    write_json(path, summary)
    payload = generate_report(workdir)
    assert payload["workflow"]["status"] == "failed"
    assert "Workflow: Failed" in (workdir / "qc/report.html").read_text()


@pytest.mark.parametrize("status", ["missing", "pending", "unknown", "running"])
def test_incomplete_cells_never_mean_success(workdir, status):
    path = workdir / "qc/cell_qc.tsv"
    path.write_text(path.read_text().replace("empty", status))
    assert generate_report(workdir)["workflow"]["status"] == "incomplete"


def test_dry_run_with_old_qc_does_not_mean_success(workdir):
    path = workdir / "workflow_summary.json"
    value = json.loads(path.read_text())
    value["dry_run"] = True
    write_json(path, value)
    assert generate_report(workdir)["workflow"]["status"] == "dry_run"


def test_zero_candidates_has_zero_counts_but_no_host_fraction(workdir):
    (workdir / "qc/cell_qc.tsv").write_text("cell_id\tcircRNA_count\ttotal_circRNA_support\tdetector_status\talignment_status\ncellA\t0\t0\tempty\taligned\n")
    (workdir / "qc/circ_qc.tsv").write_text("circ_id\tn_cells_detected\ttotal_support\thost_gene\n")
    payload = generate_report(workdir)
    assert payload["metrics"]["total_support"]["value"] == 0
    assert payload["metrics"]["candidates"]["value"] == 0
    assert payload["metrics"]["zero_candidate_fraction"]["value"] == 1
    assert payload["metrics"]["host_gene_coverage"]["value"] is None


def test_zero_cells_never_produces_nan_or_division_by_zero(workdir):
    (workdir / "qc/cell_qc.tsv").write_text("cell_id\tcircRNA_count\ttotal_circRNA_support\n")
    payload = generate_report(workdir)
    assert payload["metrics"]["cells"]["value"] == 0
    assert payload["metrics"]["zero_candidate_fraction"]["value"] is None
    assert all(row["prevalence"] is None for row in payload["per_candidate"])
    assert "NaN" not in (workdir / "qc/metrics.json").read_text()


@pytest.mark.parametrize("body", ["{broken", "[]", "null"])
def test_invalid_json_still_gives_diagnostic_report(workdir, body):
    (workdir / "workflow_summary.json").write_text(body)
    assert generate_report(workdir)["workflow"]["status"] == "unknown"


@pytest.mark.parametrize("body", ["cell_id\tcircRNA_count\nA\t1\nA\t2\n", "cell_id\tcircRNA_count\nA\n", "cell_id\tcell_id\nA\tB\n", "wrong\nA\n"])
def test_malformed_tsv_never_silently_changes_population(workdir, body):
    (workdir / "qc/cell_qc.tsv").write_text(body)
    payload = generate_report(workdir)
    assert payload["metrics"]["cells"]["value"] is None
    assert payload["workflow"]["status"] != "completed"


def test_inconsistency_keeps_source_values_and_disables_prevalence(workdir):
    path = workdir / "qc/circ_qc.tsv"
    path.write_text(path.read_text().replace("\t2\t8", "\t4\t99"))
    payload = generate_report(workdir)
    assert payload["metrics"]["total_support"]["value"] == 10
    assert payload["metrics"]["candidate_total_support"]["value"] == 101
    assert payload["per_candidate"][0]["n_cells_detected"] == 4
    assert payload["per_candidate"][0]["prevalence"] is None
    assert payload["workflow"]["status"] == "incomplete"
    assert any(row["ok"] is False for row in payload["consistency_checks"])


class Tags(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags = []
        self.attrs = []

    def handle_starttag(self, tag, attrs):
        self.tags.append(tag)
        self.attrs.extend(attrs)


def test_untrusted_text_escaped_and_html_is_offline(workdir):
    attack = '<script>alert("XSS")</script><img src=x onerror=alert(1)>'
    path = workdir / "qc/cell_qc.tsv"
    path.write_text(path.read_text().replace("cellA", attack))
    path = workdir / "qc/circ_qc.tsv"
    path.write_text(path.read_text().replace("GENE1", attack))
    path = workdir / "workflow_summary.json"
    summary = json.loads(path.read_text())
    summary.update(protocol=attack, warnings=[attack, "Log: /private/alice/logs/run.log", "See https://internal.example/?token=secret"])
    write_json(path, summary)
    generate_report(workdir)
    page = (workdir / "qc/report.html").read_text()
    parser = Tags()
    parser.feed(page)
    assert "script" not in parser.tags and "img" not in parser.tags and "iframe" not in parser.tags
    assert "svg" in parser.tags
    assert not any(key.startswith("on") or key == "src" for key, value in parser.attrs)
    assert all(value.startswith("#") for key, value in parser.attrs if key == "href")
    assert "&lt;script&gt;" in page and "&quot;XSS&quot;" in page
    assert "/private/alice" not in page and "internal.example" not in page
    assert "Content-Security-Policy" in page


@pytest.mark.parametrize("name", ["report.html", "metrics.json"])
def test_symlink_destinations_cannot_overwrite_source(workdir, name):
    before = snapshot(workdir)
    (workdir / "qc" / name).symlink_to(workdir / "workflow_summary.json")
    with pytest.raises(ValueError, match="unsafe report destination"):
        generate_report(workdir)
    assert snapshot(workdir) == before


def test_hardlinked_destination_does_not_modify_source(workdir):
    before = snapshot(workdir)
    (workdir / "qc/report.html").hardlink_to(workdir / "workflow_summary.json")
    generate_report(workdir)
    assert snapshot(workdir) == before


def test_symlink_qc_directory_is_refused(tmp_path):
    root = tmp_path / "run"
    other = tmp_path / "other"
    root.mkdir()
    other.mkdir()
    (root / "qc").symlink_to(other, target_is_directory=True)
    with pytest.raises(ValueError, match="real directory"):
        generate_report(root)
    assert not list(other.iterdir())


def test_outside_source_is_not_read(workdir, tmp_path):
    secret = tmp_path / "secret.json"
    secret.write_text('{"protocol": "SENSITIVE"}')
    (workdir / "workflow_summary.json").unlink()
    (workdir / "workflow_summary.json").symlink_to(secret)
    payload = generate_report(workdir)
    assert payload["workflow"]["status"] == "unknown"
    assert "SENSITIVE" not in (workdir / "qc/report.html").read_text()


def test_write_failure_surfaces_and_preserves_previous_report(workdir, monkeypatch):
    import circyto.pipeline.qc_report as report
    generate_report(workdir)
    before = {path: path.read_bytes() for path in workdir.rglob("*") if path.is_file()}
    def fail(*args, **kwargs):
        raise OSError("simulated disk error")
    monkeypatch.setattr(report.tempfile, "NamedTemporaryFile", fail)
    result = CliRunner().invoke(app, ["report", "--workdir", str(workdir)])
    assert result.exit_code == 1 and "simulated disk error" in result.output
    assert all(path.read_bytes() == value for path, value in before.items())
    assert not list((workdir / "qc").glob(".report-*"))


def test_cli_help_and_path_validation(tmp_path):
    runner = CliRunner()
    top = runner.invoke(app, ["--help"])
    assert top.exit_code == 0
    for value in ("Start here", "report", "doctor", "full-length-circrna", "smartseq3-ciri3"):
        assert value in top.output
    report = runner.invoke(app, ["report", "--help"])
    assert report.exit_code == 0 and "--workdir" in report.output
    for workflow in ("full-length-circrna", "smartseq3-ciri3"):
        help_result = runner.invoke(app, ["workflow", workflow, "--help"])
        assert "--no-report" in help_result.output
    legacy = runner.invoke(app, ["run", "--help"])
    assert "--fastq-dir" in legacy.output and "LEGACY" in legacy.output
    assert runner.invoke(app, ["report", "--workdir", str(tmp_path / "missing")]).exit_code != 0
    assert not (tmp_path / "missing").exists()


def test_collector_does_not_write(workdir):
    before = snapshot(workdir)
    collect_report_metrics(workdir)
    assert snapshot(workdir) == before
    assert not (workdir / "qc/report.html").exists()


def test_original_identifiers_preserved(workdir):
    path = workdir / "qc/cell_qc.tsv"
    path.write_text(path.read_text().replace("cellA", "/plate/A01"))
    assert generate_report(workdir)["per_cell"][0]["cell_id"] == "/plate/A01"


def test_optional_cell_failure_overrides_stale_status_counts(workdir):
    write_json(workdir / "ciri3/detector_run_summary.json", {
        "status_counts": {"success": 2, "empty": 1},
        "cells": [{"cell_id": "cellA", "status": "failed"}],
    })
    assert generate_report(workdir)["workflow"]["status"] == "failed"


@pytest.mark.parametrize("status_counts", [{"failed": "NaN", "success": 3}, ["success"]])
def test_invalid_status_counts_never_mean_success(workdir, status_counts):
    write_json(workdir / "ciri3/detector_run_summary.json", {"status_counts": status_counts})
    assert generate_report(workdir)["workflow"]["status"] == "incomplete"


def test_current_running_status_overrides_old_completion_timestamp(workdir):
    path = workdir / "workflow_summary.json"
    value = json.loads(path.read_text())
    value["status"] = "running"
    write_json(path, value)
    assert generate_report(workdir)["workflow"]["status"] == "incomplete"


def test_large_integer_support_is_exact(workdir):
    value = 2**55 + 7
    (workdir / "qc/cell_qc.tsv").write_text(f"cell_id\tcircRNA_count\ttotal_circRNA_support\nA\t1\t{value}\n")
    (workdir / "qc/circ_qc.tsv").write_text(f"circ_id\tn_cells_detected\ttotal_support\nC\t1\t{value}\n")
    payload = generate_report(workdir)
    assert payload["metrics"]["total_support"]["value"] == value
    assert payload["metrics"]["candidate_total_support"]["value"] == value


def test_source_alias_to_output_is_refused(workdir):
    target = workdir / "qc/report.html"
    source = workdir / "qc/cell_qc.tsv"
    target.write_bytes(source.read_bytes())
    source.unlink()
    source.symlink_to(target)
    before = source.read_bytes()
    with pytest.raises(ValueError, match="source artifact resolves"):
        generate_report(workdir)
    assert source.read_bytes() == before


def test_documentation_demo_matches_matrix_and_checked_in_metrics(tmp_path):
    import runpy
    import numpy as np
    from scipy.io import mmread

    repo = Path(__file__).resolve().parents[1]
    demo = runpy.run_path(str(repo / "examples/qc_report_demo.py"))
    root = tmp_path / "demo"
    demo["create_demo"](root)
    with pytest.raises(ValueError, match="never overwrites"):
        demo["create_demo"](root)
    payload = generate_report(root)
    matrix = mmread(root / "matrix/circ_counts.mtx").toarray()
    assert matrix.shape == (6, 12)
    assert np.sum(matrix) == payload["metrics"]["total_support"]["value"] == 77
    assert np.sum(matrix, axis=0).tolist() == [row["total_support"] for row in payload["per_cell"]]
    assert np.sum(matrix, axis=1).tolist() == [row["total_support"] for row in payload["per_candidate"]]
    assert np.count_nonzero(matrix, axis=0).tolist() == [row["candidate_count"] for row in payload["per_cell"]]
    assert np.count_nonzero(matrix, axis=1).tolist() == [row["n_cells_detected"] for row in payload["per_candidate"]]
    assert payload == json.loads((repo / "docs/examples/qc_report/qc/metrics.json").read_text())
