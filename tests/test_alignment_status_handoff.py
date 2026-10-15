"""Exercise the per-aligner -> combined -> QC/AnnData -> HTML handoff.

The positive case uses the real per-aligner writer with a tiny fake executable.
Other cases supply child-summary evidence, never an already-correct combined
summary. Positive counts deliberately remain present in non-completed controls.
"""

from collections import Counter
import json
import os
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite

from circyto.pipeline.align_manifest import prepare_alignment_cache, read_source_manifest
from circyto.pipeline.qc_report import generate_report
from circyto.pipeline.run_ciri3 import _write_alignment_summary
from circyto.pipeline.workflow_full_length_circrna import (
    FullLengthCircRNAWorkflowParams,
    _build_execution_summary,
    _read_manifest_rows,
    _workflow_paths,
)


def _json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n")


def _case(tmp_path, layouts=("single-end",)):
    root = tmp_path / "workflow"
    paths = _workflow_paths(root)
    for name, path in paths.items():
        if name != "workflow_summary":
            path.mkdir(parents=True, exist_ok=True)
    read = tmp_path / "read.fastq"
    read.write_text("@r1\nACGT\n+\n!!!!\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "G1";\n')
    manifest = tmp_path / "cells.tsv"
    ids = [f"cell{i}" for i in range(len(layouts))]
    manifest.write_text(
        "cell_id\tread1\tread2\tprotocol\tread_layout\tstrandedness\n"
        + "".join(
            f"{cell}\t{read}\t{read if layout == 'paired-end' else ''}\tramda\t{layout}\tunstranded\n"
            for cell, layout in zip(ids, layouts)
        )
    )
    rows = read_source_manifest(manifest, validate_files=True)
    # MatrixMarket uses candidate x cell; keep a zero entry and unequal values.
    values = np.array([[4 + i, 0 if i == 0 else 2] for i in range(len(ids))], dtype=np.int32)
    mmwrite(paths["matrix"] / "circ_counts.mtx", sparse.csr_matrix(values.T))
    (paths["matrix"] / "circ_index.txt").write_text("chr1:1|4\nchr1:2|4\n")
    (paths["matrix"] / "cell_index.txt").write_text("\n".join(ids) + "\n")
    (paths["matrix"] / "circ_feature_table.tsv").write_text(
        "circ_id\tchrom\tstart\tend\tstrand\thost_gene\n"
        "chr1:1|4\tchr1\t1\t4\t+\tG1\nchr1:2|4\tchr1\t2\t4\t-\tG2\n"
    )
    _json(paths["ciri3"] / "detector_run_summary.json", {
        "cells": [{"cell_id": cell, "status": "success", "seconds": 0.1} for cell in ids],
        "status_counts": {"success": len(ids)}, "completed_cells": len(ids), "failed_cells": 0,
    })
    return paths, rows, FullLengthCircRNAWorkflowParams(
        manifest=manifest, outdir=root, protocol="ramda", genome_fasta=ref, gtf=gtf,
    ), values


def _child(paths, records, bucket="bwa_mem", **overrides):
    # Same small record shape as the A100 per-aligner summary, without real data.
    summary = paths["align"] / bucket / "alignment_prepare_summary.json"
    payload = {"cells": records, "status_counts": dict(Counter(r.get("status", "unknown") for r in records))}
    payload.update(overrides)
    _json(summary, payload)
    return summary.with_name("alignment_manifest.tsv")


def _record(cell="cell0", status="aligned", **fields):
    return {"cell_id": cell, "status": status, "seconds": 0.2,
            "alignment_path": "retained/alignment.sam", "artifact_bucket": "bwa_mem",
            "mapper_mode": "0", "read_layout": "single-end", "cache_key": "tiny",
            "sortedness": "unsorted", "reused_alignment": status.startswith("reused"), **fields}


def _handoff(paths, rows, params, manifests):
    combined_path = paths["align"] / "alignment_prepare_summary.json"
    _write_alignment_summary(combined_path, source_rows=rows, manifest_paths=manifests, dry_run=False)
    _, raw_rows = _read_manifest_rows(params.manifest)
    summary = _build_execution_summary(
        params=params, paths=paths, source_rows=rows, raw_rows=raw_rows,
        resolved_protocols={r.cell_id: "ramda" for r in rows},
        resolved_strandedness={r.cell_id: "unstranded" for r in rows},
        skip_demux_effective=True, elapsed_seconds=0.5, warnings=[], rna_import=None,
        started_at="2026-10-15T00:00:00Z", completed_at="2026-10-15T00:00:01Z",
    )
    _json(paths["workflow_summary"], summary)
    generate_report(paths["root"])
    return (
        json.loads(combined_path.read_text()),
        pd.read_csv(paths["qc"] / "cell_qc.tsv", sep="\t", keep_default_na=False),
        ad.read_h5ad(paths["anndata"] / "circ_counts.h5ad"),
        json.loads((paths["qc"] / "metrics.json").read_text()),
    )


def test_real_single_end_writer_reaches_qc_object_and_report(tmp_path, monkeypatch):
    paths, rows, params, values = _case(tmp_path)
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    bwa = fake_bin / "bwa"
    bwa.write_text(f"#!{sys.executable}\nprint('@HD\\tVN:1.6\\tSO:unsorted')\nprint('@SQ\\tSN:chr1\\tLN:4')\n")
    bwa.chmod(0o755)
    monkeypatch.setenv("PATH", str(fake_bin) + os.pathsep + os.environ["PATH"])
    child_manifest = prepare_alignment_cache(
        manifest=params.manifest, outdir=paths["align"] / "bwa_mem", aligner="bwa-mem",
        detector_hint="ciri3", ref_fa=params.genome_fasta, threads=1, parallel=1,
    )
    child = json.loads(child_manifest.with_name("alignment_prepare_summary.json").read_text())
    assert child["cells"][0]["status"] == "aligned"
    combined, qc, obj, report = _handoff(paths, rows, params, [child_manifest, paths["align"] / "alignment_manifest.tsv"])
    assert qc.alignment_status.tolist() == ["aligned"]
    assert obj.obs.alignment_status.tolist() == ["aligned"]
    assert combined["cells"] == child["cells"]
    assert combined["status_counts"] == {"aligned": 1}
    assert combined["completed_cells"] == 1
    assert report["workflow"]["status"] == "completed"
    assert report["candidate_prevalence_denominator"] == 1
    assert report["per_cell"][0]["alignment_status"] == "aligned"
    assert "aligned" in (paths["qc"] / "report.html").read_text()
    np.testing.assert_array_equal(obj.X.toarray(), values)
    assert obj.X.dtype == np.int32
    assert obj.var_names.tolist() == ["chr1:1|4", "chr1:2|4"]


@pytest.mark.parametrize("layouts", [("paired-end",), ("single-end", "paired-end")])
def test_routed_children_and_merged_manifest_are_counted_once(tmp_path, layouts):
    paths, rows, params, values = _case(tmp_path, layouts)
    manifests = [
        _child(paths, [_record(r.cell_id, read_layout=r.read_layout)], "star" if r.read_layout == "paired-end" else "bwa_mem")
        for r in reversed(rows)
    ]
    combined, qc, obj, report = _handoff(paths, rows, params, manifests + manifests + [paths["align"] / "alignment_manifest.tsv"])
    assert [r["cell_id"] for r in combined["cells"]] == [r.cell_id for r in rows]
    assert combined["status_counts"] == {"aligned": len(rows)}
    assert report["candidate_prevalence_denominator"] == len(rows)
    np.testing.assert_array_equal(obj.X.toarray(), values)
    # Re-aggregation must not treat the previously written merged summary as a child.
    again, _, _, _ = _handoff(paths, rows, params, manifests + [paths["align"] / "alignment_manifest.tsv"])
    for key in ("cells", "status_counts", "completed_cells", "alignment_summary_sources"):
        assert again[key] == combined[key]


@pytest.mark.parametrize("status", ["failed", "unknown", "missing", None, "reused_cached", "reused_input"])
def test_status_evidence_not_positive_support_controls_completion(tmp_path, status):
    paths, rows, params, values = _case(tmp_path)
    records = [] if status is None else [_record(status=status)]
    child = _child(paths, records)
    combined, qc, obj, report = _handoff(paths, rows, params, [child])
    expected = status or "missing"
    assert qc.alignment_status.tolist() == [expected]
    assert obj.obs.alignment_status.tolist() == [expected]
    assert combined["status_counts"] == dict(Counter(r["status"] for r in records))
    completed = status in {"reused_cached", "reused_input"}
    assert combined["completed_cells"] == int(completed)
    assert report["candidate_prevalence_denominator"] == (1 if completed else None)
    assert (report["workflow"]["status"] == "completed") == completed
    np.testing.assert_array_equal(obj.X.toarray(), values)


def test_dry_run_does_not_reuse_old_success_evidence(tmp_path):
    paths, rows, _, _ = _case(tmp_path)
    child = _child(paths, [_record()])
    target = paths["align"] / "alignment_prepare_summary.json"
    _write_alignment_summary(target, source_rows=rows, manifest_paths=[child], dry_run=True)
    combined = json.loads(target.read_text())
    assert combined["cells"] == [{"cell_id": "cell0", "status": "planned"}]
    assert combined["status_counts"] == {"planned": 1}
    assert combined["completed_cells"] == 0


@pytest.mark.parametrize("case", ["duplicate", "conflicting", "foreign", "missing_id", "count_conflict", "dry_success"])
def test_ambiguous_child_evidence_is_rejected_without_writing_success(tmp_path, case):
    paths, rows, _, _ = _case(tmp_path)
    records = [_record()]
    extra = {}
    if case == "duplicate":
        records.append(_record())
    elif case == "conflicting":
        records.append(_record(status="failed"))
    elif case == "foreign":
        records.append(_record("other_cell"))
    elif case == "missing_id":
        records[0].pop("cell_id")
    elif case == "count_conflict":
        extra["status_counts"] = {"failed": 1}
    else:
        extra["dry_run"] = True
    child = _child(paths, records, **extra)
    target = paths["align"] / "alignment_prepare_summary.json"
    with pytest.raises(ValueError):
        _write_alignment_summary(target, source_rows=rows, manifest_paths=[child], dry_run=False)
    assert not target.exists()


def test_duplicate_identity_across_distinct_children_is_rejected(tmp_path):
    paths, rows, _, _ = _case(tmp_path)
    manifests = [_child(paths, [_record()], bucket) for bucket in ("bwa_mem", "star")]
    with pytest.raises(ValueError, match="[Dd]uplicate"):
        _write_alignment_summary(paths["align"] / "alignment_prepare_summary.json", source_rows=rows, manifest_paths=manifests, dry_run=False)


def test_legacy_record_only_summary_retains_reuse_status(tmp_path):
    paths, rows, params, _ = _case(tmp_path)
    child = _child(paths, [_record(status="reused_input")])
    _json(child.with_name("alignment_prepare_summary.json"), {"cells": [_record(status="reused_input")]})
    combined, _, _, report = _handoff(paths, rows, params, [child])
    assert combined["status_counts"] == {"reused_input": 1}
    assert report["candidate_prevalence_denominator"] == 1


def test_legacy_counts_only_cannot_assign_completion_to_a_cell(tmp_path):
    paths, rows, params, values = _case(tmp_path)
    child = _child(paths, [])
    _json(child.with_name("alignment_prepare_summary.json"), {"status_counts": {"aligned": 1}})
    combined, qc, obj, report = _handoff(paths, rows, params, [child])
    assert combined["cells"] == []
    assert combined["status_counts"] == {}
    assert qc.alignment_status.tolist() == ["missing"]
    assert report["candidate_prevalence_denominator"] is None
    np.testing.assert_array_equal(obj.X.toarray(), values)
