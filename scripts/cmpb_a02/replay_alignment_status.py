#!/usr/bin/env python3
"""A100-only postprocessing replay; never executes an aligner or detector.

Original outputs remain read-only inputs. The independent checker is reused
without changing its assertions; regenerated QC/HDF5/report files are copies.
"""
import argparse
import copy
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True

from run_stage import digest, dump, stamp

CELL = "IMR90_A_100"


class ReplayGuard:
    def __init__(self, original, allowed):
        self.original, self.allowed = original.resolve(), allowed.resolve()
        self.blocked = []

    def __call__(self, event, args):
        if event == "subprocess.Popen" or event in {"os.system", "os.exec", "os.posix_spawn", "os.fork", "os.forkpty"}:
            self.blocked.append({"event": event, "argv": repr(args)})
            raise RuntimeError("External execution prohibited during status replay")
        candidates = []
        if event == "open":
            path, mode, flags = args
            if (isinstance(mode, str) and any(c in mode for c in "wax+")) or (
                isinstance(flags, int) and flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND)
            ):
                candidates.append(path)
        elif event in {"os.remove", "os.rmdir", "os.mkdir", "os.chmod", "os.utime", "os.truncate"}:
            candidates.append(args[0])
        elif event in {"os.rename", "os.link", "os.symlink"}:
            candidates.extend(args[:2])
        for value in candidates:
            if not isinstance(value, (str, bytes, os.PathLike)):
                continue
            path = Path(os.fsdecode(value)).resolve()
            if path.is_relative_to(self.original) and not path.is_relative_to(self.allowed):
                raise RuntimeError("Original A02 evidence is read-only during replay: " + str(path))


def source_snapshot(original):
    result, seen = {}, {}
    # Bounded named input/evidence trees, not a new archival search.
    for directory in ["inputs", "ref", "external", "workflow_firstcell", "import_firstcell", "records"]:
        for path in sorted((original / directory).rglob("*")):
            if not path.is_file():
                continue
            stat = path.stat()
            key = (stat.st_dev, stat.st_ino)
            if key not in seen:
                seen[key] = digest(path)
            result[str(path.relative_to(original))] = {
                "sha256": seen[key], "bytes": stat.st_size,
                "symlink": os.readlink(path) if path.is_symlink() else None,
            }
    return result


def copy_file(source, target, copies):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise FileExistsError(target)
    shutil.copyfile(source, target)
    assert source.stat().st_ino != target.stat().st_ino
    assert digest(source) == digest(target)
    copies.append({"source": str(source), "copy": str(target), "sha256": digest(source)})


def bound_evidence(original):
    from execution_records import alignment_output_binding

    binding = json.loads((original / "records/native_completion_binding.json").read_text())
    if binding["cell_id"] != CELL:
        raise ValueError("Completion binding belongs to a different cell")
    a = json.loads((original / binding["alignment_record"]).read_text())
    d = json.loads((original / binding["detector_record"]).read_text())
    if a["exit_code"] != 0 or d["exit_code"] != 0 or not d["completion_bound"]:
        raise ValueError("No bound successful original execution")
    for record in [a, d]:
        if not record.get("inputs_unchanged"):
            raise ValueError("Original execution inputs were not unchanged")
        for name, sha in record["input_sha256"].items():
            if digest(name) != sha:
                raise ValueError("Recorded input hash mismatch: " + name)
    if digest(d["output_path"]) != binding["native_sha256"] or digest(d["output_path"]) != d["output_sha256"]:
        raise ValueError("Native output hash binding broken")
    aligned = alignment_output_binding(original, a, d)
    if aligned != binding["alignment_output_binding"]:
        raise ValueError("Alignment evidence no longer matches the original binding")
    failures = json.loads((original / "records/acceptance/status_acceptance.json").read_text())
    for path, sha in failures["source_sha256"].items():
        if digest(original / path) != sha:
            raise ValueError("Original failed evidence changed: " + path)
    return binding, failures


def replay(original, root):
    import anndata as ad
    import circyto
    import pandas as pd
    from circyto.pipeline.align_manifest import read_source_manifest
    from circyto.pipeline.run_ciri3 import _write_alignment_summary
    from circyto.pipeline.workflow_full_length_circrna import (
        FullLengthCircRNAWorkflowParams, _build_execution_summary,
        _read_manifest_rows, _resolve_row_protocol, _resolve_row_strandedness,
        _workflow_paths,
    )
    from circyto.pipeline.qc_report import generate_report
    from independent_checks import acceptance, equal_adata, exact, Tables, write_checks

    started = time.monotonic()
    work = root / "workflow_firstcell"
    if work.exists() or (root / "import_firstcell").exists():
        raise FileExistsError("Replay output already exists; no overwrite/resume")
    if "site-packages" not in circyto.__file__ or not Path(circyto.__file__).is_relative_to(root / "installed-env"):
        raise ValueError("Replay requires the separate repaired installed wheel")
    build_record = json.loads((root / "records/build_acceptance.json").read_text())
    wheel_members = next(iter(build_record["wheel_python_members"].values()))
    for name, sha in wheel_members.items():
        if digest(Path(circyto.__file__).parent.parent / name) != sha:
            raise ValueError("Installed Python member differs from repaired wheel: " + name)
    guard = ReplayGuard(original, root)
    sys.addaudithook(guard)
    # Exercise real spawn interception before doing any postprocessing.
    for argv in [["bwa", "mem"], ["bwa", "index"], ["STAR", "--runMode", "genomeGenerate"], ["java", "-jar", "CIRI3.jar"]]:
        try:
            subprocess.run(argv, check=True)
        except RuntimeError as exc:
            if "External execution prohibited" not in str(exc):
                raise
        else:
            raise AssertionError("Replay execution guard did not block a tool")
    probes = list(guard.blocked)
    guard.blocked.clear()
    dump(root / "records/replay_state.json", {"status": "hashing_original_sources", "started_utc": stamp()})
    before = source_snapshot(original)
    dump(root / "records/source_preservation_before.json", before)
    binding, original_failures = bound_evidence(original)
    old = original / "workflow_firstcell"
    old_summary = json.loads((old / "workflow_summary.json").read_text())
    paths = _workflow_paths(work)
    for name, path in paths.items():
        if name != "workflow_summary":
            path.mkdir(parents=True, exist_ok=False) if name == "root" else path.mkdir(parents=True, exist_ok=True)
    copies = []
    for name in ["matrix", "ciri3", "rna", "manifests"]:
        for source in sorted((old / name).rglob("*")):
            if source.is_file():
                copy_file(source, work / source.relative_to(old), copies)
    for name in ["alignment_manifest.tsv"]:
        copy_file(old / "align" / name, paths["align"] / name, copies)
    for name in ["rna_qc.tsv", "rna_gene_qc.tsv"]:
        if (old / "qc" / name).exists():
            copy_file(old / "qc" / name, paths["qc"] / name, copies)
    # Reuse the independently validated native-import bundle and actual same-cell
    # RNA as byte-preserved copies; no reimport or synthetic replacement is needed.
    for name in ["import_firstcell", "inputs/rna", "inputs/original"]:
        for source in sorted((original / name).rglob("*")):
            if source.is_file():
                copy_file(source, root / source.relative_to(original), copies)
    dump(root / "records/replay_copies.json", copies)
    manifest = Path(old_summary["paths"]["manifest"])
    rows = read_source_manifest(manifest, validate_files=True)
    if [r.cell_id for r in rows] != [CELL]:
        raise ValueError("Replay is restricted to the original one-cell manifest")
    _, raw_rows = _read_manifest_rows(manifest)
    child_manifest = old / "align/bwa_mem/alignment_manifest.tsv"
    _write_alignment_summary(
        paths["align"] / "alignment_prepare_summary.json", source_rows=rows,
        manifest_paths=[child_manifest, paths["align"] / "alignment_manifest.tsv"], dry_run=False,
    )
    params = FullLengthCircRNAWorkflowParams(
        manifest=manifest, outdir=work, protocol=old_summary["protocol"],
        genome_fasta=Path(old_summary["genome_fasta"]), gtf=Path(old_summary["gtf"]),
        gene_expression_method=old_summary["planned_multimodal"]["gene_expression_method"],
        report=False,
    )
    summary = _build_execution_summary(
        params=params, paths=paths, source_rows=rows, raw_rows=raw_rows,
        resolved_protocols={r.cell_id: _resolve_row_protocol(r, params.protocol) for r in rows},
        resolved_strandedness={r.cell_id: _resolve_row_strandedness(r, None) for r in rows},
        skip_demux_effective=old_summary["skip_demux_effective"],
        elapsed_seconds=old_summary["elapsed_seconds"], warnings=old_summary["warnings"],
        rna_import=old_summary["rna_import"], started_at=old_summary["started_at"], completed_at=old_summary["completed_at"],
    )
    summary["postprocessing_replay"] = {
        "source_workflow": str(old), "source_workflow_uuid": old_summary["workflow_uuid"],
        "source_summary_sha256": digest(old / "workflow_summary.json"),
        "repair_commit": build_record["source_commit"], "regenerated_at": stamp(),
        "scope": "status/QC/AnnData/report only; original execution timestamps and timing retained",
        "alignment_detection_rerun": False,
    }
    dump(paths["workflow_summary"], summary)
    generate_report(work)
    dump(root / "records/replay_state.json", {"status": "independent_validation", "updated_utc": stamp()})
    old_obj = ad.read_h5ad(old / "anndata/circ_counts.h5ad")
    new_obj = ad.read_h5ad(work / "anndata/circ_counts.h5ad")
    changes = {}
    for column in ["alignment_status", "alignment_seconds"]:
        changes["obs." + column] = {"before": old_obj.obs[column].astype(str).tolist(), "after": new_obj.obs[column].astype(str).tolist()}
    a, b = old_obj.copy(), new_obj.copy()
    a.obs = a.obs.drop(columns=["alignment_status", "alignment_seconds"])
    b.obs = b.obs.drop(columns=["alignment_status", "alignment_seconds"])
    equal_adata(a, b, "unchanged workflow scientific object")
    old_qc = pd.read_csv(old / "qc/cell_qc.tsv", sep="\t")
    new_qc = pd.read_csv(work / "qc/cell_qc.tsv", sep="\t")
    exact(old_qc.drop(columns=["alignment_status", "alignment_seconds"]),
          new_qc.drop(columns=["alignment_status", "alignment_seconds"]), "QC non-status fields")
    assert digest(old / "qc/circ_qc.tsv") == digest(work / "qc/circ_qc.tsv")
    acceptance(root)  # unchanged independent CSV/CSR/RNA/ordinary object checks
    report = json.loads((work / "qc/metrics.json").read_text())
    combined = json.loads((work / "align/alignment_prepare_summary.json").read_text())
    child = json.loads(child_manifest.with_name("alignment_prepare_summary.json").read_text())
    assert combined["cells"] == child["cells"]
    assert combined["status_counts"] == child["status_counts"]
    recorded = child["cells"][0]["status"]
    values = [report["per_cell"][0]["alignment_status"], str(new_obj.obs.loc[CELL, "alignment_status"]),
              report["workflow"]["status"], str(report["candidate_prevalence_denominator"])]
    checks = []
    for previous, observed in zip(original_failures["checks"], values):
        expected = previous["expected"]
        checks.append({**previous, "observed": observed, "status": "PASS" if expected == observed else "FAIL"})
    assert recorded == str(new_qc.loc[0, "alignment_status"]) == values[0] == values[1]
    parsed = Tables()
    parsed.feed((work / "qc/report.html").read_text())
    assert [CELL, "142", "474", "success", recorded] in next(t for t in parsed.tables if t[0][0] == "Cell")
    assert "completed" in "".join(parsed.text)
    assert new_obj.shape == (1, 142) and int(new_obj.X.sum()) == 474
    rt = root / "records/acceptance/workflow_roundtrip.h5ad"
    new_obj.write_h5ad(rt)
    equal_adata(new_obj, ad.read_h5ad(rt), "regenerated workflow ordinary roundtrip")
    write_checks(root / "records/status_checks.tsv", checks)
    dump(root / "records/scientific_preservation.json", {
        "workflow_values_axes_annotations_layers_metadata_equal_except": changes,
        "workflow_X_dtype": str(new_obj.X.dtype), "candidates": new_obj.n_vars,
        "support": int(new_obj.X.sum()), "native_source_sha256": binding["native_sha256"],
        "before_workflow_h5ad_sha256": digest(old / "anndata/circ_counts.h5ad"),
        "after_workflow_h5ad_sha256": digest(work / "anndata/circ_counts.h5ad"),
        "circ_qc_bytes_equal": True, "workflow_roundtrip": "PASS",
    })
    after = source_snapshot(original)
    dump(root / "records/source_preservation_after.json", after)
    if before != after:
        raise AssertionError("Original source/evidence files changed")
    if guard.blocked:
        raise AssertionError("Postprocessing attempted forbidden external execution")
    if any(c["status"] != "PASS" for c in checks):
        raise AssertionError("Original semantic assertions still fail")
    modules = {name: __import__(name).__file__ for name in ["circyto", "anndata", "mudata", "numpy", "pandas", "scipy", "h5py"]}
    dump(root / "records/replay_acceptance.json", {
        "status": "PASS", "scope": "preserved A100 postprocessing; not an independent biological run",
        "repair_commit": build_record["source_commit"], "artifacts": build_record["artifacts"],
        "source_files_verified_unchanged": len(before), "source_preservation_manifest_sha256": digest(root / "records/source_preservation_before.json"),
        "execution_guard_probes": probes, "postprocessing_external_attempts": guard.blocked,
        "checks": checks, "elapsed_seconds": time.monotonic()-started, "finished_utc": stamp(),
        "module_origins": modules, "python": sys.version,
        "installed_python_members_equal_wheel": len(wheel_members),
        "resolved_distributions": {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()},
    })
    progress = copy.deepcopy(json.loads((original / "records/pilot_progress.json").read_text()))
    progress.update(stage="A100_status_revalidation_passed", first_cell_original_outcome="failed",
                    first_cell_repaired_pipeline_acceptance="passed", source_commit=build_record["source_commit"], updated_utc=stamp())
    progress["cells"][0]["status_revalidation"] = "passed; original failed acceptance retained"
    assert all(c["execution_status"] == "not_started" for c in progress["cells"][1:])
    dump(root / "records/pilot_progress.json", progress)
    dump(root / "records/replay_state.json", {"status": "passed", "updated_utc": stamp()})
    print(json.dumps({"replay": "PASS", "source_files_unchanged": len(before), "candidates": 142, "support": 474}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    replay(args.original.resolve(strict=True), args.root.resolve(strict=True))
