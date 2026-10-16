#!/usr/bin/env python3
"""A02 gates and records only; no detector/quantifier implementation."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shlex
import shutil

from run_stage import digest, dump, stamp

CELL = "IMR90_A_100"


def freeze(root):
    sources = {p.name: digest(p) for p in Path(__file__).parent.glob("*.py")}
    identity = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    target = root / "helpers" / ("execution-" + identity[:16])
    target.mkdir()
    for name in sources:
        shutil.copyfile(Path(__file__).with_name(name), target / name)
    previous = {}
    for tool, real in [
        ("bwa", root / "env-bwa/bin/bwa"),
        ("java", root / "env/bin/java"),
    ]:
        wrapper = root / "audit-bin" / tool
        previous[tool] = {"sha256": digest(wrapper), "text": wrapper.read_text()}
        argv = [
            "/usr/bin/python3",
            str(target / "tool_audit.py"),
            "--root",
            str(root),
            "--tool",
            tool,
            "--real",
            str(real),
            "--",
        ]
        wrapper.write_text("#!/bin/sh\nexec " + shlex.join(argv) + ' "$@"\n')
        wrapper.chmod(0o755)
    # Keep the exact acquisition helper used before the header-audit clarification.
    legacy = Path("/tmp/circyto-a02-requests/prepare_inputs.acquisition-v1.py")
    shutil.copyfile(legacy, root / "helpers" / legacy.name)
    dump(
        root / "records/execution_helpers.json",
        {
            "recorded_utc": stamp(),
            "directory": str(target),
            "source_sha256": sources,
            "previous_wrapper_bytes": previous,
            "wrappers": {
                tool: {
                    "path": str(root / "audit-bin" / tool),
                    "sha256": digest(root / "audit-bin" / tool),
                }
                for tool in ["bwa", "java"]
            },
            "acquisition_helper_v1_sha256": digest(root / "helpers" / legacy.name),
            "acquisition_v1_snapshot_timing": "captured while acquisition ran, before editing its source; not a launch-time hash",
        },
    )
    print(json.dumps({"frozen_helpers": str(target), "sha256": identity}))


def index_record(root):
    fa = root / "ref/GRCh38.primary_assembly.genome.fa"
    records = [
        json.loads(p.read_text())
        for p in (root / "records/tool_commands").glob("*.json")
    ]
    index = [r for r in records if r["kind"] == "index"]
    if (
        len(index) != 1
        or index[0].get("exit_code") != 0
        or not index[0].get("inputs_unchanged")
    ):
        raise ValueError("Exactly one successful, bound index command required")
    expected = [str((root / "env-bwa/bin/bwa").resolve()), "index", str(fa)]
    if index[0]["argv"] != expected:
        raise ValueError("Index command differs from registration")
    outputs = {}
    for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
        path = Path(str(fa) + suffix)
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError("Missing BWA index: " + suffix)
        outputs[str(path.relative_to(root))] = {
            "bytes": path.stat().st_size,
            "sha256": digest(path),
        }
    dump(
        root / "records/index_manifest.json",
        {
            "argv": index[0]["argv"],
            "tool_sha256": index[0]["executable_sha256"],
            "input_sha256": index[0]["input_sha256"],
            "outputs": outputs,
            "new_validation_only": True,
        },
    )
    progress = json.loads((root / "records/pilot_progress.json").read_text())
    progress.update(stage="index_validated_first_cell_ready", updated_utc=stamp())
    dump(root / "records/pilot_progress.json", progress)
    print(json.dumps({"index_files_hashed": len(outputs)}))


def mark_workflow_start(root):
    for record in [
        "index_manifest.json",
        "fastq_validation_header_v2.json",
        "environment_manifest.json",
        "rna_extraction.json",
    ]:
        if not (root / "records" / record).is_file():
            raise ValueError("Required preflight record missing: " + record)
    if (root / "workflow_firstcell").exists():
        raise FileExistsError("First-cell workflow directory must be new")
    progress = json.loads((root / "records/pilot_progress.json").read_text())
    progress.update(stage="first_cell_workflow_requested", updated_utc=stamp())
    progress["cells"][0].update(
        execution_status="workflow_requested", new_detector_status="unknown"
    )
    dump(root / "records/pilot_progress.json", progress)


def alignment_output_binding(root, alignment, detector, *, workflow_name="workflow_firstcell",
                             cell_id=CELL, group_id="SRX26321174", reference=None,
                             cohort_ids=None):
    """Bind a pipe-produced SAM using contemporaneous workflow/detector records.

    Preserve the wrapper record: it did not hash a pipe's output stream. The
    reviewed workflow copies that stream to .partial and renames only on exit 0.
    No count or file-presence-only inference is accepted here.
    """
    if alignment.get("exit_code") != 0:
        raise ValueError("Alignment completion cannot be bound after nonzero exit")
    if alignment.get("completion_bound"):
        path = Path(alignment["output_path"])
        expected = alignment["output_sha256"]
        evidence = {"method": "wrapper file-output hash"}
    else:
        if alignment.get("output_path") is not None or not alignment.get(
            "inputs_unchanged"
        ):
            raise ValueError("Alignment completion cannot be bound")
        manifest = root / workflow_name / "align/bwa_mem/alignment_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        expected_ids = cohort_ids if cohort_ids is not None else [cell_id]
        if len(rows) != len(expected_ids) or {r["cell_id"] for r in rows} != set(expected_ids):
            raise ValueError("Exactly the declared fresh alignment cohort is required")
        selected = [r for r in rows if r["cell_id"] == cell_id]
        if len(selected) != 1:
            raise ValueError("Exactly one fresh alignment manifest row required per cell")
        row = selected[0]
        expected_fields = {
            "cell_id": cell_id,
            "aligner": "bwa-mem",
            "read_layout": "single-end",
            "reused_alignment": "false",
            "group_id": group_id,
            "reference": str(reference or root / "ref/GRCh38.primary_assembly.genome.fa"),
        }
        if any(row.get(k) != v for k, v in expected_fields.items()):
            raise ValueError("Fresh alignment manifest identity/settings disagree")
        path = Path(row["sam"])
        if (
            not path.is_relative_to(root / workflow_name / "align/bwa_mem/cache")
            or path.name != "alignment.sam"
            or path.with_suffix(".sam.partial").exists()
        ):
            raise ValueError("SAM is not the finalized fresh workflow alignment")
        detector_input = detector["argv"][detector["argv"].index("-I") + 1]
        if str(path) != detector_input:
            raise ValueError("Manifest SAM does not match the emitted detector input")
        expected = detector["input_sha256"][str(path)]
        log = path.parent / (cell_id + ".align.log")
        if not log.is_file() or "[main] Real time:" not in log.read_text():
            raise ValueError("BWA completion log is missing")
        evidence = {
            "method": "fresh workflow manifest + BWA exit 0 + detector pre-input hash",
            "wrapper_direct_output_hash_available": False,
            "manifest": str(manifest.relative_to(root)),
            "manifest_sha256": digest(manifest),
            "alignment_log": str(log.relative_to(root)),
            "alignment_log_sha256": digest(log),
            "limitation": "No independent hash of BWA stdout stream; reviewed workflow byte-copy/rename contract binds the SAM path.",
        }
    if digest(path) != expected:
        raise ValueError("Alignment-to-detector hash link broken")
    evidence.update(output_path=str(path), output_sha256=expected)
    return evidence


def bind_native(root):
    paths = list((root / "records/tool_commands").glob("*.json"))
    records = [(p, json.loads(p.read_text())) for p in paths]
    align = [(p, r) for p, r in records if r["kind"] == "alignment"]
    detect = [(p, r) for p, r in records if r["kind"] == "detector"]
    if len(align) != 1 or len(detect) != 1:
        raise ValueError("Exactly one alignment and one detector command required")
    ap, a = align[0]
    dp, d = detect[0]
    if not d.get("completion_bound") or d.get("exit_code") != 0:
        raise ValueError(
            "Cell completion cannot be bound to successful output; retain unknown/failed evidence"
        )
    if digest(d["output_path"]) != d["output_sha256"]:
        raise ValueError("Bound output changed")
    aligned = alignment_output_binding(root, a, d)
    env = json.loads((root / "records/environment_manifest.json").read_text())
    fa = str(root / "ref/GRCh38.primary_assembly.genome.fa")
    fq = str(root / "inputs/fastq/SRR30918126.fastq.gz")
    gtf = str(root / "ref/gencode.v45.primary_assembly.annotation.gtf")
    raw = root / "workflow_firstcell/ciri3" / (CELL + ".ciri3_run/ciri3_raw.tsv")
    expected_bwa = [
        env["external_executables"]["bwa"]["path"],
        "mem",
        "-t",
        "8",
        "-k",
        "15",
        "-T",
        "15",
        fa,
        fq,
    ]
    if a["argv"] != expected_bwa:
        raise ValueError("Emitted BWA argv differs from registration")
    expected_java = [
        env["external_executables"]["java"]["path"],
        "-jar",
        str(root / "external/CIRI3_Java_18.0.1.jar"),
        "-I",
        aligned["output_path"],
        "-O",
        str(raw),
        "-F",
        fa,
        "-T",
        "8",
        "-G",
        str(raw.parent / "ciri3_internal.log"),
        "-A",
        gtf,
        "-S",
        "0",
        "-Ma",
        "0",
    ]
    if d["argv"] != expected_java:
        raise ValueError("Emitted CIRI3 argv differs from registration")
    if d["input_sha256"][aligned["output_path"]] != aligned["output_sha256"]:
        raise ValueError("Alignment-to-detector hash link broken")
    if (
        a["input_sha256"][fq]
        != json.loads((root / "records/fastq_validation_header_v2.json").read_text())[
            "sha256"
        ]
    ):
        raise ValueError("Aligned FASTQ differs from validated complete library")
    for tool, record in [("bwa", a), ("java", d)]:
        if record["executable_sha256"] != env["external_executables"][tool]["sha256"]:
            raise ValueError("Tool executable changed")
    if (
        d["input_sha256"][str(root / "external/CIRI3_Java_18.0.1.jar")]
        != env["external_jar_sha256"]
    ):
        raise ValueError("Executed JAR changed")
    for input_path in [fa, gtf]:
        if d["input_sha256"][input_path] != digest(input_path):
            raise ValueError("Reference changed")
    with (root / "inputs/native_firstcell.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["cell_id", "sample_id", "path", "detector_version", "result_scope"]
        )
        writer.writerow([CELL, "GSM8558852", str(raw), "3.0.1", "single-sample"])
    with (root / "inputs/cells_firstcell.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cell_id", "detector_status"])
        writer.writerow([CELL, "completed"])
    binding = {
        "cell_id": CELL,
        "gsm": "GSM8558852",
        "srr": "SRR30918126",
        "new_detector_status": "completed",
        "original_detector_status": "unknown",
        "alignment_record": str(ap.relative_to(root)),
        "detector_record": str(dp.relative_to(root)),
        "native_sha256": digest(raw),
        "alignment_output_binding": aligned,
        "basis": "executed argv + unchanged input/tool hashes + exit 0 + output hash; independent of counts",
    }
    dump(root / "records/native_completion_binding.json", binding)
    progress = json.loads((root / "records/pilot_progress.json").read_text())
    progress.update(
        stage="new_native_completed_acceptance_pending", updated_utc=stamp()
    )
    progress["cells"][0].update(
        execution_status="workflow_finished",
        new_detector_status="completed",
        native_sha256=digest(raw),
    )
    dump(root / "records/pilot_progress.json", progress)
    print(json.dumps(binding))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "action",
        choices=["freeze", "index-record", "mark-workflow-start", "bind-native"],
    )
    ap.add_argument("--root", type=Path, required=True)
    args = ap.parse_args()
    {
        "freeze": freeze,
        "index-record": index_record,
        "mark-workflow-start": mark_workflow_start,
        "bind-native": bind_native,
    }[args.action](args.root.resolve())


if __name__ == "__main__":
    main()
