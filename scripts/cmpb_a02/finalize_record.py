#!/usr/bin/env python3
"""Produce small A02 audit records after execution has stopped; no analysis runs."""

import argparse
import csv
import datetime as dt
import json
from pathlib import Path
import shutil
import time

from run_stage import digest, dump, footprint, stamp

COHORT = [
    ("IMR90_A_100", "GSM8558852", "SRR30918126"),
    ("IMR90_A_101", "GSM8558853", "SRR30918117"),
    ("IMR90_A_102", "GSM8558854", "SRR30918116"),
    ("IMR90_A_103", "GSM8558855", "SRR30918115"),
]


def load(path):
    return json.loads(path.read_text()) if path.is_file() else {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--outcome", choices=["passed", "failed", "blocked"], required=True)
    ap.add_argument("--reason", required=True)
    args = ap.parse_args()
    root = args.root.resolve()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    record_root = root / "records"
    stages = [load(p) for p in record_root.glob("*.json")]
    stages = sorted(
        [
            r
            for r in stages
            if isinstance(r, dict) and r.get("schema") == "circyto.cmpb.a02.stage.v1"
        ],
        key=lambda r: r["started_utc"],
    )
    if any(r["status"] == "launching" for r in stages):
        raise RuntimeError("A managed stage remains unfinished; do not finalize")
    by_name = {r["stage"]: r for r in stages}
    if args.outcome == "passed":
        required = [
            "workflow-firstcell",
            "bind-firstcell-native",
            "report-new-workflow",
            "import-firstcell-native-rna",
            "independent-value-report-checks-v2"
            if "independent-value-report-checks-v2" in by_name
            else "independent-value-report-checks",
            "open-offline-reports",
            "workflow-status-acceptance",
        ]
        if any(by_name.get(name, {}).get("status") != "succeeded" for name in required):
            raise ValueError("Cannot declare PASS without every required actual stage")
        if not load(record_root / "acceptance/browser_acceptance.json").get("complete"):
            raise ValueError("Both actual offline report openings must complete")
        for name in ["native_value_checks.tsv", "report_value_checks.tsv"]:
            with (record_root / "acceptance" / name).open(newline="") as handle:
                checks = list(csv.DictReader(handle, delimiter="\t"))
            if not checks or any(row["status"] != "PASS" for row in checks):
                raise ValueError(
                    "Cannot declare PASS with missing/failed independent checks"
                )
    commands = [load(p) for p in (record_root / "tool_commands").glob("*.json")]
    align = [r for r in commands if r.get("kind") == "alignment"]
    detect = [r for r in commands if r.get("kind") == "detector"]
    binding = load(record_root / "native_completion_binding.json")
    summary = load(record_root / "acceptance/native_summary.json")
    fastq = load(record_root / "fastq_validation_header_v2.json") or load(
        record_root / "fastq_validation.json"
    )
    env = load(record_root / "environment_manifest.json")
    rna = load(record_root / "rna_extraction.json")
    status = (
        "completed"
        if binding.get("new_detector_status") == "completed"
        else (
            "failed"
            if detect and detect[-1].get("exit_code", 0) != 0
            else "unknown"
            if detect
            else "unprocessed"
        )
    )
    raw = root / "workflow_firstcell/ciri3/IMR90_A_100.ciri3_run/ciri3_raw.tsv"
    if not summary and raw.is_file():
        from independent_checks import native

        try:
            observed, _ = native(raw)
            summary = {
                "native_records": len(observed),
                "support": sum(int(r["#junction_reads"]) for r in observed),
                "native_sha256": digest(raw),
                "count_scope": "retained output only; completion is separately "
                + status,
            }
            dump(record_root / "acceptance/observed_native_output.json", summary)
        except (ValueError, KeyError) as error:
            dump(
                record_root / "acceptance/unparsed_native_output.json",
                {"sha256": digest(raw), "error": str(error)},
            )
    progress = load(record_root / "pilot_progress.json")
    progress.update(
        stage="first_cell_acceptance_" + args.outcome,
        updated_utc=stamp(),
        four_cell_pilot_complete=False,
        first_cell_outcome=args.outcome,
        reason=args.reason,
    )
    progress["cells"][0].update(
        execution_status="finished"
        if by_name.get("workflow-firstcell", {}).get("status") == "succeeded"
        else args.outcome,
        new_detector_status=status,
        acceptance=args.outcome,
    )
    for row in progress["cells"][1:]:
        if row["execution_status"] != "not_started":
            raise ValueError("An unapproved additional cell started")
    dump(record_root / "pilot_progress.json", progress)
    dump(out / "pilot_progress.json", progress)
    rows = []
    columns = [
        "cell_id",
        "gsm",
        "srr",
        "selection",
        "route",
        "execution_status",
        "original_detector_status",
        "new_detector_status",
        "fastq_sha256",
        "full_library_reads",
        "layout",
        "fasta_sha256",
        "gtf_sha256",
        "bwa_exit_code",
        "detector_exit_code",
        "native_sha256",
        "native_candidates",
        "native_support",
        "same_cell_RNA_source_sha256",
        "same_cell_RNA_slice_sha256",
        "historical_RNA_reference",
        "acceptance",
    ]
    for i, (cell, gsm, srr) in enumerate(COHORT):
        row = {key: "" for key in columns}
        row.update(
            cell_id=cell,
            gsm=gsm,
            srr=srr,
            selection="fixed ascending-GSM cohort; no outcome-based replacement",
            route="NEW Route B",
            execution_status=progress["cells"][i]["execution_status"],
            original_detector_status="unknown",
            new_detector_status=status if i == 0 else "unprocessed",
            acceptance=args.outcome if i == 0 else "not_started",
        )
        if i == 0:
            row.update(
                fastq_sha256=fastq.get("sha256", ""),
                full_library_reads=fastq.get("reads", ""),
                layout=fastq.get("layout", ""),
                fasta_sha256=load(
                    record_root / "reference-GRCh38.primary_assembly.genome.fa.gz.json"
                )
                .get("uncompressed", {})
                .get("sha256", ""),
                gtf_sha256=load(
                    record_root
                    / "reference-gencode.v45.primary_assembly.annotation.gtf.gz.json"
                )
                .get("uncompressed", {})
                .get("sha256", ""),
                bwa_exit_code=align[-1].get("exit_code", "unknown")
                if align
                else "not_started",
                detector_exit_code=detect[-1].get("exit_code", "unknown")
                if detect
                else "not_started",
                native_sha256=binding.get(
                    "native_sha256", summary.get("native_sha256", "")
                ),
                native_candidates=summary.get("native_records", ""),
                native_support=summary.get("support", ""),
                same_cell_RNA_source_sha256=rna.get("source_sha256", ""),
                same_cell_RNA_slice_sha256=rna.get("output_sha256", ""),
                historical_RNA_reference="unverified; original declarations preserved",
            )
        rows.append(row)
    with (out / "pilot_cell_trace.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    for name in ["native_value_checks.tsv", "report_value_checks.tsv"]:
        source = record_root / "acceptance" / name
        if source.is_file():
            shutil.copyfile(source, out / name)
            if args.outcome == "passed":
                with source.open(newline="") as h:
                    if any(
                        row["status"] != "PASS"
                        for row in csv.DictReader(h, delimiter="\t")
                    ):
                        raise ValueError("A value check did not pass")
        else:
            (out / name).write_text(
                "check\tscope\texpected\tobserved\tstatus\nacceptance\tfirst cell\texecution and independent checks\t"
                + args.reason.replace("\t", " ").replace("\n", " ")
                + "\tNOT_RUN\n"
            )
    guards = [load(record_root / (r["stage"] + ".guard.json")) for r in stages]
    budget = load(record_root / "budget.json")
    last_stage_finished = max(
        dt.datetime.fromisoformat(r["finished_utc"]).timestamp() for r in stages
    )
    usage = {
        "recorded_utc": stamp(),
        "evaluation_root": str(root),
        "caps": budget,
        "active_stage_elapsed_seconds": sum(r["elapsed_seconds"] for r in stages),
        "shared_preparation_stage_seconds": sum(
            r["elapsed_seconds"]
            for r in stages
            if r["category"] == "shared_preparation"
        ),
        "first_cell_stage_seconds": sum(
            r["elapsed_seconds"] for r in stages if r["category"] == "firstcell"
        ),
        "control_probe_seconds": sum(
            r["elapsed_seconds"] for r in stages if r["category"] == "probe"
        ),
        "conservative_elapsed_since_budget_start_seconds": time.time()
        - budget["started_epoch"],
        "execution_envelope_through_last_stage_seconds": last_stage_finished
        - budget["started_epoch"],
        "last_managed_stage_finished_utc": dt.datetime.fromtimestamp(
            last_stage_finished, dt.timezone.utc
        ).isoformat(),
        "handoff_gap_after_execution_seconds": time.time() - last_stage_finished,
        "sampled_cgroup_memory_peak_bytes": max(
            g.get("memory_current_peak_bytes", 0) for g in guards
        ),
        "sampled_process_tree_RSS_peak_bytes_where_recorded": max(
            g.get("process_tree_rss_peak_bytes", 0) for g in guards
        ),
        "evaluation_footprint_bytes_now": footprint(root),
        "sampled_evaluation_footprint_peak_bytes": max(
            g.get("disk_peak_bytes", 0) for g in guards
        ),
        "measurement_limits": [
            "RSS and cgroup memory sampled every 0.25 s; disk each ~1 s including traversal time.",
            "Initial acquisition/preparation recorded cgroup memory; summed process RSS added before index and first-cell work.",
            "CPU quota controller unavailable; eight-CPU affinity inheritance and all-thread violation stop tested.",
            "Disk auto-stop at 79 GiB provides headroom below 80 GiB; not a filesystem quota.",
            "Storage uses max(logical,allocated) per file, counting hard-linked cache copies conservatively.",
            "Final audit writing resumed after a server restart; calendar elapsed includes the non-execution handoff gap, reported separately. No budget reset or resumed scientific stage.",
        ],
        "stages": [
            {
                "stage": r["stage"],
                "category": r["category"],
                "status": r["status"],
                "elapsed_seconds": r["elapsed_seconds"],
                "unit_result": r["systemd_properties"],
                "guard": g,
            }
            for r, g in zip(stages, guards)
        ],
    }
    dump(record_root / "resource_usage.json", usage)
    dump(out / "resource_usage.json", usage)
    dump(
        out / "commands.json",
        {
            "base_head": "17b733c8f21bb128a035e7b100e704099e2bdb0e",
            "evaluated_source": "6c280426ebf10abcbcaeab72f87c9d34060d6ee8",
            "stages": stages,
            "actual_tool_commands": commands,
        },
    )
    dump(out / "environment_manifest.json", env)
    inputs = {
        "artifacts": load(record_root / "input_artifacts.json"),
        "ENA": load(record_root / "ena_comparison.json"),
        "FASTQ": {
            k: v for k, v in fastq.items() if k not in ["first_headers", "last_headers"]
        },
        "reference_provider": load(record_root / "reference_provider_evidence.json"),
        "FASTA": load(
            record_root / "reference-GRCh38.primary_assembly.genome.fa.gz.json"
        ),
        "GTF": load(
            record_root
            / "reference-gencode.v45.primary_assembly.annotation.gtf.gz.json"
        ),
        "compatibility": load(record_root / "reference_compatibility.json"),
        "index": load(record_root / "index_manifest.json"),
        "RNA": rna,
        "native_completion_binding": binding,
        "registration_preservation": load(
            Path("/tmp/circyto-a02-requests/preservation_check.json")
        ),
        "data_provenance_class": "NEW executed technical validation; original processed RNA preserved; historical detection metadata remains missing",
    }
    inputs["execution_and_selection_manifests"] = {}
    for name in [
        "registered_cohort.tsv",
        "manifest_firstcell.tsv",
        "native_firstcell.tsv",
        "cells_firstcell.tsv",
    ]:
        path = root / "inputs" / name
        if path.is_file():
            with path.open(newline="") as handle:
                content = list(csv.DictReader(handle, delimiter="\t"))
            inputs["execution_and_selection_manifests"][name] = {
                "sha256": digest(path),
                "rows": content,
            }
    dump(out / "input_manifest.json", inputs)
    artifacts = {}
    for prefix in [
        "inputs/rna",
        "workflow_firstcell",
        "import_firstcell",
        "records/acceptance",
        "records/acceptance_attempt1",
    ]:
        directory = root / prefix
        if directory.exists():
            for path in sorted(directory.rglob("*")):
                if path.is_file() and path.suffix.lower() not in [
                    ".sam",
                    ".bam",
                    ".bai",
                ]:
                    artifacts[str(path.relative_to(root))] = {
                        "bytes": path.stat().st_size,
                        "sha256": digest(path),
                    }
    dump(
        out / "output_artifact_manifest.json",
        {
            "files": artifacts,
            "large_files_remain_outside_Git": True,
            "alignment_hashes": "commands.json actual_tool_commands; input_manifest.json native_completion_binding.alignment_output_binding",
        },
    )
    dump(
        out / "control_tests.json", load(record_root / "final-probe_control_tests.json")
    )
    dump(
        out / "executed_helper_manifest.json",
        load(record_root / "execution_helpers.json"),
    )
    print(
        json.dumps(
            {
                "outcome": args.outcome,
                "small_record_directory": str(out),
                "copied_raw_sequences": False,
            }
        )
    )


if __name__ == "__main__":
    main()
