#!/usr/bin/env python3
"""Small A02 handoff records; reads original evidence and writes only new paths."""

import argparse
import csv
import json
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET

from continuation_inputs import preserved, ROWS, SOURCE, WHEEL_SHA
from continuation_records import CELLS
from run_stage import digest, dump, footprint, stamp


def read(path, default=None):
    return json.loads(path.read_text()) if path.exists() else default


def sanitize(value, root, original, repo):
    replacements = [
        (str(root), "${CONTINUATION_ROOT}"),
        (str(original), "${A02_ROOT}"),
        (str(repo), "${ANALYSIS_WORKTREE}"),
        (str(repo.parent / "circyto-cmpb-evaluation"), "${REPAIR_WORKTREE}"),
    ]

    def clean(s):
        for a, b in replacements:
            s = s.replace(a, b)
        return s

    if isinstance(value, str):
        return clean(value)
    if isinstance(value, list):
        return [sanitize(v, root, original, repo) for v in value]
    if isinstance(value, dict):
        return {clean(k): sanitize(v, root, original, repo) for k, v in value.items()}
    return value


def preservation(root, original, repo):
    before = read(root / "records/preservation_before.json")
    after = preserved(original)
    if before != after:
        delta = [k for k in set(before) | set(after) if before.get(k) != after.get(k)]
        dump(root / "records/preservation_conflicts.json", delta)
        raise ValueError("Original/repair source preservation failed")
    protected = Path("/tmp/circyto-cmpb-first-session-audit/preserved-files.json")
    # Location is supplied explicitly by the existing audit; no archive search.
    expected = read(protected)
    main = Path(
        subprocess.check_output(
            ["git", "worktree", "list", "--porcelain"], cwd=repo, text=True
        )
        .splitlines()[0]
        .split(" ", 1)[1]
    )
    actual = {
        str(p.relative_to(main))
        for p in (main / "manuscript/results_v2_recheck").rglob("*")
        if p.is_file()
    }
    actual.update(
        p for p in expected if not p.startswith("manuscript/results_v2_recheck/")
    )
    if actual != set(expected) or any(
        digest(main / p) != v["sha256"] or (main / p).stat().st_size != v["bytes"]
        for p, v in expected.items()
    ):
        raise ValueError("Protected manuscript/recheck path set or bytes changed")
    old = read(
        repo
        / "manuscript/cmpb/a02_firstcell/status_revalidation/protected_work_preservation.json"
    )
    for ref, sha in old["original_branch_tips"].items():
        if (
            subprocess.check_output(
                ["git", "rev-parse", ref], cwd=repo, text=True
            ).strip()
            != sha
        ):
            raise ValueError("Original branch tip moved: " + ref)
    if (
        subprocess.check_output(
            ["git", "rev-parse", "fix/workflow-alignment-status"], cwd=repo, text=True
        ).strip()
        != "05e5171419bffb96a7890c492b7e4ab7531ba18a"
    ):
        raise ValueError("Repair branch tip moved")
    old_paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", "05e5171", "--", "manuscript/cmpb"],
        cwd=repo,
        text=True,
    ).splitlines()
    for p in old_paths:
        if (repo / p).read_bytes() != subprocess.check_output(
            ["git", "show", "05e5171:" + p], cwd=repo
        ):
            raise ValueError("Earlier manuscript evidence changed: " + p)
    objects = old["original_H5MU_sha256"]
    for p, sha in objects.items():
        source = Path(
            p.replace("${ORIGINAL_OBJECT_ROOT}", "/tmp/circyto-import/load_work")
        )
        if digest(source) != sha:
            raise ValueError("Historical object changed")
    result = dict(
        status="PASS",
        original_and_repair_files_compared=len(before),
        path_sets_equal=True,
        before_manifest_sha256=digest(root / "records/preservation_before.json"),
        protected_files_compared=len(actual),
        protected_manifest_sha256=digest(protected),
        earlier_cmpb_files_compared=len(old_paths),
        original_branch_tips_unchanged=True,
        repair_branch_unchanged=True,
        original_H5MU_sha256=objects,
    )
    dump(root / "records/preservation_after.json", result)
    return result


def packet(root, original, repo, audit_directory="records/acceptance"):
    out = repo / "manuscript/cmpb/a02_complete"
    preservation_result = preservation(root, original, repo)
    budget = read(root / "records/budget.json")
    prior = read(
        repo / "manuscript/cmpb/a02_firstcell/status_revalidation/resource_usage.json"
    )
    summary = read(root / audit_directory / "cohort_summary.json", {})
    browser = read(root / audit_directory / "browser_acceptance.json", {})
    bindings = read(root / "records/joint_native_sources.json", {}).get("cells", [])
    if not bindings:
        old_binding = read(original / "records/native_completion_binding.json")
        bindings = [
            dict(
                cell_id=CELLS[0],
                alignment_status="aligned",
                detector_status=old_binding["new_detector_status"],
                native_sha256=old_binding["native_sha256"],
                binding_record=str(original / "records/native_completion_binding.json"),
            )
        ]
        bindings.extend(read(root / "records/native_completion_bindings.json", []))
    progress = []
    old_fastq = read(original / "records/fastq_validation_header_v2.json")
    trace = []
    for i, cell in enumerate(CELLS):
        b = next((x for x in bindings if x["cell_id"] == cell), {})
        row = ROWS[i - 1] if i else (100, "GSM8558852", "SRR30918126")
        fq = read(root / f"records/fastq-{row[2]}.json", {}) if i else old_fastq
        trace.append(
            dict(
                cell_id=cell,
                gsm=row[1],
                run=row[2],
                raw_reads=fq.get("reads", ""),
                raw_bases=fq.get("bases", ""),
                fastq_sha256=fq.get("sha256", ""),
                native_candidates=summary.get("per_cell_candidates", {}).get(cell, ""),
                detector_support=summary.get("per_cell_support", {}).get(cell, ""),
                native_sha256=b.get("native_sha256", ""),
                bound_alignment_status=b.get("alignment_status", "unknown"),
                bound_detector_status=b.get("detector_status", "unknown"),
                scientific_source="6c280426ebf10abcbcaeab72f87c9d34060d6ee8"
                if i == 0
                else SOURCE,
                derived_status_source=SOURCE if i == 0 else "normal repaired workflow",
                original_acceptance="FAIL retained; repaired derived acceptance PASS"
                if i == 0
                else "new execution",
                historical_detector_status="unknown",
                rna_profile="historical genuine same-cell simple-overlap; historical reference unverified",
                joint_rna_match="exact biological ID"
                if summary.get("passed")
                else "pending",
                workflow_report_cohort="A100 status replay"
                if i == 0
                else "remaining three",
                workflow_report_status="completed"
                if summary.get("passed")
                else "pending",
                native_record=b.get("binding_record", b.get("detector_record", "")),
            )
        )
        progress.append(
            dict(
                cell_id=cell,
                gsm=row[1],
                run=row[2],
                execution_status="preserved original"
                if i == 0
                else "completed"
                if b
                else "not established",
                detector_status=b.get("detector_status", "unknown"),
                acceptance="PASS" if summary.get("passed") else "pending",
                scientific_source=trace[-1]["scientific_source"],
                derived_status_source=trace[-1]["derived_status_source"],
                historical_detector_status="unknown",
            )
        )
    stages = []
    for p in sorted((root / "records").glob("*.json")):
        r = read(p)
        if (
            not isinstance(r, dict)
            or r.get("schema") != "circyto.cmpb.a02.stage.v1"
            or "elapsed_seconds" not in r
        ):
            continue
        g = read(root / "records" / (r["stage"] + ".guard.json"), {})
        stages.append(
            dict(
                stage=r["stage"],
                status=r["status"],
                seconds=r["elapsed_seconds"],
                started_utc=r["started_utc"],
                finished_utc=r["finished_utc"],
                exit_code=r["launcher_exit_code"],
                sampled_memory_bytes=g.get("memory_current_peak_bytes", 0),
                sampled_RSS_bytes=g.get("process_tree_rss_peak_bytes", 0),
                sampled_disk_bytes=g.get("disk_peak_bytes", 0),
                record=str(p),
                record_sha256=digest(p),
            )
        )
    stages.sort(key=lambda r: r["started_utc"])
    active = sum(r["seconds"] for r in stages)
    resource = dict(
        prior_active_seconds=prior["combined_active_stage_elapsed_seconds"],
        original=prior["original"],
        repair_active_seconds=prior["repair_active_stage_elapsed_seconds"],
        original_deadline_epoch=prior["original_budget_deadline_epoch"],
        repair_budget=prior["repair_budget"],
        continuation_budget=budget,
        continuation_active_seconds=active,
        all_A02_active_seconds=prior["combined_active_stage_elapsed_seconds"] + active,
        calendar_seconds_at_packet=time.time() - budget["started_epoch"],
        within_deadline=time.time() <= budget["deadline_epoch"],
        current_whole_A02_footprint_bytes=footprint(original),
        stages=stages,
        limits=[
            "20 GiB cgroup memory; zero swap; affinity CPUs 0-7 and all-descendant-thread watchdog",
            "Whole A02 disk stop at 79 GiB; 80 GiB cap; sampling is not a filesystem quota",
            "0.25 s memory/RSS and approximately 1 s disk sampling; no quota controller delegated",
            "Stage elapsed excludes human review pauses; calendar window includes them",
        ],
    )
    complete = bool(
        summary.get("passed")
        and browser.get("complete")
        and resource["within_deadline"]
    )
    progress_record = dict(
        updated_utc=stamp(),
        four_cell_acceptance="PASS" if complete else "INCOMPLETE",
        cells=progress,
        mixed_execution_lineage=True,
        historical_provenance_recovered=False,
    )
    dump(root / "records/pilot_progress.json", progress_record)
    sources = {}
    for subdir in [
        "inputs",
        "workflow_remaining3",
        "import_four",
        "import_four_repeat",
        "helpers/runtime",
    ]:
        for p in sorted((root / subdir).rglob("*")):
            if p.is_file() and "__pycache__" not in p.parts:
                sources[str(p.relative_to(root))] = dict(
                    bytes=p.stat().st_size, sha256=digest(p)
                )
    test_path = root / "records/focused_helpers_v3.xml"
    if not test_path.exists():
        test_path = root / "records/focused_helpers_v2.xml"
    if not test_path.exists():
        test_path = root / "records/focused_helpers.xml"
    tests = ET.parse(test_path).getroot()
    suites = list(tests.iter("testsuite"))
    test_result = {
        k: sum(int(s.get(k, 0)) for s in suites)
        for k in ["tests", "failures", "errors", "skipped"]
    }
    checks = {}
    for name in ["native_value_checks.tsv", "report_value_checks.tsv"]:
        p = root / audit_directory / name
        if p.exists():
            with p.open(newline="") as h:
                rows = list(csv.DictReader(h, delimiter="\t"))
            checks[name] = dict(
                rows=len(rows),
                failed=sum(r["status"] != "PASS" for r in rows),
                sha256=digest(p),
            )
            (out / name).write_text(sanitize(p.read_text(), root, original, repo))
    manifest = dict(
        base_head="05e5171419bffb96a7890c492b7e4ab7531ba18a",
        production_source=SOURCE,
        wheel_sha256=WHEEL_SHA,
        helper_hashes=read(root / "records/continuation_helpers.json"),
        corrected_audit_helper_hashes=read(
            root / "records/continuation_helpers_v3.json"
        ),
        preflight_record=dict(
            path="records/preflight.json",
            sha256=digest(root / "records/preflight.json"),
        ),
        fixed_source_records=bindings,
        files=sources,
        checks=checks,
        focused_tests=test_result,
        preservation=preservation_result,
        command_records=[
            dict(
                path=str(p.relative_to(root)),
                sha256=digest(p),
                argv=read(p).get("argv"),
            )
            for p in (root / "records").glob("*.json")
            if isinstance(read(p), dict)
            and read(p).get("schema") == "circyto.cmpb.a02.stage.v1"
        ],
        environment_record="records/preflight.json",
        prior_failures_retained=True,
        old_full_suite_reused_not_rerun="523 passed, 10 skipped, 9 warnings from prior repair review only",
    )
    for filename, value in [
        ("pilot_progress.json", progress_record),
        ("resource_usage.json", resource),
        ("source_manifest.json", manifest),
    ]:
        dump(out / filename, sanitize(value, root, original, repo))
    with (out / "pilot_cell_trace.tsv").open("w", newline="") as h:
        writer = csv.DictWriter(
            h, fieldnames=list(trace[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(sanitize(trace, root, original, repo))
    dump(
        root / "records/packet_result.json",
        dict(
            complete=complete,
            summary=summary,
            checks=checks,
            tests=test_result,
            preservation=preservation_result,
        ),
    )
    print(
        json.dumps(
            dict(
                complete=complete,
                checks=checks,
                focused_tests=test_result,
                footprint=resource["current_whole_A02_footprint_bytes"],
            )
        )
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--original", type=Path, required=True)
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--audit-directory", default="records/acceptance")
    a = ap.parse_args()
    packet(a.root.resolve(), a.original.resolve(), a.repo.resolve(), a.audit_directory)
