#!/usr/bin/env python3
"""Factual packet from retained ledgers; no scientific analysis or rerun."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import os
import platform
import re
import shutil
import statistics
import subprocess
import xml.etree.ElementTree as ET

from run_stage import digest, dump, footprint, stamp

SEGMENTS = ["a04", "a04_continued", "a04_final"]


def read(path, default=None):
    return json.loads(path.read_text()) if path.exists() else default


def sanitize(value, root, repo):
    if isinstance(value, str):
        for a, b in [
            (str(root), "${EVAL_ROOT}"),
            (str(root.parent), "${A02_ROOT}"),
            (str(repo), "${ANALYSIS_WORKTREE}"),
            ("/tmp/circyto-import/load_work", "${ORIGINAL_OBJECT_ROOT}"),
        ]:
            value = value.replace(a, b)
        value = re.sub(r"/sys/fs/cgroup/user.slice/user-\d+\.slice/user@\d+\.service", "${USER_CGROUP}", value)
        return value
    if isinstance(value, dict):
        return {
            sanitize(k, root, repo): sanitize(v, root, repo) for k, v in value.items()
        }
    if isinstance(value, list):
        return [sanitize(v, root, repo) for v in value]
    return value


def write_table(path, rows, fields=None):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def prevalence_rows(path, integrity):
    """Read retained rows, explicitly flagging an interrupted gzip stream.

    This never repairs, replaces or truncates the original evidence. A missing
    final gzip footer remains a retention failure even if earlier rows decode.
    """
    integrity.update(gzip_integrity=True, error=None)
    try:
        with gzip.open(path, "rt", newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if None in row or None in row.values():
                    raise ValueError("Incomplete prevalence row")
                yield row
    except (EOFError, gzip.BadGzipFile, UnicodeError, ValueError) as exc:
        integrity.update(gzip_integrity=False, error=repr(exc))


def active_segment(root):
    return next((s for s in reversed(SEGMENTS) if (root / s).exists()), "a04")


def case_segment(root, case):
    return next((s for s in reversed(SEGMENTS)
                 if (root / s / "cases" / case).exists()
                 or (root / s / "cases" / (case + ".tar.gz")).exists()), "a04")


def grid_complete(root, observed, planned):
    segment = active_segment(root)
    stage_name = {"a04": "a04-grid", "a04_continued": "a04-continued-grid", "a04_final": "a04-final-grid"}[segment]
    record = read(root / "records" / (segment + "_complete.json"), {})
    stage = read(root / "records" / (stage_name + ".json"), {})
    return bool(
        {(r["case"], r["arm"]) for r in observed}
        == {(c["case"], arm) for c in planned["cases"] for arm in planned["arms"]}
        and len(observed) == len(planned["cases"]) * len(planned["arms"])
        and record.get("passed") is True
        and stage.get("launcher_exit_code") == 0
        and not any((root / "records" / (s + "_stop.json")).exists() for s in SEGMENTS)
    )


def pending_outcome(root, case, arm):
    directory = root / case_segment(root, case) / "cases" / case / arm
    execution = read(directory / "execution.json")
    if execution is not None:
        return "EXECUTED_CHECK_INCOMPLETE", execution
    if directory.exists():
        return "STARTED_NO_COMPLETION_RECORD", None
    return "NOT_STARTED", None


def grid_rows(root):
    """Choose the latest *case* attempt; retain every earlier administrative arm."""
    selected, prior, preferred = [], [], {}
    for segment in SEGMENTS:
        path = root / segment / "robustness_outcomes.tsv"
        if not path.exists():
            continue
        with path.open() as f:
            for r in csv.DictReader(f, delimiter="\t"):
                r["segment"] = segment
                choice = case_segment(root, r["case"])
                preferred[r["case"]] = choice
                (selected if segment == choice else prior).append(r)
    old_keys = {(r["case"], r["arm"], r["segment"]) for r in prior}
    for case, latest in preferred.items():
        for segment in SEGMENTS[:SEGMENTS.index(latest)]:
            old = root / segment / "cases" / case
            if not old.is_dir():
                continue
            for arm in ["circyto", "reference", "unsafe"]:
                if (case, arm, segment) in old_keys or not (old / arm).exists():
                    continue
                execution = read(old / arm / "execution.json")
                prior.append(dict(case=case, arm=arm, segment=segment,
                    outcome="EXECUTED_CHECK_INCOMPLETE" if execution else "STARTED_NO_COMPLETION_RECORD",
                    exit_code=execution.get("exit_code", "") if execution else ""))
    return selected, prior, preferred


def combined_prevalence(root, preferred, integrity):
    """Make a new lossless primary-attempt table, preserving all source streams."""
    folder = root / "a04_combined"
    folder.mkdir(exist_ok=True)
    target = folder / "completed_cohort_prevalence.tsv.gz"
    if target.exists():
        raise FileExistsError("Combined source table already recorded")
    integrity.update(gzip_integrity=True, error=None, sources=[])
    with gzip.open(target, "wt", newline="", compresslevel=1) as handle:
        writer = None
        for segment in SEGMENTS:
            path = root / segment / "completed_cohort_prevalence.tsv.gz"
            if not path.exists():
                continue
            status = dict(segment=segment, sha256=digest(path), selected_rows=0, prior_attempt_rows=0)
            for row in prevalence_rows(path, status):
                if preferred.get(row["case"], segment) != segment:
                    status["prior_attempt_rows"] += 1
                    continue
                if writer is None:
                    writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                writer.writerow(row)
                status["selected_rows"] += 1
                yield row
            integrity["sources"].append(status)
            integrity["gzip_integrity"] &= status["gzip_integrity"]


def preserve(root, repo):
    checks = {}
    for name in ["inputs_before", "protected_before"]:
        before = read(root / "records" / f"{name}.json")
        bad = []
        for p, r in before.items():
            expected = r["sha256"] if isinstance(r, dict) else r
            if not Path(p).exists() or digest(p) != expected:
                bad.append(p)
        checks[name] = dict(
            files=len(before),
            sha256_mismatches=bad,
            manifest_sha256=digest(root / "records" / f"{name}.json"),
        )
    before = read(root / "records/A02_inventory_before.json")
    after = {}
    for directory, dirs, files in os.walk(root.parent):
        if Path(directory) == root:
            dirs[:] = []
            continue
        for name in files:
            p = Path(directory) / name
            s = p.lstat()
            after[str(p.relative_to(root.parent))] = [
                s.st_size,
                s.st_mtime_ns,
                os.readlink(p) if p.is_symlink() else None,
            ]
    changed = [p for p in set(before) | set(after) if before.get(p) != after.get(p)]
    checks["A02_inventory"] = dict(
        before_entries=len(before),
        after_entries=len(after),
        changed_entries=changed,
        scope="path/size/mtime/symlink preservation; not a fresh full-tree SHA assertion",
    )
    previous = read(repo / "manuscript/cmpb/a02_complete/source_manifest.json")[
        "preservation"
    ]["original_H5MU_sha256"]
    objects = {}
    for logical, expected in previous.items():
        path = Path(
            logical.replace("${ORIGINAL_OBJECT_ROOT}", "/tmp/circyto-import/load_work")
        )
        objects[logical] = dict(expected=expected, observed=digest(path))
    checks["original_H5MU"] = objects
    environment = read(root / "records/environment.json")
    site = Path(environment["module_origins"]["circyto"]).parents[1]
    members = environment["installed_python_members"]
    checks["installed_python_members"] = dict(
        files=len(members),
        sha256_mismatches=[name for name, expected in members.items() if digest(site / name) != expected],
    )
    refs = read(root / "records/refs_before.json")
    actual = dict(
        x.split()
        for x in subprocess.check_output(
            ["git", "for-each-ref", "--format=%(refname) %(objectname)", "refs/heads/"],
            cwd=repo,
            text=True,
        ).splitlines()
    )
    checks["original_ref_changes"] = {
        k: [v, actual.get(k)]
        for k, v in refs.items()
        if k != "refs/heads/analysis/cmpb-a03-a04" and actual.get(k) != v
    }
    checks["passed"] = (
        not changed
        and not checks["original_ref_changes"]
        and not checks["installed_python_members"]["sha256_mismatches"]
        and all(v["observed"] == v["expected"] for v in objects.values())
        and all(
            not checks[k]["sha256_mismatches"]
            for k in ["inputs_before", "protected_before"]
        )
    )
    return checks


def main(root, repo):
    out = repo / "manuscript/cmpb/a03_a04"
    out.mkdir(exist_ok=True)
    for label, source_name in [
        ("storage_checkpoint", "circyto-a0304-storage-checkpoint.py"),
        ("time_amendment", "circyto-a0304-time-amendment.py"),
        ("time_checkpoint", "circyto-a0304-final-checkpoint.py"),
    ]:
        source = Path("/tmp") / source_name
        if source.exists():
            source_sha = digest(source)
            retained = root / "helpers" / (label + "-" + source_sha + ".py")
            if not retained.exists():
                shutil.copyfile(source, retained)
            dump(root / "records" / (label + "_source.json"), dict(
                original_argv=["python3", str(source)], retained_source=str(retained),
                sha256=source_sha, operation="Owner-authorized administrative checkpoint/amendment; source retained with bookkeeping",
            ))
    preservation = preserve(root, repo)
    dump(root / "records/preservation_after.json", preservation)
    if not preservation["passed"]:
        raise RuntimeError("Original preservation mismatch; review required")
    control = read(root / "records/a03_control.json", {})
    measure = read(root / "records/a03_measure.json", {})
    rows = []
    for r in control.get("runs", []) + measure.get("runs", []):
        rows.append(
            dict(
                phase=r["phase"],
                arm=r["arm"],
                measured=r["measured"],
                wall_seconds=r["elapsed_seconds"],
                peak_rss_kib=r["peak_rss_kib"],
                user_seconds=r["user_seconds"],
                system_seconds=r["system_seconds"],
                input_bytes=r.get("input_bytes", ""),
                output_bytes=r["output_bytes"],
                exit_code=r["exit_code"],
                equivalent_required_endpoint=r["acceptance"]["passed"],
                allowed_CPUs="0,1",
                BLAS_OpenMP_threads=1,
                cache_policy="one separate warmup per arm; warm OS cache; no privileged flush; process startup included",
                record="a03/" + r["phase"] + "-" + r["arm"] + "/execution.json",
            )
        )
    write_table(out / "reference_comparison.tsv", rows)
    operations = []
    for name, common, c, r in [
        (
            "Input mapping",
            "Same biological cell/sample/native-file mapping and explicit metadata",
            "User supplies mapping TSV",
            "User supplies same mapping TSV",
        ),
        (
            "Completion evidence",
            "Same A02 bound statuses; unknown never inferred complete",
            "User supplies cells TSV",
            "User supplies same cells TSV",
        ),
        (
            "Native parsing/validation",
            "Same 13-column format and declarations",
            "Installed import-circ implementation",
            "Bespoke reference CSV/identity/declaration validation",
        ),
        (
            "Sparse object construction",
            "Exact integers; ID-aware axes",
            "Installed exporter",
            "Bespoke construction with ordinary SciPy/AnnData",
        ),
        (
            "RNA integration",
            "Same full RNA, own axis, maps, fields",
            "Installed MuData integration",
            "Ordinary MuData API with explicit global membership",
        ),
        (
            "Evidence",
            "Hash inputs before/after; copy non-RNA source bytes",
            "Installed provenance implementation",
            "Bespoke hashing/copy checks",
        ),
        (
            "QC and HTML",
            "Per-cell/candidate QC TSV, JSON, offline first-100 HTML",
            "Installed report",
            "Bespoke report, different cosmetic format",
        ),
        (
            "Annotation computation",
            "No GTF computation in either measured arm",
            "Not performed",
            "Not performed",
        ),
        (
            "Human effort",
            "Not measured; no code-length or usability endpoint",
            "Not assessed",
            "Not assessed",
        ),
    ]:
        operations.append(
            dict(
                operation=name,
                shared_required_work=common,
                circyto=c,
                independent_reference=r,
                interpretation="Descriptive operation inventory, not a human-effort score",
            )
        )
    write_table(out / "operation_inventory.tsv", operations)
    planned = read(
        root / "a04/planned_cases.json",
        dict(cases=[], arms=["circyto", "reference", "unsafe"]),
    )
    observed, prior_attempts, preferred = grid_rows(root)
    if prior_attempts:
        prior_fields = sorted(set().union(*(r.keys() for r in prior_attempts)))
        write_table(out / "administrative_prior_attempts.tsv", sanitize(prior_attempts, root, repo), prior_fields)
    accounted = {(r["case"], r["arm"]): r for r in observed}
    fields = (
        list(observed[0])
        if observed
        else [
            "case",
            "scenario",
            "variant",
            "model",
            "nominal_fraction",
            "seed",
            "arm",
            "expected_valid",
            "outcome",
        ]
    )
    allrows = []
    for c in planned["cases"]:
        for arm in planned["arms"]:
            r = accounted.get((c["case"], arm))
            if r is None:
                pending, execution = pending_outcome(root, c["case"], arm)
                r = {k: "" for k in fields}
                r.update(
                    case=c["case"],
                    scenario=c["scenario"],
                    variant=c["variant"],
                    model=c.get("model", "deterministic"),
                    nominal_fraction=c.get("nominal_fraction", ""),
                    seed=c.get("seed", ""),
                    arm=arm,
                    expected_valid=c["valid"],
                    outcome=pending,
                )
                if execution:
                    for key in ["exit_code", "elapsed_seconds", "peak_rss_kib"]:
                        if key in fields:
                            r[key] = execution.get(key, "")
            allrows.append(r)
    if allrows:
        write_table(
            out / "robustness_outcomes.tsv", sanitize(allrows, root, repo), fields
        )
    summary = defaultdict(Counter)
    for r in allrows:
        key = (
            r["scenario"],
            r.get("model", ""),
            r.get("nominal_fraction", ""),
            r["arm"],
        )
        s = summary[key]
        s["planned"] += 1
        s[r["outcome"]] += 1
        if r.get("value_checks_passed") == "True":
            s["all_value_checks_passed"] += 1
        if r.get("value_checks_passed") == "False":
            s["value_contract_failed"] += 1
        for k in [
            "count_errors",
            "annotation_errors",
            "cell_assignment_errors",
            "state_errors",
            "report_errors",
            "source_errors",
            "contract_errors",
        ]:
            if r.get(k, "") != "":
                s[k] += int(r[k])
    figure = [
        dict(
            panel="Figure3A-B",
            scenario=k[0],
            model=k[1],
            nominal_fraction=k[2],
            arm=k[3],
            metric=metric,
            value=value,
            scope="simulation configurations; seeds are not biological replicates",
        )
        for k, s in sorted(summary.items(), key=lambda x: str(x[0]))
        for metric, value in sorted(s.items())
    ]
    for arm in ["circyto", "reference"]:
        runs = [r for r in measure.get("runs", []) if r["arm"] == arm and r["measured"]]
        for metric, values in [
            ("wall_seconds", [r["elapsed_seconds"] for r in runs]),
            ("peak_rss_kib", [r["peak_rss_kib"] for r in runs]),
        ]:
            if values:
                for statistic, value in [
                    ("minimum", min(values)),
                    ("median", statistics.median(values)),
                    ("maximum", max(values)),
                ]:
                    figure.append(
                        dict(
                            panel="Figure4B",
                            scenario="A02-real-four-cell",
                            model="warm cache/fresh process",
                            nominal_fraction="",
                            arm=arm,
                            metric=metric + "_" + statistic,
                            value=value,
                            scope="three process runs; not scaling or biological replicates",
                        )
                    )
    full = root / "a04_combined/completed_cohort_prevalence.tsv.gz"
    extract = []
    nrows = 0
    seen = set()
    prevalence_summary = defaultdict(Counter)
    prevalence_integrity = dict(gzip_integrity=None, error="file absent")
    if (root / "a04/completed_cohort_prevalence.tsv.gz").exists():
        for r in combined_prevalence(root, preferred, prevalence_integrity):
            nrows += 1
            key = (r["case"], r["arm"])
            group = (r["scenario"], r["model"], r["nominal_fraction"], r["arm"])
            if key not in seen:
                extract.append(r)
                seen.add(key)
            s = prevalence_summary[group]
            s["candidate_rows"] += 1
            if r["oracle_undefined"] == "True":
                s["oracle_undefined"] += 1
            if r["arm_undefined"] == "True":
                s["arm_undefined"] += 1
            if r["absolute_error_numerator"] != "" and int(
                r["absolute_error_numerator"]
            ):
                s["nonzero_same_cohort_error"] += 1
        if extract:
            write_table(out / "completed_cohort_prevalence.tsv", extract)
    for k, s in sorted(prevalence_summary.items(), key=lambda x: str(x[0])):
        for metric, value in sorted(s.items()):
            figure.append(
                dict(
                    panel="Figure3C",
                    scenario=k[0],
                    model=k[1],
                    nominal_fraction=k[2],
                    arm=k[3],
                    metric=metric,
                    value=value,
                    scope="candidate/configuration rows; full-cohort composition is separate",
                )
            )
    write_table(out / "figure_source_summary.tsv", figure)
    extra = read(root / "records/additional_output_audit.json", {})
    additional_rows = [
        dict(
            case=r["case"],
            arm=r["arm"],
            passed=r["passed"],
            checks=len(r["checks"]),
            failed_checks=";".join(k for k, value in r["checks"].items() if not value),
        )
        for r in extra.get("outcomes", [])
    ]
    if additional_rows:
        write_table(out / "additional_output_checks.tsv", additional_rows)
    stages = []
    for p in (root / "records").glob("*.json"):
        j = read(p)
        if (
            isinstance(j, dict)
            and j.get("schema") == "circyto.cmpb.a02.stage.v1"
            and "elapsed_seconds" in j
        ):
            g = read(root / "records" / (j["stage"] + ".guard.json"), {})
            stages.append(
                dict(
                    stage=j["stage"],
                    status=j["status"],
                    elapsed_seconds=j["elapsed_seconds"],
                    started_utc=j["started_utc"],
                    finished_utc=j["finished_utc"],
                    exit_code=j["launcher_exit_code"],
                    memory_peak=g.get("memory_current_peak_bytes"),
                    RSS_peak=g.get("process_tree_rss_peak_bytes"),
                    disk_peak=g.get("disk_peak_bytes"),
                    argv=j["argv"],
                    log=j["log"],
                    record_sha256=digest(p),
                )
            )
    stages.sort(key=lambda x: x["started_utc"])
    resource = dict(
        session_budget=read(root / "records/budget.json"),
        a04_budget=read(root / "records/a04_budget.json"),
        stages=stages,
        active_stage_seconds=sum(x["elapsed_seconds"] for x in stages),
        control_probes=[dict(
            probe=p["probe"], passed=p["passed"], exit_code=p["exit_code"],
            surviving_test_pids=p["surviving_test_pids"],
            expected_stop_reason=p["guard"].get("stop_reason"),
            record="records/a0304-probe-" + p["probe"] + ".json",
        ) for p in read(root / "records/a0304-probe_control_tests.json", [])],
        control_probe_manifest_sha256=digest(root / "records/a0304-probe_control_tests.json"),
        nonzero_stage_interpretation="Control probes intentionally test nonzero exits/termination; a04-grid and a04-continued-grid exit 130 are recorded owner-authorized administrative checkpoints. Scientific arm exits and acceptance remain separately tabulated.",
        new_storage_bytes=footprint(root),
        whole_A02_storage_bytes=footprint(root.parent),
        recorded_utc=stamp(),
        measurement_limits="Cgroup sampled peak includes cache; summed process RSS is sampled, timed child peak RSS comes from /usr/bin/time. Disk watchdog is not filesystem quota. Human pauses excluded from algorithm timings.",
    )
    dump(out / "resource_usage.json", sanitize(resource, root, repo))
    env = read(root / "records/environment.json")
    env["recorded_host_properties"] = dict(
        system=platform.system(),
        release=platform.release(),
        machine=platform.machine(),
        cpu_model=next(
            (
                line.split(":", 1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")
            ),
            "not recorded",
        ),
        process_affinity=sorted(os.sched_getaffinity(0)),
        source="same-host final bookkeeping; no hostname or credentials recorded",
    )
    dump(out / "environment.json", sanitize(env, root, repo))
    manifests = {}
    for prefix in ["a03", *SEGMENTS, "a04_combined"]:
        for p in (root / prefix).rglob("*"):
            if p.is_file():
                manifests[str(p.relative_to(root))] = dict(
                    bytes=p.stat().st_size, sha256=digest(p)
                )
    final_tests = root / "records/a04-handoff-tests.xml"
    tests = ET.parse(final_tests if final_tests.exists() else root / "records/a04-helper-tests.xml").getroot()
    suites = list(tests.iter("testsuite"))
    test_counts = {
        k: sum(int(s.get(k, 0)) for s in suites)
        for k in ["tests", "failures", "errors", "skipped"]
    }
    complete = grid_complete(root, observed, planned) and bool(
        prevalence_integrity["gzip_integrity"]
    ) and extra.get("passed") is True
    result = dict(
        A02_preservation_checks_passed=True,
        A03=control.get("passed", False) and measure.get("passed", False),
        A04_complete=complete,
        A04_attempts=len(observed),
        A04_planned=len(planned["cases"]) * 3,
        A04_stop=read(root / "records/a04_stop.json"),
        A04_stage=read(root / "records/a04-grid.json"),
        A04_guard=read(root / "records/a04-grid.guard.json"),
        A04_continuation_stage=read(root / "records/a04-continued-grid.json"),
        A04_continuation_guard=read(root / "records/a04-continued-grid.guard.json"),
        A04_continuation_stop=read(root / "records/a04_continued_stop.json"),
        A04_final_stage=read(root / "records/a04-final-grid.json"),
        A04_final_guard=read(root / "records/a04-final-grid.guard.json"),
        A04_final_stop=read(root / "records/a04_final_stop.json"),
        time_amendment=read(root / "records/time_amendment.json"),
        time_checkpoint=read(root / "records/time_checkpoint.json"),
        time_amendment_source=read(root / "records/time_amendment_source.json"),
        time_checkpoint_source=read(root / "records/time_checkpoint_source.json"),
        storage_amendment=read(root / "records/storage_amendment.json"),
        storage_checkpoint_source=read(root / "records/storage_checkpoint_source.json"),
        prior_administrative_attempt_rows=len(prior_attempts),
        focused_tests=test_counts,
        full_suite="Prior 523-test production suite is reused; not rerun",
        preservation=preservation,
        prevalence_decoded_rows=nrows,
        prevalence_integrity=prevalence_integrity,
        additional_output_audit=dict(
            performed=bool(extra),
            passed=extra.get("passed"),
            stopped=extra.get("stopped"),
            reason=extra.get("reason"),
            rows=len(additional_rows),
            source=read(root / "records/additional_output_audit_source.json"),
            full_record_sha256=digest(root / "records/additional_output_audit.json") if extra else None,
        ),
        prevalence_review_extract_rule="First candidate in each attempted valid configuration/arm; complete lossless gzip is local and hashed",
        completed_cohort_table=dict(
            path=str(full.relative_to(root)),
            sha256=digest(full) if full.exists() else None,
            bytes=full.stat().st_size if full.exists() else 0,
        ),
        files=manifests,
        input_manifest=read(root / "records/inputs_before.json"),
        A03_frozen_sources=read(root / "records/frozen_a03_sources.json"),
        A04_frozen_sources=read(root / "records/frozen_a04_sources.json"),
        A04_continuation_frozen_sources=read(root / "records/frozen_a04_continuation_sources.json"),
        A04_final_frozen_sources=read(root / "records/frozen_a04_final_sources.json"),
        new_helper_source_hashes={
            p.name: digest(p) for p in Path(__file__).parent.glob("*.py")
        },
        base_head="307758f57a39c009393d03fc0f69e4638b2498c6",
        production_source="e6843a38469d49f96d44758a4e37893ebbe2734b",
    )
    dump(out / "input_output_manifest.json", sanitize(result, root, repo))
    dump(
        root / "records/handoff_result.json",
        dict(
            A03=result["A03"],
            A04_complete=complete,
            attempts=len(observed),
            planned=result["A04_planned"],
            preservation_passed=preservation["passed"],
            focused_tests=test_counts,
        ),
    )
    print(
        json.dumps(
            dict(
                A03=result["A03"],
                A04_complete=complete,
                attempts=len(observed),
                planned=result["A04_planned"],
                new_bytes=resource["new_storage_bytes"],
            )
        )
    )


def close_resource_record(root, repo):
    """Include this packet writer's completed ledger, without another evaluation."""
    path = repo / "manuscript/cmpb/a03_a04/resource_usage.json"
    resource = read(path)
    known = {x["stage"] for x in resource["stages"]}
    for record_path in sorted((root / "records").glob("*.json")):
        stage = read(record_path)
        if not isinstance(stage, dict) or stage.get("schema") != "circyto.cmpb.a02.stage.v1":
            continue
        if "elapsed_seconds" not in stage or stage["stage"] in known:
            continue
        guard = read(root / "records" / (stage["stage"] + ".guard.json"), {})
        resource["stages"].append(dict(
            stage=stage["stage"], status=stage["status"],
            elapsed_seconds=stage["elapsed_seconds"], started_utc=stage["started_utc"],
            finished_utc=stage["finished_utc"], exit_code=stage["launcher_exit_code"],
            memory_peak=guard.get("memory_current_peak_bytes"),
            RSS_peak=guard.get("process_tree_rss_peak_bytes"), disk_peak=guard.get("disk_peak_bytes"),
            argv=stage["argv"], log=stage["log"], record_sha256=digest(record_path),
        ))
    resource["stages"].sort(key=lambda x: x["started_utc"])
    resource["active_stage_seconds"] = sum(x["elapsed_seconds"] for x in resource["stages"])
    resource["new_storage_bytes"] = footprint(root)
    resource["whole_A02_storage_bytes"] = footprint(root.parent)
    resource["recorded_utc"] = stamp()
    resource["finalization_note"] = "Includes completed packet-writing stage; storage sampled before this small metadata update. This bookkeeping is not another measured algorithm run."
    dump(path, sanitize(resource, root, repo))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--close-resource-record", action="store_true")
    a = ap.parse_args()
    if a.close_resource_record:
        close_resource_record(a.root.resolve(), a.repo.resolve())
    else:
        main(a.root.resolve(), a.repo.resolve())
