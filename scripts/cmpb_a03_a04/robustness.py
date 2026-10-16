#!/usr/bin/env python3
"""Sequential fixed A04 grid; planned cases and every arm outcome persist."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import time

from execute import run_arm
from perturb import plan, generate
from robustness_oracle import verify
from run_stage import dump, digest, stamp, footprint

ARMS = ["circyto", "reference", "unsafe"]


def archive(directory):
    """Lossless own-output retention; never accepts paths outside this case tree."""
    files = {
        str(p.relative_to(directory)): digest(p)
        for p in directory.rglob("*")
        if p.is_file()
    }
    target = directory.with_suffix(".tar.gz")
    if target.exists():
        raise FileExistsError(target)
    with tarfile.open(target, "w:gz", compresslevel=1) as t:
        for name in sorted(files):
            t.add(directory / name, arcname=name, recursive=False)
    with tarfile.open(target, "r:gz") as t:
        for name, h in files.items():
            stream = t.extractfile(name)
            actual = hashlib.sha256()
            for chunk in iter(lambda: stream.read(1048576), b""):
                actual.update(chunk)
            if actual.hexdigest() != h:
                raise ValueError("Archive member corruption: " + name)
    record = dict(
        archive=str(target),
        sha256=digest(target),
        bytes=target.stat().st_size,
        members=files,
        lossless_member_hashes_verified=True,
    )
    dump(target.with_suffix(".manifest.json"), record)
    # Remove only this just-created redundant uncompressed case, after hash equality.
    shutil.rmtree(directory)
    return record


def freeze(root, label="frozen_a04_sources"):
    files = list(Path(__file__).parent.glob("*.py")) + [
        Path(__file__).parents[2] / "manuscript/cmpb/a03_a04" / n
        for n in [
            "execution_amendment.md",
            "reference_contract.md",
            "a04_contract_details.md",
        ]
    ]
    if label in ["frozen_a04_continuation_sources", "frozen_a04_final_sources"]:
        files += [Path(__file__).parents[2] / "manuscript/cmpb/a03_a04" / n for n in ["storage_amendment.md", "output_audit_amendment.md"]]
    if label == "frozen_a04_final_sources":
        files.append(Path(__file__).parents[2] / "manuscript/cmpb/a03_a04/time_amendment.md")
    p = root / "records" / (label + ".json")
    if p.exists():
        raise FileExistsError(p)
    hashes = {str(f): digest(f) for f in files}
    for f in files:
        dest = root / "helpers" / (hashes[str(f)] + "-" + f.name)
        if not dest.exists():
            shutil.copyfile(f, dest)
    dump(p, dict(frozen_utc=stamp(), sources=hashes))


def append(path, row):
    exists = path.exists()
    with path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row), delimiter="\t", lineterminator="\n")
        if not exists:
            w.writeheader()
        w.writerow(row)
        f.flush()


def checkpoint_completed_cases(root, segment="a04"):
    if (root / "records" / (segment + "_stop.json")).exists():
        raise RuntimeError("Scientific stop prevents administrative continuation")
    """Skip only complete, retained original cases; never bypass a scientific stop."""
    if (root / "records/a04_stop.json").exists():
        raise RuntimeError("Scientific stop prevents administrative continuation")
    if not (root / "records/storage_amendment.json").exists():
        raise RuntimeError("Recorded owner storage amendment required")
    old = root / segment
    with (old / "robustness_outcomes.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    arms = {}
    for row in rows:
        arms.setdefault(row["case"], set()).add(row["arm"])
    with (old / "archives.tsv").open() as f:
        retained = list(csv.DictReader(f, delimiter="\t"))
    complete = set()
    for record in retained:
        path = root / record["archive"]
        manifest = json.loads(path.with_suffix(".manifest.json").read_text())
        if not path.is_file() or manifest["sha256"] != record["sha256"] or not manifest["lossless_member_hashes_verified"]:
            raise RuntimeError("Checkpoint retention identity mismatch")
        if arms.get(record["case"]) == set(ARMS):
            complete.add(record["case"])
    return complete


def run(root, continuation=False, final_segment=False):
    if continuation and final_segment:
        raise ValueError("Only one administrative segment may be selected")
    label = "a04_final" if final_segment else "a04_continued" if continuation else "a04"
    folder = root / label
    cases = plan()
    if continuation or final_segment:
        completed_cases = checkpoint_completed_cases(root)
        if final_segment:
            if not (root / "records/time_checkpoint.json").exists():
                raise RuntimeError("Recorded authorized time checkpoint required")
            completed_cases |= checkpoint_completed_cases(root, "a04_continued")
        cases = [case for case in cases if case["case"] not in completed_cases]
        folder.mkdir(exist_ok=False)
    freeze(root, "frozen_a04_final_sources" if final_segment else "frozen_a04_continuation_sources" if continuation else "frozen_a04_sources")
    dump(
        folder / "planned_cases.json",
        dict(
            cases=cases,
            arms=ARMS,
            total_cases=len(cases),
            total_arm_attempts=len(cases) * 3,
            seed_replicates_are_not_biological=True,
        ),
    )
    outcomes = folder / "robustness_outcomes.tsv"
    prevalence = folder / "completed_cohort_prevalence.tsv.gz"
    if outcomes.exists() or prevalence.exists():
        raise FileExistsError("New grid output required")
    completed = 0
    budget = json.loads((root / "records/a04_budget.json").read_text())
    new_cap = json.loads((root / "records/budget.json").read_text())[
        "new_disk_stop_bytes"
    ]
    for case in cases:
        if time.time() >= budget["deadline_epoch"] or footprint(root) >= new_cap:
            raise RuntimeError("A04 resource budget exhausted")
        directory = folder / "cases" / case["case"]
        config, details = generate(root / "a04/base", case, directory / "inputs")
        full_counts = dict(
            zip(
                details.get("full_cohort_candidates", []),
                details.get("full_cohort_positive_counts", []),
            )
        )
        with gzip.open(prevalence, "at", newline="", compresslevel=1) as pfile:
            for arm in ARMS:
                r, out = run_arm(directory / "inputs/config.json", arm, directory / arm)
                check = {}
                rows = []
                if r["exit_code"] == 0 and case["valid"]:
                    try:
                        check, rows = verify(config, out, arm)
                    except Exception as exc:
                        check = dict(passed=False, oracle_exception=repr(exc))
                        dump(directory / arm / "oracle_exception.json", check)
                classification = (
                    ("accepted_valid" if r["exit_code"] == 0 else "false_rejection")
                    if case["valid"]
                    else (
                        "invalid_accepted"
                        if r["exit_code"] == 0
                        else "correct_rejection"
                    )
                )
                record = dict(
                    case=case["case"],
                    scenario=case["scenario"],
                    variant=case["variant"],
                    model=case.get("model", "deterministic"),
                    nominal_fraction=case.get("nominal_fraction", ""),
                    seed=case.get("seed", ""),
                    selected_count=details.get("selected_count", ""),
                    realized_fraction=details.get("realized_fraction", ""),
                    moved_identities=details.get("moved_identities", 0),
                    arm=arm,
                    expected_valid=case["valid"],
                    outcome=classification,
                    exit_code=r["exit_code"],
                    value_checks_passed=check.get("passed", "not_applicable"),
                    count_errors=check.get("count_errors", ""),
                    annotation_errors=check.get("annotation_errors", ""),
                    cell_assignment_errors=check.get("cell_assignment_errors", ""),
                    state_errors=check.get("status_errors", ""),
                    report_errors=check.get("report_errors", ""),
                    source_errors=check.get("source_errors", ""),
                    contract_errors=check.get("contract_errors", ""),
                    unavailable_fields=json.dumps(check.get("unavailable_fields", [])),
                    undefined_candidates=check.get("undefined_candidates", ""),
                    oracle_exception=check.get("oracle_exception", ""),
                    elapsed_seconds=r["elapsed_seconds"],
                    peak_rss_kib=r["peak_rss_kib"],
                    execution_record=str(
                        (directory / arm / "execution.json").relative_to(root)
                    ),
                    retention="See archives.tsv; unexpected failed case remains unpacked",
                )
                append(outcomes, record)
                for v in rows:
                    v = dict(
                        case=case["case"],
                        scenario=case["scenario"],
                        model=case.get("model", "deterministic"),
                        nominal_fraction=case.get("nominal_fraction", ""),
                        realized_fraction=details.get("realized_fraction", ""),
                        seed=case.get("seed", ""),
                        arm=arm,
                        **v,
                        full_cohort_numerator=full_counts.get(v["candidate"], ""),
                        full_cohort_denominator=23 if full_counts else "",
                    )
                    writer = csv.DictWriter(
                        pfile, fieldnames=list(v), delimiter="\t", lineterminator="\n"
                    )
                    if completed == 0 and pfile.tell() == 0:
                        writer.writeheader()
                    writer.writerow(v)
                pfile.flush()
                completed += 1
                dump(
                    root / "records" / (label + "_progress.json"),
                    dict(
                        updated_utc=stamp(),
                        planned_cases=len(cases),
                        planned_arm_attempts=len(cases) * 3,
                        completed_arm_attempts=completed,
                        current_case=case["case"],
                        last_arm=arm,
                        last_outcome=classification,
                    ),
                )
                if arm != "unsafe" and (
                    classification in ["false_rejection", "invalid_accepted"]
                    or case["valid"]
                    and not check.get("passed")
                ):
                    dump(
                        root / "records" / (label + "_stop.json"),
                        dict(
                            record=record,
                            reason="Production defect requires owner review"
                            if arm == "circyto"
                            else "Independent reference/oracle correction requires dated review",
                            retained_directory=str(directory),
                        ),
                    )
                    raise RuntimeError(
                        "STOP: unexpected " + arm + " result in " + case["case"]
                    )
        a = archive(directory)
        append(
            folder / "archives.tsv",
            dict(
                case=case["case"],
                archive=str(Path(a["archive"]).relative_to(root)),
                sha256=a["sha256"],
                bytes=a["bytes"],
                members=len(a["members"]),
            ),
        )
        if completed % 30 == 0:
            print(
                json.dumps(
                    dict(
                        completed=completed,
                        planned=len(cases) * 3,
                        case=case["case"],
                        new_bytes=footprint(root),
                    )
                ),
                flush=True,
            )
    dump(
        root / "records" / (label + "_complete.json"),
        dict(
            passed=True,
            completed_arm_attempts=completed,
            planned_arm_attempts=len(cases) * 3,
            finished_utc=stamp(),
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--continuation", action="store_true")
    ap.add_argument("--final-segment", action="store_true")
    a = ap.parse_args()
    run(a.root.resolve(), a.continuation, a.final_segment)
