#!/usr/bin/env python3
"""A02 continuation preparation only; no aligner, indexing, or detector launch."""

import argparse
import csv
import importlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import zipfile

from prepare_inputs import check_fastq, fetch
from replay_alignment_status import source_snapshot
from run_stage import digest, dump, footprint, stamp

WHEEL_SHA = "28c7cbde3c31cf90dd7c247487ba5f410c5b82f6f97d7a2281a0ce52912f7383"
SOURCE = "e6843a38469d49f96d44758a4e37893ebbe2734b"
ROWS = [
    (
        101,
        "GSM8558853",
        "SRR30918117",
        "SRX26321183",
        "SAMN44108154",
        3935595,
        594274845,
        230818973,
        "09ede82e0f86db6147ce82809bf7fb22",
    ),
    (
        102,
        "GSM8558854",
        "SRR30918116",
        "SRX26321184",
        "SAMN44108153",
        4002138,
        604322838,
        236808421,
        "ed7668e9efa214a1bda5d3ea2bf899eb",
    ),
    (
        103,
        "GSM8558855",
        "SRR30918115",
        "SRX26321185",
        "SAMN44108152",
        3954925,
        597193675,
        233399236,
        "f90306b3462108160d7e35d4b707729b",
    ),
]


def expected(row):
    n, gsm, run, experiment, sample, reads, bases, size, md5 = row
    return dict(
        run_accession=run,
        experiment_accession=experiment,
        sample_accession=sample,
        study_accession="PRJNA1169833",
        library_layout="SINGLE",
        read_count=str(reads),
        base_count=str(bases),
        fastq_ftp=f"ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/0{run[-2:]}/{run}/{run}.fastq.gz",
        fastq_bytes=str(size),
        fastq_md5=md5,
    )


def compare_metadata(observed):
    by_run = {r["run_accession"]: r for r in observed}
    if len(by_run) != len(observed) or set(by_run) != {r[2] for r in ROWS}:
        raise ValueError("Archive returned duplicate/missing/foreign runs")
    differences = []
    for row in ROWS:
        for key, value in expected(row).items():
            if by_run[row[2]].get(key) != value:
                differences.append(
                    dict(
                        run=row[2],
                        field=key,
                        registered=value,
                        observed=by_run[row[2]].get(key),
                    )
                )
    return differences


def preserved(original):
    result = source_snapshot(original)
    # Include all retained repair evidence/artifacts; environments are verified
    # separately by their complete installed Python members and distributions.
    repair = original / "status_revalidation_20261016"
    for path in sorted(repair.rglob("*")):
        if path.is_file() and not any(
            p in {"test-env", "installed-env", "__pycache__", "tmp", "cache"}
            for p in path.relative_to(repair).parts
        ):
            result[str(path.relative_to(original))] = {
                "sha256": digest(path),
                "bytes": path.stat().st_size,
                "symlink": os.readlink(path) if path.is_symlink() else None,
            }
    return result


def preflight(root, original, repo):
    final = original / "status_revalidation_20261016/final_candidate"
    prefix = final / "installed-env"
    if (
        Path(sys.prefix) != prefix
        or os.environ.get("PYTHONPATH")
        or os.environ.get("PYTHONNOUSERSITE") != "1"
    ):
        raise ValueError("Wrong isolated installed environment")
    wheel = final / "artifacts/circyto-0.10.0-py3-none-any.whl"
    if digest(wheel) != WHEEL_SHA:
        raise ValueError("Accepted artifact changed")
    old = json.loads((original / "records/environment_manifest.json").read_text())
    reviewed = json.loads(
        (
            repo
            / "manuscript/cmpb/a02_firstcell/status_revalidation/environment_manifest.json"
        ).read_text()
    )
    versions = {d.metadata["Name"]: d.version for d in metadata.distributions()}
    if versions != reviewed["resolved_distributions"]:
        raise ValueError("Resolved dependencies changed")
    origins = {
        n: str(Path(importlib.import_module(n).__file__).resolve())
        for n in reviewed["module_origins"]
    }
    if not all(Path(v).is_relative_to(prefix) for v in origins.values()):
        raise ValueError("Module imported from outside the repaired environment")
    members = {}
    with zipfile.ZipFile(wheel) as z:
        for n in z.namelist():
            if n.startswith("circyto/") and n.endswith(".py"):
                data = z.read(n)
                if data != (
                    Path(origins["circyto"]).parent.parent / n
                ).read_bytes() or data != subprocess.check_output(
                    ["git", "show", SOURCE + ":" + n], cwd=repo
                ):
                    raise ValueError("Source/wheel/installation mismatch: " + n)
                members[n] = digest(Path(origins["circyto"]).parent.parent / n)
    before = preserved(original)
    approved = json.loads(
        (final / "records/source_preservation_after.json").read_text()
    )
    if any(before.get(p) != rec for p, rec in approved.items()):
        raise ValueError("Preserved original scientific evidence changed")
    dump(root / "records/preservation_before.json", before)
    index = json.loads((original / "records/index_manifest.json").read_text())
    for p, rec in index["outputs"].items():
        if digest(original / p) != rec["sha256"]:
            raise ValueError("Original BWA index changed: " + p)
    for tool, rec in old["external_executables"].items():
        if digest(rec["path"]) != rec["sha256"]:
            raise ValueError("External executable changed: " + tool)
    if (
        digest(original / "external/CIRI3_Java_18.0.1.jar")
        != old["external_jar_sha256"]
    ):
        raise ValueError("External JAR changed")
    previous = json.loads(
        (
            repo
            / "manuscript/cmpb/a02_firstcell/status_revalidation/resource_usage.json"
        ).read_text()
    )
    totals = []
    for ledger in [original, original / "status_revalidation_20261016"]:
        records = [
            json.loads(p.read_text()) for p in (ledger / "records").glob("*.json")
        ]
        totals.append(
            sum(
                r["elapsed_seconds"]
                for r in records
                if isinstance(r, dict)
                and r.get("schema") == "circyto.cmpb.a02.stage.v1"
            )
        )
    if abs(sum(totals) - previous["combined_active_stage_elapsed_seconds"]) > 0.001:
        raise ValueError("Prior active resource ledgers do not reconcile")
    helpers = root / "helpers/runtime"
    helpers.mkdir(parents=True)
    for name in ["run_stage.py", "tool_audit.py"]:
        shutil.copyfile(Path(__file__).with_name(name), helpers / name)
    bindir = root / "audit-bin"
    bindir.mkdir()
    for tool, rec in old["external_executables"].items():
        argv = [
            "/usr/bin/python3",
            str(helpers / "tool_audit.py"),
            "--root",
            str(root),
            "--tool",
            tool,
            "--real",
            rec["path"],
            "--",
        ]
        wrapper = bindir / tool
        wrapper.write_text("#!/bin/sh\nexec " + shlex.join(argv) + ' "$@"\n')
        wrapper.chmod(0o755)
    env = dict(
        PATH=f"{bindir}:{prefix}/bin:{original}/env-bwa/bin:/usr/bin:/bin",
        PYTHONNOUSERSITE="1",
        CIRCYTO_CIRI3_JAR=str(original / "external/CIRI3_Java_18.0.1.jar"),
        CIRCYTO_CIRI3_JAVA=str(bindir / "java"),
        CIRCYTO_CIRI3_EXTRA_ARGS="-S 0 -Ma 0",
        JAVA_TOOL_OPTIONS="-Xmx16g -XX:ActiveProcessorCount=8",
    )
    os.environ.update(env)
    dump(root / "records/execution_environment.json", env)
    commands = []
    for name, argv, codes in [
        ("pip-check", [sys.executable, "-m", "pip", "check"], [0]),
        ("pip-freeze", [sys.executable, "-m", "pip", "freeze", "--all"], [0]),
        ("bwa-version", [str(bindir / "bwa")], [1]),
        ("java-version", [str(bindir / "java"), "-version"], [0]),
        ("doctor", [str(prefix / "bin/circyto"), "doctor"], [0]),
    ]:
        result = subprocess.run(argv, capture_output=True, text=True)
        output = root / "records" / (name + ".txt")
        output.write_text(result.stdout + result.stderr)
        commands.append(
            dict(
                argv=argv,
                exit_code=result.returncode,
                output=str(output),
                sha256=digest(output),
            )
        )
        if result.returncode not in codes:
            raise ValueError("Preflight command failed: " + name)
    from circyto.detectors.ciri3 import inspect_ciri3_runtime

    runtime = inspect_ciri3_runtime()
    if not runtime["direct_ready"] or runtime["jar"] != env["CIRCYTO_CIRI3_JAR"]:
        raise ValueError("Actual CIRI3 configuration differs")
    dump(
        root / "records/preflight.json",
        dict(
            recorded_utc=stamp(),
            wheel_sha256=WHEEL_SHA,
            source_commit=SOURCE,
            installed_members=members,
            versions=versions,
            origins=origins,
            external_executables=old["external_executables"],
            external_jar_sha256=old["external_jar_sha256"],
            commands=commands,
            runtime=runtime,
            prior_active_seconds=totals,
            prior_total_active_seconds=sum(totals),
            footprint_bytes=footprint(original),
            host_disk=shutil.disk_usage(original)._asdict(),
            memory_info=Path("/proc/meminfo").read_text(),
            preserved_file_count=len(before),
        ),
    )
    print(
        json.dumps(
            {
                "preflight": "passed",
                "preserved_files": len(before),
                "prior_active_seconds": sum(totals),
            }
        )
    )


def acquire(root, original, repo):
    fields = list(expected(ROWS[0]))
    query = "%20OR%20".join("run_accession%3D%22" + r[2] + "%22" for r in ROWS)
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query="
        + query
        + "&fields="
        + ",".join(fields)
        + "&format=tsv&limit=20"
    )
    fetch(root, url, "inputs/ena_remaining3.tsv")
    with (root / "inputs/ena_remaining3.tsv").open(newline="") as h:
        observed = list(csv.DictReader(h, delimiter="\t"))
    differences = compare_metadata(observed)
    dump(
        root / "records/ena_comparison.json",
        dict(
            registration="17b733c",
            differences=differences,
            accessed_utc=stamp(),
            response_sha256=digest(root / "inputs/ena_remaining3.tsv"),
        ),
    )
    if differences:
        raise ValueError("Archive metadata differs; retained for review")
    validations = []
    for row in ROWS:
        record = expected(row)
        name = row[2] + ".fastq.gz"
        existing = original / "inputs/fastq" / name
        path = existing if existing.exists() else root / "inputs/fastq" / name
        if not path.exists():
            fetch(
                root,
                "https://" + record["fastq_ftp"],
                str(path.relative_to(root)),
                max_seconds=1800,
            )
        valid = check_fastq(path, record)
        valid.update(
            cell_id=f"IMR90_A_{row[0]}",
            gsm=row[1],
            run=row[2],
            experiment=row[3],
            path=str(path),
            reused_existing=path == existing,
            archive_record=record,
        )
        dump(root / "records" / ("fastq-" + row[2] + ".json"), valid)
        validations.append(valid)
    with (root / "inputs/manifest_remaining3.tsv").open("x", newline="") as h:
        writer = csv.writer(h, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "cell_id",
                "platform",
                "read1",
                "read2",
                "bam",
                "library_id",
                "n_input_reads",
                "protocol",
                "strandedness",
                "read_layout",
            ]
        )
        for v in validations:
            writer.writerow(
                [
                    v["cell_id"],
                    "plate",
                    v["path"],
                    "",
                    "",
                    v["experiment"],
                    v["reads"],
                    "ramda",
                    "unstranded",
                    "single",
                ]
            )
    print(
        json.dumps(
            {
                "validated_complete_libraries": len(validations),
                "reads": [v["reads"] for v in validations],
            }
        )
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["preflight", "acquire"])
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--original", type=Path, required=True)
    ap.add_argument("--repo", type=Path, required=True)
    args = ap.parse_args()
    globals()[args.action](
        args.root.resolve(), args.original.resolve(), args.repo.resolve()
    )


if __name__ == "__main__":
    main()
