#!/usr/bin/env python3
"""A02-only, fail-closed user-service runner; never imports circyto.

Every launch has a fresh service and output record. The six-file registration is
read-only. No shell command evaluation, resume, cleanup, or automatic next stage.
"""

from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

GIB = 1024**3
CPUS = set(range(8))
SAFE_ENV = {
    "PATH",
    "PYTHONNOUSERSITE",
    "TMPDIR",
    "XDG_CACHE_HOME",
    "PIP_CACHE_DIR",
    "MAMBA_ROOT_PREFIX",
    "CONDA_PKGS_DIRS",
    "MAMBA_PKGS_DIRS",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "CIRCYTO_CIRI3_JAR",
    "CIRCYTO_CIRI3_JAVA",
    "CIRCYTO_CIRI3_EXTRA_ARGS",
    "JAVA_TOOL_OPTIONS",
    "A02_ROOT",
    "A02_REAL_BWA",
    "A02_REAL_JAVA",
}


def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    tmp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def footprint(root):
    """Conservative sum of allocated-or-logical bytes, counting retained copies.

    Do not follow symlinks out of the evaluation root. Acquired data/env caches
    must live below this root; shared external system executables are identified
    separately, not claimed as new evaluation storage.
    """
    total = 0
    for directory, _, files in os.walk(root, followlinks=False):
        for name in files:
            try:
                st = (Path(directory) / name).lstat()
                total += max(st.st_size, st.st_blocks * 512)
            except FileNotFoundError:  # package installers can rename live files
                continue
    return total


def cgroup_path():
    lines = Path("/proc/self/cgroup").read_text().splitlines()
    line = next(line for line in lines if line.startswith("0::"))
    return Path("/sys/fs/cgroup") / line[3:].lstrip("/")


def affinity_violations(group):
    bad = []
    for procs in [group / "cgroup.procs", *group.glob("**/cgroup.procs")]:
        for pid in set(procs.read_text().split()):
            for task in Path("/proc", pid, "task").glob("*"):
                try:
                    mask = os.sched_getaffinity(int(task.name))
                    if not mask <= CPUS:
                        bad.append(
                            {
                                "pid": int(pid),
                                "tid": int(task.name),
                                "mask": sorted(mask),
                            }
                        )
                except ProcessLookupError:
                    continue
    return bad


def process_tree_rss(group):
    pids = set()
    for procs in [group / "cgroup.procs", *group.glob("**/cgroup.procs")]:
        pids.update(procs.read_text().split())
    rss = 0
    for pid in pids:
        try:
            for line in Path("/proc", pid, "status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    rss += int(line.split()[1]) * 1024
        except FileNotFoundError:
            continue
    return rss


def initialize(root):
    for name in [
        "records",
        "logs",
        "inputs",
        "ref",
        "tmp",
        "cache",
        "helpers",
        "specs",
    ]:
        (root / name).mkdir(parents=True, exist_ok=True)
    budget = root / "records/budget.json"
    if not budget.exists():
        now = time.time()
        dump(
            budget,
            {
                "started_utc": stamp(),
                "started_epoch": now,
                "deadline_epoch": now + 16 * 3600,
                "wall_cap_seconds": 16 * 3600,
                "disk_cap_bytes": 80 * GIB,
                "disk_stop_bytes": 79 * GIB,
                "cpus": sorted(CPUS),
                "shared_memory_max_bytes": 24 * GIB,
                "first_cell_memory_max_bytes": 20 * GIB,
                "first_cell_wall_cap_seconds": 4 * 3600,
            },
        )
    return json.loads(budget.read_text())


def launch(root, spec_path):
    budget = initialize(root)
    # One managed stage at a time: limits are not parallel allowances.
    lock = (root / "records/stage.lock").open("a")
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    spec = json.loads(spec_path.read_text())
    stage = spec["stage"]
    if not stage or any(
        c not in "abcdefghijklmnopqrstuvwxyz0123456789-_" for c in stage
    ):
        raise ValueError("stage must be a simple unique identifier")
    if not isinstance(spec["argv"], list) or not all(
        isinstance(x, str) for x in spec["argv"]
    ):
        raise ValueError("argv must be a string list")
    if set(spec.get("env", {})) - SAFE_ENV:
        raise ValueError("non-allowlisted environment key")
    record = root / "records" / (stage + ".json")
    if record.exists():
        raise FileExistsError("Stage already recorded; no automatic rerun: " + stage)
    remaining = int(budget["deadline_epoch"] - time.time())
    memory = 24 * GIB
    if spec["category"] == "firstcell":
        first = root / "records/firstcell_budget.json"
        if not first.exists():
            dump(
                first,
                {"started_utc": stamp(), "deadline_epoch": time.time() + 4 * 3600},
            )
        remaining = min(
            remaining,
            int(json.loads(first.read_text())["deadline_epoch"] - time.time()),
        )
        memory = 20 * GIB
    if spec["category"] == "probe":
        remaining = min(remaining, int(spec.get("wall_seconds", 10)))
        memory = min(memory, int(spec.get("memory_bytes", 128 * 1024**2)))
    if remaining <= 0 or footprint(root) >= budget["disk_stop_bytes"]:
        raise RuntimeError("Aggregate time/disk budget exhausted before launch")
    unit = "circyto-a02-" + stage
    copied = root / "specs" / (stage + ".json")
    spec.update(
        {
            "root": str(root),
            "memory_max_bytes": memory,
            "wall_seconds": remaining,
            "unit": unit,
            "disk_stop_bytes": min(
                int(spec.get("disk_stop_bytes", 79 * GIB)), 79 * GIB
            ),
        }
    )
    dump(copied, spec)
    helper = root / "helpers" / ("run_stage-" + digest(__file__)[:16] + ".py")
    if not helper.exists():
        shutil.copyfile(__file__, helper)
    argv = [
        "systemd-run",
        "--user",
        "--unit=" + unit,
        "--wait",
        "--pipe",
        "--property=MemoryMax=" + str(memory),
        "--property=MemorySwapMax=0",
        "--property=RuntimeMaxSec=" + str(remaining) + "s",
        "--property=TimeoutStopSec=1s",
        "--property=KillMode=control-group",
        "--property=OOMPolicy=kill",
        "--property=CPUAffinity=0-7",
        "--property=CPUQuota=800%",
        "--property=TasksMax=512",
        "/usr/bin/python3",
        str(helper),
        "guard",
        "--root",
        str(root),
        "--spec",
        str(copied),
    ]
    rec = {
        "schema": "circyto.cmpb.a02.stage.v1",
        "stage": stage,
        "category": spec["category"],
        "started_utc": stamp(),
        "status": "launching",
        "argv": spec["argv"],
        "launcher_argv": argv,
        "environment": spec.get("env", {}),
        "helper_sha256": digest(helper),
        "spec_sha256": digest(copied),
        "log": "logs/" + stage + ".log",
    }
    snapshots = root / "helpers/source_snapshots"
    snapshots.mkdir(exist_ok=True)
    rec["invoked_python_source_sha256"] = {}
    for argument in spec["argv"]:
        candidate = Path(argument)
        if candidate.suffix == ".py" and candidate.is_file():
            for source in candidate.parent.glob("*.py"):
                sha = digest(source)
                target = snapshots / (sha + "-" + source.name)
                if not target.exists():
                    shutil.copyfile(source, target)
                rec["invoked_python_source_sha256"][str(source)] = sha
    dump(record, rec)
    started = time.monotonic()
    with (root / rec["log"]).open("x") as log:
        result = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=False)
    properties = subprocess.run(
        [
            "systemctl",
            "--user",
            "show",
            unit,
            "--property=Result",
            "--property=ExecMainStatus",
            "--property=MemoryCurrent",
            "--property=CPUUsageNSec",
            "--property=ActiveState",
            "--property=ControlGroup",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    rec.update(
        {
            "finished_utc": stamp(),
            "elapsed_seconds": time.monotonic() - started,
            "launcher_exit_code": result.returncode,
            "status": "succeeded" if result.returncode == 0 else "failed",
            "systemd_properties": properties.stdout.strip().splitlines(),
            "systemd_query_exit": properties.returncode,
            "log_sha256": digest(root / rec["log"]),
            "footprint_bytes_after": footprint(root),
        }
    )
    dump(record, rec)
    print(
        json.dumps(
            {
                k: rec[k]
                for k in [
                    "stage",
                    "status",
                    "launcher_exit_code",
                    "elapsed_seconds",
                    "log",
                ]
            }
        )
    )
    lock.close()
    return result.returncode


def guard(root, spec_path):
    spec = json.loads(spec_path.read_text())
    group = cgroup_path()
    controls = {
        n: (group / n).read_text().strip() if (group / n).exists() else None
        for n in ["memory.max", "memory.swap.max", "cpu.max", "cgroup.controllers"]
    }
    path = root / "records" / (spec["stage"] + ".guard.json")
    state = {
        "started_utc": stamp(),
        "cgroup": str(group),
        "controls": controls,
        "affinity": sorted(os.sched_getaffinity(0)),
        "status": "preflight",
        "memory_current_peak_bytes": 0,
        "process_tree_rss_peak_bytes": 0,
        "disk_peak_bytes": footprint(root),
        "poll_seconds": 0.25,
        "disk_poll_seconds": 1,
        "cpu_enforcement": "inherited affinity 0-7 plus all-thread affinity watchdog; CPU quota only if cpu.max exists",
    }
    dump(path, state)
    if (
        controls["memory.max"] != str(spec["memory_max_bytes"])
        or controls["memory.swap.max"] != "0"
    ):
        raise RuntimeError("Cgroup memory/swap limits not applied")
    if set(state["affinity"]) != CPUS:
        raise RuntimeError("Eight-CPU affinity not applied")
    env = os.environ.copy()
    for key in ["PYTHONPATH", "CIRCYTO_CIRI3_CMD_TEMPLATE", "CIRCYTO_CIRI3_BIN"]:
        env.pop(key, None)
    env.update(
        {
            "PYTHONNOUSERSITE": "1",
            "TMPDIR": str(root / "tmp"),
            "XDG_CACHE_HOME": str(root / "cache"),
            "PIP_CACHE_DIR": str(root / "cache/pip"),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    env.update(spec.get("env", {}))
    child = subprocess.Popen(spec["argv"], cwd=root, env=env)
    state.update({"status": "running", "child_pid": child.pid})
    dump(path, state)
    last_disk = 0.0
    try:
        while child.poll() is None:
            current = int((group / "memory.current").read_text())
            state["memory_current_peak_bytes"] = max(
                state["memory_current_peak_bytes"], current
            )
            rss = process_tree_rss(group)
            state["process_tree_rss_peak_bytes"] = max(
                state["process_tree_rss_peak_bytes"], rss
            )
            if rss >= spec["memory_max_bytes"]:
                state.update(status="stopped", stop_reason="process-tree RSS threshold")
                dump(path, state)
                return 125
            if affinity_violations(group):
                state.update(status="stopped", stop_reason="CPU affinity widened")
                dump(path, state)
                return 125  # systemd KillMode kills every descendant on main exit
            if time.monotonic() - last_disk >= 1:
                used = footprint(root)
                state["disk_peak_bytes"] = max(state["disk_peak_bytes"], used)
                if used >= spec["disk_stop_bytes"]:
                    state.update(
                        status="stopped", stop_reason="automatic disk threshold"
                    )
                    dump(path, state)
                    return 125
                last_disk = time.monotonic()
                state["updated_utc"] = stamp()
                dump(path, state)
            time.sleep(0.25)
        state.update(
            status="finished", child_exit_code=child.returncode, finished_utc=stamp()
        )
        dump(path, state)
        return child.returncode if child.returncode >= 0 else 128 - child.returncode
    except BaseException as exc:
        state.update(status="guard_error", error=repr(exc), finished_utc=stamp())
        dump(path, state)
        raise  # main exit invokes systemd's process-tree cleanup


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["launch", "guard", "init"])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--spec", type=Path)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    if args.action == "init":
        print(json.dumps(initialize(root)))
        return 0
    return (launch if args.action == "launch" else guard)(
        root, args.spec.resolve(strict=True)
    )


if __name__ == "__main__":
    sys.exit(main())
