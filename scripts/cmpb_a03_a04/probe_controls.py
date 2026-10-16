#!/usr/bin/env python3
"""Tiny destructive-to-own-process tests only; no genomic inputs or circyto imports."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from run_stage import dump, footprint


def alive(pid):
    try:
        return Path("/proc", str(pid), "stat").read_text().split()[2] != "Z"
    except FileNotFoundError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--tag", default="probe")
    ap.add_argument("--footprint-root", type=Path)
    args = ap.parse_args()
    root = args.root.resolve()
    specs = root / "specs"
    specs.mkdir(exist_ok=True)
    common = (
        "import json,os,pathlib,subprocess,time; "
        'child=subprocess.Popen(["/bin/sleep","30"]); '
        'pathlib.Path("records/PROBE.pids.json").write_text(json.dumps([os.getpid(),child.pid])); '
    )
    cases = [
        (
            "inheritance",
            'import json,os,subprocess; print(json.dumps({"parent_affinity":sorted(os.sched_getaffinity(0))})); subprocess.run(["/usr/bin/python3","-c","import os; print(sorted(os.sched_getaffinity(0)))"],check=True)',
            0,
            10,
            128 * 1024**2,
        ),
        ("exit-failure", "raise SystemExit(7)", 7, 10, 128 * 1024**2),
        ("wall-stop", common + "time.sleep(30)", 1, 2, 128 * 1024**2),
        (
            "memory-stop",
            common + "x=bytearray(256*1024**2); time.sleep(30)",
            1,
            10,
            64 * 1024**2,
        ),
        (
            "disk-stop",
            common
            + 'pathlib.Path("tmp/PROBE-disk.bin").write_bytes(b"x"*(4*1024**2)); time.sleep(30)',
            125,
            10,
            128 * 1024**2,
        ),
        (
            "affinity-stop",
            common + "os.sched_setaffinity(0,{8}); time.sleep(30)",
            125,
            10,
            128 * 1024**2,
        ),
    ]
    results = []
    for name, code, expected, wall, memory in cases:
        stage = args.tag + "-" + name
        spec = {
            "stage": stage,
            "category": "probe",
            "wall_seconds": wall,
            "memory_bytes": memory,
            "argv": ["/usr/bin/python3", "-c", code.replace("PROBE", stage)],
        }
        if args.footprint_root:
            spec["footprint_root"] = str(args.footprint_root.resolve(strict=True))
        if name == "disk-stop":
            spec["disk_stop_bytes"] = footprint(args.footprint_root or root) + 2 * 1024**2
        path = specs / (stage + ".request.json")
        dump(path, spec)
        result = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).with_name("run_stage.py")),
                "launch",
                "--root",
                str(root),
                "--spec",
                str(path),
            ],
            check=False,
        )
        record = json.loads((root / "records" / (stage + ".json")).read_text())
        guard = json.loads((root / "records" / (stage + ".guard.json")).read_text())
        pid_file = root / "records" / (stage + ".pids.json")
        pids = json.loads(pid_file.read_text()) if pid_file.exists() else []
        time.sleep(0.1)
        survivors = [pid for pid in pids if alive(pid)]
        passed = not survivors and (
            result.returncode == expected if expected != 1 else result.returncode != 0
        )
        if name == "inheritance":
            text = (root / record["log"]).read_text()
            passed = passed and text.count("[0, 1]") == 2
        if name == "memory-stop":
            passed = passed and "Result=oom-kill" in record["systemd_properties"]
        if name == "wall-stop":
            passed = passed and "Result=timeout" in record["systemd_properties"]
        if name == "disk-stop":
            passed = passed and guard.get("stop_reason") == "automatic disk threshold"
        if name == "affinity-stop":
            passed = passed and guard.get("stop_reason") == "CPU affinity widened"
        results.append(
            {
                "probe": name,
                "passed": passed,
                "exit_code": result.returncode,
                "surviving_test_pids": survivors,
                "record": record,
                "guard": guard,
            }
        )
        dump(root / "records" / (args.tag + "_control_tests.json"), results)
        if not passed:
            raise RuntimeError("Resource probe failed: " + name)
    print(
        json.dumps(
            {
                "control_tests_passed": len(results),
                "test_memory_mib": [64, 128],
                "retained_test_payload_bytes": 4 * 1024**2,
            }
        )
    )


if __name__ == "__main__":
    main()
