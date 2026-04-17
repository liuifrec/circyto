from __future__ import annotations

import sys
import shutil
from typing import List
from circyto.cli.doctor_meta import DETECTOR_SPECS, detector_runtime_status, resolve_asset

import typer

DETECTORS = {
    "find-circ3": {
        "required_cmds": ["bowtie2", "samtools"],
        "required_assets": [],
    },
    "ciri-full": {
        "required_cmds": ["bwa", "java"],
        "required_assets": ["CIRI-full-jar"],
    },
}

doctor_app = typer.Typer(
    help="Diagnose circyto environment and external dependencies, including CIRI3 runtime requirements.",
    invoke_without_command=True,
)


def _print(status: str, name: str, msg: str):
    typer.echo(f"[{status:<4}] {name:<10} {msg}")


@doctor_app.callback()
def doctor():
    """
    Run environment diagnostics.
    """
    typer.echo("circyto doctor\n")

    fatal_missing: List[str] = []
    missing_optional: List[str] = []

    # ---- Python -------------------------------------------------------------
    ver = sys.version_info
    py_ver = f"{ver.major}.{ver.minor}.{ver.micro}"
    if (ver.major, ver.minor) >= (3, 10):
        _print("OK", "python", py_ver)
    else:
        _print("MISS", "python", f"{py_ver} (>=3.10 required)")
        fatal_missing.append("python>=3.10")
    found_cmds = {}
    found_assets = {}

    # ---- Required executables ----------------------------------------------
    def check(cmd: str, required_for: str | None = None):
        path = shutil.which(cmd)
        found_cmds[cmd] = path is not None
        if path:
            _print("OK", cmd, path)
        else:
            if required_for:
                _print("MISS", cmd, f"required for {required_for}")
                missing_optional.append(cmd)
            else:
                _print("WARN", cmd, "not found (optional)")


    check("bowtie2", "find-circ3")
    check("samtools", "find-circ3")
    check("bwa", "ciri-full")
    check("java", "ciri-full")
    check("STAR")  # optional

    # ---- Detector assets ----------------------------------------------------
    ciri_full_jar = resolve_asset("CIRI-full-jar")
    found_assets["CIRI-full-jar"] = ciri_full_jar.found

    if ciri_full_jar.found and ciri_full_jar.resolved_path is not None:
        source = f" ({ciri_full_jar.source})" if ciri_full_jar.source else ""
        _print("OK", "CIRI-full", f"{ciri_full_jar.resolved_path}{source}")
    else:
        _print("WARN", "CIRI-full", "jar not found")
        for checked in ciri_full_jar.checked_paths:
            typer.echo(f"       checked: {checked}")

    ciri3 = detector_runtime_status("ciri3")
    ciri3_details = ciri3["details"]
    ciri3_status = ciri3["status"]
    ciri3_name = "CIRI3"
    if ciri3_status == "READY":
        mode = ciri3_details.get("preferred_mode") or "unknown"
        if mode == "direct-jar":
            _print("OK", ciri3_name, f"ready via direct jar: {ciri3_details.get('jar')}")
            if ciri3_details.get("java"):
                typer.echo(f"       java: {ciri3_details['java']}")
        else:
            _print("OK", ciri3_name, "ready via configured template")
            if ciri3_details.get("bin"):
                typer.echo(f"       wrapper: {ciri3_details['bin']}")
    elif ciri3_status == "PARTIAL":
        _print("WARN", ciri3_name, ciri3["reason"])
        if ciri3_details.get("jar"):
            typer.echo(f"       jar: {ciri3_details['jar']}")
        if ciri3_details.get("bin"):
            typer.echo(f"       wrapper: {ciri3_details['bin']}")
        if ciri3_details.get("java"):
            typer.echo(f"       java: {ciri3_details['java']}")
        for key in ("jar", "bin", "java"):
            for checked in ciri3_details.get("checked_paths", {}).get(key, []):
                typer.echo(f"       checked {key}: {checked}")
        for err in ciri3_details.get("template_errors", []):
            typer.echo(f"       template: {err}")
    else:
        _print("WARN", ciri3_name, ciri3["reason"])
        for key in ("home", "jar", "bin", "java"):
            for checked in ciri3_details.get("checked_paths", {}).get(key, []):
                typer.echo(f"       checked {key}: {checked}")

    # ---- Summary ------------------------------------------------------------
    typer.echo("\nDetector readiness:")
    for spec in DETECTOR_SPECS:
        runtime = detector_runtime_status(spec.name)
        status = runtime["status"]
        msg = runtime["reason"]
        typer.echo(f"  {spec.name:<10}: {status}" + (f" {msg}" if msg else ""))


    typer.echo("\nSummary:")
    if fatal_missing:
        typer.echo("  ❌ Environment NOT ready")
        typer.echo("  Missing:")
        for m in sorted(set(fatal_missing)):
            typer.echo(f"   - {m}")
        raise typer.Exit(code=1)
    else:
        if missing_optional:
            typer.echo("  ✅ Baseline environment looks good")
            typer.echo("  Optional detector executables not found:")
            for m in sorted(set(missing_optional)):
                typer.echo(f"   - {m}")
            typer.echo("  Install workflow-specific tools only for the detectors you plan to run.")
        else:
            typer.echo("  ✅ Environment looks good")
        raise typer.Exit(code=0)
