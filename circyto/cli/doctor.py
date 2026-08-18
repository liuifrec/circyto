from __future__ import annotations

import shutil
import sys
from importlib import metadata
from typing import Iterable

import typer

from circyto.cli.doctor_meta import DETECTOR_SPECS, detector_runtime_status, resolve_asset
from circyto.detectors.ciri3 import inspect_ciri3_runtime
from circyto.pipeline.ciri_long_adapter import check_ciri_long_readiness


doctor_app = typer.Typer(
    help="Diagnose core package readiness and optional workflow dependencies without network access.",
    invoke_without_command=True,
)


CORE_PACKAGES = (
    ("circyto", "circyto"),
    ("typer", "typer"),
    ("rich", "rich"),
    ("pandas", "pandas"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("anndata", "anndata"),
)

OPTIONAL_TOOLS = (
    ("find-circ3", "the find-circ3 detector"),
    ("bowtie2", "find-circ3 alignment"),
    ("samtools", "alignment handling and Nanopore interoperability"),
    ("bwa", "CIRI2/CIRI-full, CIRI3 BWA mode, and CIRI-long"),
    ("perl", "the bundled CIRI2 implementation"),
    ("java", "CIRI-full and direct CIRI3 execution"),
    ("STAR", "CIRI3 STAR-mode workflows"),
    ("minimap2", "generic Nanopore cDNA alignment interoperability"),
    ("CIRI-long", "chemistry-gated RCRT/circRNA-enriched long-read workflows"),
)


def _print(status: str, name: str, message: str) -> None:
    typer.echo(f"[{status:<4}] {name:<12} {message}")


def _section(title: str) -> None:
    typer.echo(f"\n{title}:")


def _package_readiness() -> list[str]:
    missing: list[str] = []
    version = sys.version_info
    python_version = f"{version.major}.{version.minor}.{version.micro}"
    if (version.major, version.minor) >= (3, 10):
        _print("OK", "python", f"{python_version} (required >=3.10)")
    else:
        _print("MISS", "python", f"{python_version}; install Python >=3.10")
        missing.append("python>=3.10")

    for distribution, label in CORE_PACKAGES:
        try:
            package_version = metadata.version(distribution)
        except metadata.PackageNotFoundError:
            _print("MISS", label, f"Python package not installed; reinstall circyto to restore `{distribution}`")
            missing.append(distribution)
        else:
            _print("OK", label, package_version)
    return missing


def _asset_line(label: str, *, purpose: str) -> None:
    asset = resolve_asset(label)
    if asset.found and asset.resolved_path is not None:
        source = f" ({asset.source})" if asset.source else ""
        _print("OK", label, f"{asset.resolved_path}{source}")
        return
    _print("OPT", label, f"not found; needed for {purpose}")
    for checked in asset.checked_paths:
        typer.echo(f"       checked: {checked}")


def _external_tools() -> list[str]:
    missing: list[str] = []
    for executable, purpose in OPTIONAL_TOOLS:
        path = shutil.which(executable)
        if path:
            _print("OK", executable, path)
        else:
            _print(
                "OPT",
                executable,
                f"not found; install `{executable}` and add it to PATH to use {purpose}",
            )
            missing.append(executable)
    return missing


def _print_ciri3_details() -> None:
    runtime = inspect_ciri3_runtime()
    jar = runtime.get("jar")
    if jar:
        _print("OK", "CIRI3-jar", str(jar))
    else:
        _print(
            "OPT",
            "CIRI3-jar",
            "not found; use the packaged asset or set CIRCYTO_CIRI3_JAR/CIRCYTO_CIRI3_HOME",
        )
    java = runtime.get("java")
    if java:
        version = runtime.get("java_version")
        required = runtime.get("required_java_major")
        recommended = runtime.get("recommended_java_major")
        typer.echo(f"       java: {java}")
        typer.echo(
            f"       java_version: {version} "
            f"(required >= {required}, recommended {recommended})"
        )
        if runtime.get("java_version_ok"):
            _print("OK", "CIRI3-Java", f"{version} (required >={required}; recommended {recommended})")
        else:
            _print("OPT", "CIRI3-Java", f"{version} is too old; install Java >={required} (Java {recommended} recommended)")
    elif jar:
        _print("OPT", "CIRI3-Java", "not found; install Java >=12 (Java 17 recommended) for the CIRI3 jar")


def _detector_readiness() -> None:
    for spec in DETECTOR_SPECS:
        runtime = detector_runtime_status(spec.name)
        typer.echo(f"  {spec.name:<10}: {runtime['status']}" + (f" — {runtime['reason']}" if runtime["reason"] else ""))

    ciri_long = check_ciri_long_readiness()
    if ciri_long["ok"]:
        tools = ciri_long["tools"]
        version = tools["ciri_long"]["version"] if tools["ciri_long"] else "unknown"
        typer.echo(f"  {'ciri-long':<10}: READY — CIRI-long {version}; supply a reference to check BWA indexes")
    else:
        reason = "; ".join(str(item) for item in ciri_long["errors"])
        typer.echo(f"  {'ciri-long':<10}: NOT READY (optional) — {reason}")


def _protocol_limitations(lines: Iterable[str]) -> None:
    for line in lines:
        typer.echo(f"  - {line}")


@doctor_app.callback()
def doctor() -> None:
    """Run local, network-free environment diagnostics."""
    typer.echo("circyto doctor")

    _section("Required core Python package readiness")
    fatal_missing = _package_readiness()

    _section("Bundled/configured detector assets")
    _asset_line("CIRI2-adapter", purpose="CIRI2 and single-end CIRI-full fallback")
    _asset_line("CIRI-full-jar", purpose="paired-end CIRI-full")
    _asset_line("CIRI-full-adapter", purpose="CIRI-full orchestration")
    _print_ciri3_details()

    _section("Optional external executables")
    missing_optional = _external_tools()

    _section("Detector-specific readiness")
    _detector_readiness()

    _section("Protocol-specific limitations")
    _protocol_limitations(
        (
            "CIRI3, CIRI2, CIRI-full, and find-circ3 are short-read detector paths; generic Nanopore BAMs are unsupported inputs.",
            "`circyto nanopore` provides experimental cDNA alignment/QC/provenance interoperability and does not validate or call circRNAs.",
            "`circyto ciri-long` is optional and accepts only chemistry-confirmed RCRT/circRNA-enriched libraries, not ordinary ONT cDNA.",
            "Single-end `ciri-full` uses the bundled CIRI2 fallback and is not CIRI-full full-length reconstruction.",
        )
    )

    _section("Summary")
    if fatal_missing:
        typer.echo("  ❌ Core Python environment is NOT ready")
        for package in sorted(set(fatal_missing)):
            typer.echo(f"   - missing required package: {package}")
        raise typer.Exit(code=1)

    typer.echo("  ✅ Core Python package is ready")
    if missing_optional:
        typer.echo("  Optional workflow tools are missing; this does not indicate a circyto package failure.")
        typer.echo("  Install only the tools required by the workflow you plan to run.")
    else:
        typer.echo("  All recognized optional workflow executables were found.")
    typer.echo("  No network checks were performed.")
    raise typer.Exit(code=0)
