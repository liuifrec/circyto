from __future__ import annotations

import json
import shutil
from typing import Dict, Any

import typer

from circyto.cli.doctor_meta import DETECTOR_SPECS, detector_runtime_status, resolve_asset
from circyto.detectors import build_default_engines, get_detector_capabilities

detectors_app = typer.Typer(
    help="List available detectors, readiness status, and external runtime requirements.",
    invoke_without_command=True,
)


def _have_cmd(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _have_asset(asset: str) -> bool:
    return resolve_asset(asset).found


def _status_for(spec) -> str:
    return detector_runtime_status(spec.name)["status"]


def _runtime_missing_dependencies(spec, runtime: dict[str, Any]) -> list[str]:
    missing = list(runtime.get("details", {}).get("missing_cmds", []))
    missing.extend(runtime.get("details", {}).get("missing_assets", []))
    if spec.name != "ciri3":
        return missing

    details = runtime.get("details", {})
    if details.get("direct_ready") or (
        details.get("template_configured") and not details.get("template_errors")
    ):
        return []
    if details.get("jar") and not details.get("java"):
        missing.append("java")
    elif details.get("jar") and details.get("java") and not details.get("java_version_ok"):
        missing.append(f"java>={details.get('required_java_major')}")
    elif not details.get("jar") and not details.get("template_configured"):
        missing.append("CIRI3 jar or CIRCYTO_CIRI3_CMD_TEMPLATE")
    if details.get("template_errors"):
        missing.append("valid CIRCYTO_CIRI3_CMD_TEMPLATE")
    return sorted(set(str(item) for item in missing))


def _optional_dependencies(spec) -> list[dict[str, Any]]:
    return [
        {
            "name": cmd,
            "available": _have_cmd(cmd),
        }
        for cmd in spec.optional_cmds
    ]


@detectors_app.callback()
def detectors(
    json_out: bool = typer.Option(False, "--json", help="Output machine-readable JSON."),
    hints: bool = typer.Option(False, "--hints", help="Print install hints."),
):
    """
    Print detectors and dependency readiness.
    """
    rows = []
    engines = build_default_engines()
    for spec in DETECTOR_SPECS:
        status = _status_for(spec)
        runtime = detector_runtime_status(spec.name)
        needs = spec.required_cmds + spec.required_assets
        engine = engines.get(spec.name)
        caps = get_detector_capabilities(engine) if engine is not None else None
        missing_dependencies = _runtime_missing_dependencies(spec, runtime)
        optional_dependencies = _optional_dependencies(spec)
        input_modes = []
        if caps is not None:
            if caps.accepts_fastq:
                input_modes.append("fastq")
            if caps.accepts_alignment:
                input_modes.append("alignment")
        asset_details = {}
        for asset in spec.required_assets:
            resolution = resolve_asset(asset)
            asset_details[asset] = {
                "found": resolution.found,
                "path": str(resolution.resolved_path) if resolution.resolved_path else None,
                "source": resolution.source,
                "checked_paths": [str(path) for path in resolution.checked_paths],
            }
        rows.append(
            {
                "name": spec.name,
                "status": status,
                "available": status == "READY",
                "missing_dependencies": missing_dependencies,
                "optional_dependencies": optional_dependencies,
                "input_modes": input_modes,
                "recommended_execution_mode": caps.recommended_execution_mode if caps else None,
                "validation_status": spec.validation_status,
                "notes": [*spec.notes, runtime["reason"]],
                "reason": runtime["reason"],
                "type": spec.det_type,
                "needs": needs,
                "assets": asset_details,
                "runtime": runtime["details"],
                "capabilities": {
                    "accepts_fastq": caps.accepts_fastq if caps else None,
                    "accepts_alignment": caps.accepts_alignment if caps else None,
                    "supports_multisample_alignment": caps.supports_multisample_alignment if caps else None,
                    "requires_unsorted_sam": caps.requires_unsorted_sam if caps else None,
                    "supports_star": caps.supports_star if caps else None,
                    "supports_bwa": caps.supports_bwa if caps else None,
                    "recommended_execution_mode": caps.recommended_execution_mode if caps else None,
                    "max_parallel": caps.max_parallel if caps else None,
                },
                "hints": spec.hint_lines,
            }
        )

    if json_out:
        typer.echo(json.dumps({"detectors": rows}, indent=2))
        raise typer.Exit(code=0)

    # Pretty table (simple, no extra deps)
    typer.echo("NAME         STATUS     TYPE   MODE                NEEDS")
    for r in rows:
        needs_str = ", ".join(r["needs"]) if r["needs"] else "-"
        mode = r["capabilities"]["recommended_execution_mode"] or "-"
        typer.echo(f"{r['name']:<12} {r['status']:<10} {r['type']:<5} {mode:<19} {needs_str}")
        typer.echo(f"  reason: {r['reason']}")
        caps = r["capabilities"]
        typer.echo(
            "  caps: "
            f"fastq={caps['accepts_fastq']} alignment={caps['accepts_alignment']} "
            f"multisample={caps['supports_multisample_alignment']} "
            f"unsorted_sam={caps['requires_unsorted_sam']} "
            f"bwa={caps['supports_bwa']} star={caps['supports_star']} "
            f"max_parallel={caps['max_parallel']}"
        )
        runtime = r["runtime"]
        if r["name"] == "ciri3":
            if runtime.get("jar"):
                typer.echo(f"  jar: {runtime['jar']}")
            if runtime.get("bin"):
                typer.echo(f"  wrapper: {runtime['bin']}")
            if runtime.get("java"):
                typer.echo(f"  java: {runtime['java']}")
            typer.echo(
                f"  java_version: {runtime.get('java_version')} "
                f"(required >= {runtime.get('required_java_major')}, "
                f"recommended {runtime.get('recommended_java_major')})"
            )
            if runtime.get("java_version_error"):
                typer.echo(f"  java_check: {runtime['java_version_error']}")
            typer.echo(f"  execution: {runtime.get('preferred_mode') or 'unconfigured'}")
            for err in runtime.get("template_errors", []):
                typer.echo(f"  template: {err}")
        for asset, details in r["assets"].items():
            if details["path"]:
                typer.echo(f"  asset: {asset} -> {details['path']}")
            elif details["checked_paths"]:
                typer.echo(f"  asset: {asset} -> missing")
                for checked in details["checked_paths"]:
                    typer.echo(f"    checked: {checked}")
        if r["name"] == "ciri3":
            for key in ("home", "jar", "bin", "java"):
                checked = runtime.get("checked_paths", {}).get(key, [])
                if not checked:
                    continue
                typer.echo(f"  checked {key}:")
                for path in checked:
                    typer.echo(f"    {path}")

    if hints:
        typer.echo("\nHints:")
        for r in rows:
            typer.echo(f"  {r['name']}:")
            for line in r["hints"]:
                typer.echo(f"    - {line}")

    raise typer.Exit(code=0)
