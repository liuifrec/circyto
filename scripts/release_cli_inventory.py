#!/usr/bin/env python3
"""Write the release-audit inventory for the installed circyto CLI."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from typer.main import get_command

from circyto import __version__
from circyto.cli.circyto import app


LEGACY = {"prepare", "run", "run-manifest", "collect", "make"}
ALIASES = {"align-manifest": "prepare-alignment-cache", "run-batch": "run-detector"}
EXPERIMENTAL_PREFIXES = (
    "workflow",
    "nanopore",
    "ciri-long",
    "prepare-scomatic-input",
    "run-scomatic",
    "import-scomatic",
    "merge-scomatic",
    "normalize-scomatic-results",
    "export-scomatic-inputs",
    "join-circ-snv-summary",
    "scanpy-",
    "export-biogenesis",
    "import-scrr-rt",
    "merge-scrr-rt",
)
MANUSCRIPT_FACING = {
    "workflow smartseq3-ciri3",
    "workflow full-length-circrna",
    "annotate-circs",
    "repair-host-genes",
    "build-scrr-cell-map",
    "remap-scrr-mudata-obs",
    "import-scrr-cnv",
    "merge-scrr-cnv",
    "import-scrr-rt",
    "merge-scrr-rt",
    "prepare-scomatic-input",
    "run-scomatic",
    "import-scomatic",
    "merge-scomatic",
    "summarize-benchmark",
    "summarize-circ-host-genes",
    "export-mudata",
    "export-biogenesis",
}
PURPOSE_OVERRIDES = {
    "demux smartseq2": "Demultiplex pooled paired-end SMART-Seq2 FASTQs into per-cell files and a manifest.",
    "manifest validate": "Validate a circyto FASTQ manifest.",
    "manifest validate-alignment": "Validate a circyto alignment manifest.",
}


def _external_tools(path: str) -> list[str]:
    if path in {"doctor", "detectors"}:
        return ["workflow-dependent; reported by the command"]
    if path.startswith("nanopore align"):
        return ["minimap2", "samtools"]
    if path.startswith("nanopore inspect-bsj"):
        return ["samtools"]
    if path.startswith("ciri-long ") and not path.endswith(("import", "validate-manifest")):
        return ["CIRI-long", "bwa"]
    if path in {"run-ciri3", "validate-ciri3-template"}:
        return ["java>=12", "CIRI3 jar or wrapper", "bwa or STAR", "samtools"]
    if path.startswith("workflow "):
        return ["java>=12", "CIRI3 jar or wrapper", "bwa or STAR", "samtools"]
    if path in {"prepare-alignment-cache", "align-manifest", "plan-alignment-cache"}:
        return ["bwa or STAR", "samtools"]
    if path == "run-detector-from-alignments":
        return ["detector-dependent", "java>=12 for CIRI3"]
    if path in {"run-detector", "run-batch", "run-multidetector", "smoke", "smoke smartseq2"}:
        return ["detector-dependent; see `circyto detectors`"]
    if path in {"prepare-scomatic-input"}:
        return ["samtools"]
    if path == "run-scomatic":
        return ["external SComatic checkout/environment", "samtools"]
    if path in {"run", "run-manifest", "make"}:
        return ["java", "bwa", "samtools", "CIRI-full assets"]
    return []


def _required_parameters(command: Any) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for parameter in command.params:
        if not parameter.required:
            continue
        kind = "argument" if type(parameter).__name__.endswith("Argument") else "option"
        spelling = parameter.name
        if kind == "option" and parameter.opts:
            spelling = parameter.opts[0]
        rows.append({"name": parameter.name, "cli_spelling": spelling, "kind": kind})
    return rows


def _walk(command: Any, path: list[str], rows: list[dict[str, Any]]) -> None:
    if path:
        command_path = " ".join(path)
        help_text = (command.help or PURPOSE_OVERRIDES.get(command_path, "")).strip()
        lifecycle = "current"
        if command_path in LEGACY:
            lifecycle = "legacy"
        elif command_path in ALIASES:
            lifecycle = "compatibility-alias"
        rows.append(
            {
                "command_path": command_path,
                "one_line_purpose": help_text.splitlines()[0] if help_text else "",
                "required_arguments": _required_parameters(command),
                "optional_external_tools": _external_tools(command_path),
                "experimental": command_path.startswith(EXPERIMENTAL_PREFIXES),
                "manuscript_facing": command_path in MANUSCRIPT_FACING,
                "lifecycle": lifecycle,
                "alias_for": ALIASES.get(command_path),
            }
        )
    for name, child in sorted(getattr(command, "commands", {}).items()):
        _walk(child, [*path, name], rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("docs/cli_inventory.json"))
    args = parser.parse_args()

    commands: list[dict[str, Any]] = []
    _walk(get_command(app), [], commands)
    payload = {
        "schema_version": "circyto.release_cli_inventory.v1",
        "circyto_version": __version__,
        "command_count": len(commands),
        "notes": {
            "experimental": "True denotes an explicitly experimental or exploratory public command family.",
            "manuscript_facing": "True denotes commands used by the repository's manuscript reproducibility path.",
            "optional_external_tools": "Python extras are declared in pyproject.toml; this field focuses on non-core workflow tooling.",
        },
        "commands": commands,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(f"wrote {len(commands)} commands to {args.output}")


if __name__ == "__main__":
    main()
