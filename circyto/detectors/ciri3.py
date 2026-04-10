# circyto/detectors/ciri3.py
from __future__ import annotations

import csv
import os
import shlex
import string
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from .base import DetectorBase, DetectorCapabilities, DetectorRunInputs, DetectorResult
from circyto.paths import Ciri3Resolution, resolve_ciri3_installation

REQUIRED_TEMPLATE_FIELDS = {
    "alignment",
    "alignment_format",
    "cell_id",
    "threads",
    "raw_output",
    "outdir",
}
OPTIONAL_TEMPLATE_FIELDS = {
    "ref_fa",
    "gtf",
    "extra_args",
    "alignment_input",
    "mapper_mode",
    "sortedness",
    "artifact_bucket",
    "chimeric_junction",
    "unmapped_mate1",
    "unmapped_mate2",
    "bwa_sam",
    "read_layout",
    "group_id",
    "log_path",
}
ALLOWED_TEMPLATE_FIELDS = REQUIRED_TEMPLATE_FIELDS | OPTIONAL_TEMPLATE_FIELDS


def _shell_join(args: list[str]) -> str:
    return shlex.join(args)


def _normalize_ciri3_output(raw_path: Path, out_path: Path) -> None:
    if not raw_path.exists():
        raise FileNotFoundError(f"CIRI3 raw output not found: {raw_path}")

    with raw_path.open("r", encoding="utf-8", newline="") as src, out_path.open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        dst.write("circ_id\tchr\tstart\tend\tstrand\tsupport\n")
        reader = csv.DictReader(src, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"CIRI3 raw output has no header row: {raw_path}")
        for row in reader:
            circ_id = (
                row.get("circ_id")
                or row.get("circRNA_ID")
                or row.get("circRNA_id")
                or row.get("junction")
                or ""
            ).strip()
            chrom = (row.get("chr") or row.get("chrom") or row.get("contig") or "").strip()
            start = (row.get("start") or row.get("circRNA_start") or row.get("donor") or "").strip()
            end = (row.get("end") or row.get("circRNA_end") or row.get("acceptor") or "").strip()
            strand = (row.get("strand") or row.get("strandness") or ".").strip() or "."
            support = (
                row.get("support")
                or row.get("junction_reads")
                or row.get("#junction_reads")
                or row.get("bsj_reads")
                or row.get("count")
                or "0"
            ).strip()

            if not circ_id and chrom and start and end:
                circ_id = f"{chrom}:{start}|{end}|{strand}"
            if not circ_id or not chrom or not start or not end:
                continue
            dst.write(f"{circ_id}\t{chrom}\t{start}\t{end}\t{strand}\t{support}\n")


def _template_fields(template: str) -> set[str]:
    fields: set[str] = set()
    formatter = string.Formatter()
    for _, field_name, _, _ in formatter.parse(template):
        if field_name:
            fields.add(field_name)
    return fields


def validate_ciri3_template(template: str) -> tuple[bool, list[str], dict[str, Any]]:
    errors: list[str] = []
    fields = _template_fields(template)
    missing = sorted(REQUIRED_TEMPLATE_FIELDS - fields)
    unknown = sorted(fields - ALLOWED_TEMPLATE_FIELDS)
    if missing:
        errors.append(f"Missing required placeholders: {', '.join(missing)}")
    if unknown:
        errors.append(f"Unknown placeholders: {', '.join(unknown)}")
    summary = {
        "fields": sorted(fields),
        "required_fields": sorted(REQUIRED_TEMPLATE_FIELDS),
        "optional_fields": sorted(OPTIONAL_TEMPLATE_FIELDS),
    }
    return (len(errors) == 0), errors, summary


def _tail_text(path: Path, max_lines: int = 40) -> str:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-max_lines:])


def _count_ciri3_rows(path: Path) -> int:
    if not path.exists():
        return 0
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return len([line for line in lines[1:] if line.strip()])


def _option_value(args: list[str], *names: str) -> str | None:
    for i, token in enumerate(args):
        if token in names and i + 1 < len(args):
            return args[i + 1]
        for name in names:
            prefix = f"{name}="
            if token.startswith(prefix):
                return token[len(prefix) :]
    return None


def _has_option(args: list[str], *names: str) -> bool:
    return _option_value(args, *names) is not None or any(token in names for token in args)


def _resolve_ciri3_cli_args(
    *,
    extra_args: str,
    mapper_mode: str,
) -> list[str]:
    tokens = shlex.split(extra_args.strip()) if extra_args.strip() else []
    if not _has_option(tokens, "-S", "--stringency"):
        tokens.extend(["-S", "0"])
    if not _has_option(tokens, "-Ma", "--mapper"):
        tokens.extend(["-Ma", mapper_mode])
    return tokens


def _infer_mapper_mode(
    *,
    aligner: str | None,
    row_meta: dict[str, Any],
) -> str:
    explicit = str(row_meta.get("mapper_mode", "")).strip()
    if explicit in {"0", "1"}:
        return explicit
    aligner_value = str(aligner or row_meta.get("aligner", "")).strip().lower()
    if aligner_value == "star":
        return "1"
    return "0"


def _build_ciri3_alignment_input(
    *,
    alignment: Path,
    row_meta: dict[str, Any],
    mapper_mode: str,
) -> tuple[str, dict[str, Any]]:
    sortedness = str(row_meta.get("sortedness", "")).strip() or "unknown"
    if mapper_mode == "1":
        if alignment.suffix.lower() != ".sam":
            raise RuntimeError(
                f"CIRI3 STAR mode requires STAR-generated aligned SAM input; got {alignment.suffix or 'unknown'} for {alignment}"
            )
        chimeric = str(row_meta.get("chimeric_junction", "")).strip()
        if not chimeric:
            raise RuntimeError("CIRI3 STAR mode requires chimeric_junction in the alignment manifest row.")
        bwa_sam = str(row_meta.get("bwa_sam", "")).strip()
        if not bwa_sam:
            raise RuntimeError("CIRI3 STAR mode requires bwa_sam in the alignment manifest row.")
        input_parts = [chimeric, str(alignment)]
        input_kind = "star-hybrid"
        input_parts.append(bwa_sam)
        return (
            ",".join(input_parts),
            {
                "input_file_type": "STAR tuple",
                "input_sortedness": sortedness,
                "mapper_mode": mapper_mode,
                "chimeric_junction": chimeric,
                "bwa_sam": bwa_sam or None,
                "ciri3_input_kind": input_kind,
            },
        )

    if alignment.suffix.lower() != ".sam":
        raise RuntimeError(
            f"CIRI3 direct execution requires unsorted SAM input; got {alignment.suffix or 'unknown'} for {alignment}"
        )
    return (
        str(alignment),
        {
            "input_file_type": "SAM",
            "input_sortedness": sortedness,
            "mapper_mode": mapper_mode,
            "ciri3_input_kind": "bwa-sam",
        },
    )


def inspect_ciri3_runtime(
    *,
    command_template: str | None = None,
) -> dict[str, Any]:
    runtime_template = command_template or os.environ.get("CIRCYTO_CIRI3_CMD_TEMPLATE")
    resolution = resolve_ciri3_installation()
    details: dict[str, Any] = {
        "template_configured": runtime_template is not None,
        "template_source": "argument" if command_template is not None else ("env:CIRCYTO_CIRI3_CMD_TEMPLATE" if runtime_template else None),
        "home": str(resolution.home.resolved_path) if resolution.home.resolved_path else None,
        "home_source": resolution.home.source,
        "jar": str(resolution.jar.resolved_path) if resolution.jar.resolved_path else None,
        "jar_source": resolution.jar.source,
        "bin": str(resolution.bin.resolved_path) if resolution.bin.resolved_path else None,
        "bin_source": resolution.bin.source,
        "java": str(resolution.java.resolved_path) if resolution.java.resolved_path else None,
        "java_source": resolution.java.source,
        "checked_paths": {
            "home": [str(path) for path in resolution.home.checked_paths],
            "jar": [str(path) for path in resolution.jar.checked_paths],
            "bin": [str(path) for path in resolution.bin.checked_paths],
            "java": [str(path) for path in resolution.java.checked_paths],
        },
    }

    direct_ready = resolution.jar.resolved_path is not None and resolution.java.resolved_path is not None
    wrapper_ready = resolution.bin.resolved_path is not None and runtime_template is not None
    template_valid = False
    template_errors: list[str] = []
    if runtime_template is not None:
        template_valid, template_errors, template_summary = validate_ciri3_template(runtime_template)
        details["template_summary"] = template_summary
    else:
        details["template_summary"] = None

    mode = None
    if direct_ready:
        mode = "direct-jar"
    elif runtime_template is not None and template_valid:
        mode = "template"
    elif resolution.bin.resolved_path is not None:
        mode = "wrapper-partial"
    elif resolution.jar.resolved_path is not None:
        mode = "jar-partial"

    details["direct_ready"] = direct_ready
    details["wrapper_ready"] = wrapper_ready
    details["preferred_mode"] = mode
    details["template_errors"] = template_errors
    return details


def _runtime_errors(details: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    if details["direct_ready"]:
        return errors

    if details["template_configured"]:
        errors.extend(details.get("template_errors", []))
        if not details.get("template_errors"):
            return errors

    if details.get("jar") and not details.get("java"):
        errors.append("CIRI3 jar found, but Java is missing.")
    elif details.get("jar"):
        errors.append("CIRI3 jar found, but direct execution is not usable.")
    elif details.get("bin") and not details["template_configured"]:
        errors.append("CIRI3 wrapper found, but no CIRCYTO_CIRI3_CMD_TEMPLATE is configured.")
    elif not details.get("jar") and not details.get("bin"):
        errors.append("No CIRI3 jar or wrapper detected.")

    if not details["template_configured"] and not details["direct_ready"]:
        errors.append("No direct jar contract available and no CIRCYTO_CIRI3_CMD_TEMPLATE configured.")
    return errors


def _build_direct_ciri3_command(
    *,
    resolution: Ciri3Resolution,
    alignment_input: str,
    raw_output: Path,
    ref_fa: Path,
    gtf: Path | None,
    threads: int,
    cli_args: list[str],
    internal_log: Path,
) -> list[str]:
    if resolution.java.resolved_path is None or resolution.jar.resolved_path is None:
        raise RuntimeError("Direct CIRI3 execution requires both a java executable and a CIRI3 jar.")

    cmd = [
        str(resolution.java.resolved_path),
        "-jar",
        str(resolution.jar.resolved_path),
        "-I",
        alignment_input,
        "-O",
        str(raw_output),
        "-F",
        str(ref_fa),
        "-T",
        str(max(1, threads)),
        "-G",
        str(internal_log),
    ]
    if gtf is not None:
        cmd.extend(["-A", str(gtf)])
    if cli_args:
        cmd.extend(cli_args)
    return cmd


@dataclass
class Ciri3Detector(DetectorBase):
    """
    Alignment-native CIRI3 backend.

    This backend is intentionally strict about its execution contract:
    - it consumes BAM/SAM via DetectorRunInputs.input_mode="alignment"
    - it can execute directly via java -jar when the local runtime contract is ready
    - it can also run through CIRCYTO_CIRI3_CMD_TEMPLATE or the instance command_template

    Template placeholders:
      {alignment}, {alignment_format}, {cell_id}, {threads}, {ref_fa}, {gtf},
      {raw_output}, {outdir}, {extra_args}
    """

    name: str = "ciri3"
    input_type: str = "alignment"
    supports_paired_end: bool = True
    max_parallel: int = 8
    capabilities: DetectorCapabilities = DetectorCapabilities(
        accepts_fastq=False,
        accepts_alignment=True,
        prefers_paired=True,
        supports_single_end=True,
        supports_multisample_alignment=True,
        requires_unsorted_sam=True,
        supports_star=True,
        supports_bwa=True,
        max_parallel=8,
        recommended_execution_mode="alignment-first",
    )

    command_template: str | None = None
    extra_args: str = ""

    def is_available(self) -> bool:
        details = inspect_ciri3_runtime(command_template=self.command_template)
        return bool(details["direct_ready"] or (details["template_configured"] and not details.get("template_errors")))

    def version(self) -> str | None:
        details = inspect_ciri3_runtime(command_template=self.command_template)
        jar = details.get("jar")
        java = details.get("java")
        if not jar or not java:
            return "CIRI3 (jar/java not found)"
        for args in (["-H"], ["DE_BSJ", "-H"]):
            try:
                result = subprocess.run(
                    [java, "-jar", jar, *args],
                    check=False,
                    capture_output=True,
                    text=True,
                )
            except OSError:
                break
            text = (result.stdout or result.stderr or "").strip()
            if text:
                return text.splitlines()[0][:200]
        return "CIRI3 (version unknown)"

    def expected_inputs(self) -> dict[str, str]:
        return {
            "input_mode": "alignment",
            "accepted_extensions": ".sam",
            "supports_multisample_alignment": "true",
            "requires_unsorted_sam": "true",
            "supports_bwa": "true",
            "supports_star": "true",
            "normalization_schema": "circ_id chr start end strand support",
            "required_template_placeholders": ",".join(sorted(REQUIRED_TEMPLATE_FIELDS)),
        }

    def resolve_extra_args(self) -> str:
        if self.extra_args.strip():
            return self.extra_args
        return os.environ.get("CIRCYTO_CIRI3_EXTRA_ARGS", "-S 0")

    def resolve_command_template(self) -> str | None:
        return self.command_template or os.environ.get("CIRCYTO_CIRI3_CMD_TEMPLATE")

    def validate_runtime(self, *, template: str | None = None) -> tuple[bool, list[str], dict[str, Any]]:
        details = inspect_ciri3_runtime(command_template=template or self.command_template)
        errors = _runtime_errors(details)
        if details["direct_ready"]:
            return True, [], details
        if details["template_configured"] and not details.get("template_errors"):
            return True, [], details
        return False, errors, details

    def build_template_context(self, inputs: DetectorRunInputs, *, raw_output: Path, run_dir: Path, log_path: Path) -> dict[str, Any]:
        alignment = inputs.alignment_path
        assert alignment is not None
        row_meta = inputs.extra.get("alignment_manifest_row", {}) if inputs.extra else {}
        mapper_mode = _infer_mapper_mode(aligner=row_meta.get("aligner"), row_meta=row_meta)
        alignment_input = str(alignment)
        if mapper_mode == "1":
            alignment_input, _ = _build_ciri3_alignment_input(
                alignment=alignment,
                row_meta=row_meta,
                mapper_mode=mapper_mode,
            )
        return {
            "alignment": alignment,
            "alignment_input": alignment_input,
            "alignment_format": alignment.suffix.lstrip(".") or "bam",
            "cell_id": inputs.cell_id,
            "threads": inputs.threads,
            "ref_fa": inputs.ref_fa or "",
            "gtf": inputs.gtf or "",
            "raw_output": raw_output,
            "outdir": run_dir,
            "extra_args": self.resolve_extra_args(),
            "mapper_mode": mapper_mode,
            "sortedness": row_meta.get("sortedness", ""),
            "artifact_bucket": row_meta.get("artifact_bucket", ""),
            "chimeric_junction": row_meta.get("chimeric_junction", ""),
            "unmapped_mate1": row_meta.get("unmapped_mate1", ""),
            "unmapped_mate2": row_meta.get("unmapped_mate2", ""),
            "bwa_sam": row_meta.get("bwa_sam", ""),
            "read_layout": inputs.effective_read_layout(),
            "group_id": inputs.alignment_group or "",
            "log_path": log_path,
        }

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        if inputs.input_mode != "alignment":
            raise ValueError("CIRI3 requires alignment input mode")
        return self.run_from_alignment(inputs)

    def run_from_alignment(self, inputs: DetectorRunInputs) -> DetectorResult:
        alignment = inputs.alignment_path
        if alignment is None:
            raise ValueError("CIRI3 requires BAM or SAM input")

        outdir = inputs.outdir
        outdir.mkdir(parents=True, exist_ok=True)
        run_dir = outdir / f"{inputs.cell_id}.ciri3_run"
        run_dir.mkdir(parents=True, exist_ok=True)
        raw_output = run_dir / "ciri3_raw.tsv"
        normalized_output = outdir / f"{inputs.cell_id}.tsv"
        log_path = outdir / f"{inputs.cell_id}.ciri3.log"

        runtime_template = self.resolve_command_template()
        ok, errors, details = self.validate_runtime(template=runtime_template)
        if not ok:
            raise RuntimeError(
                "Invalid CIRI3 runtime configuration: "
                + "; ".join(errors)
                + f". Details: {details}"
            )

        context = self.build_template_context(inputs, raw_output=raw_output, run_dir=run_dir, log_path=log_path)
        execution_mode = str(details.get("preferred_mode") or "unknown")
        internal_log = run_dir / "ciri3_internal.log"
        started_at = time.perf_counter()
        row_meta = inputs.extra.get("alignment_manifest_row", {}) if inputs.extra else {}
        mapper_mode = _infer_mapper_mode(aligner=row_meta.get("aligner"), row_meta=row_meta)
        cli_args = _resolve_ciri3_cli_args(extra_args=self.resolve_extra_args(), mapper_mode=mapper_mode)
        stringency = _option_value(cli_args, "-S", "--stringency")
        with log_path.open("w", encoding="utf-8") as log_handle:
            if self.command_template is not None:
                prep_meta = {
                    "input_file_type": "template",
                    "input_sortedness": str(row_meta.get("sortedness", "")).strip() or "unknown",
                    "mapper_mode": mapper_mode,
                }
                assert runtime_template is not None
                try:
                    cmd = runtime_template.format(**context)
                except KeyError as exc:
                    raise RuntimeError(f"CIRI3 template references unsupported placeholder: {exc}") from exc
                log_handle.write(f"$ {cmd}\n")
                log_handle.flush()
                result = subprocess.run(
                    cmd,
                    shell=True,
                    executable="/bin/bash",
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
                command_string = cmd
                execution_mode = "template"
            elif details["direct_ready"]:
                resolution = resolve_ciri3_installation()
                if inputs.ref_fa is None:
                    raise RuntimeError("CIRI3 direct execution requires ref_fa.")
                alignment_input, prep_meta = _build_ciri3_alignment_input(
                    alignment=alignment,
                    row_meta=row_meta,
                    mapper_mode=mapper_mode,
                )
                log_handle.write(
                    f"# input_file_type={prep_meta.get('input_file_type')} "
                    f"sortedness={prep_meta.get('input_sortedness')} "
                    f"mapper_mode={mapper_mode} stringency={stringency or 'default'}\n"
                )
                cmd_list = _build_direct_ciri3_command(
                    resolution=resolution,
                    alignment_input=alignment_input,
                    raw_output=raw_output,
                    ref_fa=inputs.ref_fa,
                    gtf=inputs.gtf,
                    threads=inputs.threads,
                    cli_args=cli_args,
                    internal_log=internal_log,
                )
                log_handle.write(f"$ {_shell_join(cmd_list)}\n")
                log_handle.flush()
                result = subprocess.run(
                    cmd_list,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    text=True,
                )
                command_string = _shell_join(cmd_list)
            else:
                raise RuntimeError(
                    "Invalid CIRI3 runtime configuration: no explicit command_template and no direct jar contract."
                )
        if result.returncode != 0:
            tail = _tail_text(log_path)
            suffix = f"\n--- log tail ---\n{tail}" if tail else ""
            raise RuntimeError(f"CIRI3 command failed for {inputs.cell_id}; see {log_path}.{suffix}")

        raw_row_count = _count_ciri3_rows(raw_output)
        with log_path.open("a", encoding="utf-8") as log_handle:
            log_handle.write(f"# raw_row_count={raw_row_count}\n")
        _normalize_ciri3_output(raw_output, normalized_output)
        normalized_row_count = _count_ciri3_rows(normalized_output)
        return DetectorResult(
            detector=self.name,
            cell_id=inputs.cell_id,
            outdir=outdir,
            tsv_path=normalized_output,
            run_dir=run_dir,
            log_path=log_path,
            meta={
                "threads": inputs.threads,
                "elapsed_seconds": round(time.perf_counter() - started_at, 3),
                "input_mode": "alignment",
                "read_layout": inputs.effective_read_layout(),
                "alignment_path": str(alignment),
                "command": command_string,
                "command_template_configured": runtime_template is not None,
                "command_template_fields": sorted(_template_fields(runtime_template)) if runtime_template else [],
                "execution_mode": execution_mode,
                "input_file_type": prep_meta.get("input_file_type"),
                "input_sortedness": prep_meta.get("input_sortedness"),
                "mapper_mode": mapper_mode,
                "stringency": stringency,
                "ciri3_jar": details.get("jar"),
                "ciri3_wrapper": details.get("bin"),
                "java_executable": details.get("java"),
                "ciri3_home": details.get("home"),
                "ciri3_internal_log": str(internal_log) if internal_log.exists() else None,
                "ciri3_extra_args": " ".join(cli_args),
                "raw_output_path": str(raw_output),
                "input_mode": "alignment",
                "detector_backend": self.name,
                "raw_row_count": raw_row_count,
                "normalized_row_count": normalized_row_count,
                **prep_meta,
                "expected_inputs": self.expected_inputs(),
            },
        )
