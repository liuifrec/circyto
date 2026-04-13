from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from circyto.detectors import build_default_engines
from circyto.detectors.base import get_detector_capabilities
from circyto.detectors.ciri3 import Ciri3Detector
from circyto.manifest.v1 import validate_manifest_tsv
from circyto.pipeline.run_detector import read_manifest, run_detector_manifest
from circyto.pipeline.collect import collect_matrix
from circyto.pipeline.collect_find_circ3 import collect_find_circ3_matrix
from circyto.pipeline.align_manifest import prepare_alignment_cache, run_detector_alignment_manifest
from circyto.writers.convert import convert_matrix_files
from circyto.demux.smartseq2 import SmartSeq2DemuxParams, demux_smartseq2_pooled
from circyto.paths import get_packaged_smoke_demo_dir

smoke_app = typer.Typer(
    help="One-command smoke tests (public data / dev).",
    invoke_without_command=True,
    no_args_is_help=False,
)
console = Console()


def _detector_output_has_calls(path: Path, detector_name: str) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False

    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if detector_name in {"ciri-full", "ciri2"}:
        data_lines = [line for line in text[1:] if line.strip()]
        return len(data_lines) > 0

    data_lines = [line for line in text if line.strip() and not line.startswith("#")]
    return len(data_lines) > 0


def _assert_detector_outputs_nonempty(run_dir: Path, detector_name: str) -> None:
    if detector_name == "find-circ3":
        output_paths = sorted(run_dir.glob("*/*_splice_sites.bed"))
    else:
        output_paths = sorted(run_dir.glob("*.tsv"))

    if not output_paths:
        raise RuntimeError(
            f"Detector stage completed without producing per-cell outputs under {run_dir}."
        )

    nonempty = [p for p in output_paths if _detector_output_has_calls(p, detector_name)]
    if nonempty:
        console.print(
            f"[green]Detector outputs with calls[/green]: {len(nonempty)}/{len(output_paths)}"
        )
        return

    raise RuntimeError(
        f"Detector stage produced only empty outputs under {run_dir}. "
        "Refusing to continue to an all-zero matrix."
    )


def _is_nonempty_file(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _should_reuse_demux(manifest_path: Path) -> bool:
    ok, errors, summary = validate_manifest_tsv(manifest_path, strict=True)
    if not ok:
        console.print(
            f"[yellow]Existing demux manifest invalid; rerunning demux[/yellow]: {'; '.join(errors[:3])}"
        )
        return False
    if summary["cells"] == 0:
        console.print("[yellow]Existing demux manifest has 0 cells; rerunning demux[/yellow]")
        return False
    return True


def _log_stage_done(stage: str, started_at: float, extra: str = "") -> None:
    elapsed = time.perf_counter() - started_at
    suffix = f" {extra}" if extra else ""
    console.print(f"[green]{stage} complete[/green] in {elapsed:.2f}s{suffix}")


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_matrix_market_dims(matrix_path: Path) -> tuple[int, int, int]:
    if not matrix_path.exists():
        return 0, 0, 0
    with matrix_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("%"):
                continue
            parts = text.split()
            if len(parts) >= 3:
                try:
                    return int(parts[0]), int(parts[1]), int(parts[2])
                except ValueError:
                    return 0, 0, 0
            return 0, 0, 0
    return 0, 0, 0


def _count_detector_rows(tsv_path: Path) -> int:
    if not tsv_path.exists():
        return 0
    lines = tsv_path.read_text(encoding="utf-8", errors="replace").splitlines()
    return len([line for line in lines[1:] if line.strip()])


def _pick_manifest_value(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = (row.get(key) or "").strip()
        if value:
            return value
    return ""


def _read_smoke_manifest_rows(manifest: Path) -> tuple[list[str], list[dict[str, str]]]:
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "cell_id" not in reader.fieldnames:
            raise typer.BadParameter(f"Manifest missing required column 'cell_id': {manifest}")
        header = list(reader.fieldnames)
        rows = [{k: "" if v is None else str(v) for k, v in row.items()} for row in reader]
    if not rows:
        raise typer.BadParameter(f"Manifest contains 0 data rows: {manifest}")
    return header, rows


def _validate_smoke_manifest_selection(
    manifest: Path,
    rows: list[dict[str, str]],
    *,
    fixture_kind: str,
) -> dict[str, str]:
    layouts: set[str] = set()
    selected_fastqs: list[str] = []
    selected_cells: list[str] = []
    for index, row in enumerate(rows, start=1):
        cell_id = (row.get("cell_id") or "").strip() or f"row{index}"
        read1_raw = _pick_manifest_value(row, "read1", "r1")
        read2_raw = _pick_manifest_value(row, "read2", "r2")
        declared_layout = (row.get("read_layout") or "").strip()
        if declared_layout and declared_layout not in {"single-end", "paired-end"}:
            raise typer.BadParameter(
                f"Invalid read_layout for smoke fixture row cell_id={cell_id}: {declared_layout}"
            )
        inferred_layout = declared_layout or ("paired-end" if read2_raw else "single-end")
        if not read1_raw:
            raise typer.BadParameter(
                f"Smoke fixture row cell_id={cell_id} is missing read1/r1. "
                f"Fixture kind={fixture_kind}. Pass --manifest explicitly if needed."
            )
        if inferred_layout == "single-end" and read2_raw:
            raise typer.BadParameter(
                f"Smoke fixture row cell_id={cell_id} declares single-end but also provides read2/r2={read2_raw}. "
                "Pass --manifest explicitly if needed."
            )
        if inferred_layout == "paired-end" and not read2_raw:
            raise typer.BadParameter(
                f"Smoke fixture row cell_id={cell_id} requires paired-end input but is missing read2/r2. "
                "Pass --manifest explicitly if needed."
            )
        read1 = (manifest.parent / read1_raw) if not Path(read1_raw).is_absolute() else Path(read1_raw)
        if not read1.exists():
            raise typer.BadParameter(
                f"Smoke fixture read1 missing for cell_id={cell_id}: {read1} "
                f"(mode={inferred_layout}). Pass --manifest explicitly if needed."
            )
        selected_fastqs.append(str(read1.resolve()))
        if inferred_layout == "paired-end":
            read2 = (manifest.parent / read2_raw) if not Path(read2_raw).is_absolute() else Path(read2_raw)
            if not read2.exists():
                raise typer.BadParameter(
                    f"Smoke fixture read2 missing for cell_id={cell_id}: {read2} "
                    f"(mode=paired-end). Pass --manifest explicitly if needed."
                )
            selected_fastqs.append(str(read2.resolve()))
        layouts.add(inferred_layout)
        selected_cells.append(cell_id)

    if len(layouts) != 1:
        raise typer.BadParameter(
            f"Smoke fixture selection mixes read layouts: {', '.join(sorted(layouts))}. "
            "Choose a single-layout fixture set or pass --manifest explicitly."
        )

    return {
        "read_layout": next(iter(layouts)),
        "fixture_cells": ",".join(selected_cells),
        "selected_fastqs": ",".join(selected_fastqs),
    }


def _copy_smoke_resource(src: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return dest


def _stage_packaged_smoke_demo(
    *,
    outdir: Path,
    read_layout: str,
) -> tuple[Path, Path, Path, dict[str, str]]:
    demo_dir = get_packaged_smoke_demo_dir()
    if not demo_dir.exists():
        raise typer.BadParameter(f"Packaged smoke demo directory not found: {demo_dir}")
    if read_layout not in {"single-end", "paired-end"}:
        raise typer.BadParameter(f"Unsupported smoke read_layout: {read_layout}")

    stage_dir = outdir / "demo"
    stage_dir.mkdir(parents=True, exist_ok=True)
    staged_ref = _copy_smoke_resource(demo_dir / "reference.fa", stage_dir / "reference.fa")
    staged_gtf = _copy_smoke_resource(demo_dir / "reference.gtf", stage_dir / "reference.gtf")

    if read_layout == "single-end":
        staged_r1 = _copy_smoke_resource(demo_dir / "single_end.fastq", stage_dir / "single_end.fastq")
        staged_r2: Path | None = None
        cell_id = "smoke_se"
    else:
        staged_r1 = _copy_smoke_resource(demo_dir / "paired_end_R1.fastq", stage_dir / "paired_end_R1.fastq")
        staged_r2 = _copy_smoke_resource(demo_dir / "paired_end_R2.fastq", stage_dir / "paired_end_R2.fastq")
        cell_id = "smoke_pe"

    manifest = outdir / "smoke_manifest.tsv"
    header = [
        "cell_id",
        "platform",
        "read1",
        "read2",
        "bam",
        "library_id",
        "n_input_reads",
        "group_id",
        "read_layout",
    ]
    row = {
        "cell_id": cell_id,
        "platform": "smartseq2",
        "read1": str(staged_r1.resolve()),
        "read2": str(staged_r2.resolve()) if staged_r2 is not None else "",
        "bam": "",
        "library_id": "packaged_smoke_demo",
        "n_input_reads": "2",
        "group_id": "smoke",
        "read_layout": read_layout,
    }
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    return (
        manifest,
        staged_ref,
        staged_gtf,
        {
            "fixture_kind": "packaged-demo",
            "fixture_manifest": str(manifest),
            "fixture_cells": cell_id,
            "read_layout": read_layout,
            "selected_fastqs": ",".join(
                [str(staged_r1.resolve())] + ([str(staged_r2.resolve())] if staged_r2 is not None else [])
            ),
            "packaged_demo_dir": str(demo_dir),
            "staged_demo_dir": str(stage_dir),
        },
    )


def _build_smoke_detector(detector_name: str, detector_engine):
    if detector_name == "ciri3" and isinstance(detector_engine, Ciri3Detector):
        script = get_packaged_smoke_demo_dir() / "ciri3_smoke_template.sh"
        if not script.exists():
            raise typer.BadParameter(f"Packaged smoke template not found: {script}")
        return Ciri3Detector(
            command_template=f"bash {script} {{alignment}} {{raw_output}} {{cell_id}}",
        )
    return detector_engine


def _select_alignment_first_fixture(
    *,
    manifest: Path | None,
    outdir: Path,
    subset_cells: int,
    read_layout: str,
) -> tuple[Path, Path | None, Path | None, dict[str, str]]:
    if subset_cells < 1:
        raise typer.BadParameter("--subset-cells must be >= 1")

    if manifest is not None:
        if not manifest.exists():
            raise typer.BadParameter(f"Smoke manifest not found: {manifest}")
        header, rows = _read_smoke_manifest_rows(manifest)
        selected = rows[:subset_cells]
        selected_manifest = outdir / "smoke_manifest.tsv"
        with selected_manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            writer.writeheader()
            for row in selected:
                writer.writerow({key: row.get(key, "") for key in header})
        meta = _validate_smoke_manifest_selection(
            selected_manifest,
            selected,
            fixture_kind="manifest-override",
        )
        return selected_manifest, None, None, {
            "fixture_kind": "manifest-override",
            "fixture_manifest": str(manifest),
            **meta,
        }

    if subset_cells != 1:
        raise typer.BadParameter(
            "The packaged smoke demo provides one stable demo sample per layout. "
            "Use --manifest to smoke test more than one cell."
        )
    manifest_path, staged_ref, staged_gtf, meta = _stage_packaged_smoke_demo(
        outdir=outdir,
        read_layout=read_layout,
    )
    return manifest_path, staged_ref, staged_gtf, meta


def _resolve_smoke_references(
    *,
    ref_fa: Path | None,
    gtf: Path | None,
    genome_dir: Path | None,
    aligner: str,
    fixture_ref: Path | None,
    fixture_gtf: Path | None,
) -> tuple[Path, Path | None, Path | None]:
    resolved_ref = ref_fa or fixture_ref
    resolved_gtf = gtf or fixture_gtf
    resolved_genome_dir = genome_dir
    if resolved_ref is None or not resolved_ref.exists():
        raise typer.BadParameter(f"Reference FASTA not found: {resolved_ref}")
    if resolved_gtf is not None and not resolved_gtf.exists():
        raise typer.BadParameter(f"GTF not found: {resolved_gtf}")
    return resolved_ref, resolved_gtf, resolved_genome_dir


def _require_smoke_command(name: str, *, purpose: str) -> None:
    if shutil.which(name) is None:
        raise typer.BadParameter(f"Smoke requires '{name}' for {purpose}, but it was not found in PATH.")


def _run_smoke_command(cmd: list[str], *, cwd: Path, label: str) -> None:
    result = subprocess.run(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        tail = "\n".join((result.stdout or "").splitlines()[-40:])
        raise RuntimeError(
            f"Smoke {label} failed with exit={result.returncode}: {' '.join(cmd)}"
            + (f"\n--- output tail ---\n{tail}" if tail else "")
        )


def _ensure_bwa_index(ref_fa: Path) -> list[str]:
    suffixes = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    if all((ref_fa.parent / f"{ref_fa.name}{suffix}").exists() for suffix in suffixes):
        return []
    _require_smoke_command("bwa", purpose="BWA smoke alignment")
    _run_smoke_command(["bwa", "index", str(ref_fa)], cwd=ref_fa.parent, label="bwa index")
    return [str(ref_fa.parent / f"{ref_fa.name}{suffix}") for suffix in suffixes]


def _ensure_star_genome_dir(ref_fa: Path, gtf: Path | None, genome_dir: Path) -> Path:
    if gtf is None:
        raise typer.BadParameter("STAR smoke requires a GTF annotation.")
    genome_dir.mkdir(parents=True, exist_ok=True)
    if (genome_dir / "genomeParameters.txt").exists():
        return genome_dir
    _require_smoke_command("STAR", purpose="STAR smoke alignment")
    _run_smoke_command(
        [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--runThreadN",
            "1",
            "--genomeDir",
            str(genome_dir),
            "--genomeFastaFiles",
            str(ref_fa),
            "--sjdbGTFfile",
            str(gtf),
            "--genomeSAindexNbases",
            "3",
        ],
        cwd=genome_dir,
        label="STAR genomeGenerate",
    )
    return genome_dir


def _summarize_alignment_first_smoke(
    *,
    outdir: Path,
    detector: str,
    aligner: str,
    manifest: Path,
    alignment_manifest: Path,
    detector_dir: Path,
    matrix_dir: Path,
    allow_empty: bool,
) -> dict:
    detector_summary_path = detector_dir / "detector_run_summary.json"
    prepare_summary_path = outdir / "align" / "alignment_prepare_summary.json"
    detector_summary = json.loads(detector_summary_path.read_text(encoding="utf-8"))
    prepare_summary = json.loads(prepare_summary_path.read_text(encoding="utf-8"))
    matrix_path = matrix_dir / "circ_counts.mtx"
    circ_index = matrix_dir / "circ_index.txt"
    cell_index = matrix_dir / "cell_index.txt"
    n_rows, n_cols, nnz = _parse_matrix_market_dims(matrix_path)
    cell_record = detector_summary["cells"][0] if detector_summary.get("cells") else {}
    tsv_path = Path(cell_record["tsv_path"]) if cell_record.get("tsv_path") else detector_dir / "unknown.tsv"
    normalized_rows = int(cell_record.get("normalized_row_count") or 0)
    workflow_success = all(
        path.exists()
        for path in (manifest, alignment_manifest, detector_summary_path, matrix_path, circ_index, cell_index)
    )
    biological_nonempty = normalized_rows > 0 or nnz > 0
    smoke_pass = workflow_success and (allow_empty or biological_nonempty)
    return {
        "detector": detector,
        "aligner": aligner,
        "mode": "alignment-first",
        "manifest": str(manifest),
        "alignment_manifest": str(alignment_manifest),
        "prepare_summary": str(prepare_summary_path),
        "detector_summary": str(detector_summary_path),
        "tsv_path": str(tsv_path),
        "matrix_path": str(matrix_path),
        "circ_index_path": str(circ_index),
        "cell_index_path": str(cell_index),
        "cells_used": len(prepare_summary.get("cells", [])),
        "workflow_stages": {
            "prepare_alignment_cache": True,
            "run_detector_from_alignments": True,
            "collect_matrix": True,
        },
        "input_fixture": prepare_summary.get("manifest"),
        "input_file_type": cell_record.get("input_file_type"),
        "mapper_mode": cell_record.get("mapper_mode"),
        "normalized_row_count": normalized_rows,
        "matrix_rows": n_rows,
        "matrix_cols": n_cols,
        "matrix_nnz": nnz,
        "allow_empty": allow_empty,
        "biological_nonempty": biological_nonempty,
        "workflow_success": workflow_success,
        "smoke_pass": smoke_pass,
    }


@smoke_app.callback(invoke_without_command=True)
def smoke_root(
    ctx: typer.Context,
    detector: str = typer.Option("ciri3", "--detector", "-d"),
    aligner: str = typer.Option("bwa-mem", "--aligner"),
    mode: str = typer.Option("alignment-first", "--mode"),
    read_layout: str = typer.Option("paired-end", "--read-layout"),
    outdir: Path = typer.Option(Path("work/smoke"), "--outdir", "-o"),
    ref_fa: Optional[Path] = typer.Option(None, "--ref-fa"),
    gtf: Optional[Path] = typer.Option(None, "--gtf"),
    manifest: Optional[Path] = typer.Option(None, "--manifest"),
    genome_dir: Optional[Path] = typer.Option(None, "--genome-dir"),
    threads: int = typer.Option(2, "--threads"),
    parallel: int = typer.Option(1, "--parallel"),
    subset_cells: int = typer.Option(1, "--subset-cells"),
    allow_empty: bool = typer.Option(True, "--allow-empty/--require-nonempty"),
    tmpdir: Optional[Path] = typer.Option(None, "--tmpdir"),
    use_testdata: bool = typer.Option(True, "--use-testdata/--use-local-manifest"),
    keep_workdir: bool = typer.Option(True, "--keep-workdir/--cleanup"),
    json_summary: bool = typer.Option(False, "--json-summary"),
) -> None:
    """
    Run a fast installation and workflow smoke test.

    The default smoke path is alignment-first CIRI3 with packaged demo inputs.
    Success means the workflow completed cleanly. Empty outputs are still a PASS
    unless `--require-nonempty` is requested.
    """
    if ctx.invoked_subcommand is not None:
        return

    detector = detector.lower().strip()
    aligner = aligner.lower().strip()
    mode = mode.lower().strip()
    if mode != "alignment-first":
        raise typer.BadParameter(
            "circyto smoke currently supports alignment-first smoke only. "
            "Use 'circyto smoke smartseq2 ...' for the older per-cell-fastq smoke path."
        )
    if detector != "ciri3":
        raise typer.BadParameter(
            "circyto smoke currently supports detector=ciri3 only. "
            "Other detectors should be smoke-tested via their existing workflows."
        )
    if read_layout not in {"single-end", "paired-end"}:
        raise typer.BadParameter("--read-layout must be single-end or paired-end")
    if aligner not in {"bwa-mem", "star"}:
        raise typer.BadParameter("Smoke aligner must be one of: bwa-mem, star")
    if not use_testdata and manifest is None:
        raise typer.BadParameter("--use-local-manifest requires --manifest")
    if aligner == "star" and read_layout != "paired-end":
        raise typer.BadParameter("STAR smoke currently requires --read-layout paired-end")

    engines = build_default_engines()
    if detector not in engines:
        avail = ", ".join(sorted(engines.keys()))
        raise typer.BadParameter(f"Unknown detector '{detector}'. Available: {avail}")
    det_engine = engines[detector]
    det_engine = _build_smoke_detector(detector, det_engine)
    caps = get_detector_capabilities(det_engine)
    if not caps.accepts_alignment:
        raise typer.BadParameter(f"Detector '{detector}' does not support alignment-first smoke.")
    if aligner == "star" and not caps.supports_star:
        raise typer.BadParameter(f"Detector '{detector}' does not support STAR alignment-first smoke.")
    if aligner == "bwa-mem" and not caps.supports_bwa:
        raise typer.BadParameter(f"Detector '{detector}' does not support BWA alignment-first smoke.")

    if hasattr(det_engine, "validate_runtime"):
        ok, errors, _ = det_engine.validate_runtime()
        if not ok:
            raise typer.BadParameter("CIRI3 runtime not ready: " + "; ".join(errors))

    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_override = manifest.resolve() if manifest is not None else None
    smoke_manifest, fixture_ref, fixture_gtf, fixture_meta = _select_alignment_first_fixture(
        manifest=manifest_override if not use_testdata else None,
        outdir=outdir,
        subset_cells=subset_cells,
        read_layout=read_layout,
    )
    if use_testdata and manifest_override is not None:
        smoke_manifest, fixture_ref, fixture_gtf, fixture_meta = _select_alignment_first_fixture(
            manifest=manifest_override,
            outdir=outdir,
            subset_cells=subset_cells,
            read_layout=read_layout,
        )
    resolved_ref, resolved_gtf, resolved_genome_dir = _resolve_smoke_references(
        ref_fa=ref_fa,
        gtf=gtf,
        genome_dir=genome_dir,
        aligner=aligner,
        fixture_ref=fixture_ref,
        fixture_gtf=fixture_gtf,
    )

    generated_index_paths = _ensure_bwa_index(resolved_ref)
    if aligner == "star":
        resolved_genome_dir = _ensure_star_genome_dir(
            resolved_ref,
            resolved_gtf,
            resolved_genome_dir or (outdir / "demo" / "star_index"),
        )

    extra_flags = ""
    previous_star_tmpdir = os.environ.get("CIRCYTO_STAR_TMPDIR")
    if aligner == "star" and resolved_genome_dir is not None:
        extra_flags = f"--genomeDir {resolved_genome_dir}"
    if tmpdir is not None:
        os.environ["CIRCYTO_STAR_TMPDIR"] = str(tmpdir)

    prepare_dir = outdir / "align"
    detector_dir = outdir / detector
    matrix_dir = outdir / "matrix"
    summary_path = outdir / "smoke_summary.json"
    prepare_dir.mkdir(parents=True, exist_ok=True)
    detector_dir.mkdir(parents=True, exist_ok=True)
    matrix_dir.mkdir(parents=True, exist_ok=True)
    selected_fastqs = [item for item in fixture_meta.get("selected_fastqs", "").split(",") if item]
    if not json_summary:
        console.print(
            f"fixture={fixture_meta['fixture_kind']} manifest={smoke_manifest} read_layout={fixture_meta.get('read_layout', 'unknown')}"
        )
        if selected_fastqs:
            console.print("selected_fastqs=" + ", ".join(selected_fastqs))
        console.print(f"reference={resolved_ref}")
        if aligner == "star":
            console.print(f"genome_dir={resolved_genome_dir}")
    started = time.perf_counter()
    try:
        alignment_manifest = prepare_alignment_cache(
            manifest=smoke_manifest,
            outdir=prepare_dir,
            aligner=aligner,
            ref_fa=resolved_ref,
            detector_hint=detector,
            threads=threads,
            parallel=parallel,
            chunk_size=1,
            extra_flags=extra_flags,
        )
        run_detector_alignment_manifest(
            detector=det_engine,
            manifest=alignment_manifest,
            outdir=detector_dir,
            ref_fa=resolved_ref,
            gtf=resolved_gtf,
            threads=threads,
            parallel=parallel,
        )
        collect_matrix(
            str(detector_dir),
            str(matrix_dir / "circ_counts.mtx"),
            str(matrix_dir / "circ_index.txt"),
            str(matrix_dir / "cell_index.txt"),
        )
        summary = _summarize_alignment_first_smoke(
            outdir=outdir,
            detector=detector,
            aligner=aligner,
            manifest=smoke_manifest,
            alignment_manifest=alignment_manifest,
            detector_dir=detector_dir,
            matrix_dir=matrix_dir,
            allow_empty=allow_empty,
        )
        summary.update(
            {
                "elapsed_seconds": round(time.perf_counter() - started, 3),
                "fixture": fixture_meta,
                "ref_fa": str(resolved_ref),
                "gtf": str(resolved_gtf) if resolved_gtf else None,
                "genome_dir": str(resolved_genome_dir) if resolved_genome_dir else None,
                "tmpdir": str(tmpdir) if tmpdir else None,
                "read_layout": read_layout,
                "generated_index_paths": generated_index_paths,
                "summary_path": str(summary_path),
            }
        )
        _write_json(summary_path, summary)
        if json_summary:
            typer.echo(json.dumps(summary, indent=2, sort_keys=True))
        else:
            console.print(
                "[bold]SMOKE[/bold]",
                f"detector={detector}",
                f"aligner={aligner}",
                f"read_layout={read_layout}",
                f"cells={summary['cells_used']}",
                f"matrix={summary['matrix_rows']}x{summary['matrix_cols']}",
                f"nnz={summary['matrix_nnz']}",
                f"allow_empty={allow_empty}",
                f"pass={summary['smoke_pass']}",
            )
            console.print(f"fixture={summary['fixture']['fixture_kind']} manifest={smoke_manifest}")
            console.print(f"artifacts: align={prepare_dir} detector={detector_dir} matrix={matrix_dir}")
            console.print(f"summary={summary_path}")
        if not summary["smoke_pass"]:
            raise typer.Exit(code=1)
        if not keep_workdir:
            for path in (prepare_dir, detector_dir, matrix_dir):
                if path.exists():
                    shutil.rmtree(path)
        return
    finally:
        if tmpdir is not None:
            if previous_star_tmpdir is None:
                os.environ.pop("CIRCYTO_STAR_TMPDIR", None)
            else:
                os.environ["CIRCYTO_STAR_TMPDIR"] = previous_star_tmpdir


@smoke_app.command("smartseq2")
def smoke_smartseq2(
    r1: Path = typer.Option(..., "--r1", help="Pooled R1 fastq(.gz)"),
    r2: Path = typer.Option(..., "--r2", help="Pooled R2 fastq(.gz)"),
    barcodes: Path = typer.Option(..., "--barcodes", help="TSV: cell_id<TAB>barcode"),
    outdir: Path = typer.Option(Path("work/smoke_smartseq2"), "-o", "--outdir"),
    library_id: str = typer.Option("smartseq2_smoke", "--library-id"),
    detector: str = typer.Option(
        "find-circ3", "-d", "--detector",
        help="Detector name (e.g., find-circ3, ciri2, ciri-full)"
    ),
    ref_fa: Path = typer.Option(..., "--ref-fa"),
    gtf: Optional[Path] = typer.Option(None, "--gtf"),
    threads: int = typer.Option(8, "--threads"),
    parallel: int = typer.Option(1, "--parallel"),
    limit_reads: int = typer.Option(100000, "--limit-reads"),
    # Smart-seq2 defaults (locked for pooled public data like GSE145700)
    barcode_read: str = typer.Option("R2", "--barcode-read"),
    barcode_start: int = typer.Option(0, "--barcode-start"),
    barcode_length: int = typer.Option(8, "--barcode-length"),
    overwrite: bool = typer.Option(False, "--overwrite"),
):
    """
    Demux Smart-seq2 pooled FASTQs -> run detector -> collect matrix -> convert to h5ad.
    Strong defaults for GSE145700-style pooled reads.
    """
    outdir = outdir.resolve()

    engines = build_default_engines()
    if detector not in engines:
        avail = ", ".join(sorted(engines.keys()))
        raise typer.BadParameter(f"Unknown detector '{detector}'. Available: {avail}")
    det_engine = engines[detector]

    demux_dir = outdir / "demux"
    run_dir = outdir / "run" / det_engine.name
    matrix_dir = outdir / "matrix" / det_engine.name
    h5ad_dir = outdir / "h5ad"
    h5ad_path = h5ad_dir / f"{det_engine.name}.h5ad"
    manifest_path = demux_dir / "manifest.tsv"

    if overwrite:
        for p in (demux_dir, run_dir, matrix_dir, h5ad_dir):
            if p.exists():
                shutil.rmtree(p)

    # create dirs
    demux_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)
    matrix_dir.mkdir(parents=True, exist_ok=True)
    h5ad_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"[bold]SMOKE smartseq2[/bold] outdir={outdir}")
    console.print(f"demux -> {demux_dir}")
    console.print(f"run   -> {run_dir}")
    console.print(f"matrix-> {matrix_dir}")
    console.print(f"h5ad  -> {h5ad_path}")

    # 1) DEMUX
    demux_started = time.perf_counter()
    if overwrite or not _should_reuse_demux(manifest_path):
        params = SmartSeq2DemuxParams(
            r1=r1,
            r2=r2,
            barcodes_tsv=barcodes,
            outdir=demux_dir,
            manifest_path=manifest_path,
            library_id=library_id,
            barcode_read=barcode_read,
            barcode_start=barcode_start,
            barcode_length=barcode_length,
            umi_start=None,
            umi_length=None,
            trim_barcode=True,
            trim_umi=False,
            max_mismatch=0,
            limit_reads=limit_reads,
            overwrite=overwrite,
        )
        demux_stats = demux_smartseq2_pooled(params)
        if not manifest_path.exists():
            raise FileNotFoundError(f"Manifest not found: {manifest_path}")
        _log_stage_done("Demux", demux_started, f"cells={len(demux_stats)} manifest={manifest_path}")
    else:
        manifest_rows = read_manifest(manifest_path, validate_files=True)
        _log_stage_done("Demux", demux_started, f"reused manifest={manifest_path} cells={len(manifest_rows)}")

    # 2) RUN DETECTOR
    detector_started = time.perf_counter()
    detector_results = run_detector_manifest(
        detector=det_engine,
        manifest=manifest_path,
        outdir=run_dir,
        ref_fa=ref_fa,
        gtf=gtf,
        threads=threads,
        parallel=parallel,
    ) or []
    _assert_detector_outputs_nonempty(run_dir, det_engine.name)
    detector_status_counts: dict[str, int] = {}
    for result in detector_results:
        status = "skipped_existing" if result.meta.get("skipped_existing") else "success"
        detector_status_counts[status] = detector_status_counts.get(status, 0) + 1
    status_text = ", ".join(f"{k}={v}" for k, v in sorted(detector_status_counts.items()))
    _log_stage_done("Detector", detector_started, f"detector={det_engine.name} {status_text}")

    # 3) COLLECT MATRIX
    collect_started = time.perf_counter()
    matrix_path = matrix_dir / "circ_counts.mtx"
    circ_index_path = matrix_dir / "circ_index.txt"
    cell_index_path = matrix_dir / "cell_index.txt"

    if (
        not overwrite
        and _is_nonempty_file(matrix_path)
        and circ_index_path.exists()
        and cell_index_path.exists()
    ):
        _log_stage_done("Collect", collect_started, f"reused matrix={matrix_path}")
    else:
        if det_engine.name == "find-circ3":
            collect_find_circ3_matrix(
                str(run_dir), str(matrix_path), str(circ_index_path), str(cell_index_path)
            )
        else:
            collect_matrix(str(run_dir), str(matrix_path), str(circ_index_path), str(cell_index_path))
        _log_stage_done("Collect", collect_started, f"matrix={matrix_path}")

    # 4) CONVERT -> H5AD
    convert_started = time.perf_counter()
    if not overwrite and _is_nonempty_file(h5ad_path):
        _log_stage_done("Convert", convert_started, f"reused h5ad={h5ad_path}")
    else:
        convert_matrix_files(
            matrix_path=matrix_path,
            circ_index_path=circ_index_path,
            cell_index_path=cell_index_path,
            h5ad=h5ad_path,
        )
        _log_stage_done("Convert", convert_started, f"h5ad={h5ad_path}")
    console.print(f"[bold green]DONE[/bold green] {h5ad_path}")
