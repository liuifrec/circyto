from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import time
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from circyto.detectors import build_default_engines
from circyto.detectors.base import get_detector_capabilities
from circyto.manifest.v1 import validate_manifest_tsv
from circyto.pipeline.run_detector import read_manifest, run_detector_manifest
from circyto.pipeline.collect import collect_matrix
from circyto.pipeline.collect_find_circ3 import collect_find_circ3_matrix
from circyto.pipeline.align_manifest import prepare_alignment_cache, run_detector_alignment_manifest
from circyto.writers.convert import convert_matrix_files
from circyto.demux.smartseq2 import SmartSeq2DemuxParams, demux_smartseq2_pooled
from circyto.paths import get_repo_root, get_tools_dir

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


def _select_alignment_first_fixture(
    *,
    manifest: Path | None,
    outdir: Path,
    subset_cells: int,
) -> tuple[Path, dict[str, str]]:
    if subset_cells < 1:
        raise typer.BadParameter("--subset-cells must be >= 1")

    if manifest is not None:
        rows: list[dict[str, str]] = []
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or "cell_id" not in reader.fieldnames:
                raise typer.BadParameter(f"Manifest missing required column 'cell_id': {manifest}")
            header = list(reader.fieldnames)
            for row in reader:
                rows.append({k: "" if v is None else str(v) for k, v in row.items()})
        if not rows:
            raise typer.BadParameter(f"Manifest contains 0 data rows: {manifest}")
        selected = rows[:subset_cells]
        selected_manifest = outdir / "smoke_manifest.tsv"
        with selected_manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            writer.writeheader()
            for row in selected:
                writer.writerow({key: row.get(key, "") for key in header})
        return selected_manifest, {
            "fixture_kind": "manifest-override",
            "fixture_manifest": str(manifest),
            "fixture_cells": ",".join(row["cell_id"] for row in selected),
        }

    repo_root = get_repo_root()
    tools_dir = get_tools_dir()
    if repo_root is None or tools_dir is None:
        raise typer.BadParameter("Could not resolve repo-local smoke fixtures.")
    if subset_cells != 1:
        raise typer.BadParameter(
            "Built-in alignment-first smoke fixtures only provide one stable paired-end sample. "
            "Use --manifest to smoke test more than one cell."
        )

    default_manifest = repo_root / "manifest_2.tsv"
    if default_manifest.exists():
        rows: list[dict[str, str]] = []
        with default_manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames:
                rows = [{k: "" if v is None else str(v) for k, v in row.items()} for row in reader]
        if rows:
            first = rows[0]
            source_r1 = (repo_root / first.get("r1", "")).resolve()
            source_r2 = (repo_root / first.get("r2", "")).resolve()
            if source_r1.exists() and source_r2.exists():
                subset_dir = outdir / "fixture_fastq"
                subset_dir.mkdir(parents=True, exist_ok=True)
                r1 = subset_dir / f"{first['cell_id']}_R1_1k.fastq"
                r2 = subset_dir / f"{first['cell_id']}_R2_1k.fastq"

                def _subset_fastq(src: Path, dst: Path, n_reads: int = 1000) -> None:
                    with gzip.open(src, "rt") as in_handle, dst.open("w", encoding="utf-8") as out_handle:
                        for _ in range(n_reads * 4):
                            line = in_handle.readline()
                            if not line:
                                break
                            out_handle.write(line)

                if not r1.exists():
                    _subset_fastq(source_r1, r1)
                if not r2.exists():
                    _subset_fastq(source_r2, r2)

                selected_manifest = outdir / "smoke_manifest.tsv"
                selected_manifest.write_text(
                    "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id\n"
                    f"{first['cell_id']}_1k\tsmartseq2\t{r1}\t{r2}\t\tciri3-smoke\t0\tsmoke\n",
                    encoding="utf-8",
                )
                return selected_manifest, {
                    "fixture_kind": "repo-chr21-subset",
                    "fixture_manifest": str(default_manifest),
                    "fixture_source_cell": first["cell_id"],
                    "fixture_r1": str(r1),
                    "fixture_r2": str(r2),
                    "fixture_cells": f"{first['cell_id']}_1k",
                }

    r1 = tools_dir / "CIRI-full_v2.0" / "CIRI-full_test" / "test_1.fq.gz"
    r2 = tools_dir / "CIRI-full_v2.0" / "CIRI-full_test" / "test_2.fq.gz"
    if not r1.exists() or not r2.exists():
        raise typer.BadParameter(f"Smoke FASTQs not found: {r1} {r2}")

    selected_manifest = outdir / "smoke_manifest.tsv"
    selected_manifest.write_text(
        "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id\n"
        f"smoke1\tsmartseq2\t{r1}\t{r2}\t\tciri3-smoke\t0\tsmoke\n",
        encoding="utf-8",
    )
    return selected_manifest, {
        "fixture_kind": "repo-bundled",
        "fixture_manifest": str(selected_manifest),
        "fixture_r1": str(r1),
        "fixture_r2": str(r2),
        "fixture_cells": "smoke1",
    }


def _resolve_smoke_references(
    *,
    ref_fa: Path | None,
    gtf: Path | None,
    genome_dir: Path | None,
    aligner: str,
) -> tuple[Path, Path | None, Path | None]:
    repo_root = get_repo_root()
    default_ref = repo_root / "ref" / "chr21.fa" if repo_root is not None else None
    default_gtf = repo_root / "ref" / "chr21.gtf" if repo_root is not None else None
    default_genome_dir = repo_root / "ref" / "star_index_chr21" if repo_root is not None else None

    resolved_ref = ref_fa or default_ref
    resolved_gtf = gtf or default_gtf
    resolved_genome_dir = genome_dir or default_genome_dir
    if resolved_ref is None or not resolved_ref.exists():
        raise typer.BadParameter(f"Reference FASTA not found: {resolved_ref}")
    if aligner == "star" and (resolved_genome_dir is None or not resolved_genome_dir.exists()):
        raise typer.BadParameter(
            "STAR smoke requires --genome-dir or a repo-local ref/star_index_chr21 fixture."
        )
    return resolved_ref, resolved_gtf, resolved_genome_dir


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

    The default smoke path is alignment-first CIRI3 with a tiny local fixture.
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
    if aligner not in {"bwa-mem", "star"}:
        raise typer.BadParameter("Smoke aligner must be one of: bwa-mem, star")
    if not use_testdata and manifest is None:
        raise typer.BadParameter("--use-local-manifest requires --manifest")

    engines = build_default_engines()
    if detector not in engines:
        avail = ", ".join(sorted(engines.keys()))
        raise typer.BadParameter(f"Unknown detector '{detector}'. Available: {avail}")
    det_engine = engines[detector]
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
    smoke_manifest, fixture_meta = _select_alignment_first_fixture(
        manifest=manifest_override if not use_testdata else None,
        outdir=outdir,
        subset_cells=subset_cells,
    )
    if use_testdata and manifest_override is not None:
        smoke_manifest, fixture_meta = _select_alignment_first_fixture(
            manifest=manifest_override,
            outdir=outdir,
            subset_cells=subset_cells,
        )
    resolved_ref, resolved_gtf, resolved_genome_dir = _resolve_smoke_references(
        ref_fa=ref_fa,
        gtf=gtf,
        genome_dir=genome_dir,
        aligner=aligner,
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
