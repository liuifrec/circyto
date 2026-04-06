from __future__ import annotations

import shutil
import time
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from circyto.detectors import build_default_engines
from circyto.manifest.v1 import validate_manifest_tsv
from circyto.pipeline.run_detector import read_manifest, run_detector_manifest
from circyto.pipeline.collect import collect_matrix
from circyto.pipeline.collect_find_circ3 import collect_find_circ3_matrix
from circyto.writers.convert import convert_matrix_files
from circyto.demux.smartseq2 import SmartSeq2DemuxParams, demux_smartseq2_pooled

smoke_app = typer.Typer(help="One-command smoke tests (public data / dev).")
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
