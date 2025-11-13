from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Iterable, Tuple


def _ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _run_shell(cmd: str) -> None:
    # Run through bash for compatibility with templates that include shell syntax
    subprocess.run(cmd, shell=True, check=True)


def _format_cmd(
    template: str,
    *,
    ref_fa: Path,
    gtf: Path,
    r1: Path,
    r2: Path | None,
    out_tsv: Path,
) -> str:
    return template.format(
        ref_fa=str(ref_fa),
        gtf=str(gtf),
        r1=str(r1),
        r2=str(r2) if r2 else "",
        out_tsv=str(out_tsv),
    )


def _iter_fastq_batches(fastq_dir: Path) -> Iterable[Tuple[str, Path, Path | None]]:
    """
    Expect layout:
      fastq_dir/
        Batch1/
          R1.fastq.gz
          R2.fastq.gz (optional)
        Batch2/...
    """
    for d in sorted(p for p in fastq_dir.iterdir() if p.is_dir()):
        r1 = d / "R1.fastq.gz"
        r2 = d / "R2.fastq.gz"
        if r1.exists():
            yield d.name, r1, (r2 if r2.exists() else None)


def run_cirifull_over_fastqs(
    fastq_dir: Path,
    outdir: Path,
    cmd_template: str,
    ref_fa: Path,
    gtf: Path,
    threads: int = 8,
) -> None:
    """
    Run the command template over each batch directory in fastq_dir.
    For smoke tests, the template can simply write a header into {out_tsv}.
    """
    _ensure_outdir(outdir)
    batches = list(_iter_fastq_batches(Path(fastq_dir)))
    print(f"[circyto] ciri-full(batched) jobs: {len(batches)} → {outdir}")

    for name, r1, maybe_r2 in batches:
        out_tsv = Path(outdir) / f"{name}.tsv"
        cmd = _format_cmd(
            cmd_template,
            ref_fa=Path(ref_fa),
            gtf=Path(gtf),
            r1=r1,
            r2=maybe_r2,
            out_tsv=out_tsv,
        )
        _run_shell(cmd)

    print("[circyto] run_cirifull_over_fastqs done.")


def _read_manifest_rows(manifest: Path) -> Iterable[Tuple[str, Path, Path | None]]:
    """
    Manifest TSV with headers: cell_id  r1  r2
    r2 can be missing/blank for single-end.
    """
    with open(manifest, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            cell_id = row.get("cell_id") or row.get("cell") or row.get("id") or "cell"
            r1 = Path(row.get("r1", "").strip())
            r2_raw = row.get("r2", "").strip()
            r2 = Path(r2_raw) if r2_raw else None
            if not r1:
                continue
            yield cell_id, r1, r2


def run_cirifull_with_manifest(
    manifest: Path,
    outdir: Path,
    cmd_template: str,
    ref_fa: Path,
    gtf: Path,
    threads: int = 8,
) -> None:
    """
    Run the command template for each row in a per-cell manifest TSV.
    """
    _ensure_outdir(outdir)
    rows = list(_read_manifest_rows(Path(manifest)))
    print(f"[circyto] ciri-full(manifest) jobs: {len(rows)} → {outdir}")

    for cell_id, r1, maybe_r2 in rows:
        out_tsv = Path(outdir) / f"{cell_id}.tsv"
        cmd = _format_cmd(
            cmd_template,
            ref_fa=Path(ref_fa),
            gtf=Path(gtf),
            r1=r1,
            r2=maybe_r2,
            out_tsv=out_tsv,
        )
        _run_shell(cmd)

    print("[circyto] run_cirifull_with_manifest done.")
