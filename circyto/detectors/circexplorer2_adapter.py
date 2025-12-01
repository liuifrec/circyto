# circyto/detectors/circexplorer2_adapter.py
from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Mapping, Sequence


def _run(cmd: Sequence[str], cwd: str | None = None) -> None:
    printable = " ".join(cmd)
    print(f"[circexplorer2] {printable}")
    subprocess.run(cmd, check=True, cwd=cwd)


def _clean_parse_output(parse_out: Path) -> None:
    # (unchanged: whatever you currently have here)
    lines = parse_out.read_text().splitlines()
    if not lines:
        return
    parse_out.write_text("\n".join(line for line in lines if line.strip()))


def run_circexplorer2(
    *,
    sample: Mapping[str, str],
    outdir_root: str,
    reference_fa: str,
    star_index: str | None,
    ref_flat: str | None,
    threads: int = 1,
    extra_args: Sequence[str] | None = None,
    skip_star: bool = False,
    star_tmp_base: str | None = None,
) -> None:
    """
    Run the CIRCexplorer2 pipeline for a single sample.

    Parameters
    ----------
    sample:
        Mapping with keys: "cell_id", "r1", "r2".
    outdir_root:
        Output root directory for this detector.
    reference_fa:
        Reference FASTA used for STAR alignment.
    star_index:
        STAR genome index directory. If None/empty and skip_star is False,
        a RuntimeError is raised.
    ref_flat:
        refFlat annotation file required by CIRCexplorer2.
    """
    cell_id = sample["cell_id"]
    fq1 = sample["r1"]
    fq2 = sample.get("r2", "")

    outdir_root_path = Path(outdir_root)
    cell_dir = outdir_root_path / cell_id
    cell_dir.mkdir(parents=True, exist_ok=True)

    # Normalise inputs
    star_index = (star_index or "").strip()
    ref_flat = (ref_flat or "").strip()

    if not skip_star and not star_index:
        raise RuntimeError(
            "run_circexplorer2: STAR index directory is empty.\n"
            "Check CIRCYTO_CIRCEXPLORER2_STAR_INDEX or the detector configuration."
        )

    if not skip_star and not ref_flat:
        raise RuntimeError(
            "run_circexplorer2: refFlat annotation is empty.\n"
            "Check CIRCYTO_CIRCEXPLORER2_REF_FLAT or the detector configuration."
        )

    chimeric_junction = cell_dir / f"{cell_id}_Chimeric.out.junction"
    star_prefix = cell_dir / f"{cell_id}_"

    if not skip_star:
        cmd_star = [
            "STAR",
            "--runThreadN",
            str(threads),
            "--genomeDir",
            star_index,
            "--readFilesIn",
            fq1,
            fq2,
            "--readFilesCommand",
            "zcat",
            "--outFileNamePrefix",
            str(star_prefix),
            "--outSAMtype",
            "BAM",
            "Unsorted",
            "--chimSegmentMin",
            "10",
            "--chimJunctionOverhangMin",
            "10",
            "--chimOutType",
            "Junctions",
        ]

        if star_tmp_base:
            tmp_dir = Path(star_tmp_base) / cell_id
            tmp_dir.parent.mkdir(parents=True, exist_ok=True)
            cmd_star.extend(["--outTmpDir", str(tmp_dir)])

        _run(cmd_star)
    else:
        print(
            f"[circexplorer2] Skipping STAR for {cell_id} "
            f"(skip_star={skip_star}), expecting existing junction file at {chimeric_junction}"
        )

    # ---- CIRCexplorer2 parsing & annotation ----
    parse_out = cell_dir / f"{cell_id}_CIRCexplorer2_circ.txt"

    cmd_parse = [
        "CIRCexplorer2",
        "parse",
        "-t",
        "STAR",
        "-b",
        str(chimeric_junction),
        "-o",
        str(parse_out),
    ]
    if extra_args:
        cmd_parse.extend(extra_args)

    _run(cmd_parse)

    _clean_parse_output(parse_out)

    # Annotation step
    cmd_annot = [
        "CIRCexplorer2",
        "annotate",
        "-r",
        ref_flat,
        "-g",
        reference_fa,
        "-b",
        str(parse_out),
        "-o",
        str(parse_out),
    ]
    _run(cmd_annot)
