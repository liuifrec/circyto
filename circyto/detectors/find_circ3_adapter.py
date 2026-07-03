from __future__ import annotations

import subprocess
from pathlib import Path
from typing import IO, Mapping, Optional, Sequence


def _run(cmd: Sequence[str], *, stdout: IO[bytes] | None = None, stderr: IO[str] | None = None) -> None:
    """
    Small wrapper for subprocess.run with basic logging and no shell interpolation.
    """
    printable = " ".join(map(str, cmd))
    print(f"[find_circ3] {printable}")
    subprocess.run(list(map(str, cmd)), check=True, stdout=stdout, stderr=stderr)


def _run_pipeline(
    left_cmd: Sequence[str],
    right_cmd: Sequence[str],
    *,
    stdout_path: Path,
    log_path: Path,
) -> None:
    """
    Run left_cmd | right_cmd > stdout_path without shell=True.
    """
    print(f"[find_circ3] {' '.join(map(str, left_cmd))} | {' '.join(map(str, right_cmd))} > {stdout_path}")
    with log_path.open("w", encoding="utf-8") as log, stdout_path.open("wb") as out:
        left = subprocess.Popen(
            list(map(str, left_cmd)),
            stdout=subprocess.PIPE,
            stderr=log,
        )
        assert left.stdout is not None
        right = subprocess.run(
            list(map(str, right_cmd)),
            stdin=left.stdout,
            stdout=out,
            stderr=log,
            check=False,
        )
        left.stdout.close()
        left_returncode = left.wait()
        if left_returncode != 0:
            raise subprocess.CalledProcessError(left_returncode, list(map(str, left_cmd)))
        if right.returncode != 0:
            raise subprocess.CalledProcessError(right.returncode, list(map(str, right_cmd)))


def _require_existing_file(path: str | Path, *, label: str) -> Path:
    value = Path(path)
    if not value.exists():
        raise FileNotFoundError(f"find-circ3 {label} not found: {value}")
    if not value.is_file():
        raise ValueError(f"find-circ3 {label} is not a file: {value}")
    return value


def _bowtie2_index_from_fasta(reference_fa: str) -> str:
    """
    Derive a bowtie2 index prefix from a reference FASTA path.

    Example:
      reference_fa = 'ref/chr21.fa'  ->  bowtie2 index prefix 'ref/chr21'
    """
    p = Path(reference_fa)
    if p.suffix.lower() in {".fa", ".fasta", ".fna"}:
        return str(p.with_suffix(""))
    return str(p)


def run_find_circ3(
    sample: Mapping[str, str],
    reference_fa: str,
    outdir_root: str,
    threads: int = 4,
    extra_args: Optional[list[str]] = None,
) -> Path:
    """
    Run the full find_circ3 short-read pipeline for a single sample,
    consistent with the updated find_circ3 README:

        Step 1: bowtie2 first pass -> sample.aln.bam
        Step 1b: samtools view -b -f 4 -> sample.unmapped.bam
        Step 2: find-circ3 anchors sample.unmapped.bam > sample.anchors.fastq
        Step 3: bowtie2 second pass on anchors.fastq -> sample.anchors.sam
        Step 4: find-circ3 call sample.anchors.sam ... > sample_splice_sites.bed

    Expected manifest columns:
        cell_id, r1, [r2]

    Args
    ----
    sample : Mapping[str, str]
        Row from manifest with keys "cell_id", "r1", and optional "r2".
    reference_fa : str
        Path to genome FASTA used by bowtie2 (index prefix will be derived).
    outdir_root : str
        Root output directory for this detector.
    threads : int
        Number of threads for bowtie2.
    extra_args : list[str] | None
        Extra CLI arguments to append to `find-circ3 call` if needed.

    Returns
    -------
    Path
        Path to the final splice_sites BED file.
    """
    cell_id = str(sample["cell_id"])
    fq1 = _require_existing_file(str(sample["r1"]), label="R1 FASTQ")
    fq2 = _require_existing_file(str(sample["r2"]), label="R2 FASTQ") if sample.get("r2") else None
    reference_path = _require_existing_file(reference_fa, label="reference FASTA")

    outdir = Path(outdir_root) / cell_id
    outdir.mkdir(parents=True, exist_ok=True)

    aln_bam = outdir / f"{cell_id}.aln.bam"
    unmapped_bam = outdir / f"{cell_id}.unmapped.bam"
    anchors_fq = outdir / f"{cell_id}.anchors.fastq"
    anchors_sam = outdir / f"{cell_id}.anchors.sam"
    splice_sites_bed = outdir / f"{cell_id}_splice_sites.bed"
    firstpass_log = outdir / f"{cell_id}_firstpass.log"
    anchors_log = outdir / f"{cell_id}_anchors.log"
    secondpass_log = outdir / f"{cell_id}_secondpass.log"
    call_log = outdir / f"{cell_id}_call.log"

    bowtie2_index = _bowtie2_index_from_fasta(str(reference_path))

    # --- Step 1: first-pass mapping (FASTQ -> aln.bam) ---
    if fq2:
        # Paired-end
        cmd_align = [
            "bowtie2",
            "-x",
            bowtie2_index,
            "-1",
            str(fq1),
            "-2",
            str(fq2),
            "--very-sensitive",
            "-p",
            str(threads),
        ]
    else:
        # Single-end
        cmd_align = [
            "bowtie2",
            "-x",
            bowtie2_index,
            "-U",
            str(fq1),
            "--very-sensitive",
            "-p",
            str(threads),
        ]
    _run_pipeline(cmd_align, ["samtools", "view", "-bS", "-"], stdout_path=aln_bam, log_path=firstpass_log)

    # --- Step 1b: extract unmapped reads -> unmapped.bam ---
    with (outdir / f"{cell_id}_samtools_unmapped.log").open("w", encoding="utf-8") as log:
        _run(
            [
                "samtools",
                "view",
                "-b",
                "-f",
                "4",
                "-o",
                str(unmapped_bam),
                str(aln_bam),
            ],
            stderr=log,
        )

    # --- Step 2: generate anchors (FASTQ) using Python 3 unmapped2anchors3 ---
    # Recommended grouped CLI:
    #   find-circ3 anchors sample_unmapped.bam --anchor 20 --min-qual 5 > anchors.fastq
    with anchors_fq.open("wb") as out, anchors_log.open("w", encoding="utf-8") as log:
        _run(
            [
                "find-circ3",
                "anchors",
                str(unmapped_bam),
                "--anchor",
                "20",
                "--min-qual",
                "5",
            ],
            stdout=out,
            stderr=log,
        )

    # --- Step 3: map anchors with bowtie2 (anchors.fastq -> anchors.sam) ---
    # As per README:
    # bowtie2 -q -U sample_anchors.fastq -x genome_index --reorder --mm --very-sensitive \
    #   --score-min C,-15,0 2> sample_secondpass.log > sample_anchors.sam
    with anchors_sam.open("wb") as out, secondpass_log.open("w", encoding="utf-8") as log:
        _run(
            [
                "bowtie2",
                "-q",
                "-U",
                str(anchors_fq),
                "-x",
                bowtie2_index,
                "--reorder",
                "--mm",
                "--very-sensitive",
                "--score-min",
                "C,-15,0",
                "-p",
                str(threads),
            ],
            stdout=out,
            stderr=log,
        )

    # --- Step 4: call circRNAs from anchors.sam ---
    # README example:
    # find-circ3 call sample_anchors.sam \
    #   --genome genome.fa \
    #   --name sample1 \
    #   --prefix sample1_ \
    #   --anchor 20 \
    #   --min-uniq-qual 2 \
    #   --max-mismatches 2 \
    #   --margin 5 \
    #   --strandpref \
    #   --stats sample.log \
    #   --reads sample_reads.fa \
    #   > sample_splice_sites.bed
    # Minimal, tested-good CLI based on standalone find-circ3 tests
    call_cmd = [
        "find-circ3",
        "call",
        str(anchors_sam),
        "--genome",
        str(reference_path),
        "--name",
        cell_id,
        "--prefix",
        f"{cell_id}_",
        "--anchor",
        "20",
    ]

    # Allow the detector interface to inject extra options if needed
    if extra_args:
        call_cmd.extend(extra_args)

    with splice_sites_bed.open("wb") as out, call_log.open("w", encoding="utf-8") as log:
        _run(call_cmd, stdout=out, stderr=log)

    return splice_sites_bed


