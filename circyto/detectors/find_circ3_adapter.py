from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Mapping, Optional, Sequence


def _run(cmd: Sequence[str] | str, shell: bool = False) -> None:
    """
    Small wrapper for subprocess.run with basic logging.
    """
    if isinstance(cmd, str):
        printable = cmd
    else:
        printable = " ".join(cmd)
    print(f"[find_circ3] {printable}")
    subprocess.run(cmd, check=True, shell=shell)


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
    fq1 = str(sample["r1"])
    fq2 = sample.get("r2")

    outdir = Path(outdir_root) / cell_id
    outdir.mkdir(parents=True, exist_ok=True)

    aln_bam = outdir / f"{cell_id}.aln.bam"
    unmapped_bam = outdir / f"{cell_id}.unmapped.bam"
    anchors_fq = outdir / f"{cell_id}.anchors.fastq"
    anchors_sam = outdir / f"{cell_id}.anchors.sam"
    splice_sites_bed = outdir / f"{cell_id}_splice_sites.bed"

    bowtie2_index = _bowtie2_index_from_fasta(reference_fa)

    # --- Step 1: first-pass mapping (FASTQ -> aln.bam) ---
    if fq2:
        # Paired-end
        cmd_align = (
            f"bowtie2 "
            f"-x {bowtie2_index} "
            f"-1 {fq1} -2 {fq2} "
            f"--very-sensitive "
            f"-p {threads} "
            f"2> {outdir / (cell_id + '_firstpass.log')} "
            f"| samtools view -bS - > {aln_bam}"
        )
    else:
        # Single-end
        cmd_align = (
            f"bowtie2 "
            f"-x {bowtie2_index} "
            f"-U {fq1} "
            f"--very-sensitive "
            f"-p {threads} "
            f"2> {outdir / (cell_id + '_firstpass.log')} "
            f"| samtools view -bS - > {aln_bam}"
        )
    _run(cmd_align, shell=True)

    # --- Step 1b: extract unmapped reads -> unmapped.bam ---
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
        ]
    )

    # --- Step 2: generate anchors (FASTQ) using Python 3 unmapped2anchors3 ---
    # Recommended grouped CLI:
    #   find-circ3 anchors sample_unmapped.bam --anchor 20 --min-qual 5 > anchors.fastq
    cmd_anchors = (
        f"find-circ3 anchors {unmapped_bam} "
        f"--anchor 20 "
        f"--min-qual 5 "
        f"> {anchors_fq}"
    )
    _run(cmd_anchors, shell=True)

    # --- Step 3: map anchors with bowtie2 (anchors.fastq -> anchors.sam) ---
    # As per README:
    # bowtie2 -q -U sample_anchors.fastq -x genome_index --reorder --mm --very-sensitive \
    #   --score-min C,-15,0 2> sample_secondpass.log > sample_anchors.sam
    cmd_anchor_align = (
        f"bowtie2 "
        f"-q "
        f"-U {anchors_fq} "
        f"-x {bowtie2_index} "
        f"--reorder "
        f"--mm "
        f"--very-sensitive "
        f"--score-min C,-15,0 "
        f"-p {threads} "
        f"2> {outdir / (cell_id + '_secondpass.log')} "
        f"> {anchors_sam}"
    )
    _run(cmd_anchor_align, shell=True)

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
        str(reference_fa),
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

    call_cmd_str = " ".join(call_cmd) + f" > {splice_sites_bed}"
    _run(call_cmd_str, shell=True)

    return splice_sites_bed



