# circyto/detectors/find_circ3.py
from __future__ import annotations

from pathlib import Path
import shutil
from typing import List

from .base import (
    DetectorBase,
    DetectorCapabilities,
    DetectorRunInputs,
    DetectorResult,
)
from .find_circ3_adapter import run_find_circ3

REQUIRED_FIND_CIRC3_COMMANDS = ("bowtie2", "samtools", "find-circ3")


class FindCirc3Detector(DetectorBase):
    """
    Detector wrapper for the find_circ3 short-read pipeline.
    """

    # Required by tests
    name: str = "find-circ3"
    input_type: str = "short-read"
    supports_paired_end: bool = True
    capabilities: DetectorCapabilities = DetectorCapabilities(
        accepts_fastq=True,
        accepts_alignment=False,
        prefers_paired=True,
        supports_single_end=False,
        supports_multisample_alignment=False,
        max_parallel=2,
        recommended_execution_mode="per-cell-fastq",
    )
    max_parallel: int = 2

    def missing_dependencies(self) -> list[str]:
        return [cmd for cmd in REQUIRED_FIND_CIRC3_COMMANDS if shutil.which(cmd) is None]

    def is_available(self) -> bool:
        """
        Return True only when the external find-circ3 runtime is on PATH.
        """
        return not self.missing_dependencies()

    def version(self) -> str:
        """
        Basic version; we can refine by probing `find-circ3 --version`.
        """
        return "unknown"

    def run(self, inputs: DetectorRunInputs) -> List[DetectorResult]:
        """
        Required by tests: a single-sample detector interface.

        This mirrors the interface of CIRI-full's DetectorBase.run().
        """
        run_find_circ3(
            sample={
                "cell_id": inputs.cell_id,
                "r1": str(inputs.r1),
                "r2": str(inputs.r2) if inputs.r2 else None,
            },
            reference_fa=str(inputs.ref_fa),
            outdir_root=str(inputs.outdir),
            threads=inputs.threads,
        )

        out_bed = (
            Path(inputs.outdir)
            / inputs.cell_id
            / f"{inputs.cell_id}_splice_sites.bed"
        )

        return [
            DetectorResult(
                detector=self.name,
                cell_id=inputs.cell_id,
                outdir=inputs.outdir,
                tsv_path=out_bed,
            )
        ]

    def run_from_fastq(self, inputs: DetectorRunInputs) -> List[DetectorResult]:
        return self.run(inputs)

    def run_manifest(self, inputs: DetectorRunInputs) -> List[DetectorResult]:
        """
        Manifest-level interface for multi-sample runs.
        """
        results: List[DetectorResult] = []

        for sample in inputs.iter_samples():
            run_find_circ3(
                sample=sample,
                reference_fa=str(inputs.ref_fa),
                outdir_root=str(inputs.outdir),
                threads=inputs.threads,
            )

            out_bed = (
                Path(inputs.outdir)
                / sample["cell_id"]
                / f"{sample['cell_id']}_splice_sites.bed"
            )

            results.append(
                DetectorResult(
                    detector=self.name,
                    cell_id=sample["cell_id"],
                    outdir=inputs.outdir,
                    tsv_path=out_bed,
                )
            )

        return results
