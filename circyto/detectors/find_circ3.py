# circyto/detectors/find_circ3.py
from __future__ import annotations

from pathlib import Path
from typing import List

from .base import (
    DetectorBase,
    DetectorRunInputs,
    DetectorResult,
)
from .find_circ3_adapter import run_find_circ3


class FindCirc3Detector(DetectorBase):
    """
    Detector wrapper for the find_circ3 short-read pipeline.
    """

    # Required by tests
    name: str = "find-circ3"
    input_type: str = "short-read"
    supports_paired_end: bool = True

    def is_available(self) -> bool:
        """
        Very lightweight availability check.
        We assume bowtie2, samtools, find-circ3, and find-circ3-anchors
        are on PATH if the user wants to run it.
        """
        return True

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
