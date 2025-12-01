"""
CIRCexplorer2 detector wrapper for the unified detector API.

This class adapts the low-level `circexplorer2_adapter.run_circexplorer2`
function to the generic DetectorBase protocol used by `run_detector`
and `run_multidetector`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from shutil import which
from typing import Optional

from .base import DetectorRunInputs, DetectorResult
from .circexplorer2_adapter import run_circexplorer2


@dataclass
class CircExplorer2Detector:
    """
    High-level CIRCexplorer2 detector implementing the DetectorBase protocol.

    Parameters
    ----------
    star_index:
        Path to the STAR genome index directory used by CIRCexplorer2.
    ref_flat:
        Path to the refFlat annotation file.
    star_tmp_base:
        Optional base directory for STAR temporary folders.
    skip_star:
        If True, skip the STAR step and assume STAR outputs are present.
    max_parallel:
        Soft hint for maximum parallelism; the pipeline enforces this.
    """

    star_index: Optional[str] = None
    ref_flat: Optional[str] = None
    star_tmp_base: Optional[str] = None
    skip_star: bool = False
    max_parallel: int = 1

    # Static metadata expected by the DetectorBase protocol
    name: str = "circexplorer2"
    input_type: str = "fastq"
    supports_paired_end: bool = True

    # ------------------------------------------------------------------ #
    # Protocol methods                                                   #
    # ------------------------------------------------------------------ #

    def is_available(self) -> bool:
        """
        Return True if CIRCexplorer2 appears to be usable on this system.

        We keep this deliberately lightweight: just check whether the
        CIRCexplorer2 and STAR binaries are on PATH (when STAR is required).
        """
        if self.skip_star:
            # Only require CIRCexplorer2 itself.
            return which("CIRCexplorer2") is not None

        # Need both STAR and CIRCexplorer2 to be present.
        return which("STAR") is not None and which("CIRCexplorer2") is not None

    def version(self) -> str | None:
        """
        Optional version hook. For now we don't attempt to parse the version
        string and just return None.
        """
        return None

    # ------------------------------------------------------------------ #
    # Core run logic                                                     #
    # ------------------------------------------------------------------ #

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        """
        Execute CIRCexplorer2 for a single cell/sample described by inputs.

        This method is called by `run_detector_manifest` for each row in the
        manifest.
        """
        # Sanity checks
        if inputs.r1 is None or inputs.r2 is None:
            raise ValueError("CircExplorer2Detector requires paired-end FASTQ (r1 and r2).")

        outdir_root: Path = inputs.outdir
        outdir_root.mkdir(parents=True, exist_ok=True)

        cell_dir = outdir_root / inputs.cell_id
        cell_dir.mkdir(parents=True, exist_ok=True)

        circ_out = cell_dir / f"{inputs.cell_id}_CIRCexplorer2_circ.txt"

        # Delegate to the low-level adapter
        run_circexplorer2(
            sample={
                "cell_id": inputs.cell_id,
                "r1": str(inputs.r1),
                "r2": str(inputs.r2),
            },
            reference_fa=str(inputs.ref_fa),
            outdir_root=str(outdir_root),
            star_index=self.star_index or "",
            ref_flat=self.ref_flat or "",
            star_tmp_base=self.star_tmp_base,
            skip_star=self.skip_star,
            threads=inputs.threads,
            extra_args=inputs.extra,
        )

        return DetectorResult(
            detector=self.name,
            cell_id=inputs.cell_id,
            outdir=outdir_root,
            tsv_path=circ_out,
            run_dir=cell_dir,
            log_path=None,
            meta={},
        )
