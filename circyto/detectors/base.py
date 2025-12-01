from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Protocol, runtime_checkable, Literal

DetectorInputType = Literal["fastq", "bam"]


@dataclass
class DetectorRunInputs:
    """
    Standardized inputs for running a detector on a single cell/sample.
    This is intentionally minimal and reusable across all detectors.
    """

    # Core identifiers
    cell_id: str
    outdir: Path

    # Input files
    r1: Optional[Path] = None
    r2: Optional[Path] = None
    bam: Optional[Path] = None

    # Reference / annotation
    ref_fa: Optional[Path] = None
    gtf: Optional[Path] = None

    # Execution parameters
    threads: int = 1

    # Optional extras (detector-specific flags, env, etc.)
    # Tests expect this to exist and to default to an empty dict.
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class DetectorResult:
    """
    Standardized result container for a single detector on a single cell/sample.
    """

    detector: str
    cell_id: str
    outdir: Path
    tsv_path: Path

    # Optional paths for debugging / logging / provenance
    run_dir: Optional[Path] = None
    log_path: Optional[Path] = None

    # Arbitrary metadata (runtime, version info, custom stats, ...)
    # Tests expect this to be a dict by default, not None.
    meta: Dict[str, Any] = field(default_factory=dict)


@runtime_checkable
class DetectorBase(Protocol):
    """
    Protocol that all detector engines should satisfy.

    This keeps the detector API lightweight and testable while allowing us to
    plug in different detector implementations (CIRI-full, CIRI2, find_circ3,
    CIRCexplorer2, etc.) behind a common interface.
    """

    name: str
    input_type: DetectorInputType
    supports_paired_end: bool

    # Maximum parallel jobs this detector can realistically handle at once.
    # Used by run_detector_manifest / run_multidetector to avoid oversubscription.
    max_parallel: int

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        """
        Run the detector on a single cell/sample.

        Implementations should:
        - Create any necessary output directories under inputs.outdir
        - Write a per-cell TSV (or BED-like) circRNA call file
        - Populate DetectorResult.tsv_path with that file
        - Optionally populate run_dir, log_path, meta
        """
        ...
