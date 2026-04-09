from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Protocol, runtime_checkable, Literal

DetectorInputType = Literal["fastq", "bam", "sam", "alignment", "short-read"]
DetectorExecutionMode = Literal["per-cell-fastq", "alignment-first", "multisample-alignment"]
DetectorInputMode = Literal["fastq", "alignment"]


@dataclass(frozen=True)
class DetectorCapabilities:
    accepts_fastq: bool = True
    accepts_alignment: bool = False
    prefers_paired: bool = True
    supports_single_end: bool = True
    supports_multisample_alignment: bool = False
    requires_unsorted_sam: bool = False
    supports_star: bool = False
    supports_bwa: bool = False
    max_parallel: int = 1
    recommended_execution_mode: DetectorExecutionMode = "per-cell-fastq"


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
    sam: Optional[Path] = None
    input_mode: DetectorInputMode = "fastq"
    read_layout: Optional[str] = None
    alignment_group: Optional[str] = None

    # Reference / annotation
    ref_fa: Optional[Path] = None
    gtf: Optional[Path] = None

    # Execution parameters
    threads: int = 1

    # Optional extras (detector-specific flags, env, etc.)
    # Tests expect this to exist and to default to an empty dict.
    extra: Dict[str, Any] = field(default_factory=dict)
    provenance: Dict[str, Any] = field(default_factory=dict)

    @property
    def alignment_path(self) -> Optional[Path]:
        return self.bam or self.sam

    def effective_read_layout(self) -> str:
        if self.read_layout:
            return self.read_layout
        return "paired-end" if self.r2 is not None else "single-end"

    def iter_samples(self) -> Iterable[dict[str, Any]]:
        sample = {"cell_id": self.cell_id, "r1": self.r1, "r2": self.r2, "bam": self.bam, "sam": self.sam}
        return [sample]


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
    capabilities: DetectorCapabilities

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


def get_detector_capabilities(detector: Any) -> DetectorCapabilities:
    caps = getattr(detector, "capabilities", None)
    if isinstance(caps, DetectorCapabilities):
        return caps

    max_parallel = int(getattr(detector, "max_parallel", 1) or 1)
    input_type = str(getattr(detector, "input_type", "fastq"))
    accepts_alignment = input_type in {"bam", "sam", "alignment"}
    accepts_fastq = input_type in {"fastq", "short-read"} or not accepts_alignment
    supports_paired_end = bool(getattr(detector, "supports_paired_end", True))
    return DetectorCapabilities(
        accepts_fastq=accepts_fastq,
        accepts_alignment=accepts_alignment,
        prefers_paired=supports_paired_end,
        supports_single_end=True,
        supports_multisample_alignment=False,
        max_parallel=max_parallel,
        recommended_execution_mode="alignment-first" if accepts_alignment and not accepts_fastq else "per-cell-fastq",
    )
