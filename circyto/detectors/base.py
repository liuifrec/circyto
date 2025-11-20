# circyto/detectors/base.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, Literal, Optional, Protocol, runtime_checkable


DetectorInputType = Literal["fastq", "bam"]


@dataclass
class DetectorRunInputs:
    """
    Standardized inputs for running a circRNA detector on a single cell/sample.
    This is detector-agnostic and can be reused by any engine.
    """

    cell_id: str
    r1: Path
    r2: Optional[Path]
    outdir: Path
    ref_fa: Optional[Path] = None
    gtf: Optional[Path] = None
    threads: int = 8
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class DetectorResult:
    """
    Standardized outputs from a detector run on a single cell/sample.
    tsv_path should point to a CIRI-like circRNA table compatible with `collect`.
    """

    detector: str
    cell_id: str
    outdir: Path
    tsv_path: Path

    # Optional extras for logging / debugging / downstream use
    run_dir: Optional[Path] = None
    log_path: Optional[Path] = None
    meta: Dict[str, Any] = field(default_factory=dict)


@runtime_checkable
class DetectorBase(Protocol):
    """
    Minimal interface that all detector engines must implement.

    This is purposefully small – we can add methods later without breaking
    early adapters. The key design goal is:

      DetectorRunInputs -> DetectorResult

    so that higher-level orchestration code can be detector-agnostic.
    """

    #: short, CLI-friendly name, e.g. "ciri-full", "ciri-long", "find_circ"
    name: str

    #: type of primary input this detector expects ("fastq", "bam", ...)
    input_type: DetectorInputType

    #: whether the detector supports paired-end reads natively
    supports_paired_end: bool

    def is_available(self) -> bool:
        """
        Lightweight check – return True if the detector can be run
        in the current environment (binaries / jars are present).
        """
        ...

    def version(self) -> Optional[str]:
        """
        Return a human-readable version string, or None if unknown.
        This is mainly for logging / reproducibility.
        """
        ...

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        """
        Run the detector for a single cell and return a DetectorResult.

        Implementations are free to create subdirectories, temporary files, etc.,
        but they must guarantee that:

          - `inputs.outdir` exists
          - `DetectorResult.tsv_path` points to the final circRNA table
            compatible with `collect` (or a small adapter).
        """
        ...
