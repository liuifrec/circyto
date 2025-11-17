# circyto/detectors/__init__.py
from __future__ import annotations

from typing import Dict

from .base import DetectorBase, DetectorRunInputs, DetectorResult
from .ciri_full import CiriFullDetector


def available_detectors() -> Dict[str, DetectorBase]:
    """
    Return a mapping from detector name -> detector engine for all
    detectors that appear to be available in the current environment.

    For now this only includes CIRI-full, but v0.6.0+ will extend this
    to CIRI-long, find_circ, CIRCexplorer2, etc.
    """
    engines: Dict[str, DetectorBase] = {}

    cirifull = CiriFullDetector()
    if cirifull.is_available():
        engines[cirifull.name] = cirifull

    return engines


__all__ = [
    "DetectorBase",
    "DetectorRunInputs",
    "DetectorResult",
    "CiriFullDetector",
    "available_detectors",
]
