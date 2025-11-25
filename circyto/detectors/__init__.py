# circyto/detectors/__init__.py
from __future__ import annotations

from typing import Dict

from .base import DetectorBase, DetectorRunInputs, DetectorResult
from .ciri_full import CiriFullDetector
from .ciri2 import Ciri2Detector
from .find_circ3 import FindCirc3Detector


def available_detectors() -> Dict[str, DetectorBase]:
    """
    Return a mapping from detector name -> detector engine for all
    detectors that appear to be available in the current environment.
    """
    engines: Dict[str, DetectorBase] = {}

    cirifull = CiriFullDetector()
    if cirifull.is_available():
        engines[cirifull.name] = cirifull

    ciri2 = Ciri2Detector()
    if ciri2.is_available():
        engines[ciri2.name] = ciri2

    findcirc3 = FindCirc3Detector()
    if findcirc3.is_available():
        engines[findcirc3.name] = findcirc3

    return engines


__all__ = [
    "DetectorBase",
    "DetectorRunInputs",
    "DetectorResult",
    "CiriFullDetector",
    "Ciri2Detector",
    "FindCirc3Detector",
    "available_detectors",
    "build_default_engines",
]


def build_default_engines() -> Dict[str, DetectorBase]:
    """
    Return the default detector engines keyed by CLI name.
    """
    return {
        "ciri-full": CiriFullDetector(),
        "ciri2": Ciri2Detector(),
        "find-circ3": FindCirc3Detector(),
    }
