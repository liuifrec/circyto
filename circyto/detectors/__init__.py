"""
Detector registry and convenience helpers.

This module wires together the concrete detector implementations and exposes:

- DetectorRunInputs, DetectorResult, DetectorBase  (re-exported from .base)
- available_detectors()  → mapping of detector-name → detector instance
- build_default_engines()  → same as above, used by CLI / pipelines
"""

from __future__ import annotations

import os
from typing import Dict

from .base import DetectorRunInputs, DetectorResult, DetectorBase
from .ciri_full import CiriFullDetector
from .ciri2 import Ciri2Detector
from .find_circ3 import FindCirc3Detector

try:
    # Optional: CIRCexplorer2 support
    from .circ_explorer2 import CircExplorer2Detector  # type: ignore[attr-defined]
except Exception:  # pragma: no cover - if module missing, we just skip it
    CircExplorer2Detector = None  # type: ignore[assignment]


def _include_circexplorer2() -> bool:
    """
    Decide whether to include CIRCexplorer2 in the default engine set.

    Controlled by environment variable:

      CIRCYTO_CIRCEXPLORER2_SKIP=1  → do NOT include
      (anything else / unset)       → include if import succeeded
    """
    skip_env = os.environ.get("CIRCYTO_CIRCEXPLORER2_SKIP", "0").strip().lower()
    if skip_env in {"1", "true", "yes"}:
        return False
    return CircExplorer2Detector is not None


def build_default_engines() -> Dict[str, DetectorBase]:
    """
    Instantiate the built-in detectors and return a mapping from detector name
    to detector *instance*.

    These instances are used by:
      - CLI: `run-detector`, `run-multidetector`
      - Internal pipelines: `run_detector_manifest`, `run_multidetector_pipeline`
    """
    engines: Dict[str, DetectorBase] = {
        "ciri-full": CiriFullDetector(),
        "ciri2": Ciri2Detector(),
        "find-circ3": FindCirc3Detector(),
    }

    # CIRCexplorer2 is optional and configured via env inside its run() logic.
    if _include_circexplorer2():
        engines["circexplorer2"] = CircExplorer2Detector()  # type: ignore[call-arg]

    return engines


def available_detectors() -> Dict[str, DetectorBase]:
    """
    Public helper used in tests and CLI to list detectors.

    Returns the same mapping as build_default_engines(), but is kept as a thin
    wrapper so the test API stays stable.
    """
    return build_default_engines()


__all__ = [
    "DetectorRunInputs",
    "DetectorResult",
    "DetectorBase",
    "available_detectors",
    "build_default_engines",
]
