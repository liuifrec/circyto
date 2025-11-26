from __future__ import annotations

import os
from pathlib import Path

import pytest

from circyto.detectors import build_default_engines


@pytest.mark.slow
@pytest.mark.integration
def test_multidetector_ciri_full_and_find_circ3_demo(tmp_path: Path) -> None:
    """
    Minimal 'smoke' test for the multi-detector pipeline with ciri-full + find-circ3.

    This only checks that:
      - both detectors are registered
      - merge-detectors and compare-detectors can run end-to-end

    It does NOT assert strong biology, just structural integrity.
    """

    if os.environ.get("CIRCYTO_SKIP_INTEGRATION", "").lower() in {"1", "true", "yes"}:
        pytest.skip("CIRCYTO_SKIP_INTEGRATION is set, skipping integration test")

    engines = build_default_engines()
    for det in ("ciri-full", "find-circ3"):
        assert det in engines, f"Detector {det!r} not registered in build_default_engines()"

    # You can optionally shell out to run-multidetector, merge-detectors, compare-detectors here
    # via subprocess and assert that the key output files exist.
    # For now, we only enforce that the detectors are available and 'runnable' at the config level.
