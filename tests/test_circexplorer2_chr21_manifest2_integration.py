from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest


@pytest.mark.slow
def test_circexplorer2_cli_smoke(tmp_path: Path) -> None:
    """
    Lightweight CLI smoke test for the 'circexplorer2' detector.

    This test no longer runs a full STAR+CIRCexplorer2 pipeline (which can be
    very slow and environment-dependent). Instead, it simply verifies that:

      - The circyto CLI can be invoked.
      - The 'circexplorer2' detector is recognized as a valid detector.

    The heavy integration (running STAR and CIRCexplorer2 end-to-end) is covered
    in tests/test_circexplorer2_integration.py and is allowed to skip if the
    environment is not ready.
    """

    # Allow CI or local runs to skip this if desired
    if os.environ.get("CIRCYTO_SKIP_INTEGRATION", "").lower() in {"1", "true", "yes"}:
        pytest.skip("CIRCYTO_SKIP_INTEGRATION is set, skipping CIRCexplorer2 CLI smoke test")

    repo_root = Path(__file__).resolve().parents[1]

    # Use `--help` and check the detector name appears in the help text.
    cmd = [
        "python",
        "-m",
        "circyto.cli.circyto",
        "run-detector",
        "--help",
    ]

    proc = subprocess.run(
        cmd,
        cwd=repo_root,
        capture_output=True,
        text=True,
    )

    # CLI should at least run successfully
    assert proc.returncode == 0, f"circyto CLI --help failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"

    # Help text should mention circexplorer2 as a valid detector somewhere
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert "circexplorer2" in combined, "Expected 'circexplorer2' to appear in CLI help output"
