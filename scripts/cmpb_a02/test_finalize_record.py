"""Audit finalization must reject unperformed acceptance and active jobs."""

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import finalize_record


class FinalizationGates(unittest.TestCase):
    def test_unexecuted_acceptance_cannot_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "records").mkdir()
            (root / "records/control_tests.json").write_text("[]\n")
            with patch.object(
                sys,
                "argv",
                [
                    "finalize_record.py",
                    "--root",
                    str(root),
                    "--output",
                    str(root / "small"),
                    "--outcome",
                    "passed",
                    "--reason",
                    "synthetic unit test",
                ],
            ):
                with self.assertRaisesRegex(ValueError, "required actual stage"):
                    finalize_record.main()
            self.assertFalse((root / "records/pilot_progress.json").exists())

    def test_active_stage_cannot_be_finalized(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "records").mkdir()
            (root / "records/active.json").write_text(
                json.dumps(
                    {
                        "schema": "circyto.cmpb.a02.stage.v1",
                        "started_utc": "2026-10-15",
                        "stage": "synthetic",
                        "status": "launching",
                    }
                )
            )
            with patch.object(
                sys,
                "argv",
                [
                    "finalize_record.py",
                    "--root",
                    str(root),
                    "--output",
                    str(root / "small"),
                    "--outcome",
                    "blocked",
                    "--reason",
                    "synthetic unit test",
                ],
            ):
                with self.assertRaisesRegex(RuntimeError, "unfinished"):
                    finalize_record.main()


if __name__ == "__main__":
    unittest.main()
