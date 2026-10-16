"""Fake-browser orchestration test; no biological or browser-rendering claim."""

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import browser_reports


class BrowserOrchestration(unittest.TestCase):
    def test_offline_arguments_and_two_distinct_report_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = "<html><body><table><tr><th>synthetic</th></tr><tr><td>7</td></tr></table></body></html>"
            for name in ["workflow_remaining3", "import_four"]:
                path = root / name / "qc/report.html"
                path.parent.mkdir(parents=True)
                path.write_text(fixture)
            fake = root / "fake-chrome"
            fake.write_text(
                '#!/usr/bin/python3\nimport pathlib,sys,urllib.parse\np=pathlib.Path(urllib.parse.unquote(urllib.parse.urlparse(sys.argv[-1]).path));sys.stdout.write(p.read_text());s=next(x.split("=",1)[1] for x in sys.argv if x.startswith("--screenshot="));pathlib.Path(s).write_bytes(b"synthetic-test-image-placeholder")\n'
            )
            fake.chmod(0o755)
            with (
                patch.object(sys, "argv", ["browser_reports.py", "--root", str(root),
                    "--workflow-directory", "workflow_remaining3", "--import-directory", "import_four"]),
                patch.object(browser_reports.shutil, "which", return_value=str(fake)),
                contextlib.redirect_stdout(io.StringIO()),
            ):
                browser_reports.main()
            result = json.loads(
                (root / "records/acceptance/browser_acceptance.json").read_text()
            )
            self.assertTrue(result["complete"])
            self.assertEqual(
                [x["report"] for x in result["results"]],
                ["new_workflow", "native_import"],
            )
            for item in result["results"]:
                self.assertIn("--host-resolver-rules=MAP * ~NOTFOUND", item["argv"])
                self.assertEqual(item["tables_equal"], 1)


if __name__ == "__main__":
    unittest.main()
