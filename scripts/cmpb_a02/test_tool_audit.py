"""Fake executable audit tests; these outputs are explicitly synthetic."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from run_stage import digest


class AuditChecks(unittest.TestCase):
    def test_exit_binding_and_partial_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake = root / "fake-java"
            fake.write_text(
                '#!/usr/bin/python3\nimport pathlib,sys\np=pathlib.Path(sys.argv[sys.argv.index("-O")+1]);p.write_text("synthetic partial native evidence\\n");sys.exit(int(sys.argv[-1]))\n'
            )
            fake.chmod(0o755)
            for name in ["jar", "sam", "fa", "gtf"]:
                (root / name).write_text("synthetic input\n")
            for expected in [7, 0]:
                output = root / ("native-" + str(expected))
                args = [
                    "-jar",
                    str(root / "jar"),
                    "-I",
                    str(root / "sam"),
                    "-F",
                    str(root / "fa"),
                    "-A",
                    str(root / "gtf"),
                    "-O",
                    str(output),
                    str(expected),
                ]
                command = [
                    sys.executable,
                    str(Path(__file__).with_name("tool_audit.py")),
                    "--root",
                    str(root),
                    "--tool",
                    "java",
                    "--real",
                    str(fake),
                    "--",
                    *args,
                ]
                self.assertEqual(
                    subprocess.run(command, check=False).returncode, expected
                )
                record = json.loads(
                    sorted((root / "records/tool_commands").glob("*.json"))[
                        -1
                    ].read_text()
                )
                self.assertEqual(record["argv"], [str(fake), *args])
                self.assertEqual(record["output_sha256"], digest(output))
                self.assertEqual(record["completion_bound"], expected == 0)
                self.assertEqual(
                    record["input_sha256"][str(root / "sam")], digest(root / "sam")
                )


if __name__ == "__main__":
    unittest.main()
