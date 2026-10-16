"""Tiny safety controls for the A100-only replay, never real detector inputs."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from build_status_wheel import inspect_members
from replay_alignment_status import ReplayGuard, copy_file, source_snapshot
import run_stage


class ReplaySafety(unittest.TestCase):
    def test_all_process_launch_boundaries_blocked(self):
        guard = ReplayGuard(Path("/original"), Path("/original/derived"))
        for event in ["subprocess.Popen", "os.system", "os.exec", "os.posix_spawn", "os.fork", "os.forkpty"]:
            with self.assertRaisesRegex(RuntimeError, "execution prohibited"):
                guard(event, ("fake command",))
        self.assertEqual(len(guard.blocked), 6)

    def test_original_mutations_blocked_derived_writes_allowed(self):
        guard = ReplayGuard(Path("/original"), Path("/original/derived"))
        guard("open", ("/original/source.h5ad", "r", 0))
        guard("open", ("/original/derived/copy.h5ad", "w", 0))
        for event, args in [
            ("open", ("/original/source.h5ad", "w", 0)),
            ("os.rename", ("/original/source", "/original/derived/copy", -1, -1)),
            ("os.remove", ("/original/source", -1)),
        ]:
            with self.assertRaisesRegex(RuntimeError, "read-only"):
                guard(event, args)

    def test_copy_is_independent_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, target = root / "source", root / "copy"
            source.write_bytes(b"synthetic fixture")
            records = []
            copy_file(source, target, records)
            with self.assertRaises(FileExistsError):
                copy_file(source, target, records)
            target.write_bytes(b"changed copy")
            self.assertEqual(source.read_bytes(), b"synthetic fixture")
            self.assertEqual(len(records), 1)

    def test_snapshot_detects_changed_bytes_and_added_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "inputs").mkdir()
            source = root / "inputs/tiny"
            source.write_bytes(b"before")
            before = source_snapshot(root)
            source.write_bytes(b"after!")
            self.assertNotEqual(before, source_snapshot(root))
            source.write_bytes(b"before")
            self.assertEqual(before, source_snapshot(root))
            (root / "inputs/extra").write_bytes(b"new")
            self.assertNotEqual(before, source_snapshot(root))

    def test_archive_inspection_rejects_all_excluded_assets(self):
        for name in ["x.jar", "x.class", "CIRI-full_manual.txt"]:
            with self.assertRaisesRegex(AssertionError, "Excluded"):
                inspect_members([name])

    def test_guard_accounts_for_original_plus_derived_disk(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ledger = root / "derived"
            (ledger / "records").mkdir(parents=True)
            (root / "retained-original").write_bytes(b"x" * 10000)
            group = root / "fake_cgroup"
            group.mkdir()
            (group / "memory.max").write_text("4096")
            (group / "memory.swap.max").write_text("0")
            spec = ledger / "spec.json"
            spec.write_text(json.dumps({"stage": "fake", "footprint_root": str(root),
                "memory_max_bytes": 4096, "argv": ["fake"], "disk_stop_bytes": 999999}))
            child = SimpleNamespace(pid=123, returncode=0, poll=lambda: 0)
            with patch.object(run_stage, "cgroup_path", return_value=group), \
                 patch.object(run_stage.os, "sched_getaffinity", return_value=set(range(8))), \
                 patch.object(run_stage.subprocess, "Popen", return_value=child):
                self.assertEqual(run_stage.guard(ledger, spec), 0)
            result = json.loads((ledger / "records/fake.guard.json").read_text())
            self.assertEqual(result["footprint_root"], str(root))
            self.assertGreaterEqual(result["disk_peak_bytes"], 10000)


if __name__ == "__main__":
    unittest.main()
