"""Synthetic completion-binding records; a file alone must never imply success."""

import json
from pathlib import Path
import tempfile
import unittest

from execution_records import alignment_output_binding, bind_native
from run_stage import digest, dump


class BindingChecks(unittest.TestCase):
    def test_pipe_binding_requires_fresh_manifest_exit_and_unchanged_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sam = root / "workflow_firstcell/align/bwa_mem/cache/test/alignment.sam"
            sam.parent.mkdir(parents=True)
            sam.write_text("synthetic SAM bytes\n")
            (sam.parent / "IMR90_A_100.align.log").write_text(
                "[main] Real time: 0.001 sec\n"
            )
            manifest = sam.parents[2] / "alignment_manifest.tsv"
            header = "cell_id\taligner\tread_layout\treused_alignment\tgroup_id\treference\tsam\n"
            row = (
                "\t".join(
                    [
                        "IMR90_A_100",
                        "bwa-mem",
                        "single-end",
                        "false",
                        "SRX26321174",
                        str(root / "ref/GRCh38.primary_assembly.genome.fa"),
                        str(sam),
                    ]
                )
                + "\n"
            )
            manifest.write_text(header + row)
            alignment = {
                "exit_code": 0,
                "completion_bound": False,
                "output_path": None,
                "inputs_unchanged": True,
            }
            detector = {
                "argv": ["java", "-I", str(sam)],
                "input_sha256": {str(sam): digest(sam)},
            }
            bound = alignment_output_binding(root, alignment, detector)
            self.assertEqual(bound["output_sha256"], digest(sam))
            self.assertFalse(bound["wrapper_direct_output_hash_available"])
            with self.assertRaisesRegex(ValueError, "nonzero exit"):
                alignment_output_binding(root, dict(alignment, exit_code=7), detector)
            manifest.write_text(header + row.replace("\tfalse\t", "\ttrue\t"))
            with self.assertRaisesRegex(ValueError, "identity/settings disagree"):
                alignment_output_binding(root, alignment, detector)
            manifest.write_text(header + row)
            sam.write_text("changed bytes\n")
            with self.assertRaisesRegex(ValueError, "hash link broken"):
                alignment_output_binding(root, alignment, detector)

    def test_failed_partial_then_explicit_success_binding(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            names = [
                "env-bwa/bin/bwa",
                "env/bin/java",
                "ref/GRCh38.primary_assembly.genome.fa",
                "ref/gencode.v45.primary_assembly.annotation.gtf",
                "inputs/fastq/SRR30918126.fastq.gz",
                "external/CIRI3_Java_18.0.1.jar",
                "workflow_firstcell/align/first.sam",
                "workflow_firstcell/ciri3/IMR90_A_100.ciri3_run/ciri3_raw.tsv",
            ]
            for name in names:
                p = root / name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("synthetic unit fixture\n")
            bwa, java, fa, gtf, fq, jar, sam, raw = [str(root / n) for n in names]
            env = {
                "external_executables": {
                    "bwa": {"path": bwa, "sha256": digest(bwa)},
                    "java": {"path": java, "sha256": digest(java)},
                },
                "external_jar_sha256": digest(jar),
            }
            dump(root / "records/environment_manifest.json", env)
            dump(
                root / "records/fastq_validation_header_v2.json", {"sha256": digest(fq)}
            )
            dump(
                root / "records/pilot_progress.json",
                {"cells": [{"cell_id": "IMR90_A_100"}]},
            )
            alignment = {
                "kind": "alignment",
                "completion_bound": True,
                "exit_code": 0,
                "output_path": sam,
                "output_sha256": digest(sam),
                "input_sha256": {fq: digest(fq), fa: digest(fa)},
                "executable_sha256": digest(bwa),
                "argv": [bwa, "mem", "-t", "8", "-k", "15", "-T", "15", fa, fq],
            }
            detector = {
                "kind": "detector",
                "completion_bound": False,
                "exit_code": 7,
                "output_path": raw,
                "output_sha256": digest(raw),
                "input_sha256": {p: digest(p) for p in [sam, fa, gtf, jar]},
                "executable_sha256": digest(java),
                "argv": [
                    java,
                    "-jar",
                    jar,
                    "-I",
                    sam,
                    "-O",
                    raw,
                    "-F",
                    fa,
                    "-T",
                    "8",
                    "-G",
                    str(Path(raw).parent / "ciri3_internal.log"),
                    "-A",
                    gtf,
                    "-S",
                    "0",
                    "-Ma",
                    "0",
                ],
            }
            dump(root / "records/tool_commands/alignment.json", alignment)
            dump(root / "records/tool_commands/detector.json", detector)
            with self.assertRaisesRegex(ValueError, "completion cannot be bound"):
                bind_native(root)
            self.assertFalse((root / "inputs/cells_firstcell.tsv").exists())
            detector.update(completion_bound=True, exit_code=0)
            dump(root / "records/tool_commands/detector.json", detector)
            bind_native(root)
            binding = json.loads(
                (root / "records/native_completion_binding.json").read_text()
            )
            self.assertEqual(binding["native_sha256"], digest(raw))
            self.assertEqual(binding["original_detector_status"], "unknown")
            self.assertEqual(binding["new_detector_status"], "completed")


if __name__ == "__main__":
    unittest.main()
