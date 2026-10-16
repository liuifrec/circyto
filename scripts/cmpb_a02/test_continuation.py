"""Small synthetic scope/metadata checks, not biological observations."""

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from continuation_inputs import ROWS, compare_metadata, expected
from continuation_records import native_union
from execution_records import alignment_output_binding
from independent_checks import HEADER
from continuation_checks import oracle, Checks, display_fraction, import_source_path
from continuation_packet import sanitize
import numpy as np
from run_stage import digest
import run_stage


class ContinuationChecks(unittest.TestCase):
    def test_import_source_bytes_use_bound_RNA_input_and_copied_native(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "records").mkdir()
            rna = root / "supplied-rna.h5ad"
            rna.write_bytes(b"tiny source-byte fixture, not an H5AD")
            extraction = root / "records/rna_extraction.json"
            command = root / "records/continuation-import-four.json"
            extraction.write_text(json.dumps(dict(output=str(rna), sha256=digest(rna))))
            command.write_text(json.dumps(dict(argv=["circyto", "--rna", str(rna)])))
            self.assertEqual(import_source_path(root, "rna.h5ad"), rna)
            self.assertEqual(
                import_source_path(root, "ciri3_native_000001.tsv"),
                root / "import_four/evidence/ciri3_native_000001.tsv",
            )
            rna.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "source bytes changed"):
                import_source_path(root, "rna.h5ad")

    def test_import_RNA_source_must_match_executed_command(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "records").mkdir()
            (root / "records/rna_extraction.json").write_text(
                json.dumps(dict(output=str(root / "a.h5ad"), sha256="unused"))
            )
            command = root / "records/continuation-import-four.json"
            for argv in [["circyto"], ["circyto", "--rna", str(root / "b.h5ad")]]:
                command.write_text(json.dumps(dict(argv=argv)))
                with self.assertRaises(ValueError):
                    import_source_path(root, "rna.h5ad")

    def test_fraction_zero_padding_is_exactly_equal(self):
        for left, right in [
            ("1.0", "1.00"),
            ("1.000000", "1.00"),
            ("0.25", "0.250000"),
        ]:
            self.assertEqual(display_fraction(left), display_fraction(right))

    def test_changed_display_fraction_fails_without_tolerance(self):
        for left, right in [("0.25", "0.26"), ("1", "0.999999"), ("0.5", "0.5000001")]:
            self.assertNotEqual(display_fraction(left), display_fraction(right))

    def test_packet_paths_are_logical_in_keys_and_values(self):
        value = {
            "/private/root/continuation/input": [
                "/private/root/ref",
                "/checkout/script",
            ]
        }
        self.assertEqual(
            sanitize(
                value,
                Path("/private/root/continuation"),
                Path("/private/root"),
                Path("/checkout"),
            ),
            {
                "${CONTINUATION_ROOT}/input": [
                    "${A02_ROOT}/ref",
                    "${ANALYSIS_WORKTREE}/script",
                ]
            },
        )

    def test_oracle_matches_identity_and_feature_order_with_exact_integers(self):
        records = {
            "001": [{"circRNA_ID": "c1", "#junction_reads": str(2**53 + 1)}],
            "NA": [{"circRNA_ID": "c2", "#junction_reads": "0"}],
        }
        result = oracle(records, ["NA", "001"], ["c2", "c1"])
        self.assertEqual(result.dtype, np.int64)
        self.assertEqual(result.nnz, 1)
        self.assertEqual(int(result[1, 1]), 2**53 + 1)
        self.assertEqual(result.shape, (2, 2))

    def test_failed_check_is_persisted_before_stopping(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "checks.tsv"
            with self.assertRaises(AssertionError):
                Checks(p).check("identity", "A101", "A102", "synthetic")
            self.assertIn("FAIL", p.read_text())

    def test_union_keeps_first_seen_order_and_rejects_annotation_conflict(self):
        with tempfile.TemporaryDirectory() as d:
            a, b = Path(d) / "a.tsv", Path(d) / "b.tsv"
            row = [
                "chr1:1|4",
                "chr1",
                "1",
                "4",
                "3",
                "1",
                "0",
                "1",
                "exon",
                "G",
                "+",
                "synthetic-read",
                "1",
            ]

            def write(path, rows):
                with path.open("w", newline="") as h:
                    w = csv.writer(h, delimiter="\t")
                    w.writerow(HEADER)
                    w.writerows(rows)

            write(a, [row])
            write(b, [row[:4] + ["7"] + row[5:]])
            per, union = native_union({"a": a, "b": b})
            self.assertEqual(len(union), 1)
            self.assertEqual(per["b"][0]["#junction_reads"], "7")
            changed = row.copy()
            changed[10] = "-"
            write(b, [changed])
            with self.assertRaisesRegex(ValueError, "cross-cell conflict"):
                native_union({"a": a, "b": b})

    def test_multi_cell_pipe_binding_selects_identity_and_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            sam = root / "workflow_remaining3/align/bwa_mem/cache/tiny/alignment.sam"
            sam.parent.mkdir(parents=True)
            sam.write_text("synthetic\n")
            (sam.parent / "A102.align.log").write_text("[main] Real time: 0.001 sec\n")
            manifest = sam.parents[2] / "alignment_manifest.tsv"
            rows = [
                [c, "bwa-mem", "single-end", "false", "X" + c, "/reference", str(sam)]
                for c in ["A101", "A102", "A103"]
            ]

            def write(values):
                with manifest.open("w", newline="") as h:
                    w = csv.writer(h, delimiter="\t")
                    w.writerow(
                        [
                            "cell_id",
                            "aligner",
                            "read_layout",
                            "reused_alignment",
                            "group_id",
                            "reference",
                            "sam",
                        ]
                    )
                    w.writerows(values)

            write(rows)
            a = dict(
                exit_code=0,
                completion_bound=False,
                output_path=None,
                inputs_unchanged=True,
            )
            detector = dict(
                argv=["java", "-I", str(sam)], input_sha256={str(sam): digest(sam)}
            )
            args = dict(
                workflow_name="workflow_remaining3",
                cell_id="A102",
                group_id="XA102",
                reference="/reference",
                cohort_ids=["A101", "A102", "A103"],
            )
            self.assertEqual(
                alignment_output_binding(root, a, detector, **args)["output_sha256"],
                digest(sam),
            )
            write([rows[0], rows[0], rows[2]])
            with self.assertRaisesRegex(ValueError, "declared fresh alignment cohort"):
                alignment_output_binding(root, a, detector, **args)

    def test_registered_order_and_metadata(self):
        self.assertEqual([r[0] for r in ROWS], [101, 102, 103])
        self.assertEqual(compare_metadata([expected(r) for r in reversed(ROWS)]), [])
        self.assertTrue(
            expected(ROWS[0])["fastq_ftp"].endswith(
                "/017/SRR30918117/SRR30918117.fastq.gz"
            )
        )

    def test_layout_and_hash_changes_remain_conflicts(self):
        records = [expected(r) for r in ROWS]
        records[0].update(library_layout="PAIRED", fastq_md5="changed")
        self.assertEqual(
            {d["field"] for d in compare_metadata(records)},
            {"library_layout", "fastq_md5"},
        )

    def test_duplicate_missing_foreign_records_rejected(self):
        records = [expected(r) for r in ROWS]
        for changed in [
            records[:2],
            records + [records[0]],
            [dict(r, run_accession="other") for r in records],
        ]:
            with self.assertRaisesRegex(ValueError, "duplicate/missing/foreign"):
                compare_metadata(changed)

    def test_continuation_requires_explicit_scope(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            spec = root / "request.json"
            spec.write_text(
                json.dumps(dict(stage="tiny", category="continuation", argv=["fake"]))
            )
            with patch.object(
                run_stage, "initialize", return_value={"deadline_epoch": 9999999999}
            ):
                (root / "records").mkdir()
                with self.assertRaisesRegex(ValueError, "explicit scoped budget"):
                    run_stage.launch(root, spec)


if __name__ == "__main__":
    unittest.main()
