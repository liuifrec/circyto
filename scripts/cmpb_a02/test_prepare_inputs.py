"""Tiny synthetic structural tests, not biological acceptance."""

import gzip
from pathlib import Path
import tempfile
import unittest

from prepare_inputs import check_fastq, check_reference, hashes


class InputChecks(unittest.TestCase):
    def test_complete_single_end_and_gzip_crc(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(
                    "@SRR1.1 1 length=4\nACGT\n+\n!!!!\n@SRR1.2 2 length=4\nNNNN\n+\n####\n"
                )
            h = hashes(path)
            expected = dict(
                fastq_md5=h["md5"],
                fastq_bytes=h["bytes"],
                run_accession="SRR1",
                read_count=2,
                base_count=8,
            )
            # SRA spot ordinal 2 must not be mistaken for an explicit mate 2.
            self.assertEqual(check_fastq(path, expected)["reads"], 2)
            with gzip.open(path, "wt") as handle:
                handle.write("@SRR1.1\nACGT\n+\n!!!!\n@SRR1.2\nNNNN\n+\n####\n")
            h = hashes(path)
            expected.update(fastq_md5=h["md5"], fastq_bytes=h["bytes"])
            self.assertEqual(check_fastq(path, expected)["reads"], 2)
            expected["read_count"] = 3
            with self.assertRaisesRegex(ValueError, "totals"):
                check_fastq(path, expected)

    def test_explicit_mate_two_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("@SRR1.1/2\nACGT\n+\n!!!!\n")
            h = hashes(path)
            expected = dict(
                fastq_md5=h["md5"],
                fastq_bytes=h["bytes"],
                run_accession="SRR1",
                read_count=1,
                base_count=4,
            )
            with self.assertRaisesRegex(ValueError, "Mate-2"):
                check_fastq(path, expected)
            with gzip.open(path, "wt") as handle:
                handle.write("@SRR1.1 INSTRUMENT:1:FLOWCELL:1:1:2:3/2\nACGT\n+\n!!!!\n")
            h = hashes(path)
            expected.update(fastq_md5=h["md5"], fastq_bytes=h["bytes"])
            with self.assertRaisesRegex(ValueError, "Mate-2"):
                check_fastq(path, expected)
            with gzip.open(path, "wt") as handle:
                handle.write("@SRR1.1 INSTRUMENT:1:FLOWCELL:1:1:2:3/1\nACGT\n+\n!!!!\n")
            h = hashes(path)
            expected.update(fastq_md5=h["md5"], fastq_bytes=h["bytes"])
            self.assertEqual(
                check_fastq(path, expected)["explicit_mate_header_counts"],
                {"1": 1, "2": 0, "unspecified": 0},
            )

    def test_truncated_record_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("@SRR1.1\nACGT\n+\n!\n")
            h = hashes(path)
            expected = dict(
                fastq_md5=h["md5"],
                fastq_bytes=h["bytes"],
                run_accession="SRR1",
                read_count=1,
                base_count=4,
            )
            with self.assertRaisesRegex(ValueError, "structure"):
                check_fastq(path, expected)

    def test_reference_coordinate_compatibility(self):
        with tempfile.TemporaryDirectory() as directory:
            fa, gtf = Path(directory) / "ref.fa", Path(directory) / "ref.gtf"
            fa.write_text(">chr1\nACGT\n")
            gtf.write_text('chr1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "synthetic";\n')
            self.assertTrue(
                check_reference(fa, gtf)["all_gtf_contigs_and_coordinates_compatible"]
            )
            gtf.write_text('chr1\ttest\tgene\t1\t5\t.\t+\t.\tgene_id "synthetic";\n')
            with self.assertRaisesRegex(ValueError, "coordinate"):
                check_reference(fa, gtf)


if __name__ == "__main__":
    unittest.main()
