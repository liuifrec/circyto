import json
from pathlib import Path
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

import reference
from oracle import verify


def fixture(p, large=False):
    (p / "cells.tsv").write_text(
        "cell_id\tdetector_status\n001\tcompleted\nNA\tcompleted\n"
    )
    header = "\t".join(reference.NATIVE) + "\n"
    count = str(2**53 + 1 if large else 7)
    (p / "a.tsv").write_text(
        header
        + "\t".join(
            [
                "chr1:2|8",
                "chr1",
                "2",
                "8",
                count,
                "1",
                "0",
                "1",
                "exon",
                "G",
                "+",
                "synthetic-read",
                "1",
            ]
        )
        + "\n"
    )
    (p / "b.tsv").write_text(header)
    (p / "native.tsv").write_text(
        "cell_id\tsample_id\tpath\tdetector_version\tresult_scope\n001\tS1\ta.tsv\t3.0.1\tsingle-sample\nNA\tS2\tb.tsv\t3.0.1\tsingle-sample\n"
    )
    r = ad.AnnData(
        sparse.csr_matrix(np.array([[3, 20], [10, 2]], dtype=np.int64)),
        obs=pd.DataFrame({"source": ["rna-b", "rna-a"]}, index=["NA", "001"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    r.layers["original"] = r.X.copy()
    r.obsm["embedding"] = np.array([[2.0, 3.0], [5.0, 7.0]])
    r.uns["evidence"] = "synthetic fixture"
    r.write_h5ad(p / "rna.h5ad")
    return dict(
        mode="native",
        input=str(p / "native.tsv"),
        cells=str(p / "cells.tsv"),
        rna=str(p / "rna.h5ad"),
        genome_build="synthetic",
        coordinate_convention="1-based-closed",
        support_definition="CIRI3_junction_reads",
        detector="CIRI3",
    )


class ReferenceContract(unittest.TestCase):
    def test_exact_large_integer_literal_ID_and_all_RNA_fields(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            c = fixture(p, True)
            reference.convert(c, p / "out")
            result, prevalence = verify(c, p / "out", "reference")
            self.assertTrue(result["passed"], result)
            self.assertEqual(prevalence[0]["oracle_denominator"], 2)
            self.assertEqual(
                int(ad.read_h5ad(p / "out/anndata/circ_counts.h5ad").X[0, 0]), 2**53 + 1
            )

    def test_no_completed_denominator_is_null(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            c = fixture(p)
            (p / "cells.tsv").write_text("cell_id\n001\nNA\n")
            reference.convert(c, p / "out")
            q = json.loads((p / "out/qc/metrics.json").read_text())
            self.assertIsNone(q["detection_fraction_completed"]["value"])
            self.assertEqual(q["detection_fraction_completed"]["denominator"], 0)

    def test_failed_positive_support_is_retained(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            c = fixture(p)
            (p / "cells.tsv").write_text(
                "cell_id\tdetector_status\n001\tfailed\nNA\tunknown\n"
            )
            reference.convert(c, p / "out")
            self.assertTrue(verify(c, p / "out", "reference")[0]["passed"])

    def test_reject_duplicate_metadata_and_duplicate_RNA(self):
        for which in ["metadata", "rna"]:
            with self.subTest(which=which), tempfile.TemporaryDirectory() as d:
                p = Path(d)
                c = fixture(p)
                if which == "metadata":
                    (p / "cells.tsv").write_text("cell_id\n001\n001\n")
                else:
                    a = ad.read_h5ad(p / "rna.h5ad")
                    a.obs_names = ["001", "001"]
                    a.write_h5ad(p / "rna.h5ad")
                with self.assertRaisesRegex(ValueError, "Duplicate"):
                    reference.convert(c, p / "out")

    def test_native_missing_result_is_not_completed_zero(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            c = fixture(p)
            (p / "native.tsv").write_text(
                "\n".join((p / "native.tsv").read_text().splitlines()[:2]) + "\n"
            )
            with self.assertRaisesRegex(ValueError, "lacks explicit source"):
                reference.convert(c, p / "out")

    def test_reject_duplicated_native_candidate(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            c = fixture(p)
            t = (p / "a.tsv").read_text()
            (p / "a.tsv").write_text(t + t.splitlines()[1] + "\n")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                reference.convert(c, p / "out")

    def test_reject_conflicting_native_annotation_and_mixed_declaration(self):
        for which in ["annotation", "declaration"]:
            with self.subTest(which=which), tempfile.TemporaryDirectory() as d:
                p = Path(d)
                c = fixture(p)
                if which == "annotation":
                    (p / "b.tsv").write_text(
                        (p / "a.tsv").read_text().replace("\tG\t", "\tH\t")
                    )
                else:
                    (p / "cells.tsv").write_text(
                        "cell_id\tgenome_build\n001\twrong\nNA\tsynthetic\n"
                    )
                with self.assertRaisesRegex(ValueError, "Conflicting"):
                    reference.convert(c, p / "out")

    def test_total_overflow_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            c = fixture(p)
            (p / "a.tsv").write_text(
                (p / "a.tsv").read_text().replace("\t7\t", "\t9223372036854775807\t")
            )
            (p / "b.tsv").write_text((p / "a.tsv").read_text())
            with self.assertRaisesRegex(ValueError, "Aggregate count overflow"):
                reference.convert(c, p / "out")

    def test_explicit_sparse_zero_and_absent_column(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            c = fixture(p)
            folder = p / "matrix"
            folder.mkdir()
            c.update(mode="matrix", input=str(folder))
            (folder / "cell_index.txt").write_text("NA\n001\n")
            (folder / "circ_index.txt").write_text("opaque\n")
            (folder / "circ_feature_table.tsv").write_text(
                "circ_id\tchrom\tstart\tend\tstrand\thost_gene\nopaque\tchr1\t2\t8\t+\tGene\n"
            )
            mmwrite(
                folder / "circ_counts.mtx",
                sparse.coo_matrix(
                    (np.array([0, 7], dtype=np.int64), ([0, 0], [0, 1])), shape=(1, 2)
                ),
            )
            reference.convert(c, p / "out")
            self.assertTrue(verify(c, p / "out", "reference")[0]["passed"])
            (folder / "cell_index.txt").write_text("001\n")
            mmwrite(
                folder / "circ_counts.mtx",
                sparse.coo_matrix(np.array([[7]], dtype=np.int64)),
            )
            with self.assertRaisesRegex(ValueError, "lacks explicit source"):
                reference.convert(c, p / "bad")

    def test_oracle_detects_changed_counts(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            c = fixture(p)
            reference.convert(c, p / "out")
            a = ad.read_h5ad(p / "out/anndata/circ_counts.h5ad")
            a.X.data[0] += 1
            a.write_h5ad(p / "out/anndata/circ_counts.h5ad")
            self.assertEqual(verify(c, p / "out", "reference")[0]["count_errors"], 1)


if __name__ == "__main__":
    unittest.main()
