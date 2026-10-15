"""Small synthetic oracle tests; not an A04 robustness experiment."""

import csv
from pathlib import Path
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from independent_checks import (
    CELL,
    HEADER,
    Tables,
    equal_adata,
    exact,
    native,
    verify_native_object,
)


class IndependentChecks(unittest.TestCase):
    def test_status_semantics_cannot_pass_from_rendering_agreement(self):
        from independent_checks import status_comparisons

        payload = {
            "per_cell": [{"cell_id": CELL, "alignment_status": "aligned"}],
            "workflow": {"status": "completed"},
            "candidate_prevalence_denominator": 1,
        }
        self.assertTrue(
            all(a == b for _, a, b in status_comparisons(payload, "aligned", "aligned"))
        )
        payload["per_cell"][0]["alignment_status"] = "missing"
        payload["workflow"]["status"] = "incomplete"
        payload["candidate_prevalence_denominator"] = None
        self.assertTrue(
            all(a != b for _, a, b in status_comparisons(payload, "aligned", "missing"))
        )

    def test_integer_storage_width_is_lossless_but_changes_and_floats_fail(self):
        from independent_checks import exact_integer_csr_values

        expected = sparse.csr_matrix([[1, 2147483647]], dtype=np.int64)
        actual = expected.astype(np.int32)
        exact_integer_csr_values(expected, actual, "workflow int32")
        with self.assertRaises(AssertionError):
            exact_integer_csr_values(expected, actual.astype(float), "float")
        actual.data[0] += 1
        with self.assertRaises(AssertionError):
            exact_integer_csr_values(expected, actual, "changed integer")

    def test_displayed_tables_cards_and_svg_values(self):
        parsed = Tables()
        parsed.feed(
            '<article class="metric"><div class="label">Support</div><div class="value">7</div></article>'
            "<table><tr><th>ID</th><th>value</th></tr><tr><td>001 &amp; NA</td><td>7</td></tr></table>"
            '<svg aria-label="counts"><title>counts</title><rect><title>7: 1</title></rect><text>1</text></svg>'
        )
        self.assertEqual(parsed.cards, [{"label": "Support", "value": "7"}])
        self.assertEqual(parsed.tables, [[["ID", "value"], ["001 & NA", "7"]]])
        self.assertEqual(
            parsed.svgs,
            [{"label": "counts", "titles": ["counts", "7: 1"], "texts": ["1"]}],
        )

    def test_exact_integer_and_annotation_error_detection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.tsv"
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(HEADER)
                writer.writerow(
                    [
                        "chr1:1|4",
                        "chr1",
                        "1",
                        "4",
                        str(2**53 + 1),
                        "1",
                        "0",
                        "1",
                        "exon",
                        "TEST",
                        "+",
                        "synthetic:read",
                        "1",
                    ]
                )
            records, X = native(path)
            self.assertEqual(int(X[0, 0]), 2**53 + 1)
            var = pd.DataFrame(
                dict(
                    circ_id=["chr1:1|4"],
                    chrom=["chr1"],
                    start=[1],
                    end=[4],
                    strand=["+"],
                    ciri3_gene_id=["TEST"],
                    ciri3_circRNA_type=["exon"],
                    host_gene_id=["TEST"],
                    host_gene_id_source=["CIRI3_gene_id"],
                ),
                index=["chr1:1|4"],
            )
            obj = ad.AnnData(X.copy(), obs=pd.DataFrame(index=[CELL]), var=var)
            verify_native_object(records, X, obj)
            obj.X.data[0] -= 1
            with self.assertRaises(AssertionError):
                verify_native_object(records, X, obj)
            obj.X = X.copy()
            obj.var.loc["chr1:1|4", "strand"] = "-"
            with self.assertRaises(AssertionError):
                verify_native_object(records, X, obj)

    def test_genuine_header_only_contract_without_invented_status(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic-empty.tsv"
            path.write_text("\t".join(HEADER) + "\n")
            records, X = native(path)
            self.assertEqual(records, [])
            self.assertEqual(X.shape, (1, 0))
            # The reader provides values only, never a processing status.
            self.assertEqual(X.nnz, 0)

    def test_rna_layer_and_embedding_differences_are_detected(self):
        a = ad.AnnData(
            sparse.csr_matrix(np.array([[1, 2]], dtype=np.int64)),
            obs=pd.DataFrame(index=[CELL]),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        a.layers["counts"] = a.X.copy()
        a.obsm["position"] = np.array([[1.0, 2.0]])
        a.uns["declaration"] = {"source": "synthetic fixture", "unknown": np.nan}
        b = a.copy()
        equal_adata(a, b)
        b.layers["counts"].data[0] = 2
        with self.assertRaises(AssertionError):
            equal_adata(a, b)
        b = a.copy()
        b.obsm["position"][0, 0] = 3
        with self.assertRaises(AssertionError):
            equal_adata(a, b)

    def test_sparse_explicit_zero_does_not_require_dense_allocation(self):
        a = sparse.csr_matrix(
            ([1, 0], ([0, 0], [0, 1])), shape=(1, 1000000), dtype=np.int64
        )
        b = a.copy()
        b.eliminate_zeros()
        exact(a, b)


if __name__ == "__main__":
    unittest.main()
