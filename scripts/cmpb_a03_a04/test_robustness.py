import json
from pathlib import Path
import tempfile
import unittest
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from perturb import plan, choose, generate
from reference import convert
from robustness_oracle import verify
from robustness import archive
from unsafe import run as unsafe
from run_stage import digest


def base(p):
    p.mkdir()
    ids = [f"cell{i}" for i in range(23)]
    x = np.array([[i + 1, 0, 3, 0] for i in range(23)], dtype=np.int64)
    var = pd.DataFrame(
        dict(
            circ_id=["a", "b", "c", "d"],
            chrom=["chr1"] * 4,
            start=[1, 3, 5, 7],
            end=[2, 4, 6, 8],
            strand=["+"] * 4,
            host_gene=["A", "B", "C", "D"],
        ),
        index=["a", "b", "c", "d"],
    )
    a = ad.AnnData(sparse.csr_matrix(x), obs=pd.DataFrame(index=ids), var=var)
    a.write_h5ad(p / "circ_source_copy.h5ad")
    r = ad.AnnData(
        sparse.csr_matrix(x[::-1]),
        obs=pd.DataFrame(index=ids[::-1]),
        var=pd.DataFrame(index=["g1", "g2", "g3", "g4"]),
    )
    r.layers["source"] = r.X.copy()
    r.uns["label"] = "synthetic tiny 23-cell fixture"
    r.write_h5ad(p / "rna_source_copy.h5ad")
    (p / "base_record.json").write_text(
        json.dumps(dict(sampling_weights={c: 1 / (i + 1) for i, c in enumerate(ids)}))
    )


class PerturbationChecks(unittest.TestCase):
    def test_plan_has_exact_seeds_counts_and_deterministic_controls(self):
        cases = plan()
        self.assertEqual(len(cases), 715)
        self.assertEqual(len({c["case"] for c in cases}), 715)
        for s, n in [("R01", 100), ("R02", 200), ("R05", 200), ("R06", 200)]:
            selected = [c for c in cases if c["scenario"] == s]
            self.assertEqual(len(selected), n)
            self.assertEqual(
                {c["seed"] for c in selected}, set(range(2026101500, 2026101520))
            )

    def test_registered_rounding_and_weight_validation(self):
        ids = [str(i) for i in range(23)]
        weights = {c: 1 for c in ids}
        self.assertEqual(
            [
                len(choose(ids, weights, f, 2026101500, "uniform")[1])
                for f in [0, 0.05, 0.1, 0.25, 0.5]
            ],
            [0, 1, 2, 6, 12],
        )
        weights["0"] = float("nan")
        with self.assertRaises(ValueError):
            choose(ids, weights, 0.25, 2026101500, "profile-associated")

    def test_permutation_k1_is_identity_and_source_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            base(p / "base")
            before = digest(p / "base/rna_source_copy.h5ad")
            c = next(
                c
                for c in plan()
                if c["scenario"] == "R01" and c["nominal_fraction"] == 0.05
            )
            config, detail = generate(p / "base", c, p / "case")
            self.assertEqual(detail["selected_count"], 1)
            self.assertEqual(detail["moved_identities"], 0)
            self.assertEqual(before, digest(p / "base/rna_source_copy.h5ad"))
            convert(config, p / "out")
            self.assertTrue(verify(config, p / "out", "reference")[0]["passed"])

    def test_partial_overlap_changes_membership_not_status_or_values(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            base(p / "base")
            c = next(
                c
                for c in plan()
                if c["scenario"] == "R02" and c["nominal_fraction"] == 0.25
            )
            config, details = generate(p / "base", c, p / "case")
            self.assertEqual(
                (details["n_circ_cells"], details["n_rna_cells"]), (20, 20)
            )
            convert(config, p / "out")
            self.assertTrue(verify(config, p / "out", "reference")[0]["passed"])
            unsafe(config, p / "unsafe")
            r, _ = verify(config, p / "unsafe", "unsafe")
            self.assertFalse(r["passed"])
            self.assertGreater(r["cell_assignment_errors"] + r["status_errors"], 0)

    def test_failed_state_retains_support_and_partial_case_is_explicit(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            base(p / "base")
            for scenario in ["R05", "R06"]:
                c = next(
                    c
                    for c in plan()
                    if c["scenario"] == scenario and c["nominal_fraction"] == 0.25
                )
                config, details = generate(p / "base", c, p / scenario)
                convert(config, p / (scenario + "-out"))
                r, prev = verify(config, p / (scenario + "-out"), "reference")
                self.assertTrue(r["passed"])
                self.assertEqual(prev[0]["oracle_denominator"], 17)
                x = mmread(Path(config["input"]) / "circ_counts.mtx").T.tocsr()
                if scenario == "R05":
                    self.assertEqual(
                        int(x.sum()),
                        int(ad.read_h5ad(p / "base/circ_source_copy.h5ad").X.sum()),
                    )
                else:
                    self.assertGreater(int(x.sum()), 0)

    def test_archive_keeps_exact_bytes_and_removes_only_generated_copy(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            case = p / "case"
            case.mkdir()
            (case / "fixture.txt").write_text("generated fixture")
            outside = p / "original.txt"
            outside.write_text("protected")
            r = archive(case)
            self.assertTrue(r["lossless_member_hashes_verified"])
            self.assertFalse(case.exists())
            self.assertEqual(outside.read_text(), "protected")
            self.assertTrue(Path(r["archive"]).exists())


if __name__ == "__main__":
    unittest.main()
