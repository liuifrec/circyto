from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from circyto.multimodal.sync import MUDATA_SYNC_POLICY, mudata_from_modalities, read_h5mu, write_h5mu


ad = pytest.importorskip("anndata")
mu = pytest.importorskip("mudata")
sp = pytest.importorskip("scipy.sparse")


def _partial_modalities() -> dict[str, Any]:
    rna_obs = pd.DataFrame(
        {
            "shared": pd.Series(["B", "A"], index=["cell_b", "cell_a"], dtype=object),
            "rna_only": pd.Series([10, pd.NA], index=["cell_b", "cell_a"], dtype="Int64"),
        },
        index=["cell_b", "cell_a"],
    )
    circ_obs = pd.DataFrame(
        {
            "shared": pd.Series(["A", "C"], index=["cell_a", "cell_c"], dtype=object),
            "circ_only": [1.5, np.nan],
        },
        index=["cell_a", "cell_c"],
    )
    return {
        "rna": ad.AnnData(
            X=sp.csr_matrix([[1, 0], [0, 2]], dtype=np.int32),
            obs=rna_obs,
            var=pd.DataFrame(
                {"kind": ["gene", "gene"], "gene_symbol": ["GENE2", "GENE1"]},
                index=["gene_2", "gene_1"],
            ),
        ),
        "circ": ad.AnnData(
            X=sp.csr_matrix([[3], [4]], dtype=np.int32),
            obs=circ_obs,
            var=pd.DataFrame({"kind": ["circ"], "host_gene": ["GENE1"]}, index=["circ_1"]),
        ),
    }


def _aligned_modalities(extra: str | None, *, empty_circ: bool = False) -> dict[str, Any]:
    obs = pd.DataFrame(
        {
            "batch": pd.Series(["b2", "b1"], index=["cell_b", "cell_a"], dtype=object),
            "quality": [2.5, np.nan],
        },
        index=["cell_b", "cell_a"],
    )
    modalities: dict[str, Any] = {
        "rna": ad.AnnData(
            X=sp.csr_matrix([[1, 0], [0, 2]], dtype=np.int32),
            obs=obs.assign(rna_metric=[7, 8]),
            var=pd.DataFrame(
                {"gene_id": ["gene_2", "gene_1"], "gene_name": ["GENE2", "GENE1"]},
                index=["gene_2", "gene_1"],
            ),
        ),
        "circ": ad.AnnData(
            X=sp.csr_matrix((2, 0), dtype=np.int32)
            if empty_circ
            else sp.csr_matrix([[3, 0], [0, 4]], dtype=np.int32),
            obs=obs.assign(circ_metric=[1.0, np.nan]),
            var=pd.DataFrame(index=pd.Index([], dtype=str))
            if empty_circ
            else pd.DataFrame(
                {"circ_id": ["circ_2", "circ_1"], "host_gene": ["GENE2", "GENE1"]},
                index=["circ_2", "circ_1"],
            ),
        ),
    }
    if extra == "candidate_snv":
        modalities[extra] = ad.AnnData(
            X=np.array([[5], [0]], dtype=np.float32),
            obs=obs.assign(snv_metric=[0.2, np.nan]),
            var=pd.DataFrame({"variant_id": ["chr1:10:A>G"]}, index=["chr1:10:A>G"]),
        )
    elif extra == "rt":
        modalities[extra] = ad.AnnData(
            X=np.array([[0, 1], [1, 1]], dtype=np.int8),
            obs=obs.assign(rt_metric=["early", "late"]),
            var=pd.DataFrame(
                {
                    "seqname": pd.Categorical(["chr1", "chr1"]),
                    "start": [0, 50_000],
                    "end": [50_000, 100_000],
                },
                index=["chr1:0-50000", "chr1:50000-100000"],
            ),
        )
    elif extra == "cnv":
        modalities[extra] = ad.AnnData(
            X=np.array([[2, 3], [1, 2]], dtype=np.int16),
            obs=obs.assign(cnv_metric=[2.5, np.nan]),
            var=pd.DataFrame(
                {
                    "seqname": pd.Categorical(["chr2", "chr2"]),
                    "start": [0, 50_000],
                    "end": [50_000, 100_000],
                },
                index=["chr2:0-50000", "chr2:50000-100000"],
            ),
        )
    return modalities


def _dense(matrix: Any) -> np.ndarray:
    return matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)


def test_historical_partial_overlap_construction_semantics() -> None:
    modalities = _partial_modalities()
    mdata = mudata_from_modalities(modalities)

    assert list(mdata.obs_names) == ["cell_b", "cell_a", "cell_c"]
    assert list(mdata.mod["rna"].obs_names) == ["cell_b", "cell_a"]
    assert list(mdata.mod["circ"].obs_names) == ["cell_a", "cell_c"]
    assert set(mdata.obs.columns) == {"rna:shared", "rna:rna_only", "circ:shared", "circ:circ_only"}
    assert mdata.obs.loc["cell_a", "rna:shared"] == "A"
    assert mdata.obs.loc["cell_a", "circ:shared"] == "A"
    assert pd.isna(mdata.obs.loc["cell_c", "rna:rna_only"])
    assert str(mdata.obs["rna:rna_only"].dtype) == "Int64"

    assert list(mdata.mod["rna"].var_names) == ["gene_2", "gene_1"]
    assert list(mdata.mod["circ"].var_names) == ["circ_1"]
    assert list(mdata.var_names) == ["gene_2", "gene_1", "circ_1"]
    assert set(mdata.var.columns) == {"rna:gene_symbol", "circ:host_gene", "kind"}
    assert mdata.var.loc["gene_1", "rna:gene_symbol"] == "GENE1"
    assert mdata.var.loc["circ_1", "circ:host_gene"] == "GENE1"

    assert mdata.obsmap["rna"].tolist() == [1, 2, 0]
    assert mdata.obsmap["circ"].tolist() == [0, 1, 2]
    assert mdata.varmap["rna"].tolist() == [1, 2, 0]
    assert mdata.varmap["circ"].tolist() == [0, 0, 1]
    np.testing.assert_array_equal(_dense(mdata.mod["rna"].X), [[1, 0], [0, 2]])
    np.testing.assert_array_equal(_dense(mdata.mod["circ"].X), [[3], [4]])


@pytest.mark.parametrize(
    ("extra", "empty_circ"),
    [
        (None, False),
        ("candidate_snv", False),
        ("rt", False),
        ("cnv", False),
        (None, True),
    ],
)
def test_representative_mudata_round_trip_semantics(
    tmp_path: Path,
    extra: str | None,
    empty_circ: bool,
) -> None:
    modalities = _aligned_modalities(extra, empty_circ=empty_circ)
    expected_modalities = list(modalities)
    expected_obs = {name: adata.obs.copy() for name, adata in modalities.items()}
    expected_var = {name: adata.var.copy() for name, adata in modalities.items()}
    expected_x = {name: _dense(adata.X).copy() for name, adata in modalities.items()}

    mdata = mudata_from_modalities(modalities)
    mdata.obs = modalities["rna"].obs[["batch", "quality"]].copy()
    path = tmp_path / f"{extra or 'bimodal'}{'_empty' if empty_circ else ''}.h5mu"
    write_h5mu(mdata, path)
    restored = read_h5mu(path)

    assert list(restored.mod) == expected_modalities
    assert list(restored.obs_names) == ["cell_b", "cell_a"]
    assert restored.obs["batch"].tolist() == ["b2", "b1"]
    assert str(restored.obs["quality"].dtype) == "float64"
    assert pd.isna(restored.obs.loc["cell_a", "quality"])
    assert list(restored.var_names) == [
        feature_id
        for name in expected_modalities
        for feature_id in expected_var[name].index
    ]
    for name in expected_modalities:
        assert list(restored.mod[name].obs_names) == ["cell_b", "cell_a"]
        assert list(restored.mod[name].var_names) == list(expected_var[name].index)
        assert_frame_equal(restored.mod[name].obs, expected_obs[name], check_dtype=True)
        assert_frame_equal(restored.mod[name].var, expected_var[name], check_dtype=True)
        np.testing.assert_array_equal(_dense(restored.mod[name].X), expected_x[name])
        for column in expected_obs[name].columns:
            assert f"{name}:{column}" in restored.obs.columns
    assert list(restored.mod["rna"].var_names) == ["gene_2", "gene_1"]
    assert list(restored.mod["circ"].var_names) == ([] if empty_circ else ["circ_2", "circ_1"])


def test_partial_overlap_round_trip_keeps_union_order_and_missing_values(tmp_path: Path) -> None:
    mdata = mudata_from_modalities(_partial_modalities())
    path = tmp_path / "partial.h5mu"
    write_h5mu(mdata, path)
    restored = read_h5mu(path)

    assert list(restored.obs_names) == ["cell_b", "cell_a", "cell_c"]
    assert list(restored.mod["rna"].obs_names) == ["cell_b", "cell_a"]
    assert list(restored.mod["circ"].obs_names) == ["cell_a", "cell_c"]
    assert pd.isna(restored.obs.loc["cell_b", "circ:circ_only"])
    assert pd.isna(restored.obs.loc["cell_c", "rna:rna_only"])
    assert restored.obsmap["rna"].tolist() == [1, 2, 0]
    assert restored.obsmap["circ"].tolist() == [0, 1, 2]


def test_future_default_differs_but_explicit_policy_preserves_global_metadata() -> None:
    with mu.set_options(pull_on_update=False):
        future_default = mu.MuData(_partial_modalities())
        explicit = mudata_from_modalities(_partial_modalities())
        restored_future_default = mu.MuData(_partial_modalities())

    assert list(future_default.obs_names) == list(explicit.obs_names)
    assert list(future_default.var_names) == list(explicit.var_names)
    assert future_default.obs.empty
    assert future_default.var.empty
    assert restored_future_default.obs.empty
    assert restored_future_default.var.empty
    assert set(explicit.obs.columns) == {"rna:shared", "rna:rna_only", "circ:shared", "circ:circ_only"}
    assert set(explicit.var.columns) == {"rna:gene_symbol", "circ:host_gene", "kind"}
    assert MUDATA_SYNC_POLICY == "pull_on_update=True"
