"""A05 original stored-signal summaries. No production imports or full densification."""
import argparse
import csv
import json
from pathlib import Path
import subprocess

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from inspect_sources import SOURCES, digest, plain

EXPECTED = {"Smart-seq3": (192, 2503, 2659, 5711, 2379), "IMR90": (23, 2443, 2914, 7422, 2429)}


def table(path, rows):
    assert rows, path
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def summarize(values):
    a = np.asarray(values)
    finite = a[np.isfinite(a)]
    q = np.quantile(finite, [0, .25, .5, .75, 1], method="linear") if len(finite) else [None] * 5
    return dict(n_total=len(a), n_finite=len(finite), n_missing=len(a) - len(finite),
                **dict(zip(["minimum", "q25", "median", "q75", "maximum"], q)),
                quantile_method="numpy linear; no filtering of zeros")


def read_sparse(g):
    x = ad.io.read_elem(g["X"])
    assert sparse.isspmatrix_csr(x), "Expected original CSR"
    assert x.has_canonical_format, "Duplicated coordinates must be reviewed"
    assert np.isfinite(x.data).all() and (x.data >= 0).all()
    assert (x.data == np.floor(x.data)).all()
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--original", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    repo = Path(__file__).resolve().parents[2]
    helper = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    # Exact helper bytes are additionally bound even when replaying from a packet.
    code_sha = digest(__file__)
    cells, candidates, memberships, thresholds, threshold_cells, summaries, checks, edges, axes = ([] for _ in range(9))
    for dataset, (member, sha) in SOURCES.items():
        path = args.original / member
        assert digest(path) == sha
        provenance = dict(dataset=dataset, source_object_sha256=sha, helper_commit=helper)
        with h5py.File(path, "r") as f:
            global_obs = ad.io.read_elem(f["obs"])
            co, cv = [ad.io.read_elem(f["mod/circ"][key]) for key in ["obs", "var"]]
            ro, rv = [ad.io.read_elem(f["mod/rna"][key]) for key in ["obs", "var"]]
            assert global_obs.index.is_unique and co.index.is_unique and ro.index.is_unique
            assert cv.index.is_unique and rv.index.is_unique
            assert set(global_obs.index) == set(co.index) | set(ro.index)
            # The original metadata has membership but no cell-bound completion records.
            status_fields = [c for frame in [global_obs, co, ro] for c in frame.columns if any(t in c.lower() for t in ["status", "completion", "completed"])]
            assert not status_fields, ("Review recorded processing fields", status_fields)
            x, rna = read_sparse(f["mod/circ"]), read_sparse(f["mod/rna"])
            assert x.shape == (len(co), len(cv)) and rna.shape == (len(ro), len(rv))
            pos = x > 0
            row_n = np.asarray(pos.sum(axis=1)).ravel().astype(np.int64)
            row_s = np.asarray(x.sum(axis=1, dtype=np.int64)).ravel()
            col_n = np.asarray(pos.sum(axis=0)).ravel().astype(np.int64)
            col_s = np.asarray(x.sum(axis=0, dtype=np.int64)).ravel()
            col_max = np.zeros(x.shape[1], dtype=np.int64)
            np.maximum.at(col_max, x.indices, x.data.astype(np.int64))
            rna_n = np.asarray((rna > 0).sum(axis=1)).ravel().astype(np.int64)
            rna_s = np.asarray(rna.sum(axis=1, dtype=np.int64)).ravel()
            n, p, nnz, total, anno = EXPECTED[dataset]
            actual = (x.shape[0], x.shape[1], int(np.count_nonzero(x.data)), int(x.data.sum(dtype=np.int64)), int((cv.host_gene.astype(str).str.strip() != "").sum()))
            checks.append(dict(**provenance, check="registration_totals", expected=json.dumps(EXPECTED[dataset]), actual=json.dumps(actual), passed=actual == EXPECTED[dataset]))
            assert actual == EXPECTED[dataset], (dataset, actual)
            ci, ri = {k: i for i, k in enumerate(co.index)}, {k: i for i, k in enumerate(ro.index)}
            recorded = co.membership.isin(["shared", "circ_only"]).to_numpy()
            n_members = int(recorded.sum())
            cnv_ids = set(ad.io.read_elem(f["mod/cnv/obs"]).index) if "cnv" in f["mod"] else set()
            for order, cid in enumerate(global_obs.index):
                i, j = ci.get(cid), ri.get(cid)
                membership = str(global_obs.loc[cid, "membership"])
                if i is not None:
                    assert str(co.loc[cid, "membership"]) == membership
                if j is not None:
                    assert str(ro.loc[cid, "membership"]) == membership
                rec = dict(**provenance, source_cell_order=order, cell_id=cid, recorded_membership=membership,
                    stored_circ_row_present=i is not None, stored_rna_row_present=j is not None,
                    recorded_circ_member=membership in ["shared", "circ_only"],
                    recorded_rna_member=membership in ["shared", "rna_only"], stored_cnv_row_present=cid in cnv_ids,
                    historical_circ_processing_status="unknown", completion_evidence="no cell-bound status recovered; aggregate workflow summaries do not bind completion",
                    matched_rna_circ=membership == "shared" and i is not None and j is not None,
                    circ_source_row=i, rna_source_row=j,
                    circ_positive_candidates=int(row_n[i]) if i is not None else None,
                    circ_support=int(row_s[i]) if i is not None else None,
                    rna_overlap_signal=int(rna_s[j]) if j is not None else None,
                    rna_positive_features=int(rna_n[j]) if j is not None else None,
                    circ_zero_support_row=bool(row_s[i] == 0) if i is not None else None,
                    completed_zero="unknown", rna_units="historical simple-overlap profile counts; not TPM or independently quantified transcriptomics")
                cells.append(rec)
                memberships.append({k: v for k, v in rec.items() if k not in ["circ_positive_candidates", "circ_support", "rna_overlap_signal", "rna_positive_features", "rna_units"]})
            # Retain ALL original host fields, including blank and multi-gene annotations.
            for j, candidate in enumerate(cv.index):
                rec = dict(**provenance, source_candidate_order=j, candidate_id=candidate,
                    **{k: plain(v) for k, v in cv.iloc[j].items()}, positive_stored_cells=int(col_n[j]),
                    pooled_support=int(col_s[j]), maximum_per_cell_support=int(col_max[j]),
                    stored_cell_denominator=n, stored_cell_fraction=float(col_n[j] / n),
                    recorded_circ_member_denominator=n_members,
                    positive_recorded_circ_members=int(pos[recorded, j].sum()),
                    recorded_circ_member_fraction=float(pos[recorded, j].sum() / n_members),
                    completed_cohort_denominator=None, completed_cohort_prevalence=None,
                    completion_note="undefined: historical assay completion unknown",
                    annotation_provenance="original /mod/circ/var; host_gene_source retained; original GTF digest unknown",
                    coordinate_convention="stored original coordinates; historical executing detector provenance incomplete")
                candidates.append(rec)
            for i, cid in enumerate(co.index):
                for k in range(x.indptr[i], x.indptr[i+1]):
                    if x.data[k] > 0:
                        edges.append(dict(**provenance, source_cell_order=i, cell_id=cid, source_candidate_order=int(x.indices[k]), candidate_id=cv.index[x.indices[k]], support=int(x.data[k])))
            for mode in ["individual_entry", "pooled_candidate"]:
                for threshold in [1, 2, 3, 5]:
                    # New independent mask against original X in EVERY iteration.
                    keep = x.data >= threshold if mode == "individual_entry" else col_s[x.indices] >= threshold
                    y = sparse.csr_matrix((np.where(keep, x.data, 0), x.indices.copy(), x.indptr.copy()), shape=x.shape)
                    yn = np.asarray((y > 0).sum(axis=1)).ravel()
                    ys = np.asarray(y.sum(axis=1, dtype=np.int64)).ravel()
                    ync = np.asarray((y > 0).sum(axis=0)).ravel()
                    rec = dict(**provenance, threshold_mode=mode, threshold=threshold,
                        stored_cell_denominator=n, recorded_circ_member_denominator=n_members,
                        original_candidate_denominator=p, retained_candidates=int(np.count_nonzero(ync)),
                        positive_cell_candidate_entries=int(np.count_nonzero(y.data)), retained_support=int(y.sum(dtype=np.int64)),
                        stored_cells_with_positive_support=int(np.count_nonzero(ys)), stored_zero_rows=int((ys == 0).sum()),
                        median_positive_candidates_per_stored_cell=float(np.median(yn)), median_support_per_stored_cell=float(np.median(ys)),
                        rule="original entry >= threshold" if mode == "individual_entry" else "original pooled candidate sum >= threshold; retain all original entries of eligible candidates")
                    thresholds.append(rec)
                    for i, cid in enumerate(co.index):
                        threshold_cells.append(dict(**provenance, threshold_mode=mode, threshold=threshold, source_cell_order=i, cell_id=cid, positive_candidates=int(yn[i]), support=int(ys[i]), recorded_circ_member=bool(recorded[i])))
            for subset, mask in [("all_stored_circ_rows", np.ones(n, dtype=bool)), ("recorded_circ_members", recorded)]:
                for metric, values in [("circ_positive_candidates", row_n), ("circ_support", row_s)]:
                    summaries.append(dict(**provenance, observation_unit="cell", subset=subset, metric=metric, **summarize(values[mask])))
            for metric, values in [("rna_overlap_signal", rna_s), ("rna_positive_features", rna_n)]:
                summaries.append(dict(**provenance, observation_unit="cell", subset="all_stored_rna_rows", metric=metric, **summarize(values)))
            for metric, values in [("positive_stored_cells", col_n), ("pooled_support", col_s), ("maximum_per_cell_support", col_max), ("stored_cell_fraction", col_n/n)]:
                summaries.append(dict(**provenance, observation_unit="candidate", subset="all_stored_candidates", metric=metric, **summarize(values)))
            for m, o, v, xx in [("circ", co, cv, x), ("rna", ro, rv, rna)]:
                axes.append(dict(**provenance, modality=m, rows=xx.shape[0], columns=xx.shape[1],
                    stored_entries=xx.nnz, numerical_nonzeros=int(np.count_nonzero(xx.data)), explicit_sparse_zeros=int((xx.data == 0).sum()),
                    total_signal=int(xx.sum(dtype=np.int64)), dtype=str(xx.dtype), missing_stored_values=int((~np.isfinite(xx.data)).sum()),
                    unique_cell_ids=o.index.is_unique, unique_feature_ids=v.index.is_unique, axis_orientation="cells by features"))
        assert digest(path) == sha
    for name, rows in [("descriptive_cell_metrics", cells), ("descriptive_candidate_metrics", candidates),
                       ("modality_membership_audit", memberships), ("threshold_sensitivity", thresholds),
                       ("threshold_cell_metrics", threshold_cells), ("distribution_summary", summaries),
                       ("registration_checks", checks), ("circ_nonzero_entries", edges), ("matrix_audit", axes)]:
        table(args.out / (name + ".tsv"), rows)
    (args.out / "analysis_provenance.json").write_text(json.dumps(dict(helper_commit=helper, analyze_sha256=code_sha,
        source_hashes={k: v[1] for k,v in SOURCES.items()}, sources_modified=False,
        correlations="not calculated; descriptive scatter only; no p-values or fitted association models",
        source_order="unchanged stored global cell order and circ var order", missing_encoding="blank numeric field; explicit unknown status"), indent=2) + "\n")
    print(json.dumps(dict(cells=len(cells), candidates=len(candidates), nonzero_entries=len(edges), threshold_rows=len(thresholds))))


if __name__ == "__main__":
    main()
