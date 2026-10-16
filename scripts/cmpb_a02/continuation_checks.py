#!/usr/bin/env python3
"""A02 cohort extension of the reviewed independent CSV/integer/CSR oracle.

No circyto parser, normalizer, reporting function, or production API is imported.
Only the registered four biological identities and two distinct report cohorts
are checked. This is not an A03 conversion benchmark or A04 experiment.
"""

import argparse
from collections import Counter
import csv
from fractions import Fraction
import json
from pathlib import Path
from statistics import median

import anndata as ad
import mudata as md
import numpy as np
from scipy import sparse
from scipy.io import mmread

from continuation_records import CELLS, four_rna, native_union
from independent_checks import (
    Tables,
    equal_adata,
    exact,
    exact_integer_csr_values,
    native,
    verify_native_object,
    write_checks,
)
from run_stage import digest, dump, stamp


def oracle(per_cell, ids, features):
    lookup = {c: i for i, c in enumerate(features)}
    rr, cc, vv = [], [], []
    for i, cell in enumerate(ids):
        for row in per_cell[cell]:
            rr.append(i)
            cc.append(lookup[row["circRNA_ID"]])
            vv.append(int(row["#junction_reads"]))
    matrix = sparse.csr_matrix(
        (np.array(vv, dtype=np.int64), (rr, cc)),
        shape=(len(ids), len(features)),
        dtype=np.int64,
    )
    matrix.eliminate_zeros()
    return matrix


def table(path):
    with path.open(newline="") as h:
        reader = csv.DictReader(h, delimiter="\t")
        return reader.fieldnames, list(reader)


class Checks:
    def __init__(self, path):
        self.path, self.rows = path, []

    def check(self, name, expected, observed, scope):
        ok = expected == observed
        self.rows.append(
            dict(
                check=name,
                scope=scope,
                expected=str(expected),
                observed=str(observed),
                status="PASS" if ok else "FAIL",
            )
        )
        write_checks(self.path, self.rows)
        if not ok:
            raise AssertionError(f"{scope}: {name}: {expected!r} != {observed!r}")

    def exact(self, name, expected, observed, scope, function=exact):
        try:
            function(expected, observed, name)
        except Exception as exc:
            self.check(name, "exact equality", repr(exc), scope)
        self.check(name, "exact equality", "exact equality", scope)


def check_mudata(a, b):
    if set(a.mod) != set(b.mod):
        raise AssertionError("MuData modalities changed")
    for key in a.mod:
        equal_adata(a.mod[key], b.mod[key], "MuData " + key)
    for key in ["obs", "var", "uns"]:
        exact(getattr(a, key), getattr(b, key), "MuData global " + key)
    for key in ["obsm", "varm", "obsp", "varp", "obsmap", "varmap"]:
        exact(dict(getattr(a, key)), dict(getattr(b, key)), "MuData global " + key)


def import_checks(root, original, per_cell, union, checks, audits):
    out = root / "import_four"
    circ = ad.read_h5ad(out / "anndata/circ_counts.h5ad")
    features = list(union)
    want = oracle(per_cell, CELLS, features)
    checks.check("import_cell_order", CELLS, list(circ.obs_names), "joint import")
    checks.check(
        "import_first_seen_feature_order",
        features,
        list(circ.var_names),
        "joint import",
    )
    checks.exact("every_integer_CSR_entry", want, circ.X, "joint import")
    checks.check("native_import_dtype", "int64", str(circ.X.dtype), "joint import")
    sources = json.loads((root / "records/joint_native_sources.json").read_text())[
        "cells"
    ]
    for i, cell in enumerate(CELLS):
        records = per_cell[cell]
        _, expected = native(sources[i]["native_path"])
        selected = circ[[cell], [r["circRNA_ID"] for r in records]].copy()
        verify_native_object(records, expected, selected, cell_id=cell)
        checks.check("all_native_coordinates_strands_annotations", True, True, cell)
        for row in records:
            checks.check(
                "candidate_support:" + row["circRNA_ID"],
                int(row["#junction_reads"]),
                int(circ.X[i, features.index(row["circRNA_ID"])]),
                cell,
            )
        for path in [
            Path(sources[i]["native_path"]),
            out / f"evidence/ciri3_native_{i + 1:06d}.tsv",
        ]:
            checks.check(
                "native_bytes_including_ancillary_read_IDs",
                sources[i]["native_sha256"],
                digest(path),
                cell,
            )
    checks.exact(
        "MatrixMarket_all_values",
        want,
        mmread(out / "matrix/circ_counts.mtx").T.tocsr(),
        "joint import",
    )
    checks.check(
        "matrix_text_cell_axis",
        CELLS,
        (out / "matrix/cell_index.txt").read_text().splitlines(),
        "joint import",
    )
    checks.check(
        "matrix_text_feature_axis",
        features,
        (out / "matrix/circ_index.txt").read_text().splitlines(),
        "joint import",
    )
    rna = ad.read_h5ad(root / "inputs/rna/imr90_original_overlap_profile.four.h5ad")
    checks.exact(
        "historical_RNA_all_fields", four_rna(original), rna, "RNA slice", equal_adata
    )
    checks.check(
        "historical_RNA_original_order", CELLS[::-1], list(rna.obs_names), "RNA slice"
    )
    checks.check("historical_RNA_feature_count", 63187, rna.n_vars, "RNA slice")
    multi = md.read_h5mu(out / "mudata/circyto_multimodal.h5mu")
    checks.check("global_cell_order", CELLS, list(multi.obs_names), "joint MuData")
    checks.exact(
        "circ_modality_all_fields", circ, multi.mod["circ"], "joint MuData", equal_adata
    )
    checks.exact(
        "RNA_modality_all_fields", rna, multi.mod["rna"], "joint MuData", equal_adata
    )
    for modality in ["circ", "rna"]:
        for global_i, cell in enumerate(CELLS):
            mapped = int(np.asarray(multi.obsmap[modality]).ravel()[global_i])
            checks.check(
                "exact_ID_map:" + modality,
                cell,
                str(multi.mod[modality].obs_names[mapped - 1]) if mapped else "missing",
                cell,
            )
    checks.check(
        "shared_membership",
        ["shared"] * 4,
        multi.obs.membership.tolist(),
        "joint MuData",
    )
    for key in ["circ_present", "rna_present"]:
        checks.check(key, [True] * 4, multi.obs[key].tolist(), "joint MuData")
    rt = audits
    circ.write_h5ad(rt / "circ_roundtrip.h5ad")
    checks.exact(
        "AnnData_read_write_read",
        circ,
        ad.read_h5ad(rt / "circ_roundtrip.h5ad"),
        "joint import",
        equal_adata,
    )
    multi.write_h5mu(rt / "multimodal_roundtrip.h5mu")
    check_mudata(multi, md.read_h5mu(rt / "multimodal_roundtrip.h5mu"))
    checks.check(
        "MuData_read_write_read_all_fields_maps_missingness", True, True, "joint import"
    )
    repeat = root / "import_four_repeat"
    if repeat.exists():
        checks.exact(
            "repeat_import_AnnData_contents",
            circ,
            ad.read_h5ad(repeat / "anndata/circ_counts.h5ad"),
            "repeat postprocessing",
            equal_adata,
        )
        check_mudata(multi, md.read_h5mu(repeat / "mudata/circyto_multimodal.h5mu"))
        checks.check(
            "repeat_import_MuData_contents", True, True, "repeat postprocessing"
        )
        for subdir in ["matrix", "qc", "evidence"]:
            for p in (out / subdir).iterdir():
                if p.is_file() and p.suffix != ".h5ad":
                    checks.check(
                        "repeat_text_bytes:" + str(p.relative_to(out)),
                        digest(p),
                        digest(repeat / p.relative_to(out)),
                        "repeat postprocessing",
                    )
    return circ, want


def display_fraction(value):
    """Exact displayed numeric value; decimal zero padding is not rounding."""
    return Fraction(value)


def import_source_path(root, name):
    """RNA is hash-referenced; other supplied evidence is copied by the importer."""
    if name != "rna.h5ad":
        return root / "import_four/evidence" / name
    extraction = json.loads((root / "records/rna_extraction.json").read_text())
    command = json.loads((root / "records/continuation-import-four.json").read_text())
    argv = command["argv"]
    if argv.count("--rna") != 1:
        raise ValueError("Import command has no unambiguous RNA source")
    source = Path(extraction["output"])
    if source != Path(argv[argv.index("--rna") + 1]):
        raise ValueError("RNA extraction and executed import source differ")
    if digest(source) != extraction["sha256"]:
        raise ValueError("Extracted RNA source bytes changed")
    return source


def import_report(root, circ, want, checks):
    out = root / "import_four"
    scope = "four-cell IMPORT"
    payload = json.loads((out / "qc/metrics.json").read_text())
    per_count = want.getnnz(axis=1).tolist()
    per_support = np.asarray(want.sum(axis=1)).ravel().tolist()
    for k, value in dict(
        n_cells=4,
        n_candidates=want.shape[1],
        nnz=want.nnz,
        total_support=int(want.sum()),
    ).items():
        checks.check("JSON_" + k, value, payload[k], scope)
    checks.check(
        "exact_RNA_overlap",
        dict(n_shared=4, n_circ_only=0, n_rna_only=0, rna_genome_build_verified=False),
        payload["rna_overlap"],
        scope,
    )
    for name, numerator in [
        ("detection_fraction_completed", sum(x > 0 for x in per_count)),
        ("zero_fraction_completed", sum(x == 0 for x in per_count)),
    ]:
        checks.check(
            name,
            dict(numerator=numerator, denominator=4, value=numerator / 4),
            payload[name],
            scope,
        )
    parsed = Tables()
    parsed.feed((out / "qc/report.html").read_text())
    text = "".join(parsed.text)
    checks.check(
        "HTML_summary",
        True,
        f"4 metadata cells; {want.shape[1]} candidates; {int(want.sum())} total support."
        in text,
        scope,
    )
    checks.check(
        "IMPORT_scope_label",
        True,
        "not evidence that a detector workflow completed" in text,
        scope,
    )
    for filename in ["cell_qc.tsv", "circ_qc.tsv"]:
        header, rows = table(out / "qc" / filename)
        shown = next(t for t in parsed.tables if t and t[0] == header)
        checks.check(
            "HTML_row_limit:" + filename, min(100, len(rows)), len(shown) - 1, scope
        )
        for i, row in enumerate(rows[:100]):
            for col, actual in zip(header, shown[i + 1]):
                value = row[col]
                if col == "prevalence_completed":
                    value, actual = display_fraction(value), display_fraction(actual)
                checks.check(f"HTML:{filename}:{i}:{col}", value, actual, scope)
        if filename == "cell_qc.tsv":
            checks.check("QC_cell_order", CELLS, [r["cell_id"] for r in rows], scope)
            for i, r in enumerate(rows):
                checks.check(
                    "QC_candidate_count",
                    per_count[i],
                    int(r["circRNA_count"]),
                    r["cell_id"],
                )
                checks.check(
                    "QC_support",
                    per_support[i],
                    int(r["circRNA_total_support"]),
                    r["cell_id"],
                )
                checks.check(
                    "QC_bound_status", "completed", r["detector_status"], r["cell_id"]
                )
        else:
            checks.check(
                "QC_feature_order",
                list(circ.var_names),
                [r["circ_id"] for r in rows],
                scope,
            )
            detected = want.getnnz(axis=0).tolist()
            support = np.asarray(want.sum(axis=0)).ravel().tolist()
            for i, r in enumerate(rows):
                for col, value in dict(
                    n_cells_detected=detected[i],
                    total_support=support[i],
                    n_completed_cells_detected=detected[i],
                    completed_cell_denominator=4,
                ).items():
                    checks.check(col, value, int(r[col]), r["circ_id"])
                checks.check(
                    "prevalence_for_bound_cohort",
                    detected[i] / 4,
                    float(r["prevalence_completed"]),
                    r["circ_id"],
                )
    checks.check("first_100_disclosed", True, "First 100 rows" in text, scope)
    for name, sha in payload["source_sha256"].items():
        checks.check(
            "report_source_bytes:" + name,
            sha,
            digest(import_source_path(root, name)),
            scope,
        )


def workflow_checks(root, per_cell, checks, numerical, audits):
    work = root / "workflow_remaining3"
    ids = CELLS[1:]
    features = list(dict.fromkeys(r["circRNA_ID"] for c in ids for r in per_cell[c]))
    X = oracle(per_cell, ids, features)
    counts = X.getnnz(axis=1).tolist()
    totals = np.asarray(X.sum(axis=1)).ravel().tolist()
    circ = ad.read_h5ad(work / "anndata/circ_counts.h5ad")
    scope = "three-cell WORKFLOW"
    checks.check("AnnData_cell_order", ids, list(circ.obs_names), scope)
    checks.check("AnnData_feature_set", set(features), set(circ.var_names), scope)
    want = oracle(per_cell, ids, list(circ.var_names))
    numerical.exact(
        "every_workflow_integer_value_lossless_widening",
        want,
        circ.X,
        scope,
        exact_integer_csr_values,
    )
    numerical.check("workflow_storage_dtype", "int32", str(circ.X.dtype), scope)
    union = {r["circRNA_ID"]: r for c in ids for r in per_cell[c]}
    for identifier, r in circ.var.iterrows():
        for col, field in [
            ("chrom", "chr"),
            ("start", "circRNA_start"),
            ("end", "circRNA_end"),
            ("strand", "strand"),
        ]:
            expected = union[identifier][field]
            actual = r[col]
            if col in ["start", "end"]:
                expected, actual = int(expected), int(actual)
            numerical.check("workflow_" + col, expected, actual, identifier)
    matrix_ids = (work / "matrix/circ_index.txt").read_text().splitlines()
    matrix_cells = (work / "matrix/cell_index.txt").read_text().splitlines()
    selected = [c for i, c in enumerate(ids) if totals[i] >= 1]
    checks.check(
        "collector_registered_cell_total_filter", selected, matrix_cells, scope
    )
    numerical.exact(
        "workflow_MatrixMarket_all_values",
        oracle(per_cell, selected, matrix_ids),
        mmread(work / "matrix/circ_counts.mtx").T.tocsr(),
        scope,
        exact_integer_csr_values,
    )
    _, cells = table(work / "qc/cell_qc.tsv")
    _, candidates = table(work / "qc/circ_qc.tsv")
    payload = json.loads((work / "qc/metrics.json").read_text())
    combined = json.loads((work / "align/alignment_prepare_summary.json").read_text())
    checks.check(
        "combined_actual_counts", {"aligned": 3}, combined["status_counts"], scope
    )
    checks.check(
        "completed_workflow", "completed", payload["workflow"]["status"], scope
    )
    checks.check(
        "bound_denominator", 3, payload["candidate_prevalence_denominator"], scope
    )
    checks.check("QC_axis", ids, [r["cell_id"] for r in cells], scope)
    checks.check("JSON_axis", ids, [r["cell_id"] for r in payload["per_cell"]], scope)
    for i, cell in enumerate(ids):
        for col, value in dict(
            alignment_status="aligned", detector_status="success"
        ).items():
            for label, actual in [
                ("QC", cells[i][col]),
                ("JSON", payload["per_cell"][i][col]),
                ("AnnData", circ.obs.loc[cell, col]),
            ]:
                checks.check(label + ":" + col, value, actual, cell)
        for col, value in [
            ("circRNA_count", counts[i]),
            ("total_circRNA_support", totals[i]),
        ]:
            checks.check("QC:" + col, value, int(cells[i][col]), cell)
        for col, value in [
            ("candidate_count", counts[i]),
            ("total_support", totals[i]),
        ]:
            checks.check("JSON:" + col, value, payload["per_cell"][i][col], cell)
    support = {c: int(X[:, i].sum()) for i, c in enumerate(features)}
    detected = {c: int(X[:, i].getnnz()) for i, c in enumerate(features)}
    checks.check(
        "candidate_QC_IDs", set(features), {r["circ_id"] for r in candidates}, scope
    )
    for r in candidates:
        for col, value in [
            ("total_support", support[r["circ_id"]]),
            ("n_cells_detected", detected[r["circ_id"]]),
        ]:
            checks.check("candidate_QC:" + col, value, int(r[col]), r["circ_id"])
    checks.check(
        "candidate_JSON_IDs",
        set(features),
        {r["circ_id"] for r in payload["per_candidate"]},
        scope,
    )
    for r in payload["per_candidate"]:
        for col, value in [
            ("total_support", support[r["circ_id"]]),
            ("n_cells_detected", detected[r["circ_id"]]),
            ("prevalence", detected[r["circ_id"]] / 3),
        ]:
            checks.check("candidate_JSON:" + col, value, r[col], r["circ_id"])
    metrics = payload["metrics"]
    for key, value in dict(
        cells=3,
        candidates=len(features),
        total_support=sum(totals),
        candidate_total_support=sum(totals),
        median_candidates_per_cell=median(counts),
        zero_candidate_cells=sum(n == 0 for n in counts),
        zero_candidate_fraction=sum(n == 0 for n in counts) / 3,
    ).items():
        checks.check("JSON_metric:" + key, value, metrics[key]["value"], scope)
    missing = {
        "",
        ".",
        "-",
        "na",
        "n/a",
        "nan",
        "none",
        "null",
        "unknown",
        "unassigned",
    }
    annotated = sum(
        r.get("host_gene", "").strip().lower() not in missing for r in candidates
    )
    checks.check(
        "JSON_host_annotations",
        annotated,
        metrics["host_gene_annotated_candidates"]["value"],
        scope,
    )
    checks.check(
        "JSON_host_coverage",
        annotated / len(features) if features else None,
        metrics["host_gene_coverage"]["value"],
        scope,
    )
    html_workflow(work, payload, counts, list(detected.values()), checks)
    circ.write_h5ad(audits / "workflow_roundtrip.h5ad")
    numerical.exact(
        "workflow_AnnData_read_write_read",
        circ,
        ad.read_h5ad(audits / "workflow_roundtrip.h5ad"),
        scope,
        equal_adata,
    )
    dump(
        audits / "workflow_report_audit.json",
        dict(
            workflow_status=payload["workflow"]["status"],
            warnings=payload["warnings"],
            cohort=ids,
            denominator=3,
            numerical_values_preserved=True,
            collector_cell_columns=matrix_cells,
            native_annotation_policy="Native gene/type/read IDs retained in source and joint import; workflow host fields are separate GTF-derived annotations.",
            new_RNA_summary=payload["optional_summaries"]["rna"],
        ),
    )


def html_workflow(work, payload, counts, detected, checks):
    scope = "three-cell WORKFLOW HTML"
    parsed = Tables()
    parsed.feed((work / "qc/report.html").read_text())

    def display(v):
        return (
            "Not available"
            if v is None
            else str(int(v))
            if isinstance(v, float) and v.is_integer()
            else str(v)
        )

    def rows(header, values):
        if values:
            actual = next(t for t in parsed.tables if t and t[0] == header)
            checks.check(
                "table:" + header[0],
                [[display(v) for v in r] for r in values[:100]],
                actual[1:],
                scope,
            )

    for key, label in [
        ("cells", "QC cells"),
        ("candidates", "circRNA candidates"),
        ("total_support", "Total support"),
        ("median_candidates_per_cell", "Median candidates / cell"),
        ("zero_candidate_fraction", "Zero-count cells"),
        ("host_gene_coverage", "Recorded host-gene coverage"),
    ]:
        m = payload["metrics"][key]
        value = (
            f"{m['value']:.1%}"
            if m["unit"] == "fraction" and m["value"] is not None
            else display(m["value"])
        )
        checks.check(
            "card:" + label,
            value,
            next(c["value"] for c in parsed.cards if c["label"] == label),
            scope,
        )
    rows(
        ["Cell", "Candidates", "Support", "Detector", "Alignment"],
        [
            [
                r[k]
                for k in [
                    "cell_id",
                    "candidate_count",
                    "total_support",
                    "detector_status",
                    "alignment_status",
                ]
            ]
            for r in payload["per_cell"]
        ],
    )
    rows(
        ["circRNA candidate", "Cells detected", "Prevalence", "Support", "Host gene"],
        [
            [
                r["circ_id"],
                r["n_cells_detected"],
                f"{r['prevalence']:.2%}" if r["prevalence"] is not None else None,
                r["total_support"],
                r["host_gene"] or "Not available",
            ]
            for r in payload["per_candidate"]
        ],
    )
    rows(
        ["Consistency check", "Left value", "Right value", "Result"],
        [
            [
                r["check"],
                r.get("left"),
                r.get("right"),
                "Match"
                if r["ok"] is True
                else "Mismatch"
                if r["ok"] is False
                else "Not available",
            ]
            for r in payload["consistency_checks"]
        ],
    )
    rows(
        ["Source", "SHA-256 / availability"],
        [
            [r["path"], r.get("sha256", r.get("reason", "Not available"))]
            for r in payload["sources"]
        ],
    )
    rows(
        ["RNA summary field", "Recorded value"],
        [[k, v] for k, v in payload["optional_summaries"]["rna"].items()],
    )
    for r in payload["sources"]:
        if "sha256" in r:
            checks.check(
                "source:" + r["path"], r["sha256"], digest(work / r["path"]), scope
            )
    for i, values in enumerate([counts, detected]):
        frequency = Counter(values)
        if len(frequency) > 12:
            raise ValueError("Unexpected pilot histogram scope; review separately")
        checks.check(
            "histogram:" + str(i),
            [f"{n}: {frequency[n]}" for n in sorted(frequency)],
            parsed.svgs[i]["titles"][1:],
            scope,
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--original", type=Path, required=True)
    ap.add_argument("--audit-directory", default="records/acceptance")
    args = ap.parse_args()
    root = args.root.resolve()
    original = args.original.resolve()
    audits = root / args.audit_directory
    audits.mkdir(parents=True, exist_ok=False)
    checks = Checks(audits / "native_value_checks.tsv")
    reports = Checks(audits / "report_value_checks.tsv")
    sources = json.loads((root / "records/joint_native_sources.json").read_text())
    try:
        per_cell, union = native_union(
            {r["cell_id"]: Path(r["native_path"]) for r in sources["cells"]}
        )
        circ, want = import_checks(root, original, per_cell, union, checks, audits)
        import_report(root, circ, want, reports)
        workflow_checks(root, per_cell, reports, checks, audits)
        dump(
            audits / "cohort_summary.json",
            dict(
                recorded_utc=stamp(),
                passed=True,
                unique_candidates=len(union),
                per_cell_candidates={c: len(v) for c, v in per_cell.items()},
                per_cell_support={
                    c: sum(int(r["#junction_reads"]) for r in v)
                    for c, v in per_cell.items()
                },
                total_support=int(want.sum()),
                nnz=want.nnz,
                native_checks=len(checks.rows),
                report_checks=len(reports.rows),
                states_not_observed_in_real_pilot=[
                    "failed",
                    "unknown",
                    "unprocessed",
                    "completed-zero",
                ]
                if all(want.getnnz(axis=1) > 0)
                else ["failed", "unknown", "unprocessed"],
                mixed_execution_lineage=True,
                historical_provenance_recovered=False,
                biological_accuracy_assessed=False,
            ),
        )
    except Exception as exc:
        dump(
            audits / "acceptance_failure.json",
            dict(
                error=repr(exc),
                recorded_utc=stamp(),
                completed_native_checks=len(checks.rows),
                completed_report_checks=len(reports.rows),
            ),
        )
        raise
    print(
        json.dumps(
            {
                "native_checks": len(checks.rows),
                "report_checks": len(reports.rows),
                "passed": True,
            }
        )
    )


if __name__ == "__main__":
    main()
