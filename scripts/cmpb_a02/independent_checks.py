#!/usr/bin/env python3
"""A02 independent CSV/exact-integer/CSR checks. No circyto imports.

Only this first-cell contract is covered. This is not an alternative caller or
the A03 reference-conversion benchmark. Synthetic unit fixtures are not data.
"""

from __future__ import annotations

import argparse
import copy
from collections import Counter
import csv
from html.parser import HTMLParser
import json
from pathlib import Path
import re

import anndata as ad
import h5py
import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from run_stage import digest, dump, stamp

CELL = "IMR90_A_100"
SOURCE_SHA = "bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8"
HEADER = [
    "circRNA_ID",
    "chr",
    "circRNA_start",
    "circRNA_end",
    "#junction_reads",
    "SM_MS_SMS",
    "#non_junction_reads",
    "junction_reads_ratio",
    "circRNA_type",
    "gene_id",
    "strand",
    "junction_reads_ID",
    "Score",
]


def native(path):
    records, seen = [], set()
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != HEADER:
            raise ValueError("Unsupported native schema")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed native row")
            for key in ["circRNA_start", "circRNA_end", "#junction_reads"]:
                if not re.fullmatch(r"[0-9]+", row[key]):
                    raise ValueError("Native value is not an exact nonnegative integer")
            start, end, support = [
                int(row[key])
                for key in ["circRNA_start", "circRNA_end", "#junction_reads"]
            ]
            if not 1 <= start <= end or support > np.iinfo(np.int64).max:
                raise ValueError("Invalid native interval/count")
            if row["circRNA_ID"] != f"{row['chr']}:{start}|{end}" or row[
                "strand"
            ] not in ["+", "-"]:
                raise ValueError("Native ID/coordinate/strand mismatch")
            if row["circRNA_ID"] in seen:
                raise ValueError("Duplicate native candidate")
            seen.add(row["circRNA_ID"])
            records.append(row)
    counts = np.asarray(
        [int(row["#junction_reads"]) for row in records], dtype=np.int64
    )
    X = sparse.csr_matrix(
        (counts, (np.zeros(len(records), dtype=int), np.arange(len(records)))),
        shape=(1, len(records)),
        dtype=np.int64,
    )
    X.eliminate_zeros()
    return records, X


def exact(a, b, location="value"):
    """Exact values, sparse indices and NA semantics without dense matrix expansion."""
    if sparse.issparse(a) or sparse.issparse(b):
        if not sparse.issparse(a) or not sparse.issparse(b):
            raise AssertionError(location + ": sparse representation changed")
        a, b = a.tocsr(copy=True), b.tocsr(copy=True)
        if a.shape != b.shape or a.dtype != b.dtype:
            raise AssertionError(location + ": sparse shape/dtype changed")
        a.sum_duplicates()
        b.sum_duplicates()
        a.eliminate_zeros()
        b.eliminate_zeros()
        a.sort_indices()
        b.sort_indices()
        for key in ["data", "indices", "indptr"]:
            np.testing.assert_array_equal(
                getattr(a, key), getattr(b, key), err_msg=location + "." + key
            )
    elif isinstance(a, pd.DataFrame) and isinstance(b, pd.DataFrame):
        pd.testing.assert_frame_equal(a, b, check_exact=True)
    elif isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            raise AssertionError(location + ": keys changed")
        for key in a:
            exact(a[key], b[key], location + "." + str(key))
    elif a is None or b is None:
        if a is not b:
            raise AssertionError(location + ": None changed")
    else:
        np.testing.assert_array_equal(np.asarray(a), np.asarray(b), err_msg=location)


def equal_adata(a, b, label="AnnData"):
    for key in ["X", "obs", "var", "uns"]:
        exact(getattr(a, key), getattr(b, key), label + "." + key)
    for key in ["layers", "obsm", "varm", "obsp", "varp"]:
        exact(dict(getattr(a, key)), dict(getattr(b, key)), label + "." + key)
    if a.raw is None or b.raw is None:
        if a.raw is not b.raw:
            raise AssertionError(label + ": raw presence changed")
    else:
        exact(a.raw.X, b.raw.X, label + ".raw.X")
        exact(a.raw.var, b.raw.var, label + ".raw.var")
        exact(dict(a.raw.varm), dict(b.raw.varm), label + ".raw.varm")


def exact_integer_csr_values(expected, actual, label):
    """Compare integer values across storage widths without accepting floats.

    The legacy workflow writer explicitly stores int32; the native importer
    stores int64. Widening an integer for this read-only check is lossless.
    Retain strict dtype equality for native import and object round trips.
    """
    if not sparse.issparse(actual) or actual.dtype.kind not in "iu":
        raise AssertionError(label + ": an integer sparse matrix is required")
    if actual.data.size and (
        int(actual.data.min()) < 0 or int(actual.data.max()) > np.iinfo(np.int64).max
    ):
        raise AssertionError(label + ": invalid integer support range")
    exact(expected, actual.astype(np.int64), label)


def original_slice(root):
    source = root / "inputs/original/full_length_rna_circ_cnv.hostgene_fixed.h5mu"
    if digest(source) != SOURCE_SHA:
        raise ValueError("Original object hash mismatch")
    with h5py.File(source, "r") as handle:
        rna = ad.io.read_elem(handle["mod/rna"])
    if list(rna.obs_names).count(CELL) != 1:
        raise ValueError("First-cell identity missing or duplicated in original RNA")
    for field, value in {
        "gsm_id": "GSM8558852",
        "canonical_cell_id": CELL,
        "rna_cell_id": "RNA_IMR90_A_100",
        "dna_cell_id": "DNA_IMR90_A_100",
    }.items():
        if field not in rna.obs or str(rna.obs.loc[CELL, field]) != value:
            raise ValueError("Original RNA biological identity disagrees: " + field)
    result = rna[[CELL], :].copy()
    # AnnData slicing can remove unused categories; preserve source declarations.
    result.obs = rna.obs.loc[[CELL]].copy()
    result.var = rna.var.copy()
    result.uns = copy.deepcopy(rna.uns)
    return result


def extract_rna(root):
    original = original_slice(root)
    target = root / "inputs/rna/imr90_original_overlap_profile.firstcell.h5ad"
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    original.write_h5ad(target)
    equal_adata(original, ad.read_h5ad(target), "original RNA extraction")
    # Preserve original parent-object declarations separately. Copying them into
    # a new genome_build field could falsely suggest historical verification.
    source = root / "inputs/original/full_length_rna_circ_cnv.hostgene_fixed.h5mu"
    metadata_copy = target.with_name("original_parent_metadata.h5")
    with h5py.File(source, "r") as source_handle, h5py.File(metadata_copy, "x") as dest:
        source_handle.copy("uns", dest, name="original_mudata_uns")
        exact(
            ad.io.read_elem(source_handle["uns"]),
            ad.io.read_elem(dest["original_mudata_uns"]),
            "original parent metadata",
        )
    dump(
        root / "records/rna_extraction.json",
        {
            "recorded_utc": stamp(),
            "source_sha256": SOURCE_SHA,
            "source_modality": "/mod/rna",
            "cell_id": CELL,
            "shape": list(original.shape),
            "output": str(target.relative_to(root)),
            "output_sha256": digest(target),
            "all_available_rna_fields_equal": True,
            "profile": "historical same-cell lightweight simple-overlap",
            "historical_reference_hash": "unknown",
            "new_reference_verifies_old_profile": False,
            "metadata_policy": "RNA metadata unchanged; extraction provenance kept separately here",
            "original_parent_metadata": str(metadata_copy.relative_to(root)),
            "original_parent_metadata_sha256": digest(metadata_copy),
            "parent_metadata_evidence_class": "original processed-object declarations, not recovered original run records",
        },
    )
    print(json.dumps({"original_rna_extracted": list(original.shape), "cell_id": CELL}))


def write_checks(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["check", "scope", "expected", "observed", "status"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


class Tables(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables, self.text = [], []
        self.table = self.row = self.cell = None
        self.cards, self.svgs = [], []
        self.card = self.card_field = self.svg = self.svg_field = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "article" and attrs.get("class") == "metric":
            self.card = {}
        if (
            self.card is not None
            and tag == "div"
            and attrs.get("class") in ["label", "value"]
        ):
            self.card_field = attrs["class"]
            self.card[self.card_field] = []
        if tag == "svg":
            self.svg = {"label": attrs.get("aria-label"), "titles": [], "texts": []}
        if self.svg is not None and tag in ["title", "text"]:
            self.svg_field = (tag, [])
        if tag == "table":
            self.table = []
        if tag == "tr":
            self.row = []
        if tag in ["th", "td"]:
            self.cell = []

    def handle_data(self, data):
        self.text.append(data)
        if self.card_field is not None:
            self.card[self.card_field].append(data)
        if self.svg_field is not None:
            self.svg_field[1].append(data)
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag == "div":
            self.card_field = None
        if tag == "article" and self.card is not None:
            self.cards.append({k: "".join(v) for k, v in self.card.items()})
            self.card = None
        if self.svg_field is not None and tag == self.svg_field[0]:
            self.svg["titles" if tag == "title" else "texts"].append(
                "".join(self.svg_field[1])
            )
            self.svg_field = None
        if tag == "svg" and self.svg is not None:
            self.svgs.append(self.svg)
            self.svg = None
        if tag in ["th", "td"] and self.cell is not None:
            self.row.append("".join(self.cell).strip())
            self.cell = None
        if tag == "tr" and self.row is not None:
            self.table.append(self.row)
            self.row = None
        if tag == "table" and self.table is not None:
            self.tables.append(self.table)
            self.table = None


def verify_native_object(records, expected, circ):
    if list(circ.obs_names) != [CELL] or list(circ.var_names) != [
        r["circRNA_ID"] for r in records
    ]:
        raise AssertionError(
            "Cell/candidate axis order differs from the one-cell native contract"
        )
    exact(expected, circ.X, "native vs imported CSR")
    for row in records:
        actual = circ.var.loc[row["circRNA_ID"]]
        for target, source in [
            ("circ_id", "circRNA_ID"),
            ("chrom", "chr"),
            ("strand", "strand"),
            ("ciri3_gene_id", "gene_id"),
            ("ciri3_circRNA_type", "circRNA_type"),
        ]:
            if actual[target] != row[source]:
                raise AssertionError(
                    "Native annotation changed: " + row["circRNA_ID"] + "/" + target
                )
        for target, source in [("start", "circRNA_start"), ("end", "circRNA_end")]:
            if int(actual[target]) != int(row[source]):
                raise AssertionError("Native coordinate changed")
        gene = "" if row["gene_id"] in ["", "NA"] else row["gene_id"]
        if actual["host_gene_id"] != gene or actual["host_gene_id_source"] != (
            "CIRI3_gene_id" if gene else ""
        ):
            raise AssertionError("Native gene ID mapping changed")


def acceptance(root):
    out = root / "import_firstcell"
    audits = root / "records/acceptance"
    audits.mkdir(parents=True, exist_ok=True)
    rows = []

    def passed(check, scope, expected="exact equality", observed="exact equality"):
        rows.append(
            dict(
                check=check,
                scope=scope,
                expected=expected,
                observed=observed,
                status="PASS",
            )
        )
        write_checks(audits / "native_value_checks.tsv", rows)

    native_path = (
        root / "workflow_firstcell/ciri3" / (CELL + ".ciri3_run/ciri3_raw.tsv")
    )
    records, expected = native(native_path)
    circ = ad.read_h5ad(out / "anndata/circ_counts.h5ad")
    verify_native_object(records, expected, circ)
    passed("native_CSR_axes_all_fields", f"{len(records)} native candidates; one cell")
    for index, row in enumerate(records):
        passed(
            "candidate_value_and_annotation",
            row["circRNA_ID"],
            row["#junction_reads"],
            str(int(circ.X[0, index])),
        )
    source_copy = out / "evidence/ciri3_native_000001.tsv"
    if digest(source_copy) != digest(native_path):
        raise AssertionError("Native ancillary evidence bytes changed")
    passed(
        "native_source_bytes_including_read_IDs",
        "complete original native file",
        digest(native_path),
        digest(source_copy),
    )
    mm = mmread(out / "matrix/circ_counts.mtx").T.tocsr()
    exact(expected, mm, "MatrixMarket transpose")
    if (out / "matrix/cell_index.txt").read_text().splitlines() != [CELL] or (
        out / "matrix/circ_index.txt"
    ).read_text().splitlines() != list(circ.var_names):
        raise AssertionError("Text matrix axis differs")
    passed("matrix_market_values_and_text_axes", "all entries")
    rna = ad.read_h5ad(
        root / "inputs/rna/imr90_original_overlap_profile.firstcell.h5ad"
    )
    equal_adata(original_slice(root), rna, "original source slice")
    passed(
        "source_RNA_all_fields",
        f"{rna.n_vars} features; all X/layers/annotations/embeddings/raw",
    )
    circ_roundtrip = audits / "circ_roundtrip.h5ad"
    circ.write_h5ad(circ_roundtrip)
    equal_adata(circ, ad.read_h5ad(circ_roundtrip))
    passed("ordinary_AnnData_read_write_read", "all available fields")
    multi = md.read_h5mu(out / "mudata/circyto_multimodal.h5mu")
    if list(multi.obs_names) != [CELL] or set(multi.mod) != {"circ", "rna"}:
        raise AssertionError("MuData modality/cell axis mismatch")
    equal_adata(circ, multi.mod["circ"], "MuData circ")
    equal_adata(rna, multi.mod["rna"], "MuData RNA")
    passed("same_cell_RNA_in_MuData", "exact biological identity and all fields")
    rt = audits / "multimodal_roundtrip.h5mu"
    multi.write_h5mu(rt)
    reread = md.read_h5mu(rt)
    for key in multi.mod:
        equal_adata(multi.mod[key], reread.mod[key], "MuData roundtrip " + key)
    for key in ["obs", "var", "uns"]:
        exact(getattr(multi, key), getattr(reread, key), "MuData global " + key)
    for key in ["obsm", "varm", "obsp", "varp", "obsmap", "varmap"]:
        exact(
            dict(getattr(multi, key)),
            dict(getattr(reread, key)),
            "MuData global " + key,
        )
    passed(
        "ordinary_MuData_read_write_read",
        "modalities/global axes/maps/metadata/NA masks",
    )
    report_checks(root, records, expected, circ)
    workflow_checks(root, records, expected)
    dump(
        audits / "native_summary.json",
        {
            "native_sha256": digest(native_path),
            "native_records": len(records),
            "support": sum(int(r["#junction_reads"]) for r in records),
            "nnz": expected.nnz,
            "checks_passed": len(rows),
            "historical_circ_yield_compared": False,
            "new_reference_verifies_historical_RNA": False,
        },
    )


def report_checks(root, records, expected, circ):
    """All displayed IMPORT table values; workflow JSON/HTML has a separate audit."""
    out = root / "import_firstcell"
    metrics = json.loads((out / "qc/metrics.json").read_text())
    checks = []

    def check(name, expected_value, observed, scope="IMPORT"):
        success = expected_value == observed
        checks.append(
            dict(
                check=name,
                scope=scope,
                expected=str(expected_value),
                observed=str(observed),
                status="PASS" if success else "FAIL",
            )
        )
        write_checks(root / "records/acceptance/report_value_checks.tsv", checks)
        if not success:
            raise AssertionError("Report disagreement: " + name)

    total = sum(int(r["#junction_reads"]) for r in records)
    for key, value in [
        ("n_cells", 1),
        ("n_candidates", len(records)),
        ("nnz", expected.nnz),
        ("total_support", total),
    ]:
        check(key, value, metrics[key])
    check("exact_RNA_overlap", 1, metrics["rna_overlap"]["n_shared"])
    check(
        "historical_RNA_reference_not_verified",
        False,
        metrics["rna_overlap"]["rna_genome_build_verified"],
    )
    check(
        "completed_denominator",
        1,
        metrics["detection_fraction_completed"]["denominator"],
    )
    check(
        "positive_completed_numerator",
        int(total > 0),
        metrics["detection_fraction_completed"]["numerator"],
    )
    parsed = Tables()
    parsed.feed((out / "qc/report.html").read_text())
    check(
        "HTML_summary_present",
        True,
        f"1 metadata cells; {len(records)} candidates; {total} total support."
        in "".join(parsed.text),
    )
    check(
        "HTML_import_provenance_label",
        True,
        "not evidence that a detector workflow completed" in "".join(parsed.text),
    )
    for filename, axis_name in [("cell_qc.tsv", "cell_id"), ("circ_qc.tsv", "circ_id")]:
        with (out / "qc" / filename).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            table = list(reader)
            header = reader.fieldnames
        displayed = next(t for t in parsed.tables if t and t[0] == header)
        check(
            "HTML_" + filename + "_displayed_rows",
            min(100, len(table)),
            len(displayed) - 1,
        )
        for i, row in enumerate(table[:100]):
            for col, value in zip(header, displayed[i + 1]):
                # Pandas HTML renders floating metrics to six decimals. Compare
                # exact strings for IDs/annotations and the documented rounding
                # for numeric table columns only; original TSV remains exact.
                original = row[col]
                if value != original:
                    if col not in ["prevalence_completed"]:
                        raise AssertionError(
                            "HTML table value changed: " + filename + "/" + col
                        )
                    check("HTML_" + col + "_" + str(i), f"{float(original):.6f}", value)
        check("HTML_" + filename + "_all_displayed_fields", True, True)
    check(
        "HTML_candidate_limit_disclosed", True, "First 100 rows" in "".join(parsed.text)
    )
    dump(
        root / "records/acceptance/import_report_audit.json",
        {
            "checks": len(checks),
            "HTML_tables_checked": len(parsed.tables),
            "visible_candidate_rows": min(100, len(records)),
            "full_candidate_rows": len(records),
            "browser_render": "pending separate browser check",
            "workflow_report": "pending separate audit",
        },
    )


def workflow_checks(root, records, expected):
    """Compare genuine workflow values separately, allowing documented zero columns."""
    work = root / "workflow_firstcell"
    path = root / "records/acceptance/report_value_checks.tsv"
    with path.open(newline="") as handle:
        checks = list(csv.DictReader(handle, delimiter="\t"))

    def check(name, want, observed):
        checks.append(
            dict(
                check=name,
                scope="NEW WORKFLOW",
                expected=str(want),
                observed=str(observed),
                status="PASS" if want == observed else "FAIL",
            )
        )
        write_checks(path, checks)
        if want != observed:
            raise AssertionError("Workflow disagreement: " + name)

    def table(name):
        with (work / name).open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    total = sum(int(r["#junction_reads"]) for r in records)
    native_by_id = {r["circRNA_ID"]: r for r in records}
    cells, candidates = table("qc/cell_qc.tsv"), table("qc/circ_qc.tsv")
    check("QC_cell_axis", [CELL], [r["cell_id"] for r in cells])
    check("QC_support", total, int(cells[0]["total_circRNA_support"]))
    check("QC_positive_candidates", expected.nnz, int(cells[0]["circRNA_count"]))
    check("QC_candidate_IDs", set(native_by_id), {r["circ_id"] for r in candidates})
    for row in candidates:
        raw = native_by_id[row["circ_id"]]
        check(
            "candidate_support:" + row["circ_id"],
            int(raw["#junction_reads"]),
            int(row["total_support"]),
        )
        check(
            "candidate_detected_cells:" + row["circ_id"],
            int(int(raw["#junction_reads"]) > 0),
            int(row["n_cells_detected"]),
        )
    workflow_circ = ad.read_h5ad(work / "anndata/circ_counts.h5ad")
    check("workflow_AnnData_cell_axis", [CELL], list(workflow_circ.obs_names))
    lookup = {r["circRNA_ID"]: i for i, r in enumerate(records)}
    order = [lookup[x] for x in workflow_circ.var_names]
    exact_integer_csr_values(
        expected[:, order], workflow_circ.X, "workflow selected-cell CSR"
    )
    check("workflow_CSR_all_values", True, True)
    check("workflow_CSR_storage_dtype", "int32", str(workflow_circ.X.dtype))
    for key, source in [
        ("chrom", "chr"),
        ("start", "circRNA_start"),
        ("end", "circRNA_end"),
        ("strand", "strand"),
    ]:
        for identifier, row in workflow_circ.var.iterrows():
            wanted = native_by_id[identifier][source]
            actual = row[key]
            if key in ["start", "end"]:
                wanted, actual = int(wanted), int(actual)
            if wanted != actual:
                raise AssertionError("Workflow coordinate/strand changed")
    check("workflow_native_coordinates_strands", True, True)
    matrix_cells = (work / "matrix/cell_index.txt").read_text().splitlines()
    check(
        "collector_columns_under_registered_cell_total_filter",
        [CELL] if total >= 1 else [],
        matrix_cells,
    )
    matrix_ids = (work / "matrix/circ_index.txt").read_text().splitlines()
    matrix = mmread(work / "matrix/circ_counts.mtx").T.tocsr()
    matrix_order = [lookup[x] for x in matrix_ids]
    selected = (
        expected[:, matrix_order]
        if total >= 1
        else sparse.csr_matrix((0, len(matrix_ids)), dtype=np.int64)
    )
    exact_integer_csr_values(selected, matrix, "workflow MatrixMarket")
    check("workflow_matrix_values_and_axis", True, True)
    payload = json.loads((work / "qc/metrics.json").read_text())
    metrics = payload["metrics"]
    for key, want in [
        ("cells", 1),
        ("candidates", len(records)),
        ("total_support", total),
        ("candidate_total_support", total),
        ("median_candidates_per_cell", expected.nnz),
        ("zero_candidate_cells", int(total == 0)),
    ]:
        check("JSON_" + key, want, metrics[key]["value"])
    host_missing = {
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
        r.get("host_gene", "").strip().lower() not in host_missing for r in candidates
    )
    check(
        "JSON_recorded_host_annotations",
        annotated,
        metrics["host_gene_annotated_candidates"]["value"],
    )
    check(
        "JSON_host_coverage",
        annotated / len(records) if records else None,
        metrics["host_gene_coverage"]["value"],
    )
    check(
        "JSON_cell_candidate_values",
        int(cells[0]["circRNA_count"]),
        payload["per_cell"][0]["candidate_count"],
    )
    for item in payload["per_candidate"]:
        check(
            "JSON_candidate_support:" + item["circ_id"],
            int(native_by_id[item["circ_id"]]["#junction_reads"]),
            item["total_support"],
        )
    parsed = Tables()
    parsed.feed((work / "qc/report.html").read_text())

    def display(value):
        if value is None:
            return "Not available"
        if isinstance(value, float) and value.is_integer():
            return str(int(value))
        return str(value)

    def html_table(header, rows):
        if not rows:
            return
        actual = next(t for t in parsed.tables if t and t[0] == header)
        want = [[display(v) for v in row] for row in rows[:100]]
        check("HTML_table:" + header[0], want, actual[1:])

    for key, label in [
        ("cells", "QC cells"),
        ("candidates", "circRNA candidates"),
        ("total_support", "Total support"),
        ("median_candidates_per_cell", "Median candidates / cell"),
        ("zero_candidate_fraction", "Zero-count cells"),
        ("host_gene_coverage", "Recorded host-gene coverage"),
    ]:
        metric = metrics[key]
        want = (
            f"{metric['value']:.1%}"
            if metric["unit"] == "fraction" and metric["value"] is not None
            else display(metric["value"])
        )
        check(
            "HTML_metric:" + label,
            want,
            next(c["value"] for c in parsed.cards if c["label"] == label),
        )
    html_table(
        ["Cell", "Candidates", "Support", "Detector", "Alignment"],
        [
            [
                r["cell_id"],
                r["candidate_count"],
                r["total_support"],
                r["detector_status"],
                r["alignment_status"],
            ]
            for r in payload["per_cell"]
        ],
    )
    html_table(
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
    html_table(
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
    html_table(
        ["Source", "SHA-256 / availability"],
        [
            [r["path"], r.get("sha256", r.get("reason", "Not available"))]
            for r in payload["sources"]
        ],
    )
    for source in payload["sources"]:
        if "sha256" in source:
            check(
                "report_source_bytes:" + source["path"],
                source["sha256"],
                digest(work / source["path"]),
            )
    html_table(
        ["RNA summary field", "Recorded value"],
        [[k, v] for k, v in payload["optional_summaries"]["rna"].items()],
    )
    # First-cell histograms have at most two bins, so no bin-width inference.
    check(
        "HTML_cell_histogram_bars", [f"{expected.nnz}: 1"], parsed.svgs[0]["titles"][1:]
    )
    if payload["candidate_prevalence_denominator"] is not None and records:
        counts = Counter(int(int(r["#junction_reads"]) > 0) for r in records)
        check(
            "HTML_candidate_histogram_bars",
            [f"{k}: {counts[k]}" for k in sorted(counts)],
            parsed.svgs[1]["titles"][1:],
        )
    dump(
        root / "records/acceptance/workflow_report_audit.json",
        {
            "workflow_status": payload["workflow"]["status"],
            "report_warnings": payload["warnings"],
            "collector_cell_columns": matrix_cells,
            "selected_QC_cells": [CELL],
            "zero_cell_column_policy": "cell-total >=1; workflow selected-cell H5AD may expand omitted zero columns",
            "native_vs_workflow_annotations": "Native gene/type/read IDs preserved in original native bytes and native import; workflow host annotations have their own GTF-derived fields, not native gene-symbol equivalence.",
            "original_workflow_recovered": False,
            "new_RNA_summary": payload["optional_summaries"]["rna"],
        },
    )


def status_comparisons(payload, recorded_alignment_status, object_status):
    """Separate semantic status checks from HTML-to-JSON rendering agreement."""
    observed = payload["per_cell"]
    if len(observed) != 1 or observed[0]["cell_id"] != CELL:
        raise AssertionError("First-cell report axis changed")
    return [
        (
            "recorded_alignment_vs_report",
            recorded_alignment_status,
            observed[0]["alignment_status"],
        ),
        (
            "recorded_alignment_vs_workflow_AnnData",
            recorded_alignment_status,
            object_status,
        ),
        (
            "bound_completed_workflow_report_state",
            "completed",
            payload["workflow"]["status"],
        ),
        (
            "bound_completed_cohort_denominator",
            1,
            payload["candidate_prevalence_denominator"],
        ),
    ]


def status_acceptance(root):
    """Audit the observed A02 status defect without editing production outputs."""
    audits = root / "records/acceptance"
    binding_path = root / "records/native_completion_binding.json"
    binding = json.loads(binding_path.read_text())
    if binding.get("new_detector_status") != "completed":
        raise AssertionError("No independently bound completed new cell")
    work = root / "workflow_firstcell"
    align_path = work / "align/bwa_mem/alignment_prepare_summary.json"
    align = json.loads(align_path.read_text())
    if len(align["cells"]) != 1 or align["cells"][0]["cell_id"] != CELL:
        raise AssertionError("Aligner-specific cell record mismatch")
    status = align["cells"][0]["status"]
    if status != "aligned" or align["status_counts"] != {"aligned": 1}:
        raise AssertionError("Completed new alignment is not recorded")
    payload_path = work / "qc/metrics.json"
    payload = json.loads(payload_path.read_text())
    circ_path = work / "anndata/circ_counts.h5ad"
    object_status = str(ad.read_h5ad(circ_path).obs.loc[CELL, "alignment_status"])
    comparisons = status_comparisons(payload, status, object_status)
    path = audits / "report_value_checks.tsv"
    with path.open(newline="") as handle:
        checks = list(csv.DictReader(handle, delimiter="\t"))
    new_checks = [
        dict(
            check=name,
            scope="NEW WORKFLOW STATUS SEMANTICS",
            expected=str(want),
            observed=str(actual),
            status="PASS" if want == actual else "FAIL",
        )
        for name, want, actual in comparisons
    ]
    if any(row["check"] in {r["check"] for r in checks} for row in new_checks):
        raise ValueError("Status acceptance already recorded")
    write_checks(path, checks + new_checks)
    sources = [
        binding_path,
        align_path,
        payload_path,
        circ_path,
        work / "align/alignment_prepare_summary.json",
        work / "qc/cell_qc.tsv",
    ]
    dump(
        audits / "status_acceptance.json",
        {
            "recorded_utc": stamp(),
            "checks": new_checks,
            "source_sha256": {str(p.relative_to(root)): digest(p) for p in sources},
            "production_outputs_modified": False,
            "finding": "Aligner-specific summary records aligned; root alignment summary omits per-cell/status counts. Workflow QC/AnnData contain missing alignment status; report withholds completed-cohort prevalence.",
        },
    )
    if any(row["status"] != "PASS" for row in new_checks):
        raise AssertionError(
            "Workflow alignment-status propagation failed acceptance; outputs retained unchanged"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action", choices=["extract-rna", "acceptance", "status-acceptance"]
    )
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    {
        "extract-rna": extract_rna,
        "acceptance": acceptance,
        "status-acceptance": status_acceptance,
    }[args.action](args.root.resolve())


if __name__ == "__main__":
    main()
