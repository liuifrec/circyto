"""Independent per-record output verification. Never imports the reference converter."""

import csv
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
from scipy import sparse
from scipy.io import mmread

# Reuse reviewed object-field/HTML assertions solely as an oracle, never a converter.
_spec = importlib.util.spec_from_file_location(
    "a02_oracle", Path(__file__).parents[1] / "cmpb_a02/independent_checks.py"
)
_old = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_old)


def read_rows(p):
    with Path(p).open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def truth(config):
    cells = read_rows(config["cells"])
    ids = [r["cell_id"] for r in cells]
    status = {r["cell_id"]: r.get("detector_status", "unknown") for r in cells}
    vals = {}
    features = {}
    sources = {"cells.tsv": Path(config["cells"])}
    if config["mode"] == "native":
        sources["ciri3_manifest.tsv"] = Path(config["input"])
        for i, cell in enumerate(read_rows(config["input"]), 1):
            p = (Path(config["input"]).parent / Path(cell["path"])).resolve()
            sources[f"ciri3_native_{i:06d}.tsv"] = p
            for r in read_rows(p):
                f = r["circRNA_ID"]
                v = int(r["#junction_reads"])
                vals[cell["cell_id"], f] = v
                features.setdefault(
                    f,
                    dict(
                        circ_id=f,
                        chrom=r["chr"],
                        start=int(r["circRNA_start"]),
                        end=int(r["circRNA_end"]),
                        strand=r["strand"],
                        ciri3_gene_id=r["gene_id"],
                        ciri3_circRNA_type=r["circRNA_type"],
                    ),
                )
    elif config["mode"] == "calls":
        sources["calls.tsv"] = Path(config["input"])
        for r in read_rows(config["input"]):
            f = r["circ_id"]
            vals[r["cell_id"], f] = int(r["support"])
            features.setdefault(
                f,
                {
                    k: (int(v) if k in ["start", "end"] else v)
                    for k, v in r.items()
                    if k not in ["cell_id", "support"]
                },
            )
    else:
        folder = Path(config["input"])
        for n in [
            "circ_counts.mtx",
            "cell_index.txt",
            "circ_index.txt",
            "circ_feature_table.tsv",
        ]:
            sources[n] = folder / n
        source_ids = (folder / "cell_index.txt").read_text().splitlines()
        fids = (folder / "circ_index.txt").read_text().splitlines()
        raw = {r["circ_id"]: r for r in read_rows(folder / "circ_feature_table.tsv")}
        features = {
            f: {k: (int(v) if k in ["start", "end"] else v) for k, v in raw[f].items()}
            for f in fids
        }
        x = mmread(folder / "circ_counts.mtx").tocoo()
        for j, i, v in zip(x.row, x.col, x.data):
            vals[source_ids[i], fids[j]] = int(v)
    if config.get("rna"):
        sources["rna.h5ad"] = Path(config["rna"])
    vals = {k: v for k, v in vals.items() if v > 0}
    return ids, status, features, vals, sources


def actual_values(a):
    x = a.X.tocoo()
    return {
        (str(a.obs_names[i]), str(a.var_names[j])): int(v)
        for i, j, v in zip(x.row, x.col, x.data)
        if v != 0
    }


def verify(config, out, arm):
    ids, status, features, want, sources = truth(config)
    out = Path(out)
    result = dict(
        count_errors=0,
        annotation_errors=0,
        cell_assignment_errors=0,
        status_errors=0,
        report_errors=0,
        source_errors=0,
        contract_errors=0,
        unavailable_fields=[],
        undefined_candidates=0,
    )
    a = ad.read_h5ad(out / "anndata/circ_counts.h5ad")
    got = actual_values(a)
    result["count_errors"] = sum(
        want.get(k, 0) != got.get(k, 0) for k in set(want) | set(got)
    )
    result["contract_errors"] += (
        int(list(a.obs_names) != ids)
        + int(list(a.var_names) != list(features))
        + int(str(a.X.dtype) != "int64")
    )
    for f, record in features.items():
        for k, v in record.items():
            if k not in a.var:
                result["unavailable_fields"].append("var." + k)
                continue
            av = a.var.at[f, k] if f in a.var.index else None
            result["annotation_errors"] += int(str(av) != str(v))
    completed = {c for c in ids if status[c] == "completed"}
    den = len(completed)
    counts = {c: sum(k[0] == c for k in want) for c in ids}
    support = {c: sum(v for k, v in want.items() if k[0] == c) for c in ids}
    for c in ids:
        if c not in a.obs_names:
            result["cell_assignment_errors"] += 1
            continue
        expected_status = (
            ("completed_detected" if counts[c] else "completed_zero")
            if c in completed
            else status[c]
        )
        result["status_errors"] += int(
            a.obs.at[c, "detector_status"] != status[c]
        ) + int(a.obs.at[c, "circ_observation_status"] != expected_status)
        result["report_errors"] += int(
            int(a.obs.at[c, "circRNA_count"]) != counts[c]
        ) + int(int(a.obs.at[c, "circRNA_total_support"]) != support[c])
    report = json.loads((out / "qc/metrics.json").read_text())
    for k, v in dict(
        n_cells=len(ids),
        n_candidates=len(features),
        nnz=len(want),
        total_support=sum(want.values()),
    ).items():
        result["report_errors"] += int(report.get(k) != v)
    n = sum(counts[c] > 0 for c in completed)
    expected = dict(numerator=n, denominator=den, value=n / den if den else None)
    result["report_errors"] += int(
        report.get("detection_fraction_completed") != expected
    )
    qc = {r["circ_id"]: r for r in read_rows(out / "qc/circ_qc.tsv")}
    cellqc = {r["cell_id"]: r for r in read_rows(out / "qc/cell_qc.tsv")}
    for c in ids:
        for k, v in [
            ("circRNA_count", counts[c]),
            ("circRNA_total_support", support[c]),
        ]:
            result["report_errors"] += int(int(cellqc[c][k]) != v)
    prevalence = []
    for f in features:
        numerator = sum((c, f) in want for c in completed)
        observed_num = int(qc[f]["n_completed_cells_detected"])
        observed_den = int(qc[f]["completed_cell_denominator"])
        result["report_errors"] += int(observed_num != numerator) + int(
            observed_den != den
        )
        for k, v in [
            ("n_completed_cells_detected", numerator),
            ("completed_cell_denominator", den),
            ("n_cells_detected", sum((c, f) in want for c in ids)),
            ("total_support", sum(want.get((c, f), 0) for c in ids)),
        ]:
            result["report_errors"] += int(int(a.var.at[f, k]) != v) + int(
                int(qc[f][k]) != v
            )
        rate = qc[f]["prevalence_completed"]
        observed = None if rate == "" else float(rate)
        result["report_errors"] += int(observed != (numerator / den if den else None))
        obj = float(a.var.at[f, "prevalence_completed"])
        result["report_errors"] += int(
            not np.isnan(obj) if not den else obj != numerator / den
        )
        if not den:
            result["undefined_candidates"] += 1
        error = (
            abs(Fraction(observed_num, observed_den) - Fraction(numerator, den))
            if den and observed_den
            else None
        )
        prevalence.append(
            dict(
                candidate=f,
                oracle_numerator=numerator,
                oracle_denominator=den,
                arm_numerator=observed_num,
                arm_denominator=observed_den,
                absolute_error_numerator=error.numerator if error is not None else "",
                absolute_error_denominator=error.denominator
                if error is not None
                else "",
                oracle_undefined=not bool(den),
                arm_undefined=not bool(observed_den),
            )
        )
    parser = _old.Tables()
    parser.feed((out / "qc/report.html").read_text())
    for name in ["cell_qc.tsv", "circ_qc.tsv"]:
        rows = read_rows(out / "qc" / name)
        with (out / "qc" / name).open() as h:
            header = next(csv.reader(h, delimiter="\t"))
        candidates = [t for t in parser.tables if t and t[0] == header]
        if len(candidates) != 1:
            result["contract_errors"] += 1
            continue
        actual = candidates[0][1:]
        result["report_errors"] += int(len(actual) != min(100, len(rows)))
        for src, displayed in zip(rows[:100], actual):
            for k, value in zip(header, displayed):
                expected = src[k]
                if k == "prevalence_completed":
                    # Exact declared six-place rounding, independent of zero padding.
                    if expected == "":
                        ok = value in ["", "NaN", "nan"]
                    else:
                        ok = Fraction(value) == Fraction(f"{float(expected):.6f}")
                else:
                    ok = expected == value
                result["report_errors"] += int(not ok)
    if config.get("rna"):
        rna = ad.read_h5ad(config["rna"])
        m = md.read_h5mu(out / "mudata/circyto_multimodal.h5mu")
        try:
            _old.equal_adata(rna, m.mod["rna"], "RNA every available field")
        except AssertionError:
            result["cell_assignment_errors"] += 1
            for c in set(rna.obs_names) & set(m.mod["rna"].obs_names):
                if (rna[[c]].X != m.mod["rna"][[c]].X).nnz:
                    result["cell_assignment_errors"] += 1
        union = ids + [c for c in rna.obs_names if c not in ids]
        result["contract_errors"] += int(list(m.obs_names) != union)
        for mod in ["circ", "rna"]:
            for i, c in enumerate(m.obs_names):
                index = int(np.asarray(m.obsmap[mod]).ravel()[i])
                expected = c if c in (ids if mod == "circ" else rna.obs_names) else None
                observed = m.mod[mod].obs_names[index - 1] if index else None
                result["cell_assignment_errors"] += int(observed != expected)
        for c in union:
            membership = (
                "shared"
                if c in ids and c in rna.obs_names
                else "circ_only"
                if c in ids
                else "rna_only"
            )
            result["status_errors"] += int(m.obs.at[c, "membership"] != membership)
    for name, p in sources.items():
        expected = sha(p)
        result["source_errors"] += int(report["source_sha256"].get(name) != expected)
        if name != "rna.h5ad":
            result["source_errors"] += int(sha(out / "evidence" / name) != expected)
    result["unavailable_fields"] = sorted(set(result["unavailable_fields"]))
    result["passed"] = (
        not any(
            result[k]
            for k in [
                "count_errors",
                "annotation_errors",
                "cell_assignment_errors",
                "status_errors",
                "report_errors",
                "source_errors",
                "contract_errors",
            ]
        )
        and not result["unavailable_fields"]
    )
    return result, prevalence
