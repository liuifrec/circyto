#!/usr/bin/env python3
"""Independent ID-aware study conversion; no circyto or acceptance-oracle imports."""

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
import re
import shutil

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

MAX = 2**63 - 1
STATES = {"completed", "failed", "unknown", "unprocessed", "skipped"}
DECL = ["genome_build", "coordinate_convention", "support_definition", "detector"]
NATIVE = [
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


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def read_table(path, required):
    with Path(path).open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        names = r.fieldnames
        if (
            not names
            or len(names) != len(set(names))
            or not set(required) <= set(names)
        ):
            raise ValueError("Invalid table header")
        rows = list(r)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError("Malformed row")
    return names, rows


def identity(x):
    if (
        not isinstance(x, str)
        or not x
        or x.strip() != x
        or any(c in x for c in "\r\n\t")
    ):
        raise ValueError("Invalid literal identifier")
    return x


def unique(values):
    values = [identity(x) for x in values]
    if len(values) != len(set(values)):
        raise ValueError("Duplicate identity")
    return values


def integer(x):
    if not re.fullmatch("[0-9]+", str(x)):
        raise ValueError("Nonnegative integer required")
    v = int(x)
    if v > MAX:
        raise ValueError("Integer overflow")
    return v


def declarations(row, declared):
    for name, value in declared.items():
        if name in row and row[name] != value:
            raise ValueError("Conflicting declaration: " + name)


def feature(row, declared):
    declarations(row, declared)
    d = {k: v for k, v in row.items() if k not in ["cell_id", "support"]}
    for k in ["circ_id", "chrom"]:
        identity(d[k])
    for k in ["start", "end"]:
        d[k] = integer(d[k])
    low = 0 if declared["coordinate_convention"] == "0-based-half-open" else 1
    if (
        d["start"] < low
        or d["end"] < d["start"]
        or (low == 0 and d["end"] == d["start"])
        or d["strand"] not in ["+", "-", "."]
    ):
        raise ValueError("Invalid declared interval/strand")
    d.update(declared)
    return d


def load(config):
    declared = {k: config[k] for k in DECL}
    if declared["coordinate_convention"] not in ["1-based-closed", "0-based-half-open"]:
        raise ValueError("Unsupported coordinate convention")
    if any(not v or v.strip() != v for v in declared.values()):
        raise ValueError("Missing declaration")
    _, meta = read_table(config["cells"], ["cell_id"])
    ids = unique([r["cell_id"] for r in meta])
    if not ids:
        raise ValueError("Empty metadata")
    for r in meta:
        declarations(r, declared)
        r.setdefault("detector_status", "unknown")
        if r["detector_status"] not in STATES:
            raise ValueError("Unknown processing-state token")
    sources = {"cells.tsv": Path(config["cells"])}
    records = []
    features = {}
    available = set(ids)
    mode = config["mode"]
    if mode == "native":
        if (
            declared["detector"],
            declared["coordinate_convention"],
            declared["support_definition"],
        ) != ("CIRI3", "1-based-closed", "CIRI3_junction_reads"):
            raise ValueError("Native declarations incompatible")
        mapping_path = Path(config["input"])
        names, mapping = read_table(
            mapping_path,
            ["cell_id", "sample_id", "path", "detector_version", "result_scope"],
        )
        if (
            names
            != ["cell_id", "sample_id", "path", "detector_version", "result_scope"]
            or not mapping
        ):
            raise ValueError("Native mapping schema")
        available = set(unique([r["cell_id"] for r in mapping]))
        unique([r["sample_id"] for r in mapping])
        if len({r["detector_version"] for r in mapping}) != 1:
            raise ValueError("Mixed native versions")
        seen_files = set()
        sources["ciri3_manifest.tsv"] = mapping_path
        by_cell = {r["cell_id"]: r for r in meta}
        for i, r in enumerate(mapping, 1):
            if r["cell_id"] not in ids:
                raise ValueError("Foreign native cell")
            if (
                r["detector_version"] not in ["3.0.1", "unknown"]
                or r["result_scope"] != "single-sample"
            ):
                raise ValueError("Unsupported native producer/scope")
            p = Path(r["path"])
            p = (
                (mapping_path.parent / p).resolve()
                if not p.is_absolute()
                else p.resolve()
            )
            key = (p.stat().st_dev, p.stat().st_ino)
            if key in seen_files:
                raise ValueError("Duplicated native source")
            seen_files.add(key)
            name = f"ciri3_native_{i:06d}.tsv"
            sources[name] = p
            names, raw = read_table(p, NATIVE)
            if names != NATIVE:
                raise ValueError("Unreviewed native schema")
            unique([x["circRNA_ID"] for x in raw])
            by_cell[r["cell_id"]].update(
                ciri3_sample_id=r["sample_id"],
                ciri3_detector_version=r["detector_version"],
                ciri3_source_file=name,
            )
            for x in raw:
                start, end = integer(x["circRNA_start"]), integer(x["circRNA_end"])
                if x["circRNA_ID"] != f"{x['chr']}:{start}|{end}":
                    raise ValueError("Native ID disagrees with interval")
                records.append(
                    dict(
                        cell_id=r["cell_id"],
                        circ_id=x["circRNA_ID"],
                        chrom=x["chr"],
                        start=start,
                        end=end,
                        strand=x["strand"],
                        support=x["#junction_reads"],
                        ciri3_gene_id=x["gene_id"],
                        ciri3_circRNA_type=x["circRNA_type"],
                    )
                )
    elif mode == "calls":
        sources["calls.tsv"] = Path(config["input"])
        _, records = read_table(
            config["input"],
            ["cell_id", "circ_id", "chrom", "start", "end", "strand", "support"],
        )
    elif mode == "matrix":
        folder = Path(config["input"])
        for n in [
            "circ_counts.mtx",
            "cell_index.txt",
            "circ_index.txt",
            "circ_feature_table.tsv",
        ]:
            sources[n] = folder / n
        matrix_cells = unique((folder / "cell_index.txt").read_text().splitlines())
        available = set(matrix_cells)
        fids = unique((folder / "circ_index.txt").read_text().splitlines())
        _, raw = read_table(
            folder / "circ_feature_table.tsv",
            ["circ_id", "chrom", "start", "end", "strand"],
        )
        unique([r["circ_id"] for r in raw])
        by_id = {r["circ_id"]: r for r in raw}
        if set(by_id) != set(fids) or not available <= set(ids):
            raise ValueError("Matrix axis mismatch")
        features = {key: feature(by_id[key], declared) for key in fids}
        a = mmread(folder / "circ_counts.mtx")
        if (
            not sparse.issparse(a)
            or a.dtype.kind not in "iu"
            or a.shape != (len(fids), len(matrix_cells))
        ):
            raise ValueError("Sparse integer candidate-by-cell matrix required")
        a = a.tocoo()
        pairs = set()
        for j, i, v in zip(a.row, a.col, a.data):
            pair = (int(i), int(j))
            if pair in pairs:
                raise ValueError("Duplicate sparse entry")
            pairs.add(pair)
            records.append(
                dict(cell_id=matrix_cells[i], support=integer(v), **features[fids[j]])
            )
    else:
        raise ValueError("Unsupported mode")
    if any(
        r["detector_status"] == "completed" and r["cell_id"] not in available
        for r in meta
    ):
        raise ValueError("Completed cell lacks explicit source")
    entries = {}
    total = 0
    for r in records:
        cell = identity(r["cell_id"])
        key = identity(r["circ_id"])
        if cell not in ids:
            raise ValueError("Foreign call cell")
        if (cell, key) in entries:
            raise ValueError("Duplicate cell/candidate")
        f = feature(r, declared)
        if key in features and features[key] != f:
            raise ValueError("Conflicting candidate annotation/coordinates")
        features[key] = f
        value = integer(r["support"])
        entries[cell, key] = value
        total += value
        if total > MAX:
            raise ValueError("Aggregate count overflow")
    fids = list(features)
    rows = {c: i for i, c in enumerate(ids)}
    cols = {c: j for j, c in enumerate(fids)}
    positive = [(rows[c], cols[f], v) for (c, f), v in entries.items() if v]
    X = sparse.coo_matrix(
        (
            np.array([r[2] for r in positive], dtype=np.int64),
            ([r[0] for r in positive], [r[1] for r in positive]),
        ),
        shape=(len(ids), len(fids)),
    ).tocsr()
    obs = pd.DataFrame(meta).set_index("cell_id")
    var = pd.DataFrame([features[f] for f in fids], index=pd.Index(fids))
    if not fids:
        var = pd.DataFrame(
            columns=["circ_id", "chrom", "start", "end", "strand", *DECL],
            index=pd.Index([], dtype=str),
        )
    var.index.name = None
    rna = None
    if config.get("rna"):
        sources["rna.h5ad"] = Path(config["rna"])
        rna = ad.read_h5ad(config["rna"])
        unique(list(rna.obs_names))
        unique(list(rna.var_names))
        if not set(ids) & set(rna.obs_names):
            raise ValueError("No shared RNA identities")
        for frame in [rna.obs, rna.var]:
            if "genome_build" in frame:
                vals = frame.genome_build.dropna().astype(str)
                if any(v and v != declared["genome_build"] for v in vals):
                    raise ValueError("RNA build conflict")
        for u in [rna.uns, rna.uns.get("circyto", {})]:
            if u.get("genome_build") and u["genome_build"] != declared["genome_build"]:
                raise ValueError("RNA build conflict")
    return X, obs, var, rna, sources


def write_bundle(X, obs, var, rna, sources, config, out, before):
    if out.exists():
        raise FileExistsError("Output must be new")
    for n in ["anndata", "mudata", "matrix", "qc", "evidence"]:
        (out / n).mkdir(parents=True, exist_ok=True)
    complete = [c for c in obs.index if obs.at[c, "detector_status"] == "completed"]
    den = len(complete)
    positive = {c: 0 for c in obs.index}
    support = {c: 0 for c in obs.index}
    detected = {f: 0 for f in var.index}
    totals = {f: 0 for f in var.index}
    numer = {f: 0 for f in var.index}
    coo = X.tocoo()
    for i, j, v in zip(coo.row, coo.col, coo.data):
        if not v:
            continue
        c, f = obs.index[i], var.index[j]
        positive[c] += 1
        support[c] += int(v)
        detected[f] += 1
        totals[f] += int(v)
        numer[f] += int(c in complete)
    obs = obs.copy()
    var = var.copy()
    obs["circRNA_count"] = [positive[c] for c in obs.index]
    obs["circRNA_total_support"] = [support[c] for c in obs.index]
    obs["circ_observation_status"] = [
        ("completed_detected" if positive[c] else "completed_zero")
        if c in complete
        else obs.at[c, "detector_status"]
        for c in obs.index
    ]
    for name, d in [
        ("n_cells_detected", detected),
        ("total_support", totals),
        ("n_completed_cells_detected", numer),
    ]:
        var[name] = [d[f] for f in var.index]
    var["completed_cell_denominator"] = den
    var["prevalence_completed"] = [numer[f] / den if den else np.nan for f in var.index]
    metrics = dict(
        n_cells=len(obs),
        n_candidates=len(var),
        nnz=X.nnz,
        total_support=sum(support.values()),
        source_sha256=before,
        observation_status_counts={
            str(k): int(v)
            for k, v in obs.circ_observation_status.value_counts().items()
        },
        detection_fraction_completed=dict(
            numerator=sum(positive[c] > 0 for c in complete),
            denominator=den,
            value=sum(positive[c] > 0 for c in complete) / den if den else None,
        ),
        zero_fraction_completed=dict(
            numerator=sum(positive[c] == 0 for c in complete),
            denominator=den,
            value=sum(positive[c] == 0 for c in complete) / den if den else None,
        ),
        declarations={k: config[k] for k in DECL},
    )
    a = ad.AnnData(
        X,
        obs=obs,
        var=var,
        uns={"reference_contract": "cmpb.id_aware.v1", "source_sha256": before},
    )
    a.write_h5ad(out / "anndata/circ_counts.h5ad")
    mmwrite(out / "matrix/circ_counts.mtx", X.T)
    (out / "matrix/cell_index.txt").write_text("\n".join(obs.index) + "\n")
    (out / "matrix/circ_index.txt").write_text("\n".join(var.index) + "\n")
    obs.to_csv(out / "qc/cell_qc.tsv", sep="\t")
    var.to_csv(out / "qc/circ_qc.tsv", sep="\t", index=False)
    var.to_csv(out / "matrix/circ_feature_table.tsv", sep="\t", index=False)
    if rna is not None:
        with md.set_options(pull_on_update=False):
            m = md.MuData({"circ": a, "rna": rna})
        union = list(obs.index) + [c for c in rna.obs_names if c not in obs.index]
        if list(m.obs_names) != union:
            raise ValueError("Standard MuData ordering differs from declared contract")
        global_obs = obs.reindex(union).copy()
        for k in ["detector_status", "circ_observation_status"]:
            global_obs[k] = global_obs[k].fillna("unknown")
        global_obs["circ_present"] = [c in obs.index for c in union]
        global_obs["rna_present"] = [c in rna.obs_names for c in union]
        global_obs["membership"] = [
            "shared"
            if c in obs.index and c in rna.obs_names
            else "circ_only"
            if c in obs.index
            else "rna_only"
            for c in union
        ]
        m.obs = global_obs
        m.write_h5mu(out / "mudata/circyto_multimodal.h5mu")
        metrics["rna_overlap"] = dict(
            n_shared=len(set(obs.index) & set(rna.obs_names)),
            n_circ_only=len(set(obs.index) - set(rna.obs_names)),
            n_rna_only=len(set(rna.obs_names) - set(obs.index)),
        )
    for name, path in sources.items():
        if name != "rna.h5ad":
            shutil.copyfile(path, out / "evidence" / name)
            if digest(out / "evidence" / name) != before[name]:
                raise ValueError("Copied source bytes changed")
    (out / "qc/metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    blocks = [
        '<!doctype html><meta charset="utf-8"><title>Independent reference import</title><h1>Independent reference import</h1>',
        f"<p>{len(obs)} cells; {len(var)} candidates; {sum(support.values())} support. Completed denominator {den}.</p>",
        "<p>Technical import; not detector completion evidence or biological absence. First 100 rows; rates rounded to six decimal places.</p>",
    ]
    for title, frame in [("Cells", obs.reset_index()), ("Candidates", var)]:
        blocks += [
            "<h2>"
            + title
            + "</h2><table><tr>"
            + "".join("<th>" + html.escape(str(k)) + "</th>" for k in frame.columns)
            + "</tr>"
        ]
        for _, row in frame.head(100).iterrows():
            values = [
                (
                    ""
                    if pd.isna(v)
                    else f"{v:.6f}"
                    if k == "prevalence_completed"
                    else str(v)
                )
                for k, v in row.items()
            ]
            blocks += [
                "<tr>"
                + "".join("<td>" + html.escape(v) + "</td>" for v in values)
                + "</tr>"
            ]
        blocks += ["</table>"]
    (out / "qc/report.html").write_text("\n".join(blocks) + "\n")
    for n, p in sources.items():
        if digest(p) != before[n]:
            raise ValueError("Input mutated during conversion")


def convert(config, out):
    # Hash sources before parsing, including files explicitly named in the manifest.
    paths = {"cells.tsv": Path(config["cells"])}
    if config["mode"] == "native":
        p = Path(config["input"])
        paths["ciri3_manifest.tsv"] = p
        _, rows = read_table(p, ["path"])
        for i, r in enumerate(rows, 1):
            paths[f"ciri3_native_{i:06d}.tsv"] = (p.parent / Path(r["path"])).resolve()
    elif config["mode"] == "calls":
        paths["calls.tsv"] = Path(config["input"])
    else:
        for n in [
            "circ_counts.mtx",
            "cell_index.txt",
            "circ_index.txt",
            "circ_feature_table.tsv",
        ]:
            paths[n] = Path(config["input"]) / n
    if config.get("rna"):
        paths["rna.h5ad"] = Path(config["rna"])
    hashes = {n: digest(p) for n, p in paths.items()}
    X, obs, var, rna, sources = load(config)
    if paths != sources:
        raise ValueError("Parsed source plan differs")
    write_bundle(X, obs, var, rna, sources, config, out, hashes)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args()
    try:
        convert(json.loads(a.config.read_text()), a.outdir)
    except Exception as e:
        print(json.dumps({"error": type(e).__name__, "message": str(e)}))
        raise
