"""Installed-wheel acceptance; consumes supplied authentic native files unchanged.

python validate_native_ciri3.py --fixtures native_ciri3 --out /tmp/native-check
Optional --rna takes an EXISTING RNA H5AD containing upstream_single; this is
an ID-alignment test, not a claim of biologically matched RNA for the demo.
GNU /usr/bin/time records each child CLI's wall time and peak RSS on Linux.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import circyto


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--rna", type=Path)
    args = p.parse_args()
    args.fixtures = args.fixtures.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    provenance = json.loads((args.fixtures / "provenance.json").read_text())
    for filename, entry in provenance["files"].items():
        assert sha(args.fixtures / filename) == entry["sha256"]
    with (args.fixtures / "result.txt").open() as f:
        expected = list(csv.DictReader(f, delimiter="\t"))
    # Canonical example files are used as-is, no intermediate fake detector table.
    command = [
        "circyto",
        "import-circ",
        "--ciri3-manifest",
        str(args.fixtures / "manifest.tsv"),
        "--cells",
        str(args.fixtures / "cells.tsv"),
        "--detector",
        "CIRI3",
        "--genome-build",
        "upstream_demo_chr1_unverified",
        "--coordinate-convention",
        "1-based-closed",
        "--support-definition",
        "CIRI3_junction_reads",
    ]
    if args.rna:
        command += ["--rna", str(args.rna.resolve())]
    for name in ("first", "repeat"):
        subprocess.run(
            [
                "/usr/bin/time",
                "-f",
                '{"elapsed_seconds":%e,"peak_rss_kib":%M}',
                "-o",
                str(args.out / f"{name}.time.json"),
                *command,
                "--outdir",
                str(args.out / name),
            ],
            check=True,
        )
    output = args.out / "first"
    a = ad.read_h5ad(output / "anndata/circ_counts.h5ad")
    assert sparse.isspmatrix_csr(a.X) and a.X.dtype == np.int64
    assert a.shape == (1, 6) and a.obs_names.tolist() == ["upstream_single"]
    assert a.var_names.tolist() == [row["circRNA_ID"] for row in expected]
    assert a.X.data.tolist() == [int(row["#junction_reads"]) for row in expected]
    for dst, src in [("start", "circRNA_start"), ("end", "circRNA_end")]:
        assert a.var[dst].tolist() == [int(row[src]) for row in expected]
    for dst, src in [
        ("chrom", "chr"),
        ("strand", "strand"),
        ("ciri3_gene_id", "gene_id"),
        ("ciri3_circRNA_type", "circRNA_type"),
    ]:
        assert a.var[dst].tolist() == [row[src] for row in expected]
    assert a.obs.circ_observation_status.tolist() == ["unknown"]
    assert a.var.host_gene.eq("").all()
    metrics = json.loads((output / "qc/metrics.json").read_text())
    assert (
        metrics["total_support"] == 50
        and metrics["detection_fraction_completed"]["value"] is None
    )
    assert (
        sha(output / "evidence/ciri3_native_000001.tsv")
        == provenance["files"]["result.txt"]["sha256"]
    )
    assert json.loads(a.uns["circyto"]["metrics_json"]) == metrics
    for name in [
        "qc/metrics.json",
        "qc/report.html",
        "qc/cell_qc.tsv",
        "qc/circ_qc.tsv",
        "matrix/circ_counts.mtx",
    ]:
        assert (output / name).read_bytes() == (args.out / "repeat" / name).read_bytes()
    html = (output / "qc/report.html").read_text()
    assert "unavailable" in html and "unknown" in html and "<script" not in html
    rna_check = None
    if args.rna:
        import mudata

        source = ad.read_h5ad(args.rna)
        m = mudata.read_h5mu(output / "mudata/circyto_multimodal.h5mu")
        pd.testing.assert_frame_equal(m["rna"].obs, source.obs)
        pd.testing.assert_frame_equal(m["rna"].var, source.var)
        assert (m["rna"].X != source.X).nnz == 0
        assert m["circ"].X.sum() == 50
        rna_check = {
            "scope": "technical ID alignment only; no biologically matched upstream RNA is available",
            "input_sha256": sha(args.rna),
            "n_shared": metrics["rna_overlap"]["n_shared"],
        }
    result = {
        "result": "NATIVE CIRI3 IMPORT VALIDATED",
        "scope": "Official 13-column single-sample demonstration, not a biological benchmark",
        "python": sys.version,
        "platform": platform.platform(),
        "circyto_import_path": circyto.__file__,
        "source_sha256": provenance["files"]["result.txt"]["sha256"],
        "source_bytes": 1913,
        "schema": "ciri3.native-result-13.v1",
        "shape": list(a.shape),
        "nnz": a.X.nnz,
        "support": int(a.X.sum()),
        "processing_status": "unknown",
        "native_producer_version": "unknown",
        "rna_check": rna_check,
        "measurements": {
            name: json.loads((args.out / f"{name}.time.json").read_text())
            for name in ["first", "repeat"]
        },
        "reproducibility": "identical text/Matrix Market outputs; no HDF5 byte identity asserted",
    }
    (args.out / "acceptance.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
