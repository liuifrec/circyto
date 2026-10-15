"""Outside-checkout acceptance using the authentic native example from the sdist.

--multimodal creates explicitly SYNTHETIC RNA for an identity/value test only.
No Java, aligner, network request, source-tree import or detector execution.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import circyto
from circyto.paths import get_repo_root


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--multimodal", action="store_true")
    args = parser.parse_args()
    fixtures = args.fixtures.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert "site-packages" in circyto.__file__, circyto.__file__
    assert get_repo_root() is None, "Run outside the source checkout"
    assert circyto.__version__ == metadata.version("circyto")
    package = Path(circyto.__file__).parent
    assert not any(p.suffix.lower() in {".jar", ".class"}
                   or p.name == "CIRI-full_manual.txt" for p in package.rglob("*"))
    # Child processes use this environment's console script and no PYTHONPATH.
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["PYTHONNOUSERSITE"] = "1"
    executable = str(Path(sys.executable).parent / "circyto")

    def run(arguments, name, success=True):
        result = subprocess.run([executable, *arguments], env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (out / f"{name}.log").write_text(result.stdout)
        assert (result.returncode == 0) == success, result.stdout
        return result.stdout

    for name, arguments in [("help", ["--help"]), ("doctor", ["doctor"]),
                            ("import-help", ["import-circ", "--help"]),
                            ("report-help", ["report", "--help"])]:
        run(arguments, name)
    provenance = json.loads((fixtures / "provenance.json").read_text())
    for filename, record in provenance["files"].items():
        assert hashlib.sha256((fixtures / filename).read_bytes()).hexdigest() == record["sha256"]
    with (fixtures / "result.txt").open() as stream:
        expected = list(csv.DictReader(stream, delimiter="\t"))
    command = ["import-circ", "--ciri3-manifest", str(fixtures / "manifest.tsv"),
               "--cells", str(fixtures / "cells.tsv"), "--detector", "CIRI3",
               "--genome-build", "upstream_demo_chr1_unverified",
               "--coordinate-convention", "1-based-closed",
               "--support-definition", "CIRI3_junction_reads"]
    run([*command, "--outdir", str(out / "native")], "native")
    result = ad.read_h5ad(out / "native/anndata/circ_counts.h5ad")
    assert result.shape == (1, 6)
    assert result.obs_names.tolist() == ["upstream_single"]
    assert result.var_names.tolist() == [row["circRNA_ID"] for row in expected]
    assert sparse.isspmatrix_csr(result.X) and result.X.dtype == np.int64
    assert result.X.data.tolist() == [int(row["#junction_reads"]) for row in expected]
    assert int(result.X.sum()) == 50
    for target, source in [("start", "circRNA_start"), ("end", "circRNA_end")]:
        assert result.var[target].tolist() == [int(row[source]) for row in expected]
    for target, source in [("chrom", "chr"), ("strand", "strand"),
                           ("ciri3_gene_id", "gene_id"),
                           ("ciri3_circRNA_type", "circRNA_type")]:
        assert result.var[target].tolist() == [row[source] for row in expected]
    assert result.obs.circ_observation_status.tolist() == ["unknown"]
    assert result.var.host_gene.eq("").all()
    metrics = json.loads((out / "native/qc/metrics.json").read_text())
    assert metrics["detection_fraction_completed"]["value"] is None
    html = (out / "native/qc/report.html").read_text()
    assert "unknown" in html and "<script" not in html and "<html" in html
    # An import report must not be mistaken for a recovered workflow report.
    failure = run(["report", "--workdir", str(out / "native")], "workflow-report-refusal", False)
    assert "import" in failure.lower()
    rna = ad.AnnData(
        sparse.csr_matrix(np.array([[3], [4], [9]], dtype=np.int64)),
        obs=pd.DataFrame({"fixture": ["synthetic"] * 3},
                         index=["upstream_single", "synthetic_only", "rna_only"]),
        var=pd.DataFrame(index=["SYNTHETIC_GENE"]),
    )
    rna.uns["genome_build"] = "upstream_demo_chr1_unverified"
    rna.uns["scope"] = "Synthetic identity fixture; no biological match to upstream CIRI3 demo"
    rna.write_h5ad(out / "synthetic_rna.h5ad")
    multimodal = None
    if args.multimodal:
        import mudata

        run([*command, "--rna", str(out / "synthetic_rna.h5ad"),
             "--outdir", str(out / "matched")], "matched")
        m = mudata.read_h5mu(out / "matched/mudata/circyto_multimodal.h5mu")
        assert set(m.mod) == {"rna", "circ"}
        pd.testing.assert_frame_equal(m["rna"].obs, rna.obs)
        pd.testing.assert_frame_equal(m["rna"].var, rna.var)
        assert (m["rna"].X != rna.X).nnz == 0
        assert m["rna"].uns["scope"] == rna.uns["scope"]
        pd.testing.assert_frame_equal(m["circ"].var, result.var)
        assert m["circ"].obs_names.tolist() == result.obs_names.tolist()
        assert (m["circ"].X != result.X).nnz == 0
        matched = json.loads((out / "matched/qc/metrics.json").read_text())
        assert matched["rna_overlap"]["n_shared"] == 1
        assert (out / "matched/qc/report.html").is_file()
        multimodal = {"scope": rna.uns["scope"], "shared_ids": 1,
                      "rna_shape": list(m["rna"].shape), "values_preserved": True}
    elif importlib.util.find_spec("mudata") is None:
        failure = run([*command, "--rna", str(out / "synthetic_rna.h5ad"),
                       "--outdir", str(out / "missing-mudata")], "missing-mudata", False)
        assert "mudata" in failure.lower()
    report = {
        "result": "PASS", "python": platform.python_version(),
        "platform": platform.platform(), "module_origin": circyto.__file__,
        "versions": {name: metadata.version(name) for name in
                     ["circyto", "anndata", "numpy", "scipy", "pandas", "typer", "rich"]
                     + (["mudata"] if args.multimodal else [])},
        "external_on_path": {name: shutil.which(name) for name in ["java", "bwa", "STAR"]},
        "native": {"shape": [1, 6], "support": 50, "processing_status": "unknown",
                   "scope": "Authentic format example, not an established biological single cell"},
        "multimodal": multimodal, "offline_html": str(out / "native/qc/report.html"),
        "original_workflow_report": "not tested; import directory correctly refused",
    }
    (out / "acceptance.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
