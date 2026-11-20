"""Synthetic interoperability example; no sequencing analysis or detector runs.

Run with an installed circyto[mudata]: python import_circ_demo.py /tmp/import-demo
"""

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

root = Path(sys.argv[1])
root.mkdir(parents=True, exist_ok=False)
(root / "calls.tsv").write_text(
    "cell_id\tcirc_id\tchrom\tstart\tend\tstrand\tsupport\thost_gene\tevidence_note\n"
    "cell_A\tcirc_1\tchr1\t100\t200\t+\t7\tGENE1\tsynthetic example\n"
    "cell_C\tcirc_1\tchr1\t100\t200\t+\t2\tGENE1\tpartial failed output\n"
)
(root / "cells.tsv").write_text(
    "cell_id\tdetector_status\tbatch\n"
    "cell_B\tcompleted\tdemo\n"
    "cell_A\tcompleted\tdemo\n"
    "cell_C\tfailed\tdemo\n"
    "cell_D\tskipped\tdemo\n"
    "cell_E\tunprocessed\tdemo\n"
    "cell_F\tunknown\tdemo\n"
)
rna = ad.AnnData(
    sparse.csr_matrix(np.array([[3], [4], [9]], dtype=np.int64)),
    obs=pd.DataFrame(index=["cell_A", "cell_B", "rna_only"]),
    var=pd.DataFrame(index=["GENE1"]),
)
rna.uns["genome_build"] = "GRCh38"
rna.write_h5ad(root / "rna.h5ad")
subprocess.run(
    [
        "circyto",
        "import-circ",
        "--calls",
        str(root / "calls.tsv"),
        "--cells",
        str(root / "cells.tsv"),
        "--rna",
        str(root / "rna.h5ad"),
        "--outdir",
        str(root / "result"),
        "--genome-build",
        "GRCh38",
        "--coordinate-convention",
        "0-based-half-open",
        "--support-definition",
        "synthetic_BSJ_reads",
        "--detector",
        "synthetic_example",
    ],
    check=True,
)
result = ad.read_h5ad(root / "result/anndata/circ_counts.h5ad")
assert result.obs_names.tolist() == [
    "cell_B",
    "cell_A",
    "cell_C",
    "cell_D",
    "cell_E",
    "cell_F",
]
assert result.X.toarray().ravel().tolist() == [0, 7, 2, 0, 0, 0]
print("Synthetic example passed. Open", root / "result/qc/report.html")
