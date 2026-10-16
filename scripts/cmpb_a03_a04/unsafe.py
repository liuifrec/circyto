"""Explicit unsafe illustration, sharing reference I/O validation/output only.

After strict parsing it overwrites observation states and assigns RNA IDs by
position. This is not a named tool or a typical-scientist proxy. Its invalid-input
rejections are inherited parser behaviour, not an independent comparison.
"""

import argparse
import json
from pathlib import Path
from reference import load, write_bundle, digest


def run(config, out):
    X, obs, var, rna, sources = load(config)
    hashes = {n: digest(p) for n, p in sources.items()}
    obs["detector_status"] = "completed"
    if rna is not None:
        rna = rna.copy()
        rna.obs_names = list(obs.index[: rna.n_obs]) + [
            f"unsafe_extra_{i}" for i in range(max(0, rna.n_obs - len(obs)))
        ]
    write_bundle(X, obs, var, rna, sources, config, out, hashes)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args()
    run(json.loads(a.config.read_text()), a.outdir)
