"""Read-only original-object inventory; no scientific matrix transformations."""
import hashlib
import importlib.metadata as md
import json
from pathlib import Path
import sys

import anndata as ad
import h5py
import numpy as np

SOURCES = {
    "Smart-seq3": ("emtab8735_smartseq3/full_length.hostgene_fixed.h5mu", "0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339"),
    "IMR90": ("scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu", "bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8"),
}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def plain(v):
    if isinstance(v, dict):
        return {str(k): plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [plain(x) for x in v]
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, bytes):
        return v.decode()
    return v


def main():
    root = Path(sys.argv[1])
    original = Path(sys.argv[2])
    result = {"python": sys.version, "executable": sys.executable,
              "packages": {p: md.version(p) for p in ["numpy", "scipy", "pandas", "anndata", "mudata", "h5py"]},
              "read_mode": "h5py.File(..., 'r'); AnnData read_elem for obs/var/uns; matrices only described here",
              "sources": {}}
    for name, (member, expected) in SOURCES.items():
        path = original / member
        sha = digest(path)
        assert sha == expected, (name, sha)
        with h5py.File(path, "r") as f:
            obs = ad.io.read_elem(f["obs"])
            assert obs.index.is_unique
            record = {"archive_member": "load_work/" + member, "sha256_before": sha,
                      "bytes": path.stat().st_size, "global_obs": obs.reset_index().to_dict("records"),
                      "uns": plain(ad.io.read_elem(f["uns"])), "modalities": {}}
            for mod in f["mod"]:
                g = f["mod"][mod]
                o, v = ad.io.read_elem(g["obs"]), ad.io.read_elem(g["var"])
                assert o.index.is_unique and v.index.is_unique
                x = g["X"]
                shape = list(x.shape if isinstance(x, h5py.Dataset) else x.attrs["shape"])
                assert shape == [len(o), len(v)]
                fields = []
                def visit(key, item):
                    fields.append({"path": item.name, "kind": "dataset" if isinstance(item, h5py.Dataset) else "group",
                                   "shape": list(item.shape) if isinstance(item, h5py.Dataset) else None,
                                   "dtype": str(item.dtype) if isinstance(item, h5py.Dataset) else None})
                g.visititems(visit)
                record["modalities"][mod] = {"shape": shape, "encoding": plain(x.attrs.get("encoding-type")),
                    "obs_columns": o.columns.tolist(), "var_columns": v.columns.tolist(),
                    "obs_ids": o.index.tolist(), "var_head": v.head(3).to_dict("index"),
                    "var_tail": v.tail(3).to_dict("index"),
                    "membership_counts": o["membership"].value_counts().to_dict() if "membership" in o else None,
                    "uns": plain(ad.io.read_elem(g["uns"])), "hdf5_fields": fields}
            result["sources"][name] = record
    (root / "results/source_object_inventory.json").write_text(json.dumps(plain(result), indent=2, allow_nan=False) + "\n")
    print(json.dumps(plain({k: {"sha256": v["sha256_before"], "modalities": {m: r["shape"] for m, r in v["modalities"].items()}} for k, v in result["sources"].items()})))


if __name__ == "__main__":
    main()
