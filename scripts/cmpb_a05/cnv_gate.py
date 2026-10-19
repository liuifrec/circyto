"""Inspect CNV encodings and public source equality. No mask inference or display aggregation."""
import collections
import json
from pathlib import Path
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from analyze import table
from inspect_sources import SOURCES, digest


def main():
    root, original, out = map(Path, sys.argv[1:])
    member, sha = SOURCES["IMR90"]
    path = original / member
    assert digest(path) == sha
    inventory, encodings, cell_audit, coordinates = [], [], [], []
    conflicts = []
    def item(field, value, source="original /mod/cnv", evidence="direct read-only inspection"):
        inventory.append(dict(field=field, value=json.dumps(value, ensure_ascii=False), source=source,
                              evidence=evidence, source_object_sha256=sha))
    item("callability_mask", "NOT VERIFIED; no mask manufactured")
    acquisitions = json.loads((root / "records/public_acquisition.json").read_text())
    for r in acquisitions:
        item("public_input:" + r["name"], r, r["url"], "new bounded metadata/source-table retrieval; no raw reads")
    with h5py.File(path, "r") as f:
        g = f["mod/cnv"]
        obs, var = ad.io.read_elem(g["obs"]), ad.io.read_elem(g["var"])
        item("shape", list(g["X"].shape))
        item("X_dtype", str(g["X"].dtype))
        item("layers", list(g["layers"].keys()))
        item("uns", ad.io.read_elem(g["uns"]))
        fields = []
        g.visit(lambda p: fields.append(p))
        item("all_hdf5_fields", fields)
        item("mask_or_missing_named_fields", [p for p in fields if any(s in p.lower() for s in ["mask", "missing", "callab", "valid", "blacklist"])])
        item("native_bin_size_counts", {str(k): int(v) for k,v in var.bin_size.value_counts().items()})
        item("bin_axis_unique", bool(var.index.is_unique))
        assert len(obs) == 23 and len(var) == 60607 and var.index.is_unique
        for chrom, df in var.groupby("seqname", sort=False):
            starts, ends = df.start.to_numpy(), df.end.to_numpy()
            coordinates.append(dict(chromosome=chrom, bins=len(df), first_start=int(starts[0]), last_end=int(ends[-1]),
                minimum_bin_bases=int((ends-starts).min()), maximum_bin_bases=int((ends-starts).max()),
                out_of_order=int((np.diff(starts)<0).sum()), overlapping_adjacent=int((starts[1:] < ends[:-1]).sum()),
                gaps_between_adjacent=int((starts[1:] > ends[:-1]).sum()),
                coordinate_convention="stored start/end with bin_size=end-start; public equality checked", source_object_sha256=sha))
        for i, cid in enumerate(obs.index):
            # Read one original dense row at a time; never convert genomic matrix to dense.
            state = g["X"][i, :]
            norm = g["layers/mappabilitynorm"][i, :]
            cell_audit.append(dict(cell_id=cid, source_cell_order=i, bins=len(state),
                state_nan=int(np.isnan(state).sum()), state_inf=int(np.isinf(state).sum()), state_zero=int((state==0).sum()),
                normalized_nan=int(np.isnan(norm).sum()), normalized_inf=int(np.isinf(norm).sum()), normalized_zero=int((norm==0).sum()),
                state_zero_normalized_nonzero=int(((state==0)&(norm!=0)).sum()),
                state_nonzero_normalized_zero=int(((state!=0)&(norm==0)).sum()),
                callability="unknown; numerical zeros are not interpreted", source_object_sha256=sha))
        for key, field in [("X", "state"), ("layers/mappabilitynorm", "mappabilitynorm")]:
            values = collections.Counter()
            minimum, maximum, nonfinite, zeros = float("inf"), float("-inf"), 0, 0
            for i in range(len(obs)):
                a = g[key][i, :]
                finite = a[np.isfinite(a)]
                minimum, maximum = min(minimum, float(finite.min())), max(maximum, float(finite.max()))
                nonfinite += int((~np.isfinite(a)).sum()); zeros += int((a==0).sum())
                if field == "state": values.update(map(int, a))
            item(field + "_numeric_inventory", dict(dtype=str(g[key].dtype), finite_min=minimum, finite_max=maximum,
                nonfinite=nonfinite, zeros=zeros, total_entries=len(obs)*len(var)))
            for value, count in sorted(values.items()):
                encodings.append(dict(source="original /mod/cnv/X", encoding=str(value), count=count, interpretation="stored integer only; no callability assigned", source_object_sha256=sha))
        for kind, key in [("states", "X"), ("mappabilitynorm", "layers/mappabilitynorm")]:
            public = root / "inputs" / ("GSE278958_06_scRR-seq-DNA_IMR90_G1_summary_CNV_"+kind+"_bin50kb.txt.gz")
            if not public.exists():
                item("public_comparison:"+kind, "unavailable"); continue
            offset, mismatches, coord_mismatches, missing_tokens = 0, 0, 0, 0
            symbols = collections.Counter()
            for chunk in pd.read_csv(public, sep="\t", dtype=str, keep_default_na=False, chunksize=2048):
                if offset == 0:
                    item("public_columns:"+kind, chunk.columns.tolist(), public.name)
                    expected = set(obs.dna_cell_id)
                    if set(chunk.columns[3:]) != expected:
                        conflicts.append(kind+": public cell IDs differ"); break
                    order = obs.dna_cell_id.tolist()
                end = offset + len(chunk)
                ids = chunk.seqname+":"+chunk.start+"-"+chunk.end
                coord_mismatches += int((ids.to_numpy() != var.index[offset:end].to_numpy()).sum())
                data = chunk[order].to_numpy()
                missing_tokens += int(np.isin(data, ["", "NA", "NaN", "nan", "N/A", "."]).sum())
                if kind == "states":
                    symbols.update(data.ravel())
                    valid = np.array([s.endswith("-somy") and s[:-5].isdigit() for s in data.ravel()]).all()
                    if not valid:
                        conflicts.append("public states contain unparsed symbols"); break
                    converted = np.array([int(s[:-5]) for s in data.ravel()], dtype=np.int16).reshape(data.shape).T
                else:
                    converted = data.astype(np.float32).T
                expected_array = g[key][:, offset:end]
                mismatches += int((converted != expected_array).sum())
                offset = end
            result = dict(bins_read=offset, compared_entries=offset*len(obs), coordinate_mismatches=coord_mismatches,
                value_mismatches=mismatches, missing_tokens=missing_tokens,
                comparison="exact integer; mappabilitynorm exact after documented float32 import conversion",
                public_file_sha256=digest(public), public_file_bytes=public.stat().st_size)
            item("public_comparison:"+kind, result, public.name)
            if offset != len(var) or mismatches or coord_mismatches:
                conflicts.append(kind+": public/archive axis or value discrepancy")
            for value, count in sorted(symbols.items()):
                encodings.append(dict(source=public.name, encoding=value, count=count,
                    interpretation="public categorical token; numeric parser does not supply mask semantics", source_object_sha256=sha))
    decision = "FAIL" if conflicts else "NOT VERIFIED"
    item("decision", decision)
    item("source_conflicts", conflicts)
    table(out / "cnv_source_inventory.tsv", inventory)
    table(out / "cnv_numeric_encoding_audit.tsv", encodings)
    table(out / "cnv_cell_missingness_audit.tsv", cell_audit)
    table(out / "cnv_coordinate_audit.tsv", coordinates)
    (out / "cnv_gate.json").write_text(json.dumps(dict(decision=decision, conflicts=conflicts,
        reason="Exact callability/missing-state definition for these exported bins not recovered; no mask in H5MU; no documented export mapping from filtered/model bins to 0-somy/zero.",
        source_object_sha256=sha, no_aggregation=True, no_heatmap=True), indent=2)+"\n")
    assert not (out / "cnv_display_bins.tsv").exists()
    assert digest(path) == sha
    print(json.dumps(dict(decision=decision, conflicts=conflicts)))


if __name__ == "__main__":
    main()
