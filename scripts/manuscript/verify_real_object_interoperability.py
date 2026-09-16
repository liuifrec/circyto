#!/usr/bin/env python
"""Audit the two archived objects; no analysis, migration, or production patch.

Ordinary MuData reads establish interoperability. read_elem supplies an
independent reference for the stored semantics, never a replacement reader.
Large scratch outputs must stay outside the checkout. Only audit.json is small
and portable enough to retain with the manuscript.
"""
from __future__ import annotations

import argparse
from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import traceback
import warnings

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

import h5py
import mudata as mu
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse

from circyto.multimodal.sync import read_h5mu as circyto_read_h5mu

ARCHIVE = "c99cddaae2314d3af7d5232150c5b86bb6b52e25"
BASELINE = "44697355bcab1c525ca7ef9b130e2ad0094d9e1b"
INPUTS = {
    "smartseq3": (
        "load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu",
        "0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339",
    ),
    "imr90": (
        "load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu",
        "bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8",
    ),
}
SLOTS = ("obsm", "varm", "obsp", "varp", "uns")
CANDIDATE = "chr1:117402186|117420649"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_inputs(paths: dict[str, Path]) -> None:
    # Check BOTH identities before any HDF5 read or output creation.
    for label, path in paths.items():
        if sha256(path) != INPUTS[label][1]:
            raise ValueError(f"{label}: SHA-256 mismatch; refusing to read")


def portable(text: str) -> str:
    return text.replace(str(REPO_ROOT), "<checkout>").replace(sys.prefix, "<environment>")


def operation(function):
    """Record every warning and show it; expected helper failures are evidence."""
    result = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = function()
            record = {"status": "PASS"}
        except Exception as exc:
            record = {
                "status": "FAIL", "exception": f"{type(exc).__name__}: {exc}",
                "traceback": portable(traceback.format_exc()),
            }
    record["warnings"] = [
        {"category": item.category.__name__, "message": str(item.message)}
        for item in caught
    ]
    for item in caught:
        warnings.showwarning(item.message, item.category, portable(item.filename), item.lineno)
    return result, record


def assert_equal(left, right) -> None:
    """Exact values, dtypes, missingness, categorical metadata and sparse buffers."""
    if isinstance(left, pd.DataFrame):
        pd.testing.assert_frame_equal(left, right, check_exact=True)
    elif isinstance(left, pd.Series):
        pd.testing.assert_series_equal(left, right, check_exact=True)
    elif isinstance(left, Mapping):
        assert set(left) == set(right), "mapping inventory changed"
        for key in left:
            assert_equal(left[key], right[key])
    elif sparse.issparse(left):
        assert sparse.issparse(right)
        assert left.format == right.format and left.shape == right.shape
        assert left.dtype == right.dtype
        # These objects use CSR; retaining stored zeros and support positions is
        # stronger than comparing a total or densifying a large matrix.
        assert left.format in {"csr", "csc"}
        for name in ("indptr", "indices", "data"):
            assert_equal(getattr(left, name), getattr(right, name))
    elif isinstance(left, (np.ndarray, np.generic)):
        assert left.dtype == right.dtype, "array dtype changed"
        np.testing.assert_equal(left, right)
    else:
        np.testing.assert_equal(left, right)


def assert_modality(left, right) -> None:
    # Strict order checks precede value comparisons, so positions are tied to
    # stable cell/feature identifiers rather than merely equal dimensions.
    for axis in ("obs", "var"):
        old, new = getattr(left, axis), getattr(right, axis)
        assert old.index.is_unique and new.index.is_unique
        pd.testing.assert_index_equal(old.index, new.index, exact=True)
        assert_equal(old, new)
    assert_equal(left.X, right.X)
    for slot in ("layers", *SLOTS):
        assert_equal(dict(getattr(left, slot)), dict(getattr(right, slot)))
    assert (left.raw is None) == (right.raw is None)
    if left.raw is not None:
        assert_equal(left.raw.X, right.raw.X)
        assert_equal(left.raw.var, right.raw.var)
        assert_equal(dict(left.raw.varm), dict(right.raw.varm))


def stored_reference(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        order = list(handle["mod"].attrs["mod-order"])
        return {
            "axis": int(handle.attrs["axis"]),
            "mod": {name: read_elem(handle["mod"][name]) for name in order},
            **{key: read_elem(handle[key]) for key in ("obs", "var", "obsmap", "varmap", *SLOTS)},
        }


def assert_maps(reference: dict, obj) -> None:
    for axis in ("obs", "var"):
        ids = getattr(obj, axis).index
        assert_equal(reference[axis + "map"], dict(getattr(obj, axis + "map")))
        for name, modality in obj.mod.items():
            positions = np.asarray(getattr(obj, axis + "map")[name]).ravel()
            present = positions > 0
            local_ids = getattr(modality, axis).index
            assert len(np.unique(positions[present])) == len(local_ids)
            np.testing.assert_array_equal(ids[present], local_ids[positions[present] - 1])


def compare_stored_to_native(reference: dict, obj) -> dict:
    assert obj.axis == reference["axis"]
    assert list(reference["mod"]) == list(obj.mod)
    for name, modality in reference["mod"].items():
        assert_modality(modality, obj.mod[name])
    assert_maps(reference, obj)
    for slot in SLOTS:
        assert_equal(reference[slot], dict(getattr(obj, slot)))
    changes = {}
    for axis in ("obs", "var"):
        old, new = reference[axis], getattr(obj, axis)
        pd.testing.assert_index_equal(old.index, new.index, exact=True)
        assert set(old) <= set(new), "stored global column lost"
        # Default MuData pulls modality columns into its in-memory tables.
        # Original column labels, values and missingness must all survive.
        if len(old.columns):
            pd.testing.assert_frame_equal(
                old, new.loc[:, old.columns], check_exact=True,
                check_dtype=False, check_categorical=False,
            )
        dtype_changes = {}
        for column in old:
            try:
                assert_equal(old[column], new[column])
            except AssertionError:
                # Describe representational changes; on-disk round-trip
                # comparison below still requires every original dtype/attr.
                detail = {"before": str(old[column].dtype), "after": str(new[column].dtype)}
                for label, series in (("before", old[column]), ("after", new[column])):
                    if isinstance(series.dtype, pd.CategoricalDtype):
                        detail[label + "_categories"] = series.cat.categories.tolist()
                        detail[label + "_ordered"] = series.cat.ordered
                dtype_changes[column] = detail
        changes[axis] = {
            "stored_shape": list(old.shape), "native_shape": list(new.shape),
            "added_columns": [column for column in new if column not in old],
            "removed_columns": [], "stored_values_and_missingness_equal": True,
            "in_memory_dtype_changes": dtype_changes,
        }
    return changes


def assert_native_equal(left, right) -> None:
    assert left.axis == right.axis and list(left.mod) == list(right.mod)
    for name in left.mod:
        assert_modality(left.mod[name], right.mod[name])
    for slot in ("obs", "var", "obsmap", "varmap", *SLOTS):
        old, new = getattr(left, slot), getattr(right, slot)
        assert_equal(old if isinstance(old, pd.DataFrame) else dict(old),
                     new if isinstance(new, pd.DataFrame) else dict(new))


def compare_hdf5(source: Path, output: Path) -> dict:
    """Reject every unlisted attribute/data change, including unused attributes."""
    result = {"nodes_checked": 0, "datasets_checked": 0, "attributes_checked": 0,
              "allowed_encoder_changes": []}
    with h5py.File(source, "r") as left, h5py.File(output, "r") as right:
        old_names, new_names = [], []
        left.visit(old_names.append)
        right.visit(new_names.append)
        assert old_names == new_names, "HDF5 inventory changed"
        for name in ["/", *old_names]:
            old, new = left[name], right[name]
            result["nodes_checked"] += 1
            assert type(old) is type(new)
            assert set(old.attrs) == set(new.attrs), f"HDF5 attributes changed: {name}"
            for key in old.attrs:
                result["attributes_checked"] += 1
                if key == "encoder-version" and (name == "/" or name in [f"mod/{m}" for m in left["mod"]]):
                    assert new.attrs[key] == importlib.metadata.version("mudata")
                    if old.attrs[key] != new.attrs[key]:
                        result["allowed_encoder_changes"].append(
                            {"node": name, "attribute": key, "before": old.attrs[key], "after": new.attrs[key]}
                        )
                else:
                    assert_equal(old.attrs[key], new.attrs[key])
            if isinstance(old, h5py.Dataset):
                result["datasets_checked"] += 1
                assert old.shape == new.shape and old.dtype == new.dtype, name
                assert_equal(old[()], new[()])
    return result


def metrics(obj) -> dict:
    circ = obj.mod["circ"]
    hosts = circ.var["host_gene"].astype(object).fillna("").astype(str).str.strip()
    result = {
        "modalities": list(obj.mod),
        "shapes": {name: list(a.shape) for name, a in obj.mod.items()},
        "matrices": {
            name: {"dtype": str(a.X.dtype), "format": a.X.format if sparse.issparse(a.X) else "dense",
                   "nonzero_values": int(a.X.count_nonzero() if sparse.issparse(a.X) else np.count_nonzero(a.X)),
                   "total": float(a.X.sum())}
            for name, a in obj.mod.items()
        },
        "overlap": len(set.intersection(*(set(a.obs_names) for a in obj.mod.values()))),
        "annotated_candidates": int((~hosts.str.lower().isin({"", "nan", "none", "na"})).sum()),
        "median_detected": float(np.median(np.asarray((circ.X != 0).sum(axis=1)).ravel())),
        "median_support": float(np.median(np.asarray(circ.X.sum(axis=1)).ravel())),
        "slots": {name: {slot: list(getattr(a, slot)) for slot in ("layers", *SLOTS)}
                  for name, a in obj.mod.items()},
    }
    if CANDIDATE in circ.var_names:
        values = circ[:, [CANDIDATE]].X
        result["selected_candidate"] = {
            "id": CANDIDATE, "cells": int((values != 0).sum()), "support": float(values.sum()),
        }
    return result


def check_expected(label: str, observed: dict) -> None:
    expected_shapes = {"rna": [192, 63187], "circ": [192, 2503]} if label == "smartseq3" else {
        "rna": [23, 63187], "circ": [23, 2443], "cnv": [23, 60607],
    }
    assert observed["shapes"] == expected_shapes
    assert observed["overlap"] == (192 if label == "smartseq3" else 23)
    assert observed["annotated_candidates"] == (2379 if label == "smartseq3" else 2429)
    if label == "smartseq3":
        assert observed["matrices"]["circ"]["nonzero_values"] == 2659
        assert observed["median_detected"] == 12 and observed["median_support"] == 22.5
        assert observed["selected_candidate"] == {"id": CANDIDATE, "cells": 6, "support": 15}


def audit_object(source: Path, output: Path) -> dict:
    reference = stored_reference(source)
    native, read_record = operation(lambda: mu.read_h5mu(source))
    _, helper_record = operation(lambda: circyto_read_h5mu(source))
    result = {
        "native_read": read_record, "circyto_helper_read": helper_record,
        "stored_obs_indices": {
            name: {"name": frame.index.name, "dtype": str(frame.index.dtype),
                   "unique": frame.index.is_unique}
            for name, frame in [("global", reference["obs"]),
                                *((name, a.obs) for name, a in reference["mod"].items())]
        },
    }
    assert read_record["status"] == "PASS", read_record
    result["global_read_changes"] = compare_stored_to_native(reference, native)
    result["metrics"] = metrics(native)
    before_tables = {axis: getattr(native, axis).copy(deep=True) for axis in ("obs", "var")}
    assert not output.exists(), "refusing to overwrite round-trip output"
    _, result["native_write"] = operation(lambda: native.write_h5mu(output))
    assert result["native_write"]["status"] == "PASS"
    restored, result["native_reread"] = operation(lambda: mu.read_h5mu(output))
    assert result["native_reread"]["status"] == "PASS"
    for axis, frame in before_tables.items():
        assert_equal(frame, getattr(restored, axis))
    assert_native_equal(native, restored)
    compare_stored_to_native(reference, restored)
    result["hdf5_roundtrip"] = compare_hdf5(source, output)
    result["roundtrip"] = {"filename": output.name, "sha256": sha256(output),
                           "size_bytes": output.stat().st_size, "purpose": "verification scratch only"}
    result["equivalence"] = "PASS: modality axes/metadata/values/slots; global values/maps/slots; all HDF5 datasets and non-encoder attributes"
    if "cannot join with no overlapping index names" in helper_record.get("exception", ""):
        result["name_only_diagnostic"] = diagnose_index_name(source, reference, output)
    return result


def diagnose_index_name(source: Path, reference: dict, output: Path) -> dict:
    """Single-field experiment on a scratch copy, never a public conversion."""
    names = {a.obs.index.name for a in reference["mod"].values()}
    assert len(names) == 1 and None not in names
    name = names.pop()
    diagnostic = output.with_name(output.stem + ".diagnostic-index-name.h5mu")
    with source.open("rb") as src, diagnostic.open("xb") as dst:
        shutil.copyfileobj(src, dst)
    from anndata.io import write_elem
    with h5py.File(diagnostic, "r+") as handle:
        frame = read_elem(handle["obs"])
        frame.index.name = name
        write_elem(handle, "obs", frame)
    after = stored_reference(diagnostic)
    for key in reference:
        if key == "mod":
            for modality in reference[key]:
                assert_modality(reference[key][modality], after[key][modality])
        else:
            expected = reference[key].rename_axis(name) if key == "obs" else reference[key]
            assert_equal(expected, after[key])
    _, read_record = operation(lambda: circyto_read_h5mu(diagnostic))
    assert read_record["status"] == "PASS"
    return {"change": {"global_obs_index_name_before": reference["obs"].index.name,
                       "global_obs_index_name_after": name},
            "all_other_stored_semantics_equal": True, "helper_read": read_record,
            "filename": diagnostic.name, "sha256": sha256(diagnostic),
            "purpose": "causal diagnostic scratch only; not a proposed public file"}


def environment() -> dict:
    requirements = REPO_ROOT / "manuscript/requirements-results.txt"
    pins = dict(line.split("==") for line in requirements.read_text().splitlines()
                if line and not line.startswith("#"))
    installed = {name: importlib.metadata.version(name) for name in pins}
    assert platform.python_version() == "3.10.20", "Use the documented Python version"
    assert pins == installed, "Use manuscript/requirements-results.txt"
    # Read-only check: ordinary API calls must start with the unmodified 0.3
    # default, not a prior notebook's process-wide override.
    from mudata._core.config import OPTIONS
    assert OPTIONS["pull_on_update"] is None, "Run in a fresh Python process"
    return {"python": sys.version, "platform": platform.platform(), "pins": installed,
            "requirements_sha256": sha256(requirements), "pull_on_update": None,
            "distributions": sorted(f"{d.metadata['Name']}=={d.version}" for d in importlib.metadata.distributions())}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smartseq3", type=Path, required=True)
    parser.add_argument("--imr90", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    paths = {name: getattr(args, name).resolve() for name in INPUTS}
    verify_inputs(paths)
    output = args.output.resolve()
    if output.exists() or output.is_relative_to(REPO_ROOT):
        raise ValueError("Output must be a NEW directory outside the checkout")
    report = {"status": "INCOMPLETE", "remedy": "A: original public files; ordinary MuData API",
              "time_utc": datetime.now(timezone.utc).isoformat(), "environment": environment(),
              "script_sha256": sha256(Path(__file__)), "archive_commit": ARCHIVE,
              "baseline": BASELINE, "objects": {}}
    report["git_head"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip()
    report["production_diff_from_baseline"] = subprocess.check_output(
        ["git", "diff", BASELINE, "--", "circyto", "pyproject.toml"], cwd=REPO_ROOT, text=True,
    )
    assert not report["production_diff_from_baseline"], "Frozen production files changed"
    output.mkdir(parents=True)
    try:
        for label, path in paths.items():
            result = audit_object(path, output / f"{label}.ordinary-roundtrip.h5mu")
            report["objects"][label] = result
            result["original"] = {"archive_member": INPUTS[label][0], "sha256_before": INPUTS[label][1],
                                  "size_bytes": path.stat().st_size}
            check_expected(label, result["metrics"])
        verify_inputs(paths)
        for label, path in paths.items():
            report["objects"][label]["original"]["sha256_after"] = sha256(path)
        report["status"] = "PASS"
    finally:
        (output / "audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("PASS: ordinary MuData real-object round trips; see audit.json for the separate helper results")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
