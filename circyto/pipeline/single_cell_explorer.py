"""Bounded, read-only single-cell reporting. No normalization or embedding fit."""
from __future__ import annotations

from contextlib import ExitStack
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Any
from urllib.parse import urlsplit

import numpy as np
import pandas as pd
from scipy import sparse


SCHEMA = "circyto.explorer.v1"
COMPLETE = {"completed", "completed_detected", "completed_zero", "success", "empty", "skipped_existing"}
FAILED = {"failed", "error", "aborted", "partial_failure"}
UNPROCESSED = {"unprocessed", "not_processed", "pending", "planned", "running", "skipped"}
SAFE_INTEGER = 2**53 - 1


@dataclass(frozen=True)
class ExplorerOptions:
    title: str = "Single-cell circRNA explorer"
    umap_key: str = "X_umap"
    rna_layer: str = "X"
    expression_unit: str = "stored RNA values (units not declared)"
    support_unit: str | None = None
    cell_id_key: str | None = None
    genes: tuple[str, ...] = ()
    circs: tuple[str, ...] = ()
    obs_columns: tuple[str, ...] = ()
    cells: tuple[str, ...] = ()
    max_genes: int = 32
    max_circs: int = 32
    max_cells: int = 50_000
    max_values: int = 2_000_000
    max_bytes: int = 32_000_000


def _scalar(value: Any) -> Any:
    if value is None or not np.isscalar(value) or pd.isna(value):
        return None
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, int) and abs(value) > SAFE_INTEGER:
        return str(value)  # metadata only; measurement precision is checked separately
    return value if isinstance(value, (str, int, float, bool)) else str(value)


def _ids(frame, key, label):
    if key is not None and key not in frame:
        raise ValueError(f"{label}: missing cell identity column {key!r}")
    values = list(frame[key] if key else frame.index)
    if any(not isinstance(v, str) or not v.strip() for v in values):
        raise ValueError(f"{label}: identities must be nonempty strings; no coercion is applied")
    ids = pd.Index(values)
    if not ids.is_unique:
        raise ValueError(f"{label}: duplicate identities; supply unique biological cell IDs")
    return ids


def _state(obs, index, *, rna=False):
    keys = ("rna_status",) if rna else ("circ_observation_status", "detector_status", "processing_status")
    row = obs.iloc[index]
    declared = {k: str(row[k]).lower() if pd.notna(row[k]) else "unknown" for k in keys if k in obs}
    values = list(declared.values())
    if any(v in FAILED for v in values):
        return "failed"
    if any(v in UNPROCESSED for v in values):
        return "unprocessed"
    if "missing" in values:
        return "missing"
    # Both explicit detector/observation declarations must agree on completion.
    # Generic processing metadata supplies completion only when those are absent.
    specific = [v for k, v in declared.items() if k != "processing_status"] or values
    state = ("completed" if all(v in COMPLETE for v in specific) else "unknown") if specific else ("available" if rna else "unknown")
    if not rna and "alignment_status" in obs:
        alignment = str(obs.iloc[index]["alignment_status"]).lower()
        if alignment in FAILED:
            return "failed"
        if alignment in UNPROCESSED:
            return "unprocessed"
        if alignment not in {"aligned", "reused_input", "reused_cached", "completed"}:
            return "unknown"
    return state


def _metadata(frame, index, *, circ=False):
    fields = {}
    for key, value in frame.iloc[index].items():
        key = str(key)
        # Exclude inherited cohort summaries: all displayed rates are recomputed for this view.
        allowed = (re.search(r"gene|symbol|alias|annotation|catalog|circatlas|circbase|circsc|uniform_id|match_class|evidence|doi|source$", key, re.I)
                   or key in {"circ_id", "chrom", "chr", "chromosome", "seqname", "start", "end", "strand",
                              "genome_build", "assembly", "coordinate_convention", "support_definition", "detector"})
        if allowed and not re.search(r"path|command|file|token|password", key, re.I):
            scalar = _scalar(value)
            if scalar is not None and scalar != "":
                if len(str(scalar)) > 8192:
                    raise ValueError(f"Feature annotation {key!r} exceeds 8192 characters")
                fields[key] = scalar
    return fields


def _aliases(identity, metadata, *, circ):
    aliases = {identity}
    for key, value in metadata.items():
        if (re.search(r"gene|symbol|alias|atlas.*id|base.*id|sc.*id|database_ids|uniform_id|circ_id", key, re.I)
                and not re.search(r"source|type|status|count|length|annotation", key, re.I)):
            # Uniform circRNA names can contain exon lists: circGENE(2,3,4).1.
            aliases.add(str(value))
            aliases.update(v.strip() for v in re.split(r";|,(?![^()]*\))", str(value)) if v.strip())
    aliases = {a for a in aliases if a not in {"NA", "nan", "None", "not_in_primary_annotation_scope"}}
    if circ:
        chrom = next((metadata[k] for k in ("chrom", "chr", "chromosome", "seqname") if k in metadata), None)
        if chrom is not None and "start" in metadata and "end" in metadata:
            for sep in ("|", "-"):
                coord = f"{chrom}:{metadata['start']}{sep}{metadata['end']}"
                aliases.add(coord)
                if "strand" in metadata:
                    aliases.add(f"{coord}:{metadata['strand']}")
    aliases.update(re.sub(r"\.\d+$", "", a) for a in list(aliases) if re.match(r"ENS\w*G\d+", a))
    return sorted(aliases)


def _catalog(adata, *, circ, sidecar=None):
    if adata is None:
        return []
    identity_key = None
    if not circ and not adata.var_names.is_unique and "gene_id" in adata.var:
        identity_key = "gene_id"
    identities = _ids(adata.var, identity_key, "circRNA feature" if circ else "RNA feature")
    result = []
    for i, identity in enumerate(identities):
        meta = _metadata(adata.var, i, circ=circ)
        if sidecar is not None and identity in sidecar.index:
            extra = _metadata(sidecar, sidecar.index.get_loc(identity), circ=True)
            for key, value in extra.items():
                if key in meta and str(meta[key]) != str(value):
                    raise ValueError(f"Candidate metadata conflicts with object: {identity}, {key}")
            meta.update(extra)
        aliases = _aliases(identity, meta, circ=circ)
        if str(adata.var_names[i]) not in aliases:
            aliases.append(str(adata.var_names[i]))
        result.append(dict(id=identity, metadata=meta, aliases=aliases))
    return result


def _choose(catalog, queries, limit, label, preferred=()):
    selected = []
    if queries:
        for query in queries:
            matches = [i for i, row in enumerate(catalog) if query in row["aliases"]]
            if not matches:
                raise ValueError(f"{label} selector {query!r} has no exact ID, coordinate or alias match")
            selected.extend(i for i in matches if i not in selected)
        if len(selected) > limit:
            raise ValueError(f"{label} selectors resolve to {len(selected)} features, above limit {limit}; duplicate symbols are not collapsed")
    else:
        preferred = set(preferred)
        order = sorted(range(len(catalog)), key=lambda i: (not bool(preferred.intersection(catalog[i]["aliases"])), catalog[i]["id"]))
        selected = order[:limit]
    return selected


def _columns(matrix, indices):
    """Read only selected columns, including backed CSR/CSC, without full densification."""
    if matrix is None:
        raise ValueError("Selected measurement matrix is unavailable")
    for j in indices:
        column = matrix[:, j:j + 1]
        if sparse.issparse(column):
            column = column.toarray()
        yield np.asarray(column).reshape(-1)


def _features(adata, catalog, indices, ids, cells, base_states, *, circ, layer="X", unit=None):
    result = []
    if adata is None or not indices:
        return result
    if layer != "X" and layer not in adata.layers:
        raise ValueError(f"RNA layer {layer!r} is unavailable")
    matrix = adata.X if layer == "X" else adata.layers[layer]
    positions = ids.get_indexer(cells)
    for j, column in zip(indices, _columns(matrix, indices)):
        if column.dtype.kind not in "biuf":
            raise ValueError("Selected measurement matrix must contain real numeric values")
        finite = np.isfinite(column)
        if np.any((column[finite] > SAFE_INTEGER) | (column[finite] < -SAFE_INTEGER)):
            raise ValueError("Selected measurement exceeds JavaScript's exact numeric range; report refused")
        if circ and np.any((column[finite] < 0) | (column[finite] != np.floor(column[finite]))):
            raise ValueError("circRNA X must contain nonnegative integer native junction support")
        if circ and sum(int(v) for v in column[finite]) > SAFE_INTEGER:
            raise ValueError("Selected circRNA total exceeds JavaScript's exact numeric range; report refused")
        values, states = [], []
        for i, p in enumerate(positions):
            value = _scalar(column[p]) if p >= 0 else None
            values.append(value)
            state = base_states[i]
            if value is None and state in {"completed", "available"}:
                state = "missing"
            if circ and state == "completed":
                state = "completed_detected" if value > 0 else "completed_zero"
            states.append(state)
        feature = dict(catalog[j])
        feature.update(values=values, states=states, unit=unit or feature["metadata"].get("support_definition")
                       or adata.uns.get("circyto", {}).get("support_definition") or "native junction support (units not declared)")
        result.append(feature)
    return result


def _obs_metadata(objects, cells, options, warnings):
    candidates = {}
    for prefix, obj, ids in objects:
        if obj is None:
            continue
        for column in obj.obs:
            key = f"{prefix}:{column}"
            if re.search(r"path|command|file|token|password", str(column), re.I) or column == options.cell_id_key:
                continue
            series = obj.obs[column].copy()
            series.index = ids
            values = [_scalar(v) for v in series.reindex(cells)]
            if not any(v is not None for v in values):
                continue
            numeric = (pd.api.types.is_numeric_dtype(series.dtype) and not pd.api.types.is_bool_dtype(series.dtype)
                       and all(v is None or isinstance(v, (int, float)) for v in values))
            candidates[key] = dict(key=key, kind="numeric" if numeric else "category", values=values)
    if options.obs_columns:
        missing = set(options.obs_columns) - candidates.keys()
        if missing:
            raise ValueError(f"Unknown or empty observation columns: {sorted(missing)}")
        keys = list(dict.fromkeys(options.obs_columns))
    else:
        keys = sorted(candidates, key=lambda k: (not bool(re.search(r"cluster|leiden|annotation|cell_type|condition|exposure", k, re.I)), k))[:48]
        if len(candidates) > len(keys):
            warnings.append(f"Observation fields bounded to {len(keys)} of {len(candidates)}; use --obs-column to choose fields.")
    if len(keys) > 64:
        raise ValueError("At most 64 observation fields can be embedded")
    for key in keys:
        if any(len(str(v)) > 2048 for v in candidates[key]["values"] if v is not None):
            raise ValueError(f"Observation field {key!r} exceeds the 2048-character value limit")
    return [candidates[k] for k in keys]


def build_explorer(data=None, rna=None, *, options=None, candidate_metadata=None):
    """Create a JSON-safe payload from circ AnnData or MuData and optional RNA AnnData.

    A separately supplied RNA object supplies both expression and its embedding.
    All joins use exact unique identities; input objects are never modified.
    """
    options = options or ExplorerOptions()
    if not options.expression_unit.strip() or (options.support_unit is not None and not options.support_unit.strip()):
        raise ValueError("Measurement unit labels must not be empty")
    for key in ("max_genes", "max_circs", "max_cells", "max_values", "max_bytes"):
        if getattr(options, key) < (0 if key in {"max_genes", "max_circs"} else 1):
            raise ValueError(f"Invalid {key}")
    global_obj = data if hasattr(data, "mod") else None
    circ = data.mod.get("circ") if global_obj is not None else data
    if rna is None and global_obj is not None:
        rna = data.mod.get("rna")
    if circ is None and rna is None:
        raise ValueError("No circ or RNA modality is available")
    circ_ids = _ids(circ.obs, options.cell_id_key, "circ cells") if circ is not None else pd.Index([])
    rna_ids = _ids(rna.obs, options.cell_id_key, "RNA cells") if rna is not None else pd.Index([])
    cells = rna_ids.union(circ_ids, sort=False)
    total_cells = len(cells)
    if options.cells:
        requested = pd.Index(options.cells)
        if not requested.is_unique or any(not isinstance(c, str) or not c.strip() for c in requested):
            raise ValueError("Selected cells must have unique nonempty exact identities")
        if len(requested.difference(cells)):
            raise ValueError("Selected cell IDs are absent from both modalities")
        cells = requested
    if len(cells) > options.max_cells:
        raise ValueError(f"{len(cells)} cells exceed --max-cells {options.max_cells}; supply an explicitly scoped object")
    warnings = []
    if len(circ_ids) and len(rna_ids) and not len(circ_ids.intersection(rna_ids)):
        raise ValueError("RNA and circ have no exact biological cell ID overlap")
    if candidate_metadata is not None:
        _ids(candidate_metadata, None, "candidate metadata")
        unknown = candidate_metadata.index.difference(circ.var_names if circ is not None else [])
        if len(unknown):
            raise ValueError(f"Candidate metadata has {len(unknown)} IDs absent from circ")
    circ_catalog = _catalog(circ, circ=True, sidecar=candidate_metadata)
    gene_catalog = _catalog(rna, circ=False)
    ci = _choose(circ_catalog, options.circs, options.max_circs, "circRNA")
    hosts = [alias for j in ci for key, value in circ_catalog[j]["metadata"].items() if "gene" in key
             for alias in re.split(r"[;,]", str(value))]
    gi = _choose(gene_catalog, options.genes, options.max_genes, "Gene", preferred=hosts) if rna is not None else []
    if len(cells) * (len(ci) + len(gi)) > options.max_values:
        raise ValueError("Selected cell × feature values exceed --max-values; lower --max-genes/--max-circs or scope the input")
    base = {}
    for name, obj, ids in (("circ", circ, circ_ids), ("rna", rna, rna_ids)):
        base[name] = [_state(obj.obs, p, rna=name == "rna") if p >= 0 else "missing" for p in ids.get_indexer(cells)]
    coordinates = [None] * len(cells)
    if rna is not None and options.umap_key in rna.obsm:
        xy = np.asarray(rna.obsm[options.umap_key])
        if xy.ndim != 2 or xy.shape != (rna.n_obs, 2):
            raise ValueError("Selected RNA embedding must have exactly two columns")
        for i, p in enumerate(rna_ids.get_indexer(cells)):
            if p >= 0 and np.isfinite(xy[p]).all():
                coordinates[i] = xy[p].tolist()
    else:
        warnings.append(f"RNA embedding {options.umap_key!r} unavailable; UMAP panels omitted.")
    if rna is None:
        warnings.append("RNA unavailable; gene expression panels omitted.")
    plotted = sum(x is not None for x in coordinates)
    if plotted < len(cells):
        warnings.append(f"{len(cells) - plotted} cells lack finite RNA coordinates; retained in tables and denominators.")
    objects = [("rna", rna, rna_ids), ("circ", circ, circ_ids)]
    if global_obj is not None:
        objects.append(("global", global_obj, _ids(global_obj.obs, options.cell_id_key, "MuData cells")))
    metadata = _obs_metadata(objects, cells, options, warnings)
    circs = _features(circ, circ_catalog, ci, circ_ids, cells, base["circ"], circ=True, unit=options.support_unit)
    genes = _features(rna, gene_catalog, gi, rna_ids, cells, base["rna"], circ=False, layer=options.rna_layer, unit=options.expression_unit)
    for feature in circs:
        feature["evidence_links"] = []
        for key, value in feature["metadata"].items():
            if re.search(r"evidence.*(?:url|link)|doi", key, re.I):
                for link in str(value).split(";"):
                    link = link.strip()
                    if "doi" in key.lower() and link.startswith("10."):
                        link = "https://doi.org/" + link
                    parsed = urlsplit(link)
                    if parsed.scheme in {"https", "http"} and parsed.netloc and not parsed.username and not any(c.isspace() for c in link):
                        feature["evidence_links"].append(dict(label=key, url=link))
    payload = dict(schema_version=SCHEMA, title=options.title, cells=cells.tolist(), coordinates=coordinates,
                   metadata=metadata, genes=genes, circs=circs, warnings=warnings,
                   embedding=dict(key=options.umap_key, source="RNA obsm; unchanged", plotted_cells=plotted),
                   selection=dict(cells_total=total_cells, cells_embedded=len(cells), genes_total=len(gene_catalog), genes_embedded=len(genes), circs_total=len(circ_catalog), circs_embedded=len(circs),
                                  policy="Explicit exact matches (all duplicate symbols); otherwise stable ID order, RNA host genes first. No expression ranking or normalization.",
                                  rna_layer=options.rna_layer, max_values=options.max_values, max_bytes=options.max_bytes),
                   join=dict(policy="Exact unique biological cell identity; union retained unless an explicit cell list scopes the view; no positional join", rna_cells=len(cells.intersection(rna_ids)), circ_cells=len(cells.intersection(circ_ids)), overlap=len(cells.intersection(rna_ids).intersection(circ_ids))),
                   sources=[])
    encode_payload(payload, max_bytes=options.max_bytes)
    return payload


def encode_payload(payload, *, max_bytes=32_000_000):
    text = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    if len(text.encode("utf-8")) > max_bytes:
        raise ValueError("Explorer JSON exceeds --max-bytes; reduce embedded features or metadata")
    return text


def load_explorer(path: Path, rna_path: Path | None = None, *, options=None, candidate_metadata_path=None, cells_path=None):
    """Read backed matrices where supported. Close all read-only file handles."""
    import anndata as ad

    options = options or ExplorerOptions()
    sources = []
    def record(p):
        digest = hashlib.sha256()
        with p.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        sources.append(dict(name=p.name, bytes=p.stat().st_size, sha256=digest.hexdigest()))
    with ExitStack() as stack:
        if cells_path is not None:
            from dataclasses import replace

            record(Path(cells_path))
            options = replace(options, cells=tuple(Path(cells_path).read_text(encoding="utf-8").splitlines()))
            if not options.cells:
                raise ValueError("Selected cells file is empty")
        def read(p):
            record(p)
            if p.suffix.lower() == ".h5mu":
                try:
                    from circyto.multimodal.sync import read_h5mu
                    obj = read_h5mu(p, backed="r")
                except ImportError as exc:
                    raise ValueError("MuData input requires circyto[mudata]") from exc
            elif p.suffix.lower() == ".h5ad":
                obj = ad.read_h5ad(p, backed="r")
            else:
                raise ValueError("Explorer input must be .h5ad or .h5mu")
            stack.callback(obj.file.close)
            return obj
        data = read(Path(path))
        rna = read(Path(rna_path)) if rna_path is not None else None
        if rna is not None and hasattr(rna, "mod"):
            raise ValueError("--rna must be an RNA AnnData .h5ad")
        sidecar = None
        if candidate_metadata_path is not None:
            record(Path(candidate_metadata_path))
            sidecar = pd.read_csv(candidate_metadata_path, sep="\t", dtype=str, keep_default_na=False)
            if "circ_id" not in sidecar:
                raise ValueError("Candidate metadata requires circ_id")
            sidecar = sidecar.set_index("circ_id")
        payload = build_explorer(data, rna, options=options, candidate_metadata=sidecar)
    payload["sources"] = sources
    encode_payload(payload, max_bytes=options.max_bytes)
    return payload
