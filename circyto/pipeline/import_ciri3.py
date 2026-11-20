"""Evidence-backed adapter for CIRI3's 13-column single-sample result table.

Conversion delegates to the existing CIRI3 normalizer and import-circ sparse
builder. Native files and mappings remain separate, immutable source evidence.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd

from circyto.detectors.ciri3 import _normalize_ciri3_output

NATIVE_SCHEMA = "ciri3.native-result-13.v1"
NATIVE_HEADER = [
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
MAPPING_HEADER = ["cell_id", "sample_id", "path", "detector_version", "result_scope"]
SUPPORT = "CIRI3_junction_reads"


def plan_native(manifest: Path):
    from circyto.pipeline.import_circ import _ids, _table

    mapping = _table(manifest, MAPPING_HEADER)
    if list(mapping.columns) != MAPPING_HEADER:
        raise ValueError(
            "CIRI3 manifest must contain exactly cell_id,sample_id,path,detector_version,result_scope in that order; declarations belong in command options"
        )
    if mapping.empty:
        raise ValueError("CIRI3 manifest must map at least one single-sample result")
    _ids(mapping.cell_id, "CIRI3 mapping cell")
    _ids(mapping.sample_id, "CIRI3 mapping sample")
    _ids(mapping.path, "CIRI3 source path", unique=False)
    _ids(mapping.detector_version, "CIRI3 detector version", unique=False)
    if not mapping.result_scope.eq("single-sample").all():
        raise ValueError(
            "CIRI3 result_scope must be single-sample; pooled/multi-sample result.txt is not cell-resolved"
        )
    if not mapping.detector_version.isin(["3.0.1", "unknown"]).all():
        raise ValueError(
            "Native CIRI3 version not reviewed: use 3.0.1, or unknown only when the producing version is unavailable; other versions require a reviewed long-table conversion"
        )
    if mapping.detector_version.nunique() != 1:
        raise ValueError("Mixed CIRI3 detector versions; import separately")
    sources = {"ciri3_manifest.tsv": manifest}
    samples = []
    seen = set()
    for number, row in enumerate(mapping.to_dict("records"), 1):
        path = Path(row["path"])
        if not path.is_absolute():
            path = manifest.parent / path
        path = path.resolve(strict=True)
        identity = (path.stat().st_dev, path.stat().st_ino)
        if identity in seen:
            raise ValueError(
                "One native CIRI3 file is mapped more than once; do not duplicate a sample across cells"
            )
        seen.add(identity)
        key = f"ciri3_native_{number:06d}.tsv"
        sources[key] = path
        samples.append(
            {
                "cell_id": row["cell_id"],
                "sample_id": row["sample_id"],
                "detector_version": row["detector_version"],
                "source_file": key,
                "result_scope": row["result_scope"],
            }
        )
    return sources, {
        "schema": NATIVE_SCHEMA,
        "samples": samples,
        "completion_policy": "No completion inferred from result files or logs; use explicit cell metadata",
    }


def load_native(sources, info, obs, declarations):
    from circyto.pipeline.import_circ import _ids, _integers, _load_calls, _table

    if declarations["detector"] != "CIRI3":
        raise ValueError("Native CIRI3 input requires --detector CIRI3")
    if declarations["coordinate_convention"] != "1-based-closed":
        raise ValueError(
            "Native CIRI3 requires --coordinate-convention 1-based-closed; no coordinate shifts are applied"
        )
    if declarations["support_definition"] != SUPPORT:
        raise ValueError(
            f"Native CIRI3 requires --support-definition {SUPPORT}; Score, FSJ, ratios and UMIs are not BSJ counts"
        )
    cells = list(obs.index)
    mapped = {row["cell_id"] for row in info["samples"]}
    if not mapped <= set(cells):
        raise ValueError(
            "CIRI3 manifest contains cell IDs absent from --cells metadata"
        )
    if any(
        obs.loc[cell, "detector_status"] == "completed" for cell in set(cells) - mapped
    ):
        raise ValueError(
            "Completed cells require a mapped native result, including a header-only result for zero candidates"
        )
    frames = []
    annotations = []
    with tempfile.TemporaryDirectory(prefix="circyto-ciri3-") as tmp:
        normalized = Path(tmp) / "normalized.tsv"
        for sample in info["samples"]:
            source = sources[sample["source_file"]]
            raw = _table(source, NATIVE_HEADER)
            if list(raw.columns) != NATIVE_HEADER:
                raise ValueError(
                    f"{source.name}: unsupported native CIRI3 layout; expected exactly the documented 13 columns including Score, in order"
                )
            _ids(raw.circRNA_ID, "native CIRI3 candidate")
            _ids(raw.chr, "native CIRI3 chromosome", unique=False)
            _integers(raw["#junction_reads"], "native #junction_reads")
            start = _integers(raw.circRNA_start, "native start")
            end = _integers(raw.circRNA_end, "native end")
            expected = [f"{chrom}:{s}|{e}" for chrom, s, e in zip(raw.chr, start, end)]
            if raw.circRNA_ID.tolist() != expected:
                raise ValueError(
                    f"{source.name}: circRNA_ID disagrees with chr/start/end; expected chr:start|end without a strand suffix"
                )
            if not raw.strand.isin(["+", "-"]).all():
                raise ValueError(f"{source.name}: native CIRI3 strand must be + or -")
            _normalize_ciri3_output(source, normalized)
            frame = _table(
                normalized, ["circ_id", "chr", "start", "end", "strand", "support"]
            )
            if len(frame) != len(raw):
                raise ValueError("CIRI3 normalization changed the number of candidates")
            frame = frame.rename(columns={"chr": "chrom"})
            frame["cell_id"] = sample["cell_id"]
            frames.append(frame)
            annotations.append(
                raw[["circRNA_ID", "gene_id", "circRNA_type"]].rename(
                    columns={
                        "circRNA_ID": "circ_id",
                        "gene_id": "ciri3_gene_id",
                        "circRNA_type": "ciri3_circRNA_type",
                    }
                )
            )
        calls = pd.concat(frames, ignore_index=True)
        annotated = pd.concat(annotations, ignore_index=True)
        if (
            annotated.groupby("circ_id", sort=False)[
                ["ciri3_gene_id", "ciri3_circRNA_type"]
            ]
            .nunique()
            .gt(1)
            .any()
            .any()
        ):
            raise ValueError(
                "Conflicting native CIRI3 annotations across samples; import separately or review annotations explicitly"
            )
        path = Path(tmp) / "calls.tsv"
        calls.to_csv(path, sep="\t", index=False)
        X, var, _ = _load_calls(path, cells, declarations)
    annotated = (
        annotated.drop_duplicates("circ_id").set_index("circ_id").reindex(var.index)
    )
    for col in annotated:
        var[col] = annotated[col]
    # CIRI3 gene_id is not a gene symbol; preserve it without fabricating host_gene.
    var["host_gene_id"] = var.ciri3_gene_id.where(
        ~var.ciri3_gene_id.isin(["NA", ""]), ""
    )
    var["host_gene_id_source"] = var.host_gene_id.map(
        lambda x: "CIRI3_gene_id" if x else ""
    )
    for key in ("sample_id", "detector_version", "source_file"):
        field = "ciri3_" + key
        if field in obs:
            raise ValueError(
                f"Reserved native metadata field {field}; rename the input column"
            )
        values = {row["cell_id"]: row[key] for row in info["samples"]}
        obs[field] = [values.get(cell, "") for cell in cells]
    return X, var
