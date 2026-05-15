from __future__ import annotations

import csv
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


DEFAULT_FASTERQ_THREADS = 8
CURATED_RUN_COLUMNS = [
    "dataset_id",
    "run_id",
    "sample_id",
    "protocol",
    "source",
    "fastq_1",
    "fastq_2",
    "notes",
]
SELECTED_RUN_COLUMNS = [
    "dataset_id",
    "run_id",
    "sample_id",
    "protocol",
    "organism",
    "expected_read_layout",
    "expected_reference",
    "recommended_route",
    "source",
    "fastq_1",
    "fastq_2",
    "notes",
]


@dataclass(frozen=True)
class PublicDatasetRow:
    run_id: str
    sample_id: str
    source: str
    notes: str
    fastq_1: str = ""
    fastq_2: str = ""


@dataclass(frozen=True)
class PublicDatasetSpec:
    dataset_id: str
    protocol: str
    organism: str
    expected_read_layout: str
    expected_reference: str
    recommended_route: str
    title: str
    notes: str
    plan_kind: str
    placeholder_rows: tuple[PublicDatasetRow, ...]


_DATASETS: dict[str, PublicDatasetSpec] = {
    "e-mtab-8735": PublicDatasetSpec(
        dataset_id="E-MTAB-8735",
        protocol="smartseq3",
        organism="Homo sapiens",
        expected_read_layout="paired-end",
        expected_reference="hg38",
        recommended_route="Smart-seq3 paired-end workflow",
        title="Smart-seq3 192-cell benchmark planning target",
        notes=(
            "Existing validated Smart-seq3 benchmark used in prior IFReC workflow work. "
            "This fixture includes only a tiny subset of placeholder rows; full accession expansion can be added later."
        ),
        plan_kind="ena-fastq",
        placeholder_rows=(
            PublicDatasetRow(
                run_id="TODO_REAL_ACCESSION_EMTAB8735_001",
                sample_id="EMTAB8735_cell_001",
                source="ArrayExpress / ENA planning placeholder",
                fastq_1="TODO_REAL_ACCESSION_FASTQ_1",
                fastq_2="TODO_REAL_ACCESSION_FASTQ_2",
                notes="TODO_REAL_ACCESSION placeholder Smart-seq3 run for plan generation only; full accession expansion can be added later.",
            ),
            PublicDatasetRow(
                run_id="TODO_REAL_ACCESSION_EMTAB8735_002",
                sample_id="EMTAB8735_cell_002",
                source="ArrayExpress / ENA planning placeholder",
                fastq_1="TODO_REAL_ACCESSION_FASTQ_1",
                fastq_2="TODO_REAL_ACCESSION_FASTQ_2",
                notes="TODO_REAL_ACCESSION placeholder Smart-seq3 run for plan generation only; full accession expansion can be added later.",
            ),
            PublicDatasetRow(
                run_id="TODO_REAL_ACCESSION_EMTAB8735_003",
                sample_id="EMTAB8735_cell_003",
                source="ArrayExpress / ENA planning placeholder",
                fastq_1="TODO_REAL_ACCESSION_FASTQ_1",
                fastq_2="TODO_REAL_ACCESSION_FASTQ_2",
                notes="TODO_REAL_ACCESSION placeholder Smart-seq3 run for plan generation only; full accession expansion can be added later.",
            ),
        ),
    ),
    "gse98664": PublicDatasetSpec(
        dataset_id="GSE98664",
        protocol="ramda",
        organism="Mus musculus",
        expected_read_layout="single-end",
        expected_reference="mm10/mm39",
        recommended_route="BWA+CIRI3 single-end",
        title="Original RamDA-seq RNA-only full-length total RNA benchmark",
        notes=(
            "Mouse RamDA engineering fixture only. Do not use this dataset for human hg38 biological validation. "
            "Fixture rows are placeholder-style run records to keep tests network-free and reproducible."
        ),
        plan_kind="sra-prefetch",
        placeholder_rows=(
            PublicDatasetRow(
                run_id="TODO_REAL_ACCESSION_GSE98664_001",
                sample_id="GSE98664_cell_001",
                source="GEO / SRA planning placeholder",
                notes="TODO_REAL_ACCESSION placeholder run ID for planning-only tests; replace with expanded public metadata later.",
            ),
            PublicDatasetRow(
                run_id="TODO_REAL_ACCESSION_GSE98664_002",
                sample_id="GSE98664_cell_002",
                source="GEO / SRA planning placeholder",
                notes="TODO_REAL_ACCESSION placeholder run ID for planning-only tests; replace with expanded public metadata later.",
            ),
            PublicDatasetRow(
                run_id="TODO_REAL_ACCESSION_GSE98664_003",
                sample_id="GSE98664_cell_003",
                source="GEO / SRA planning placeholder",
                notes="TODO_REAL_ACCESSION placeholder run ID for planning-only tests; replace with expanded public metadata later.",
            ),
        ),
    ),
    "gse278944": PublicDatasetSpec(
        dataset_id="GSE278944",
        protocol="scrr",
        organism="Homo sapiens / Mus musculus",
        expected_read_layout="paired-end",
        expected_reference="hg38 and mm10/mm39, dataset-dependent",
        recommended_route="Future scRepli-RamDA-seq workflow; not first CIRI3 validation target",
        title="Future scRepli-RamDA-seq DNA replication plus full-length RNA benchmark",
        notes="Fixture rows are placeholder records only; no live accession expansion is required in tests.",
        plan_kind="sra-prefetch",
        placeholder_rows=(
            PublicDatasetRow(
                run_id="TODO_REAL_ACCESSION_GSE278944_001",
                sample_id="GSE278944_cell_001",
                source="GEO / SRA planning placeholder",
                notes="TODO_REAL_ACCESSION placeholder run ID for future scRR planning only.",
            ),
            PublicDatasetRow(
                run_id="TODO_REAL_ACCESSION_GSE278944_002",
                sample_id="GSE278944_cell_002",
                source="GEO / SRA planning placeholder",
                notes="TODO_REAL_ACCESSION placeholder run ID for future scRR planning only.",
            ),
        ),
    ),
    "shin-ramda-riken": PublicDatasetSpec(
        dataset_id="shin-ramda-riken",
        protocol="shin-ramda",
        organism="Homo sapiens",
        expected_read_layout="likely single-end for primary Shin-RamDA RNA use; verify per accession",
        expected_reference="hg38",
        recommended_route="BWA+CIRI3 single-end once accession-level FASTQs are confirmed",
        title="Shin-RamDA workflow reference from rikenbit/shin-ramda-seq-paper",
        notes="README-oriented planning target only. Do not attempt FASTQ download planning for this dataset alias.",
        plan_kind="readme-only",
        placeholder_rows=(
            PublicDatasetRow(
                run_id="",
                sample_id="workflow_reference",
                source="https://github.com/rikenbit/shin-ramda-seq-paper",
                notes="Workflow reference only; inspect the upstream repository for protocol details and future accession expansion.",
            ),
        ),
    ),
}

CURATED_TABLE_FILENAMES = {
    "E-MTAB-8735": "emtab8735_runs.tsv",
    "GSE98664": "gse98664_runs.tsv",
}


def list_supported_dataset_ids() -> list[str]:
    return sorted(spec.dataset_id for spec in _DATASETS.values())


def get_public_dataset_spec(dataset_id: str) -> PublicDatasetSpec:
    key = dataset_id.strip().lower()
    if key not in _DATASETS:
        supported = ", ".join(list_supported_dataset_ids())
        raise ValueError(f"Unknown dataset id '{dataset_id}'. Supported dataset ids: {supported}")
    return _DATASETS[key]


def _normalize_protocol(protocol: str) -> str:
    return protocol.strip().lower()


def _validate_protocol(spec: PublicDatasetSpec, requested_protocol: str) -> None:
    requested = _normalize_protocol(requested_protocol)
    actual = _normalize_protocol(spec.protocol)
    if requested != actual:
        raise ValueError(
            f"Dataset {spec.dataset_id} uses protocol '{spec.protocol}', not '{requested_protocol}'."
        )


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _curated_run_table_path(spec: PublicDatasetSpec) -> Path | None:
    filename = CURATED_TABLE_FILENAMES.get(spec.dataset_id)
    if filename is None:
        return None
    return _repo_root() / "testdata" / "public_datasets" / filename


def _read_rows_from_tsv(path: Path) -> list[PublicDatasetRow]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        missing = [column for column in CURATED_RUN_COLUMNS if column not in header]
        if missing:
            raise ValueError(f"Curated run table {path} is missing required columns: {', '.join(missing)}")
        rows: list[PublicDatasetRow] = []
        for raw in reader:
            rows.append(
                PublicDatasetRow(
                    run_id=(raw.get("run_id") or "").strip(),
                    sample_id=(raw.get("sample_id") or "").strip(),
                    source=(raw.get("source") or "").strip(),
                    fastq_1=(raw.get("fastq_1") or "").strip(),
                    fastq_2=(raw.get("fastq_2") or "").strip(),
                    notes=(raw.get("notes") or "").strip(),
                )
            )
    return rows


def _resolve_rows(spec: PublicDatasetSpec) -> tuple[list[PublicDatasetRow], str, Path | None]:
    curated_path = _curated_run_table_path(spec)
    if curated_path is not None and curated_path.exists():
        return _read_rows_from_tsv(curated_path), "curated", curated_path
    return list(spec.placeholder_rows), "placeholder", curated_path


def _selected_rows(rows: Sequence[PublicDatasetRow], max_runs: int | None) -> list[PublicDatasetRow]:
    rows = list(rows)
    if max_runs is None:
        return rows
    if max_runs < 1:
        raise ValueError("--max-runs must be >= 1 when provided.")
    return rows[:max_runs]


def _shell_quote(path: Path) -> str:
    return shlex.quote(str(path))


def _warning_lines(spec: PublicDatasetSpec, rows: Sequence[PublicDatasetRow]) -> list[str]:
    warnings: list[str] = []
    if spec.dataset_id == "GSE98664":
        warnings.extend(
            [
                f"Dataset {spec.dataset_id} is annotated as {spec.organism}.",
                f"Recommended references: {spec.expected_reference}.",
                "Do not use hg38 for biological validation.",
            ]
        )

    for row in rows:
        has_r1 = bool(row.fastq_1.strip())
        has_r2 = bool(row.fastq_2.strip())
        if spec.expected_read_layout == "single-end" and has_r2:
            warnings.append(
                f"Row {row.sample_id or row.run_id} includes fastq_2 but dataset metadata expects single-end input."
            )
        if spec.expected_read_layout == "paired-end" and has_r1 and not has_r2:
            warnings.append(
                f"Row {row.sample_id or row.run_id} appears single-end but dataset metadata expects paired-end input."
            )
    return warnings


def _render_download_plan(
    spec: PublicDatasetSpec,
    rows: Sequence[PublicDatasetRow],
    outdir: Path,
    download_method: str,
    row_mode: str,
    warnings: Sequence[str],
) -> str:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f"# Dataset: {spec.dataset_id}",
        f"# Protocol: {spec.protocol}",
        f"# Organism: {spec.organism}",
        f"# Expected read layout: {spec.expected_read_layout}",
        f"# Expected reference: {spec.expected_reference}",
        f"# Recommended route: {spec.recommended_route}",
        f"# Plan kind: {spec.plan_kind}",
        f"# Row mode: {row_mode}",
        f"# Download method: {download_method}",
        "",
    ]
    if warnings:
        lines.append("# WARNING:")
        for warning in warnings:
            lines.append(f"# {warning}")
        lines.append("")

    if spec.plan_kind == "sra-prefetch":
        lines.extend(
            [
                f"mkdir -p {_shell_quote(outdir / 'sra')}",
                f"mkdir -p {_shell_quote(outdir / 'fastq')}",
                "",
            ]
        )
        if download_method == "none":
            lines.append(f"# Download disabled. This plan records the selected {row_mode} SRA runs only.")
        elif download_method in {"sra", "ena"}:
            for row in rows:
                lines.append(f"prefetch {row.run_id}")
                lines.append(
                    "fasterq-dump "
                    f"{row.run_id} --split-files --threads {DEFAULT_FASTERQ_THREADS} -O {_shell_quote(outdir / 'fastq')}"
                )
        else:
            raise ValueError(f"Unsupported download method: {download_method}")
    elif spec.plan_kind == "ena-fastq":
        if download_method == "none":
            lines.append(f"# Download disabled. Review selected_runs.tsv for {row_mode} ENA/ArrayExpress rows.")
        elif download_method in {"ena", "sra"}:
            lines.append(f"mkdir -p {_shell_quote(outdir / 'fastq')}")
            lines.append(f"# {row_mode.capitalize()} ENA/ArrayExpress-style URLs for planning only.")
            for row in rows:
                fastq_1 = row.fastq_1 or "FASTQ_1_URL_PLACEHOLDER"
                fastq_2 = row.fastq_2 or "FASTQ_2_URL_PLACEHOLDER"
                lines.append(f"# {row.run_id}\t{fastq_1}\t{fastq_2}")
        else:
            raise ValueError(f"Unsupported download method: {download_method}")
    elif spec.plan_kind == "readme-only":
        lines.append("# README-only plan. No download commands are generated for this dataset alias.")
    else:
        raise ValueError(f"Unsupported plan kind: {spec.plan_kind}")

    lines.append("")
    return "\n".join(lines)


def _write_selected_runs_tsv(
    path: Path,
    spec: PublicDatasetSpec,
    rows: Sequence[PublicDatasetRow],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(SELECTED_RUN_COLUMNS)
        for row in rows:
            writer.writerow(
                [
                    spec.dataset_id,
                    row.run_id,
                    row.sample_id,
                    spec.protocol,
                    spec.organism,
                    spec.expected_read_layout,
                    spec.expected_reference,
                    spec.recommended_route,
                    row.source,
                    row.fastq_1,
                    row.fastq_2,
                    row.notes,
                ]
            )


def _write_readme(
    path: Path,
    spec: PublicDatasetSpec,
    rows: Sequence[PublicDatasetRow],
    download_method: str,
    dry_run: bool,
    row_mode: str,
    curated_path: Path | None,
    warnings: Sequence[str],
) -> None:
    body = [
        f"# Next steps for {spec.dataset_id}",
        "",
        f"- Title: {spec.title}",
        f"- Protocol: {spec.protocol}",
        f"- Organism: {spec.organism}",
        f"- Expected read layout: {spec.expected_read_layout}",
        f"- Expected reference: {spec.expected_reference}",
        f"- Recommended route: {spec.recommended_route}",
        f"- Plan kind: {spec.plan_kind}",
        f"- Row mode: {row_mode}",
        f"- Selected rows: {len(rows)}",
        f"- Download method: {download_method}",
        f"- Dry run: {'yes' if dry_run else 'no'}",
        "",
        "Files created by this planner:",
        "- selected_runs.tsv",
        "- download_plan.sh",
        "- README_next_steps.md",
        "",
        spec.notes,
        "",
    ]

    if warnings:
        body.extend(
            [
                "## Warning",
                "",
            ]
        )
        for warning in warnings:
            body.append(f"- {warning}")
        body.append("")

    if curated_path is not None:
        body.append(f"- Curated table path: {curated_path}")
    if row_mode == "placeholder":
        body.extend(
            [
                "",
                "Curated table not found. Embedded placeholder rows were used instead.",
                "",
            ]
        )

    if spec.dataset_id == "shin-ramda-riken":
        body.extend(
            [
                "Reference workflow:",
                "- https://github.com/rikenbit/shin-ramda-seq-paper",
                "",
                "This dataset alias is documentation-only for now. Expand accession-level metadata later if needed.",
                "",
            ]
        )
    elif spec.dataset_id == "E-MTAB-8735":
        body.extend(
            [
                "This Smart-seq3 fixture is intentionally tiny and placeholder-based.",
                "Full accession expansion for the 192-cell benchmark can be added later without changing the tests.",
                "",
            ]
        )
    else:
        body.extend(
            [
                "The selected rows are planning records only.",
                "Replace or extend the curated TSV with accession-expanded metadata later if a real download workflow is needed.",
                "",
            ]
        )

    path.write_text("\n".join(body), encoding="utf-8")


def prepare_public_dataset(
    *,
    dataset_id: str,
    outdir: Path,
    max_runs: int | None,
    dry_run: bool,
    protocol: str,
    download_method: str,
) -> dict[str, object]:
    spec = get_public_dataset_spec(dataset_id)
    _validate_protocol(spec, protocol)
    resolved_rows, row_mode, curated_path = _resolve_rows(spec)
    rows = _selected_rows(resolved_rows, max_runs)
    warnings = _warning_lines(spec, rows)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    selected_runs_path = outdir / "selected_runs.tsv"
    download_plan_path = outdir / "download_plan.sh"
    readme_path = outdir / "README_next_steps.md"

    _write_selected_runs_tsv(selected_runs_path, spec, rows)
    download_plan_path.write_text(
        _render_download_plan(spec, rows, outdir, download_method, row_mode, warnings),
        encoding="utf-8",
    )
    _write_readme(readme_path, spec, rows, download_method, dry_run, row_mode, curated_path, warnings)

    return {
        "dataset_id": spec.dataset_id,
        "protocol": spec.protocol,
        "organism": spec.organism,
        "expected_read_layout": spec.expected_read_layout,
        "expected_reference": spec.expected_reference,
        "recommended_route": spec.recommended_route,
        "download_method": download_method,
        "dry_run": dry_run,
        "row_mode": row_mode,
        "curated_table_path": curated_path,
        "selected_run_count": len(rows),
        "warnings": warnings,
        "selected_runs_path": selected_runs_path,
        "download_plan_path": download_plan_path,
        "readme_path": readme_path,
    }
