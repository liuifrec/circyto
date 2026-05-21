from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _run_builder(script_name: str, csv_path: Path, outdir: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", str(ROOT / "scripts" / script_name), str(csv_path), str(outdir)],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


def test_build_imr90_rna_manifest_from_runinfo_includes_only_rnaseq_rows(tmp_path: Path) -> None:
    csv_path = tmp_path / "SraRunTable_imr90.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Run",
                "Assay Type",
                "LibrarySource",
                "LibrarySelection",
                "LibraryLayout",
                "GEO_Accession",
                "Library Name",
                "BioSample",
                "Experiment",
                "Sample Name",
                "Title",
                "Treatment",
            ]
        )
        writer.writerow(["SRR1", "RNA-Seq", "TRANSCRIPTOMIC SINGLE CELL", "cDNA", "SINGLE", "GSM1", "RNA_IMR90_A_100", "SAMN1", "SRX1", "RNA IMR90 cell 1", "IMR90 aphidicolin treated G1", "aphidicolin"])
        writer.writerow(["SRR2", "RNA-Seq", "TRANSCRIPTOMIC SINGLE CELL", "cDNA", "SINGLE", "GSM2", "RNA_IMR90_A_101", "SAMN2", "SRX2", "RNA IMR90 cell 2", "IMR90 aphidicolin treated G1", "aphidicolin"])
        writer.writerow(["SRR3", "OTHER", "GENOMIC", "RANDOM", "PAIRED", "GSM3", "DNA_IMR90_1", "SAMN3", "SRX3", "DNA IMR90", "DNA library", ""])

    outdir = tmp_path / "out"
    result = _run_builder("build_imr90_rna_manifest_from_runinfo.sh", csv_path, outdir)

    manifest = pd.read_csv(outdir / "manifest_imr90_rna_all.tsv", sep="\t")
    inventory = pd.read_csv(outdir / "imr90_rna_run_inventory.tsv", sep="\t")
    assert list(manifest["srr"]) == ["SRR1", "SRR2"]
    assert manifest["fastq_1"].tolist() == ["raw/SRR1.fastq.gz", "raw/SRR2.fastq.gz"]
    assert manifest["fastq_2"].fillna("").tolist() == ["", ""]
    assert set(manifest["read_layout"]) == {"single"}
    assert "SRR3" not in set(inventory["srr"])
    assert "RNA-side rows retained: 2" in result.stdout
    assert "DNA/genomic/OTHER rows excluded: 1" in result.stdout


def test_build_hap1_rna_manifest_from_runinfo_excludes_genomic_and_uses_paired_naming(tmp_path: Path) -> None:
    csv_path = tmp_path / "SraRunTable_hap1.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Run",
                "Assay Type",
                "LibrarySource",
                "LibrarySelection",
                "LibraryLayout",
                "GEO_Accession",
                "Library Name",
                "BioSample",
                "Experiment",
                "Sample Name",
                "Title",
            ]
        )
        writer.writerow(["SRR10", "RNA-Seq", "TRANSCRIPTOMIC SINGLE CELL", "cDNA", "PAIRED", "GSM10", "RNA_HAP1_scRR1_MidS_103", "SAMN10", "SRX10", "RNA HAP1 MidS", "HAP1 mid-S RNA"])
        writer.writerow(["SRR11", "RNA-Seq", "TRANSCRIPTOMIC SINGLE CELL", "cDNA", "PAIRED", "GSM11", "RNA_HAP1_scRR1_MidS_104", "SAMN11", "SRX11", "RNA HAP1 MidS", "HAP1 mid-S RNA"])
        writer.writerow(["SRR12", "OTHER", "GENOMIC", "RANDOM", "PAIRED", "GSM12", "DNA_HAP1", "SAMN12", "SRX12", "DNA HAP1", "HAP1 DNA"])

    outdir = tmp_path / "out"
    result = _run_builder("build_hap1_rna_manifest_from_runinfo.sh", csv_path, outdir)

    manifest = pd.read_csv(outdir / "manifest_hap1_rna_all.tsv", sep="\t")
    inventory = pd.read_csv(outdir / "hap1_rna_run_inventory.tsv", sep="\t")
    assert list(manifest["srr"]) == ["SRR10", "SRR11"]
    assert manifest["fastq_1"].tolist() == ["raw/SRR10_1.fastq.gz", "raw/SRR11_1.fastq.gz"]
    assert manifest["fastq_2"].tolist() == ["raw/SRR10_2.fastq.gz", "raw/SRR11_2.fastq.gz"]
    assert set(manifest["read_layout"]) == {"paired"}
    assert "SRR12" not in set(inventory["srr"])
    assert "RNA-side rows retained: 2" in result.stdout
    assert "DNA/genomic/OTHER rows excluded: 1" in result.stdout


def test_builders_require_rnaseq_transcriptomic_cdna_rules(tmp_path: Path) -> None:
    csv_path = tmp_path / "SraRunTable_misc.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Run", "Assay Type", "LibrarySource", "LibrarySelection", "LibraryLayout", "Library Name"])
        writer.writerow(["SRR20", "RNA-Seq", "TRANSCRIPTOMIC SINGLE CELL", "polyA", "SINGLE", "bad_selection"])
        writer.writerow(["SRR21", "RNA-Seq", "GENOMIC", "cDNA", "SINGLE", "bad_source"])
        writer.writerow(["SRR22", "OTHER", "TRANSCRIPTOMIC SINGLE CELL", "cDNA", "SINGLE", "bad_assay"])

    outdir = tmp_path / "out"
    proc = subprocess.run(
        ["bash", str(ROOT / "scripts" / "build_imr90_rna_manifest_from_runinfo.sh"), str(csv_path), str(outdir)],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    assert proc.returncode != 0
    assert "no IMR90 RNA-side rows were found" in proc.stderr
