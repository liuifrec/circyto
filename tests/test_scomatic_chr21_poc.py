from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.scomatic_chr21_poc import normalize_scomatic_candidate_table


runner = CliRunner()
ROOT = Path(__file__).resolve().parents[1]


def _write_base_workdir(tmp_path: Path) -> Path:
    root = tmp_path / "workflow"
    (root / "rna").mkdir(parents=True)
    (root / "matrix").mkdir()
    (root / "qc").mkdir()
    (root / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\tcellA\tcellB\nG1\tGENE1\t5\t2\n",
        encoding="utf-8",
    )
    (root / "rna" / "gene_feature_table.tsv").write_text("gene_id\tgene_name\nG1\tGENE1\n", encoding="utf-8")
    (root / "qc" / "rna_qc.tsv").write_text(
        "cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\ncellA\t5\t1\t0.0\t0.0\t2\ncellB\t2\t1\t0.0\t0.0\t0\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 1\n1 1 2\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_index.txt").write_text("circ1\ncirc2\n", encoding="utf-8")
    (root / "matrix" / "cell_index.txt").write_text("cellA\n", encoding="utf-8")
    return root


def test_normalize_scomatic_candidate_table_schema() -> None:
    raw = pd.DataFrame(
        [
            {
                "cell_id": "cellA",
                "chrom": "chr21",
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "gene": "GENE1",
                "read_support": 4,
                "vaf": 0.2,
            }
        ]
    )
    normalized = normalize_scomatic_candidate_table(raw)
    assert list(normalized.columns) == [
        "variant_id",
        "cell_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene",
        "filter_status",
        "candidate_variant_class",
        "read_support",
        "vaf",
        "caller",
    ]
    assert normalized.loc[0, "candidate_variant_class"] == "RNA-derived candidate variant signal"


def test_run_scomatic_chr21_poc_synthetic_and_import_join(tmp_path: Path) -> None:
    workdir = _write_base_workdir(tmp_path)
    ref = tmp_path / "chr21.fa"
    gtf = tmp_path / "chr21.gtf"
    ref.write_text(">chr21\nACGT\n", encoding="utf-8")
    gtf.write_text('chr21\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "G1"; gene_name "GENE1";\n', encoding="utf-8")
    outdir = tmp_path / "scomatic_poc"

    proc = subprocess.run(
        [
            "bash",
            str(ROOT / "scripts" / "run_scomatic_chr21_poc.sh"),
            "--workdir",
            str(workdir),
            "--reference",
            str(ref),
            "--gtf",
            str(gtf),
            "--outdir",
            str(outdir),
            "--synthetic",
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    payload = json.loads(proc.stdout)
    assert payload["n_candidates"] == 1
    candidate_path = outdir / "scomatic_candidate_summary.tsv"
    assert candidate_path.exists()

    dna_cell = tmp_path / "dna_cell_summary.tsv"
    dna_cell.write_text(
        "cell_id\tdna_library_id\tcnv_burden\treplication_score\tcell_cycle_phase\tdna_variant_count\tnotes\n"
        "cellA\tlibA\t1.5\t0.7\tS\t3\tok\n",
        encoding="utf-8",
    )
    import_result = runner.invoke(
        app,
        [
            "import-dna-snv-summary",
            "--workdir",
            str(workdir),
            "--dna-cell-summary",
            str(dna_cell),
            "--scomatic-candidate-summary",
            str(candidate_path),
        ],
    )
    assert import_result.exit_code == 0, import_result.output
    join_result = runner.invoke(app, ["summarize-dna-rna-circ", "--workdir", str(workdir), "--write-summary"])
    assert join_result.exit_code == 0, join_result.output
    join_payload = json.loads(join_result.output)
    assert join_payload["n_shared_cells"] == 1
    table_text = (workdir / "qc" / "dna_rna_circ_cell_summary.tsv").read_text(encoding="utf-8")
    assert "scomatic_candidate_count" in table_text


def test_local_scomatic_doc_keeps_candidate_signal_terminology() -> None:
    text = Path("docs/local_scomatic_chr21_poc.md").read_text(encoding="utf-8")
    assert "RNA-derived candidate variant signals" in text
