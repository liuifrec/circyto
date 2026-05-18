from __future__ import annotations

from pathlib import Path


def test_scomatic_study_design_doc_mentions_conservative_variant_terminology() -> None:
    path = Path("docs/scomatic_circrna_study_design.md")
    assert path.exists(), f"Missing design doc: {path}"
    text = path.read_text(encoding="utf-8")
    assert "RNA-derived candidate somatic variants" in text
    assert "candidate variant signals" in text
    assert "not yet sufficient" in text


def test_scomatic_study_design_doc_lists_proposed_output_tables() -> None:
    text = Path("docs/scomatic_circrna_study_design.md").read_text(encoding="utf-8")
    for name in [
        "circ_snv_cell_summary.tsv",
        "circ_snv_host_gene_summary.tsv",
        "circ_snv_candidate_variant_summary.tsv",
        "circ_snv_recurrent_circRNA_summary.tsv",
        "circ_snv_provenance.json",
    ]:
        assert name in text, f"Expected output table not documented: {name}"
