from pathlib import Path

import pandas as pd

from circyto.pipeline.annotate_host_gene import annotate_host_genes


def test_annotate_host_genes_simple(tmp_path: Path):
    # --- 1) Fake circ_feature_table.tsv ---
    circ_tsv = tmp_path / "circ_feature_table.tsv"
    df_circ = pd.DataFrame(
        [
            {
                "circ_id": "chr1:100|200|+",
                "chrom": "chr1",
                "start": 120,
                "end": 180,
                "strand": "+",
            },
            {
                "circ_id": "chr1:500|600|-",
                "chrom": "chr1",
                "start": 510,
                "end": 590,
                "strand": "-",
            },
        ]
    )
    df_circ.to_csv(circ_tsv, sep="\t", index=False)

    # --- 2) Fake GTF with two genes ---
    gtf_path = tmp_path / "genes.gtf"
    gtf_path.write_text(
        "\n".join(
            [
                'chr1\tsource\tgene\t100\t300\t.\t+\t.\tgene_id "GID1"; gene_name "GENE1";',
                'chr1\tsource\tgene\t500\t700\t.\t-\t.\tgene_id "GID2"; gene_name "GENE2";',
                "",
            ]
        )
    )

    # --- 3) Run annotation (overwrite in place) ---
    annotate_host_genes(
        circ_feature_table=circ_tsv,
        gtf_path=gtf_path,
        out=None,
        max_genes_per_circ=5,
    )

    # --- 4) Check results ---
    df_out = pd.read_csv(circ_tsv, sep="\t")

    assert "host_gene" in df_out.columns
    assert "host_gene_id" in df_out.columns
    assert "host_genes_multi" in df_out.columns
    assert "host_gene_ids_multi" in df_out.columns
    assert "host_gene_n" in df_out.columns

    row1 = df_out[df_out["circ_id"] == "chr1:100|200|+"].iloc[0]
    row2 = df_out[df_out["circ_id"] == "chr1:500|600|-"].iloc[0]

    assert row1["host_gene"] == "GENE1"
    assert row1["host_gene_id"] == "GID1"
    assert row1["host_gene_n"] >= 1

    assert row2["host_gene"] == "GENE2"
    assert row2["host_gene_id"] == "GID2"
    assert row2["host_gene_n"] >= 1
