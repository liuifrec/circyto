from __future__ import annotations

from pathlib import Path


MANUSCRIPT_OBJECTS: dict[str, dict[str, object]] = {
    "Smart-seq3": {
        "accession": "E-MTAB-8735",
        "default_path": Path("load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu"),
        "expected_shapes": {"rna": [192, 63187], "circ": [192, 2503]},
        "expected_shared_cells": 192,
        "expected_host_gene_annotated": 2379,
        "expected_host_gene_total": 2503,
        "expected_median_circRNA_count": 12.0,
        "expected_median_circRNA_total_support": 22.5,
        "source_data_type": "Public Smart-seq3 RNA sequencing data with CIRI3 circRNA calls and RNA counts.",
        "intended_role": "Manuscript-scale RNA+circRNA reproducibility object.",
        "manuscript_release_status": "Candidate frozen manuscript object; prefer release/Zenodo archive over ordinary Git.",
    },
    "HAP1": {
        "accession": "GSE278952",
        "default_path": Path("load_work/scrr_hap1/full_length_rna_circ_rt.hostgene_fixed.h5mu"),
        "expected_shapes": {"rna": [63, 63187], "circ": [63, 3209], "rt": [56, 56881]},
        "expected_shared_cells": 56,
        "expected_host_gene_annotated": 3117,
        "expected_host_gene_total": 3209,
        "source_data_type": "Public HAP1 scRR RNA data plus processed public RT/state summaries.",
        "intended_role": "Manuscript-scale RNA+circRNA+RT reproducibility object.",
        "manuscript_release_status": "Candidate frozen manuscript object; prefer release/Zenodo archive over ordinary Git.",
    },
    "IMR90": {
        "accession": "GSE278958",
        "default_path": Path("load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu"),
        "expected_shapes": {"rna": [23, 63187], "circ": [23, 2443], "cnv": [23, 60607]},
        "expected_shared_cells": 23,
        "expected_host_gene_annotated": 2429,
        "expected_host_gene_total": 2443,
        "source_data_type": "Public IMR90 scRR RNA data plus processed public CNV summary tables.",
        "intended_role": "Manuscript-scale RNA+circRNA+CNV reproducibility object.",
        "manuscript_release_status": "Candidate frozen manuscript object; prefer release/Zenodo archive over ordinary Git.",
    },
}


EXPECTED_HOST_GENE_FIELDS = [
    "host_gene",
    "host_gene_source",
    "host_gene_from_gtf",
    "host_gene_from_circatlas",
    "host_gene_from_circatlas_id",
    "host_gene_id",
    "host_genes_multi",
    "host_gene_ids_multi",
    "host_gene_n",
]
