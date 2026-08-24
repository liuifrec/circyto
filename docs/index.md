# circyto Documentation

Welcome to the `circyto` documentation.

`circyto` is a CLI/scverse-compatible framework for single-cell circular RNA detection, annotation, and multimodal integration from full-length single-cell RNA-seq and full-length single-cell multi-omic data.

## Current Capabilities

| Modality | Status |
| --- | --- |
| `rna` | validated RNA profile / MuData modality for completed full-length workflows |
| `circ` | validated circRNA matrix, QC, `h5ad`, and MuData modality |
| `cnv` | validated processed scRR GEO CNV import and IMR90 full23 tri-modal MuData |
| `rt` | implemented processed scRR replication timing/state import for HAP1-style DNA tables; real-file validation pending |
| `candidate_snv` | exploratory optional SComatic interoperability via `merge-scomatic`; RNA-derived candidate variant signals only |

| Dataset / run | Status |
| --- | --- |
| IMR90 full23 | RNA + circ + CNV tri-modal MuData validated |
| HAP1 batch10 | RNA + circ workflow plus SComatic technical smoke validated |
| HAP1 processed RT | `rt` importer and RNA+circ+RT merge implemented with synthetic tests; real GSE278952 processed-file import pending |
| HAP1 full | pending full FASTQ download and full workflow run |

## Sections

- [Getting started](getting_started.md)
- [Validated workflows summary](validated_workflows_summary.md)
- [Release notes](release_notes.md)
- [Full-length circRNA workflow](full_length_workflow.md)
- [Gene expression and velocity integration](gene_expression_velocity_integration.md)
- [RNA and velocity contract](rna_velocity_contract.md)
- [Post-hoc RNA profile](posthoc_rna_profile.md)
- [RNA circ integration](rna_circ_integration.md)
- [MuData export](mudata_export.md)
- [MuData schema](mudata_schema.md)
- [MuData downstream](mudata_downstream.md)
- [Modality schema](modality_schema.md)
- [circyto-biogenesis schema foundation](biogenesis_schema.md)
- [Host-gene provenance](host_gene_provenance.md)
- [Manuscript reproducibility](manuscript_reproducibility.md)
- [Manuscript command reference](manuscript_command_reference.md)
- [scRR multimodal roadmap](scrr_multimodal_roadmap.md)
- [scRR DNA architecture review](scrr_dna_architecture_review.md)
- [scRR CNV modality design](scrr_cnv_modality_design.md)
- [scRR cell pairing strategy](scrr_cell_pairing_strategy.md)
- [scRR cell mapping](scrr_cell_mapping.md)
- [scRR tri-modal MuData](scrr_trimodal_mudata.md)
- [scRR replication timing modality](scrr_replication_timing_modality.md)
- [scRR DNA roadmap](scrr_dna_roadmap.md)
- [DNA/SNV integration scaffold](dna_snv_integration.md)
- [Manuscript benchmark plan](manuscript_benchmark_plan.md)
- [Current project status](current_project_status.md)
- [Manuscript figure skeleton](manuscript_figure_skeleton.md)
- [Results skeleton](results_skeleton.md)
- [scverse interoperability](scverse_interoperability.md)
- [Reproducibility checklist](reproducibility_checklist.md)
- [Scanpy downstream](scanpy_downstream.md)
- [Intermediate cleanup policy](intermediate_cleanup_policy.md)
- [Workflow integrity and provenance](workflow_integrity_and_provenance.md)
- [RamDA / Shin-RamDA protocol notes](protocols_ramda_shin_ramda.md)
- [Server install and RamDA hg38 run](server_install_and_ramda_hg38.md)
- [Detectors](detectors.md)
- [Optional CIRI-long adapter](ciri_long_adapter.md)
- [Nanopore feasibility and two-path architecture](nanopore_feasibility.md)
- [Public dataset preparation](public_dataset_prepare.md)
- [Human RamDA candidate datasets](human_ramda_candidate_datasets.md)
- [Human scRR server test plan](human_scrr_server_test_plan.md)
- [SComatic interoperability](scomatic_interop.md)
- [Full-length RNA SComatic adapter](scomatic_full_length_adapter.md)
- [SComatic-integrated circRNA study design](scomatic_circrna_study_design.md)
- [Local SComatic chr21 POC](local_scomatic_chr21_poc.md)
- [Local SComatic environment setup](scomatic_environment_setup.md)
- [CIRI-full chr21 example](ciri_full_chr21_example.md)
- SMART-Seq3 workflow: see [Getting started](getting_started.md#workflow-1b-experimental-smart-seq3-to-ciri3)

Use the links above to navigate through the documentation.

