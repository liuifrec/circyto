# Validated Workflows Summary

This page summarizes the current `circyto` workflow and utility states in a compact benchmark-oriented format.

The intent is to separate:

- validated workflows
- lightweight QC and post-processing utilities
- experimental or future directions

## Detector Status

| Detector / adapter | Status | Evidence | Caveats |
| --- | --- | --- | --- |
| CIRI3 | validated primary detector | Primary detector for the validated Smart-seq3, single-end RamDA/scRR, and paired-end scRR/RamDA-style manuscript workflows | External Java/CIRI3 runtime is still reported by `circyto detectors`; STAR mode needs STAR inputs |
| CIRI-full | experimental adapter | Parser/unit coverage plus gated chr21 Smart-seq2 integration example (`tests/test_cirifull_chr21_integration.py`, `docs/ciri_full_chr21_example.md`) | Not manuscript-scale validated; single-end `ciri-full` uses a CIRI2 fallback |
| CIRI2 | experimental adapter | Parser/unit coverage plus gated single-end chr21 known-positive regression (`tests/test_ciri2_se_chr21_integration.py`) | Not manuscript-scale validated |
| find-circ3 | experimental adapter | Fixture parser/collector tests and runtime readiness reporting | Not manuscript-scale validated |
| CIRCexplorer2 | optional experimental adapter | Optional wrapper plus fixture collector tests | Not manuscript-scale validated |

No concrete circSC dataset validation artifact for CIRI2 or CIRI-full was found in the current repository, tests, docs, committed summaries, or branch-history search. TODO: if a reproducible circSC validation is recovered, document the dataset, exact commands, outputs, and caveats before upgrading either adapter status.

## Validated

| Dataset | Protocol | Read layout | Route | Reference | Output type | Validation status | Caveats |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `E-MTAB-8735` diySpike / pooled Smart-seq3 | `smartseq3` | paired-end transcript reads plus index reads | pooled FASTQs -> demux -> manifest selection -> STAR -> BWA rescue -> CIRI3 STAR tuple mode -> matrix -> QC -> `h5ad` | `hg38` | circRNA matrix, QC, `anndata/circ_counts.h5ad` | validated end-to-end on real pooled Smart-seq3 workflow input | still exposed under experimental CLI naming; MuData remains optional and not the core validated output |
| `GSE278958` IMR90 human scRR / RamDA-like | `ramda` / scRR-style full-length | single-end | FASTQ -> BWA-MEM -> direct SAM -> CIRI3 -> matrix -> QC -> `h5ad` | `hg38` | circRNA matrix, QC, `anndata/circ_counts.h5ad` | validated on real human single-end full-length scRR workflow | biological validation target is human RNA-side scRR; this does not imply velocity or orthogonally confirmed DNA variant support |
| `GSE278958` IMR90 human scRR / RamDA-like | `ramda` / scRR-style full-length + processed scRR DNA | single-end RNA plus processed GEO CNV summaries | RNA/circ MuData -> GSM remap -> CNV AnnData merge | `hg38` / GEO processed CNV bins | tri-modal RNA+circ+CNV MuData | validated on 23 cells with trimodal overlap 23 | CNV is imported from processed GEO summary tables, not raw DNA FASTQ reprocessing |
| `GSE278952` HAP1 human scRR | `ramda` / scRR-style full-length | paired-end | FASTQ pair -> STAR -> BWA rescue -> CIRI3 STAR tuple mode -> matrix -> QC -> `h5ad` | `hg38` | circRNA matrix, QC, `anndata/circ_counts.h5ad` | validated as executable on real human paired-end scRR workflow path | paired-end RamDA route still uses explicit `--allow-paired-ramda`; biological benchmarking remains narrower than single-end RamDA history |
| `GSE278952` HAP1 human scRR + processed scRR DNA RT/state | `ramda` / scRR-style full-length + processed RT/state | paired-end RNA plus processed GEO RT/state summaries | RNA/circ MuData -> GSM remap -> RT AnnData merge -> manuscript RT/circ summaries | `hg38` / GEO processed RT/state features | tri-modal RNA+circ+RT MuData | validated on committed manuscript summaries with RNA/circ/RT overlap 56 | RT/state is imported from processed GEO-style tables, not raw DNA FASTQ reprocessing |

## Lightweight QC And Post-processing

| Utility | Input | Route | Output type | Validation status | Caveats |
| --- | --- | --- | --- | --- | --- |
| `circyto add-rna-profile --method simple-overlap` | completed workflow directory with alignment manifest | reuse existing SAM/BAM alignments -> simple overlap counting against GTF gene or exon-derived intervals | `rna/gene_counts.tsv`, `rna/gene_feature_table.tsv`, `rna/rna_import_summary.json`, optional `workflow_summary.json` update | validated as a lightweight post-hoc RNA sanity profile | not a replacement for `featureCounts`; no velocity layers; not merged into the circ-only `h5ad` |
| cleanup / provenance / integrity utilities | completed workflow directory | read-only inspection or scoped cleanup over workflow-owned intermediates | summaries, integrity reports, optional post-success cleanup metadata | validated as lightweight workflow support utilities | cleanup is explicit opt-in only; never applies to user inputs |
| `scripts/check_full_length_workflow_outputs.py`, `scripts/check_rna_profile_outputs.py`, `circyto check-workflow` | completed workflow directory | read-only inspection and validation | text or JSON summaries | validated as read-only QC and integrity helpers | intended for inspection, not biological interpretation |

## Experimental Or Future

| Area | Current state | Notes |
| --- | --- | --- |
| MuData multimodal export for full-length workflows | validated for RNA+circ and IMR90 RNA+circ+CNV | `export-mudata` writes RNA+circ; `merge-scrr-cnv` writes tri-modal scRR MuData after GSM remapping |
| HAP1 replication timing/state import | validated for processed HAP1 RT/state summaries | `import-scrr-rt` writes `rt.h5ad`; `merge-scrr-rt` writes RNA+circ+RT MuData; committed manuscript result tables summarize 56 shared RNA/circ/RT cells |
| RNA velocity-compatible layers | future | no `velocyto` dependency required today; spliced / unspliced / ambiguous layers remain planned |
| SComatic interoperability | exploratory | use conservative language such as `RNA-derived candidate variant signals`; no orthogonally confirmed DNA variant claims |
| automatic integrated circRNA + RNA + velocity pipeline | future | current implementation favors explicit contracts and staged utilities rather than a monolithic production workflow |

## Reading This Page

- `validated` means the route has been executed successfully on the stated real dataset class and is part of the current supported biological direction.
- `lightweight QC` means the utility is intentionally narrow, read-only, or post-hoc, and is meant to improve inspectability or reproducibility rather than define a new biological workflow class.
- `experimental or future` means the design may exist, but it should not yet be presented as a production or benchmarked workflow.
