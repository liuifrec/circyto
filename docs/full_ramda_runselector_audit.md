# Full RamDA/scRR RunSelector Audit

This note defines the public RunSelector-based manifest generation strategy for the current human scRR / RamDA-like datasets used by `circyto`.

Audited datasets:

- IMR90:
  - `GSE278958`
  - `SRP537297`
  - `PRJNA1169833`
- HAP1:
  - `GSE278952`
  - `SRP537170`
  - GEO series text in repo notes also maps the HAP1 study to `PRJNA1169834`

## Current validated and expected RNA-side run structure

### IMR90

Validated RNA-side single-end SRRs captured in the repo metadata:

- `SRR30918126`
- `SRR30918117`

Estimated RNA-side cell count from the currently confirmed repo snapshot:

- 2 cell libraries

Expected layout:

- single-end

Expected `circyto` route:

- `single-end FASTQ -> BWA-MEM -> direct SAM -> CIRI3 -> matrix -> h5ad`

### HAP1

Validated RNA-side paired-end SRRs captured in the repo metadata:

- `SRR30911454`
- `SRR30911453`
- `SRR30911559`

Estimated RNA-side cell count from the currently confirmed repo snapshot:

- 3 cell libraries

Expected layout:

- paired-end

Expected `circyto` route:

- `paired-end FASTQ -> STAR -> BWA rescue -> CIRI3 STAR tuple mode -> matrix -> h5ad`

## Why demux is not required

For both public human datasets, the archived sequencing unit used by `circyto` is already the isolated per-cell library run:

- one SRR corresponds to one RNA-side single-cell library
- IMR90 libraries appear as one FASTQ object per SRR
- HAP1 libraries appear as one R1/R2 pair per SRR

That means the correct abstraction is:

- one manifest row per SRR-backed cell library

No pooled barcode demultiplexing step is required before the `full-length-circrna` workflow.

## Why Smart-seq3 differs structurally

The validated Smart-seq3 `E-MTAB-8735` workflow is structurally different because it starts from pooled FASTQs that contain many cells. That path requires:

- pooled FASTQ input
- explicit demultiplexing
- per-cell manifest generation
- then alignment and circRNA detection

The RamDA/scRR public datasets instead start directly from per-cell FASTQ objects, so the current `circyto` manifest abstraction is already correct.

## RunSelector builder scripts

Two scripts are provided:

- `scripts/build_imr90_rna_manifest_from_runinfo.sh`
- `scripts/build_hap1_rna_manifest_from_runinfo.sh`

Each script expects an NCBI RunSelector `RunInfo` CSV export and writes:

- a `circyto` manifest
- an RNA-side run inventory table

### IMR90

```bash
bash scripts/build_imr90_rna_manifest_from_runinfo.sh SraRunTable_imr90.csv outdir/
```

Outputs:

- `outdir/manifest_imr90_rna_all.tsv`
- `outdir/imr90_rna_run_inventory.tsv`

### HAP1

```bash
bash scripts/build_hap1_rna_manifest_from_runinfo.sh SraRunTable_hap1.csv outdir/
```

Outputs:

- `outdir/manifest_hap1_rna_all.tsv`
- `outdir/hap1_rna_run_inventory.tsv`

## Filtering strategy

The scripts are intentionally conservative, but they no longer hardcode only the pilot SRRs.

They:

- retain RunSelector rows that match the RNA-side library structure
- require:
  - `Assay Type = RNA-Seq`
  - `LibrarySource` containing `TRANSCRIPTOMIC`
  - `LibrarySelection = cDNA`
- exclude:
  - `Assay Type = OTHER`
  - `LibrarySource = GENOMIC`
  - obvious DNA / exome rows from metadata text
- enforce the dataset-expected layout:
  - IMR90 single-end
  - HAP1 paired-end
- extract metadata fields from the RunInfo CSV where available

This stays aligned with the currently validated public structure while allowing future manuscript-scale manifests to include all RNA-side rows exposed by the official RunSelector export.

## Manifest fields

Generated manifest columns:

- `sample_id`
- `fastq_1`
- `fastq_2`
- `protocol`
- `strandedness`
- `read_layout`
- `srr`
- `gsm`
- `condition`
- `cell_cycle_phase`

Current defaults:

- `protocol=ramda`
- `strandedness=unstranded`
- IMR90 `condition=IMR90_aphidicolin_treated`, `cell_cycle_phase=G1`
- HAP1 `condition=HAP1_scRR_MidS`, `cell_cycle_phase=mid-S`

## Inventory fields

The inventory tables keep a slightly wider public-metadata view for auditing:

- `sample_id`
- `srr`
- `gsm`
- `condition`
- `cell_cycle_phase`
- `read_layout`
- `library_strategy`
- `library_source`
- `organism`
- `sample_name`
- `title`
- `spots`
- `bases`

## Storage planning notes

Single-end IMR90 is materially cheaper than paired-end HAP1:

- fewer FASTQ objects
- smaller alignment intermediates
- simpler BWA direct-SAM route

Paired-end HAP1 has higher storage risk because:

- each run has two FASTQ objects
- STAR alignment intermediates can be large
- BWA rescue adds extra alignment artifacts

Practical recommendation:

- use IMR90 first for environment and reference validation
- stage HAP1 only after disk headroom is confirmed
- run `circyto cleanup-workflow --scope alignments --dry-run` after successful completion
- review cleanup JSON before any real deletion

## Remaining manual steps before manuscript-scale runs

1. Export the official RunSelector CSV for each study.
2. Run the corresponding manifest builder script.
3. Inspect the generated inventory table to confirm that only RNA-side SRRs were retained.
4. Download only the FASTQs referenced by the generated manifest.
5. Run `circyto workflow full-length-circrna` with the expected single-end or paired-end route.
6. Run the post-workflow stack:
   - `refresh-rna-qc`
   - `summarize-rna-circ --write-summary`
   - `export-mudata`
   - optional downstream QC

## Source basis

This audit combines:

- repo-tracked dataset notes in `docs/human_ramda_candidate_datasets.md`
- repo-tracked tabular run candidates in `testdata/public_datasets/human_ramda_candidate_runs.tsv`
- direct dataset-structure notes in `docs/ramda_dataset_structure.md`
- NCBI SRA Run Selector as the intended external metadata export source
