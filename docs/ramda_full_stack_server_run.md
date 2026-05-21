# RamDA Full Stack Server Run

This page describes the staged server-side plan for the current human scRR / RamDA-like datasets.

## Current server state

### IMR90 single-end pilot

Directory:

- `/user/ifrec/liuyuchen/circyto_redo/scrr_imr90`

Current files:

- `raw/SRR30918117.fastq.gz`
- `raw/SRR30918126.fastq.gz`
- `manifest.tsv`
- `work_hg38`

These correspond to the current completed IMR90 pilot:

- `SRR30918117`
- `SRR30918126`

Current validated public IMR90 RNA-side set in the repo metadata:

- `SRR30918117`
- `SRR30918126`

Practical implication:

- the current IMR90 pilot already spans the full validated public IMR90 RNA-side candidate set captured in `circyto`
- `manifest_all.tsv` for IMR90 is still useful as a stable full-dataset manifest, even if it currently contains the same two runs as the pilot manifest

### HAP1 paired-end pilot

Directory:

- `/user/ifrec/liuyuchen/circyto_redo/scrr_hap1`

Current files:

- `raw/SRR30911454_1.fastq.gz`
- `raw/SRR30911454_2.fastq.gz`
- `manifest.tsv`
- `work_hg38`

This corresponds to the current completed HAP1 pilot:

- `SRR30911454`

## Remaining HAP1 candidate runs

Candidate RNA-side paired-end HAP1 runs:

- `SRR30911454`
- `SRR30911453`
- `SRR30911559`

Only the latter two remain to be downloaded if a larger HAP1 run is desired:

- `SRR30911453`
- `SRR30911559`

## Remaining IMR90 candidate runs

Candidate RNA-side single-end IMR90 runs from the current repo metadata:

- `SRR30918117`
- `SRR30918126`

Additional remaining IMR90 downloads from the current validated set:

- none

## References

- genome FASTA:
  `/user/ifrec/liuyuchen/ref_clean/hg38/hg38.fa`
- GTF:
  `/user/ifrec/liuyuchen/ref_clean/hg38/gtf/gencode.v45.annotation.gtf`
- STAR index:
  `/user/ifrec/liuyuchen/ref_clean/hg38/star_index`

## Staged scripts

### IMR90 full stack

```bash
bash scripts/prepare_scrr_imr90_remaining_downloads.sh
bash scripts/run_scrr_imr90_full_stack.sh
```

Behavior:

- writes `manifest_all.tsv` for the complete validated IMR90 candidate set in the repo metadata
- uses `manifest_all.tsv` if present, else falls back to `manifest.tsv`
- expects single-end rows
- uses the validated BWA single-end route
- writes to:
  `/user/ifrec/liuyuchen/circyto_redo/scrr_imr90/work_hg38_fullstack`

### HAP1 remaining downloads and manifest expansion

```bash
bash scripts/prepare_scrr_hap1_remaining_downloads.sh
```

Behavior:

- checks free disk space before any download
- downloads `SRR30911453` only if missing
- downloads `SRR30911559` only if missing
- does not redownload `SRR30911454`
- preserves paired FASTQ naming:
  - `SRR30911453_1.fastq.gz`
  - `SRR30911453_2.fastq.gz`
  - `SRR30911559_1.fastq.gz`
  - `SRR30911559_2.fastq.gz`
- writes:
  - `/user/ifrec/liuyuchen/circyto_redo/scrr_hap1/manifest_all.tsv`

### HAP1 full stack

```bash
bash scripts/run_scrr_hap1_full_stack.sh
```

Behavior:

- uses `manifest_all.tsv` if present, else falls back to `manifest.tsv`
- expects paired-end rows
- uses:
  - `--star-index`
  - `--allow-paired-ramda`
- writes to:
  `/user/ifrec/liuyuchen/circyto_redo/scrr_hap1/work_hg38_fullstack`

## Full stack stages

Both full-stack scripts run:

1. `circyto workflow full-length-circrna`
2. `circyto refresh-rna-qc`
3. `circyto summarize-rna-circ --write-summary`
4. `circyto export-mudata --overwrite`
5. `circyto inspect-mudata --json`
6. `circyto summarize-mudata-qc --json`
7. `circyto cleanup-workflow --scope alignments --dry-run`

They also include an optional commented `scanpy-qc-report` line for later exploratory downstream work.

## Expected route by dataset

### IMR90

Expected route:

`single-end FASTQ -> BWA-MEM -> direct SAM -> CIRI3 -> matrix -> h5ad`

### HAP1

Expected route:

`paired-end FASTQ -> STAR -> BWA rescue -> CIRI3 STAR tuple mode -> matrix -> h5ad`

## Expected outputs

Under each `work_hg38_fullstack` directory:

- `align/`
- `ciri3/`
- `matrix/`
- `rna/`
- `qc/`
- `anndata/circ_counts.h5ad`
- `mudata/full_length.h5mu`
- `workflow_summary.json`

Additional outputs written by the wrappers:

- `qc/mudata_inspect.json`
- `qc/mudata_qc_summary.json`
- `qc/cleanup_alignments_dryrun.json`

Logs:

- `logs/run_scrr_imr90_full_stack_YYYYMMDD.log`
- `logs/prepare_scrr_imr90_remaining_downloads_YYYYMMDD.log`
- `logs/run_scrr_hap1_full_stack_YYYYMMDD.log`
- `logs/prepare_scrr_hap1_remaining_downloads_YYYYMMDD.log`

## Storage risks

Main risks:

- HAP1 paired-end alignments can generate very large STAR/BWA rescue intermediates
- even IMR90 single-end alignments can become much larger than input FASTQs
- `mudata/full_length.h5mu` is small compared with `align/`
- `qc/` and `rna/` summaries are small

Conservative recommendation:

- keep IMR90 as the lower-risk full-stack rerun
- treat HAP1 full-stack execution as the storage-heavy step
- IMR90 manifest preparation is low risk and mainly normalizes the full validated IMR90 set into `manifest_all.tsv`
- run HAP1 download preparation first, confirm free space, then execute the full stack

## Cleanup plan

Current wrapper behavior:

- run `circyto cleanup-workflow --scope alignments --dry-run`
- do not delete anything automatically

The real cleanup command is included only as commented shell lines inside the scripts. Review the dry-run JSON first and confirm that:

- `anndata/circ_counts.h5ad` exists
- `mudata/full_length.h5mu` exists
- `qc/` summaries exist
- `workflow_summary.json` is complete
