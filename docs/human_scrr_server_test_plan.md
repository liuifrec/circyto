# Human scRR Server Test Plan

This is the first real human `hg38` scRR RNA pilot plan for `circyto`.

Constraints assumed:

- run sequentially
- only 2 cells initially
- RNA-side only
- no DNA-side runs
- avoid unnecessary downloads
- current free disk is about `290 GB`

## Recommended first pilot

Use these two IMR-90 RNA runs from `GSE278958`:

- `GSM8558852 / SRX26321174 / SRR30918126 / RNA_IMR90_A_100`
- `GSM8558853 / SRX26321183 / SRR30918117 / RNA_IMR90_A_101`

Why these first:

- public and downloadable now
- confirmed single-end in SRA XML
- directly match the validated single-end RamDA/Shin-RamDA route:
  `BWA-MEM -> direct SAM -> CIRI3`

Do not use HAP1 for the first executable pilot.

Reason:

- `GSE278952` RNA-side HAP1 runs are paired-end
- the paired-end route now has local chr21 subset validation but still needs `--allow-paired-ramda` for explicit opt-in
- IMR-90 remains the safer first server check because it stays on the already validated single-end BWA route

## Recommended references

Match the human metadata with:

- `hg38` primary assembly FASTA
- human GENCODE v38 GTF

Example placeholders:

```bash
export HG38_FA=/refs/hg38/GRCh38.primary_assembly.genome.fa
export GTF=/refs/hg38/gencode.v38.annotation.gtf
```

## Server layout

```bash
mkdir -p work/human_scrr_imr90_2cell/{sra,fastq,qc,manifests,run,tmp}
```

## Prefetch

Download only the two selected RNA runs:

```bash
prefetch \
  --output-directory work/human_scrr_imr90_2cell/sra \
  SRR30918126 \
  SRR30918117
```

If your SRA toolkit requires explicit size allowance:

```bash
prefetch \
  --max-size 50G \
  --output-directory work/human_scrr_imr90_2cell/sra \
  SRR30918126 \
  SRR30918117
```

## fasterq-dump

Run sequentially and keep temp files off the main work area:

```bash
mkdir -p work/human_scrr_imr90_2cell/tmp/fasterq

fasterq-dump \
  --threads 8 \
  --split-files \
  --temp work/human_scrr_imr90_2cell/tmp/fasterq \
  --outdir work/human_scrr_imr90_2cell/fastq \
  work/human_scrr_imr90_2cell/sra/SRR30918126/SRR30918126.sra

fasterq-dump \
  --threads 8 \
  --split-files \
  --temp work/human_scrr_imr90_2cell/tmp/fasterq \
  --outdir work/human_scrr_imr90_2cell/fastq \
  work/human_scrr_imr90_2cell/sra/SRR30918117/SRR30918117.sra
```

Expected output shape for these runs:

- `work/human_scrr_imr90_2cell/fastq/RNA_IMR90_A_100_1.fastq`
- `work/human_scrr_imr90_2cell/fastq/RNA_IMR90_A_101_1.fastq`

Even though `--split-files` is used, these are single-end runs, so only `_1.fastq` is expected.

## seqkit QC

Basic stats:

```bash
seqkit stats -a \
  work/human_scrr_imr90_2cell/fastq/RNA_IMR90_A_100_1.fastq \
  work/human_scrr_imr90_2cell/fastq/RNA_IMR90_A_101_1.fastq \
  > work/human_scrr_imr90_2cell/qc/seqkit_stats.tsv
```

Quick content check:

```bash
seqkit fqchk \
  -q 20 \
  -n 100000 \
  work/human_scrr_imr90_2cell/fastq/RNA_IMR90_A_100_1.fastq \
  > work/human_scrr_imr90_2cell/qc/RNA_IMR90_A_100.fqchk.txt

seqkit fqchk \
  -q 20 \
  -n 100000 \
  work/human_scrr_imr90_2cell/fastq/RNA_IMR90_A_101_1.fastq \
  > work/human_scrr_imr90_2cell/qc/RNA_IMR90_A_101.fqchk.txt
```

## Example manifest creation

These rows should stay explicitly single-end:

```bash
cat > work/human_scrr_imr90_2cell/manifests/imr90_2cell.tsv <<'EOF'
sample_id	fastq_1	fastq_2	protocol	strandedness	read_layout
RNA_IMR90_A_100	work/human_scrr_imr90_2cell/fastq/RNA_IMR90_A_100_1.fastq		ramda	unstranded	single
RNA_IMR90_A_101	work/human_scrr_imr90_2cell/fastq/RNA_IMR90_A_101_1.fastq		ramda	unstranded	single
EOF
```

## Recommended workflow command

```bash
circyto workflow full-length-circrna \
  --manifest work/human_scrr_imr90_2cell/manifests/imr90_2cell.tsv \
  --outdir work/human_scrr_imr90_2cell/run/full_length_ciri3_imr90_2cell \
  --protocol ramda \
  --genome-fasta "$HG38_FA" \
  --gtf "$GTF" \
  --detector ciri3 \
  --threads 8 \
  --export-h5ad
```

Expected route:

```text
single-end FASTQ -> BWA-MEM -> direct SAM -> CIRI3 -> matrix + h5ad
```

## Optional preflight dry-run

```bash
circyto workflow full-length-circrna \
  --manifest work/human_scrr_imr90_2cell/manifests/imr90_2cell.tsv \
  --outdir work/human_scrr_imr90_2cell/run/full_length_ciri3_imr90_2cell_dryrun \
  --protocol ramda \
  --genome-fasta "$HG38_FA" \
  --gtf "$GTF" \
  --detector ciri3 \
  --threads 8 \
  --dry-run
```

## Storage cautions

- `fasterq-dump` can temporarily use much more space than the final FASTQ size
- keep `--temp` on a filesystem with comfortable headroom
- do not prefetch extra SRRs until the 2-cell pilot succeeds
- delete temporary `fasterq` content after successful conversion
- expect SAM outputs under `align/` to be much larger than the input FASTQs

Conservative rule:

- never download more than these two IMR-90 SRRs in the first pass

## HAP1 follow-up, not first execution

Confirmed HAP1 RNA-side paired-end candidates:

- `SRR30911454` from `GSM8558630 / SRX26315002`
- `SRR30911453` from `GSM8558631 / SRX26315003`
- `SRR30911559` from `GSM8558632 / SRX26315028`

Use these only for:

- accession confirmation
- storage planning
- paired-end hg38 follow-up after the IMR-90 2-cell pilot is clean

Example HAP1 paired-end dry-run:

```bash
circyto workflow full-length-circrna \
  --manifest work/human_scrr_hap1_2cell/manifests/hap1_2cell.tsv \
  --outdir work/human_scrr_hap1_2cell/run/full_length_ciri3_hap1_2cell_dryrun \
  --protocol ramda \
  --genome-fasta "$HG38_FA" \
  --gtf "$GTF" \
  --star-index /refs/hg38/star_gencode_v38 \
  --dry-run
```

Example HAP1 paired-end execution:

```bash
circyto workflow full-length-circrna \
  --manifest work/human_scrr_hap1_2cell/manifests/hap1_2cell.tsv \
  --outdir work/human_scrr_hap1_2cell/run/full_length_ciri3_hap1_2cell \
  --protocol ramda \
  --genome-fasta "$HG38_FA" \
  --gtf "$GTF" \
  --star-index /refs/hg38/star_gencode_v38 \
  --allow-paired-ramda \
  --export-h5ad
```

Expected HAP1 route:

```text
paired-end FASTQ -> STAR alignment-first cache -> CIRI3 STAR tuple mode -> matrix + h5ad
```

Local validation note:

- `SRR30911454` completed a real `chr21` subset run in `circyto` with STAR alignment, BWA rescue, CIRI3, matrix export, and `h5ad` export
- the validated subset run used `10k` read pairs because the initial `100k` subset was slower than desirable on the reduced `chr21` reference
- the validated subset produced zero circRNA rows, so the route is executable but not yet biologically benchmarked

## Post-run inspection

For any completed `full-length-circrna` run, use the read-only inspector:

IMR90 single-end BWA route:

```bash
python scripts/check_full_length_workflow_outputs.py \
  work/human_scrr_imr90_2cell/run/full_length_ciri3_imr90_2cell
```

HAP1 paired-end STAR route:

```bash
python scripts/check_full_length_workflow_outputs.py \
  work/human_scrr_hap1_2cell/run/full_length_ciri3_hap1_2cell
```

JSON output for machine-friendly inspection:

```bash
python scripts/check_full_length_workflow_outputs.py \
  work/human_scrr_hap1_2cell/run/full_length_ciri3_hap1_2cell \
  --json
```

The inspector is read-only and reports:

- stage graph
- planned / completed / failed cells
- circRNAs per cell
- matrix dimensions
- top recurrent circRNAs
- `h5ad` existence
- largest files
- obvious error lines in logs
