# RamDA-seq / Shin-RamDA-seq with CIRI3

`circyto run-ciri3` now supports protocol presets for:

- `--protocol ramda`
- `--protocol shin-ramda`
- `--protocol smartseq3`

For RamDA and Shin-RamDA, circyto treats the inputs as full-length total-RNA libraries:

- no UMI assumption
- no 3-prime counting assumption
- `read_layout` is taken from the manifest
- `strandedness` stays explicit metadata
- paired-end rows use STAR with chimeric output enabled for CIRI3
- single-end rows fall back to BWA-based CIRI3-compatible alignment prep

Minimal manifest columns accepted by `run-ciri3`:

```tsv
sample_id	fastq_1	fastq_2	protocol	strandedness	read_layout
```

`fastq_2` may be empty for single-end runs.

Example:

```bash
circyto run-ciri3 \
  --manifest testdata/ramda/manifest_chr21_smoke.tsv \
  --outdir runs/ramda_chr21_ciri3 \
  --genome-fasta /path/to/chr21.fa \
  --gtf /path/to/chr21.gtf \
  --star-index /path/to/star_chr21_index \
  --protocol ramda \
  --threads 8
```

Dry-run mode prints the planned aligner and CIRI3 commands:

```bash
circyto run-ciri3 \
  --manifest testdata/shin_ramda/manifest_chr21_smoke.tsv \
  --outdir runs/shin_ramda_chr21_ciri3 \
  --genome-fasta /path/to/chr21.fa \
  --gtf /path/to/chr21.gtf \
  --star-index /path/to/star_chr21_index \
  --protocol shin-ramda \
  --dry-run
```
