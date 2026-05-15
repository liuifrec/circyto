# Server Install and RamDA hg38 Run

This note is the minimal server-side path for installing `circyto` and running a first full hg38 RamDA-seq job.

It is intentionally conservative:

- do a 2-cell first run before scaling out
- expect single-end RamDA to use `bwa` plus direct `CIRI3` SAM mode
- expect large SAM and cache files

## 1. Clone and create an environment

```bash
git clone <YOUR_CIRCYTO_REPO_URL>
cd circyto

conda create -n circyto python=3.12 -y
conda activate circyto

pip install -e .
```

Verify the CLI is available:

```bash
circyto --version
```

## 2. Verify external tools

`circyto` does not bundle every runtime dependency. Confirm these are available on the server:

- `bwa`
- `samtools`
- `java`
- the CIRI3 jar used by your local/server install

Recommended check:

```bash
circyto doctor
```

For single-end RamDA with `run-ciri3`, `bwa` plus direct CIRI3 execution is the expected route. `STAR` is not required for the single-end RamDA path.

## 3. Prepare a full-run manifest

Use the placeholder template at:

- `testdata/ramda/manifest_full_hg38_template.tsv`

Copy it to a working manifest and replace the placeholder FASTQ paths and sample IDs with your real server paths.

Required columns:

- `sample_id`
- `fastq_1`
- `fastq_2`
- `protocol`
- `strandedness`
- `read_layout`

For RamDA single-end input:

- `protocol=ramda`
- `strandedness=unstranded`
- `read_layout=single`
- leave `fastq_2` empty

## 4. First server run

Run a small 2-cell job first:

```bash
circyto run-ciri3 \
  --manifest path/to/manifest_full_hg38.tsv \
  --outdir path/to/run_ramda_hg38_2cell \
  --genome-fasta /path/to/hg38.fa \
  --gtf /path/to/gencode.gtf \
  --protocol ramda \
  --threads N
```

Expected single-end RamDA execution path:

```text
FASTQ -> BWA-MEM -> CIRI3 direct SAM mode (-Ma 0) -> matrix/log outputs
```

This is different from the paired-end Smart-seq3 STAR hybrid path.

## 5. Output expectations

Expect these directories:

- `align/`
- `ciri3/`
- `matrix/`
- `logs/`
- `qc/`

Important practical note:

- `align/` can contain large SAM files
- cache directories can grow quickly on hg38
- disk pressure is a real constraint for full-length total-RNA runs

## 6. Interpreting a small pilot run

A 2-cell validation run is successful if:

- alignment completes
- CIRI3 exits successfully
- expected output files exist
- the matrix is readable

Zero circRNA calls do not automatically mean pipeline failure, especially for a small pilot or a reduced test scope. Check the logs before treating an empty result as an error.
