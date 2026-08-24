# Local SComatic Environment Setup

This page documents the current local/server setup boundary for running tiny SComatic smoke tests alongside `circyto`.

Scope:

- local proof-of-concept only
- no hg38-scale execution here
- no modification of the upstream SComatic checkout supplied with `--scomatic-dir`

## Current WSL status

WSL local real-SComatic setup is currently blocked on this machine.

Observed failure mode:

- repeated `Bus error (core dumped)` during minimal `conda install` transactions
- most likely within the native dependency stack needed by `R`, `bedtools`, `pybedtools`, or related compiled packages

Because that is a native environment failure, not a `circyto` logic error, do not continue local WSL conda install attempts after the first repeated bus error.

WSL remains supported for:

- synthetic SComatic integration POC
- `circyto import-dna-snv-summary`
- `circyto summarize-dna-rna-circ`

## Why Python 3.10 is required

The upstream `SComatic` `requirements.txt` pins:

- `numpy==1.21.6`
- `pandas==1.3.5`
- `scipy==1.7.3`
- `pysam==0.16.0.1`

`numpy==1.21.6` is published for Python `<3.11`, so `pip install -r requirements.txt` fails on Python 3.11+.

For local reproducibility, use a dedicated Python 3.10 environment.

## Recommended safer routes

### A. HPC / server conda environment

Preferred for the first real smoke test:

```bash
conda create -n scomatic-smoke python=3.10 -y
conda activate scomatic-smoke
conda install -c conda-forge -c bioconda \
  samtools bedtools r-base r-vgam pandas pysam pybedtools -y
python /path/to/SComatic/scripts/BaseCellCounter/BaseCellCounter.py --help
```

If that works, then run the `circyto` chr21 real-smoke wrapper from the same environment.

### B. Container / mamba / micromamba environment

If WSL conda transactions remain unstable:

- use an HPC login node or batch job environment
- or use `mamba` / `micromamba`
- or use a containerized environment with the same Python 3.10 + native toolchain

## Historical local recipe

This is still the intended dependency shape, but on this WSL machine it should currently be treated as documentation, not a recommended live install path:

```bash
conda create -n scomatic python=3.10 -y
conda activate scomatic
conda install -c conda-forge -c bioconda \
  samtools bedtools datamash r-base r-vgam -y
pip install -r /path/to/SComatic/requirements.txt
Rscript -e 'library(VGAM); cat("VGAM OK\n")'
```

## Why `r-vgam` via conda fixes the R permission error

The upstream helper:

```r
install.packages(pkgs = c("VGAM"), repos = "http://cran.us.r-project.org")
```

tries to install into the default R site library, which commonly resolves to a non-writable system path such as:

- `/usr/local/lib/R/site-library`

Installing `r-vgam` through conda places the package inside the active conda environment instead, which avoids the permission failure and keeps the SComatic stack isolated.

## Environment setup helper

Print the recommended commands:

```bash
bash scripts/setup_scomatic_local_env.sh
```

Run the setup directly only on a stable HPC/server or containerized environment:

```bash
bash scripts/setup_scomatic_local_env.sh --run --env-name scomatic --scomatic-dir /path/to/SComatic
```

## Environment validation

Run:

```bash
bash scripts/check_scomatic_environment.sh --scomatic-dir /path/to/SComatic
```

This checks:

- current Python version
- `numpy` install
- `samtools`
- `bedtools`
- `datamash` if present
- `Rscript`
- `VGAM`
- required SComatic script paths

## chr21 local smoke-test boundary

### WSL note

Do not attempt further local WSL real-SComatic conda installs after a repeated `Bus error`.

At that point, stop and switch to:

- HPC/server conda
- or container / mamba / micromamba

### Server-side checklist

After the environment passes:

```bash
conda create -n scomatic-smoke python=3.10 -y
conda activate scomatic-smoke
conda install -c conda-forge -c bioconda \
  samtools bedtools r-base r-vgam pandas pysam pybedtools -y
python /path/to/SComatic/scripts/BaseCellCounter/BaseCellCounter.py --help

bash scripts/run_scomatic_chr21_poc.sh \
  --workdir work/local_chr21_simple_overlap_validation_20260519/real_se10k_rerun \
  --reference ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --outdir work/scomatic_chr21_real_smoke \
  --scomatic-dir /path/to/SComatic \
  --real-smoke
```

This is still only an environment/protocol smoke test:

- it prepares a tiny local BAM if needed
- it verifies the SComatic environment
- it runs a minimal `BaseCellCounter.py --help` launch test
- it does not perform genome-scale calling

## Troubleshooting

If `Bus error (core dumped)` appears during conda package installation:

- stop immediately
- classify it as a native environment instability
- do not keep retrying the same WSL install path
- switch to an HPC/server or containerized environment
