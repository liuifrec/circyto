<p align="center">
  <img src="assets/circCyto_logo.png" width="440" alt="circyto logo">
</p>

<p align="center">
  <a href="https://github.com/liuifrec/circyto/releases">
  <img alt="Release" src="https://img.shields.io/github/v/tag/liuifrec/circyto?sort=semver">
  </a>
  <a href="https://github.com/liuifrec/circyto/blob/main/LICENSE">
  <img alt="License" src="https://img.shields.io/github/license/liuifrec/circyto?cacheSeconds=3600&v=2">
  </a>
  <a href="https://pypi.org/project/circyto/"><img alt="Python" src="https://img.shields.io/badge/python-3.10%2B-blue"></a>
</p>

## circyto


A single-cell circRNA detection CLI to orchestrate detectors, build circRNA×cell matrices, and export multimodal AnnData.

## Current status

`circyto` is currently on the **v0.8.x** release line (`pyproject.toml` currently reports `0.8.3`).

- The public detector workflow is manifest-driven and stable enough for routine use.
- `circyto doctor` and `circyto detectors` are live commands, not planned features.
- Bundled detector asset resolution is intended to be **cwd-independent**.
- Detector integrations are still heterogeneous; detector-specific limits are documented rather than hidden.
- CIRI3 is a real alignment-first backend with direct `java -jar` support when its local runtime contract is available.

---

## What it is

`circyto` provides a reproducible command-line interface over multiple circRNA detectors (currently **CIRI-full**, **CIRI2**, **find-circ3**, and an alignment-first **CIRI3** backend) and standardizes outputs for downstream single-cell analysis (Scanpy / Seurat / scVI).

If you're new, your fastest path is:

1) install `circyto` + one detector  
2) run the bundled smoke test  
3) collect a `circ_counts.mtx` matrix  
4) move on to your own reference/manifest

---

## Install

### System prerequisites (high level)

You will need:

- **Python 3.10+**
- For `find-circ3` workflows: **bowtie2** and **samtools**
- For `ciri-full` workflows: **bwa** and **Java (JRE)**
- For `ciri3` workflows: **Java** plus a usable **CIRI3 jar**
- For `ciri3` BWA alignment-first workflows: **bwa** and **samtools**
- For `ciri3` STAR alignment-first workflows: **STAR**, **bwa**, and **samtools**

> Tip: if something fails with “command not found”, jump to `docs/getting_started.md#external-dependencies`.

### Install `circyto` from source

```bash
git clone https://github.com/liuifrec/circyto
cd circyto

# supported baseline on Linux / WSL / HPC
conda env create -f environment.yml
conda activate circyto

# installs the repo itself into that environment
pip install -e .
```

The repo-root [environment.yml](/mnt/d/circyto/environment.yml) is the supported source-install baseline. It targets **Python 3.10** and keeps the default solve limited to the core Python runtime stack on `conda-forge`: `numpy`, `scipy`, `pandas`, `anndata`, `h5py`, `typer`, `rich`, and `tqdm`.

It intentionally does **not** install detector executables such as `bwa`, `samtools`, `STAR`, `bowtie2`, `perl`, or Java. Those tools are workflow-specific, much more solver-fragile on mixed Linux/HPC systems, and are better installed separately through site modules, `apt`, `brew`, or a detector-specific conda environment. The default environment is meant to support:

- `pip install -e .`
- `circyto --help`
- `circyto doctor`
- `circyto detectors --json`

### Install at least one detector

#### Option A: find-circ3 (recommended for your first smoke test)

```bash
git clone https://github.com/liuifrec/find_circ3.git
pip install -e ./find_circ3

find-circ3 --help
```

#### Option B: CIRI-full

`circyto` expects the **CIRI-full JAR** to be available under `tools/` (see repository `tools/` layout).
You also need `bwa` and a Java runtime.

See: `docs/detectors.md#ciri-full`.

#### Option C: CIRI3

CIRI3 requires:

- `java` (mandatory)
- a local CIRI3 jar, either vendored under `tools/CIRI3/` or configured explicitly
- `bwa` for BWA-based alignment-first runs
- `STAR` only for STAR-based runs
- `samtools` for alignment handling and inspection

Supported environment variables:

- `CIRCYTO_CIRI3_HOME`
- `CIRCYTO_CIRI3_JAR`
- `CIRCYTO_CIRI3_JAVA`
- `CIRCYTO_CIRI3_EXTRA_ARGS`
- `CIRCYTO_CIRI3_CMD_TEMPLATE`

Expected vendored layout:

```text
tools/
  CIRI3/
    CIRI3_Java_*.jar
```

Minimal setup example:

```bash
conda activate circyto
export CIRCYTO_CIRI3_HOME=/path/to/CIRI3
circyto doctor
circyto detectors
```

STAR alignment-first note:

- `circyto` treats `CIRI3 + STAR` as the official upstream hybrid workflow, not as STAR-only.
- The STAR prepare stage must emit `Aligned.out.sam`, `Chimeric.out.junction`, `Unmapped.out.mate1`, and `Unmapped.out.mate2`.
- `circyto` then runs `bwa mem -T 19` on the unmapped mates and passes `Chimeric.out.junction,Aligned.out.sam,bwa_rescue.sam` into CIRI3 with `-Ma 1`.
- STAR alignment prep now runs STAR inside a local Linux temp directory and then copies the required artifacts back into the cache. Override the temp base with `CIRCYTO_STAR_TMPDIR` if your cluster requires a specific local scratch path.
- This STAR alignment-first contract is currently specific to `ciri3`. Other bundled detectors should be considered BWA-oriented or FASTQ-native unless documented otherwise.

---
## Quick check (recommended)

After installation, verify your environment and detector availability:

```bash
circyto doctor
circyto detectors
```

- **circyto doctor** checks required external dependencies.
- **circyto detectors** shows which detectors are ready to run.
- In the default baseline environment, missing detector executables are reported as optional warnings rather than a failed installation.

For CIRI3 specifically:

- `READY` means circyto found a usable CIRI3 jar and `java`, so direct `java -jar` execution is available.
- `NOT READY` means circyto could not find a usable CIRI3 jar or execution path.
- `PARTIAL` means local CIRI3 assets were found but the runtime contract is incomplete.

These commands should work regardless of your current working directory as long as the installed package can see the bundled assets you configured.

## Minimal example (bundled chr21 smoke test)

The recommended installation smoke path is:

```bash
circyto smoke --detector ciri3 --aligner bwa-mem
```

If STAR prerequisites are present, you can also exercise the official hybrid route:

```bash
circyto smoke --detector ciri3 --aligner star
```

The smoke command builds a tiny local chr21 subset, runs the current alignment-first CIRI3 workflow, and writes a compact `smoke_summary.json` under the chosen `--outdir`.

What this smoke test proves:

- runtime resolution is usable
- the alignment-first workflow can complete
- per-cell outputs can be collected into a sparse matrix
- the normalized TSV and matrix plumbing are working

What it does **not** prove:

- broad biological validity
- performance on your full dataset
- equivalence between different detector algorithms
- guaranteed non-empty biological calls on a tiny subset

By default, smoke treats an empty detector output or empty matrix as a PASS if the full workflow completed correctly. Use `--require-nonempty` when you want smoke to fail on all-zero outputs.

## Detector behavior note

`ciri-full` does **not** use one identical internal execution path for every layout:

- **paired-end manifest rows**: `circyto` runs the upstream bundled **CIRI-full Java Pipeline**
- **single-end manifest rows**: upstream CIRI-full Pipeline is paired-end only, so `circyto` falls back to the bundled **CIRI2-based detection path** and normalizes the output to the same TSV schema

This preserves manifest compatibility and downstream collection, but **single-end `ciri-full` runs are not identical to upstream CIRI-full full-length reconstruction**.

---

## Alignment-first workflow

For large demultiplexed datasets, circyto now includes an alignment-first execution track so alignments can be prepared once and reused across downstream detector runs.

### CIRI3 setup

Expected repo-local layout:

```text
tools/
  CIRI3/
    CIRI3_Java_*.jar
```

Supported environment variables:

- `CIRCYTO_CIRI3_HOME`
- `CIRCYTO_CIRI3_JAR`
- `CIRCYTO_CIRI3_JAVA`
- `CIRCYTO_CIRI3_EXTRA_ARGS`
- `CIRCYTO_CIRI3_CMD_TEMPLATE`

Verify the runtime contract before a real alignment-first CIRI3 run:

```bash
export CIRCYTO_CIRI3_HOME=/path/to/CIRI3
circyto doctor
circyto detectors
```

Mode-specific minimum tools:

- BWA mode: `bwa`, `samtools`, `java`
- STAR mode: `STAR`, `samtools`, `java`

Readiness semantics:

- `READY`: CIRI3 jar and Java were detected and direct `java -jar` execution is usable.
- `NOT READY`: circyto could not find a usable CIRI3 jar or execution path.
- `PARTIAL`: local CIRI3 assets were detected, but the direct or template execution contract is incomplete.

Direct CIRI3 execution requires a detected CIRI3 jar plus Java. Java is mandatory for direct mode. `STAR` is not required unless STAR mode is used. Template execution remains available through `--command-template` or `CIRCYTO_CIRI3_CMD_TEMPLATE`; template mode does not require `--ref-fa` unless the template itself uses `{ref_fa}`.

Validated local BWA + CIRI3 example:

```bash
circyto plan-alignment-cache \
  --manifest manifest.tsv \
  --aligner bwa-mem \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/alignment_cache

circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner bwa-mem \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/alignment_cache \
  --sentinel-cells 8 \
  --chunk-size 48

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --outdir work/ciri3_run \
  --ref-fa ref/genome.fa \
  --chunk-size 64
```

Assumptions:

- the `bwa` index matches `ref/genome.fa`
- the alignment manifest rows resolve to unsorted SAM for direct CIRI3 execution
- the same reference build is used for alignment preparation and detector execution

Validated local settings:

- aligner: `bwa mem -k 15 -T 15`
- CIRI3: `-S 0 -Ma 0`

STAR + CIRI3 example:

```bash
circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner star \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/alignment_cache_star

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache_star/alignment_manifest.tsv \
  --outdir work/ciri3_star_run \
  --ref-fa ref/genome.fa
```

STAR assumptions:

- a matching STAR genome index is available to the alignment step
- STAR chimeric outputs are present in the alignment manifest
- the same reference build is used throughout

BWA + CIRI3 is the validated local production path. STAR + CIRI3 is supported in code for alignment-first workflows, but is not yet fully validated end-to-end in this release.

Recover after a partial failure:

```bash
circyto summarize-run-state \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --run-dir work/ciri3_run \
  --mode detector

circyto export-run-subset \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --run-dir work/ciri3_run \
  --subset incomplete \
  --out work/ciri3_incomplete.tsv

# For stale outputs only:
circyto export-run-subset \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --run-dir work/ciri3_run \
  --subset stale \
  --out work/ciri3_stale.tsv
```

Validate the local plumbing with repo-shipped assets:

```bash
circyto smoke --detector ciri3 --aligner bwa-mem --outdir work/smoke_bwa
circyto smoke --detector ciri3 --aligner star --outdir work/smoke_star
```

See `docs/alignment_first_execution.md` for cache keys, resume behavior, chunk recovery, and CIRI3 execution-mode details.

---

## Next steps

- **Getting started guide** (full workflows, references/manifests): `docs/getting_started.md`
- **Detector catalog** (what each detector needs & produces): `docs/detectors.md`
- **Alignment-first execution** (shared preprocessing and reusable BAM/SAM manifests): `docs/alignment_first_execution.md`
- **CLI contract** (“where do detector/outdir go?”): `docs/cli_policy.md`
- **Project strategy & milestones**: `ROADMAP.md`

---

## Citation

A methods manuscript is under preparation. In the meantime, please cite this repository:

> Liu, Y.-C. et al. “circyto: a unified CLI for single-cell circRNA detection and multimodal matrices.” GitHub repository: https://github.com/liuifrec/circyto
