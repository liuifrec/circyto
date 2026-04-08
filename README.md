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

---

## What it is

`circyto` provides a reproducible command-line interface over multiple circRNA detectors (currently **CIRI-full**, **CIRI2**, and **find-circ3**) and standardizes outputs for downstream single-cell analysis (Scanpy / Seurat / scVI).

If you're new, your fastest path is:

1) install `circyto` + one detector  
2) run the bundled chr21 smoke test  
3) collect a `circ_counts.mtx` matrix  
4) move on to your own reference/manifest

---

## Install

### System prerequisites (high level)

You will need:

- **Python 3.10+**
- For `find-circ3` workflows: **bowtie2** and **samtools**
- For `ciri-full` workflows: **bwa** and **Java (JRE)**

> Tip: if something fails with “command not found”, jump to `docs/getting_started.md#external-dependencies`.

### Install `circyto` from source

```bash
git clone https://github.com/liuifrec/circyto
cd circyto

# (recommended) use a virtualenv
python -m venv .venv
source .venv/bin/activate

pip install -U pip
pip install -e .
```

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

---
## Quick check (recommended)

After installation, verify your environment and detector availability:

```bash
circyto doctor
circyto detectors
```

- **circyto doctor** checks required external dependencies.
- **circyto detectors** shows which detectors are ready to run.

These commands should work regardless of your current working directory as long as the installed package can see the bundled assets you configured.

## Minimal example (bundled chr21 smoke test)

This repository includes a small chr21 reference (`ref/chr21.fa`) and a small manifest (`manifest_2.tsv`) intended for smoke testing.

`manifest_2.tsv` is a **paired-end** example manifest.

What this smoke test proves:

- the detector command line runs with the bundled reference/manifest
- per-cell outputs can be collected into a sparse matrix
- the normalized TSV and matrix plumbing are working

What it does **not** prove:

- broad biological validity
- performance on your full dataset
- equivalence between different detector algorithms

### 1) Run a detector on the bundled manifest

```bash
circyto run-detector find-circ3 \
  --manifest manifest_2.tsv \
  --outdir work/find_circ3_chr21 \
  --ref-fa ref/chr21.fa \
  --threads 4 \
  --parallel 2
```

### 2) Collect a circRNA × cell MatrixMarket matrix

```bash
mkdir -p work/find_circ3_chr21_matrix

circyto collect-matrix \
  --detector find-circ3 \
  --indir work/find_circ3_chr21 \
  --matrix work/find_circ3_chr21_matrix/circ_counts.mtx \
  --circ-index work/find_circ3_chr21_matrix/circ_index.txt \
  --cell-index work/find_circ3_chr21_matrix/cell_index.txt
```

**Expected outcome**: `circ_counts.mtx` exists and has non-zero entries.

## Detector behavior note

`ciri-full` does **not** use one identical internal execution path for every layout:

- **paired-end manifest rows**: `circyto` runs the upstream bundled **CIRI-full Java Pipeline**
- **single-end manifest rows**: upstream CIRI-full Pipeline is paired-end only, so `circyto` falls back to the bundled **CIRI2-based detection path** and normalizes the output to the same TSV schema

This preserves manifest compatibility and downstream collection, but **single-end `ciri-full` runs are not identical to upstream CIRI-full full-length reconstruction**.

---

## Next steps

- **Getting started guide** (full workflows, references/manifests): `docs/getting_started.md`
- **Detector catalog** (what each detector needs & produces): `docs/detectors.md`
- **CLI contract** (“where do detector/outdir go?”): `docs/cli_policy.md`
- **Project strategy & milestones**: `ROADMAP.md`

---

## Citation

A methods manuscript is under preparation. In the meantime, please cite this repository:

> Liu, Y.-C. et al. “circyto: a unified CLI for single-cell circRNA detection and multimodal matrices.” GitHub repository: https://github.com/liuifrec/circyto
