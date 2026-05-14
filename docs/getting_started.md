# Getting started

This guide is the canonical “how do I run this without surprises?” workflow reference.

If you only want a minimal smoke test, see the README.

---

## Quick checklist

Before you start:

- You have **Python 3.10+**
- You installed `circyto` (`pip install -e .` from the repo)
- You installed **at least one detector**:
  - `find-circ3` (recommended for a first run), or
  - `ciri-full` / `ciri2` (requires Java + bwa + the bundled CIRI-full assets)
  - `ciri3` (requires Java + a CIRI3 jar, with BWA or STAR depending on alignment mode)

If your first run fails with `command not found`, skip ahead to **External dependencies**.

For Linux and HPC deployments, prefer the repo-root environment file instead of building the toolchain piecemeal:

```bash
conda env create -f environment.yml
conda activate circyto
pip install -e .
```

That environment is the supported runtime baseline for source installs. It targets **Python 3.10** and keeps the solve focused on the Python scientific stack on `conda-forge`: `numpy`, `scipy`, `pandas`, `anndata`, `h5py`, `typer`, `rich`, and `tqdm`.

It intentionally does **not** include detector executables such as `bwa`, `samtools`, `STAR`, `bowtie2`, `find-circ3`, `perl`, or Java. Install those separately only for the workflows you actually plan to run. This keeps the default environment reproducible across Linux, WSL, and HPC conda installs.

---

## Quick smoke

For a fast installation and workflow check, use:

```bash
circyto smoke --detector ciri3 --aligner bwa-mem
```

If STAR mode is configured:

```bash
circyto smoke --detector ciri3 --aligner star
```

Smoke validates runtime resolution, workflow mechanics, and matrix collection. It is not full biological validation. Tiny smoke subsets may legitimately produce empty outputs; that still counts as success unless you add `--require-nonempty`.

---

## Concepts

### What is a “detector”?

A detector is the underlying circRNA calling tool (e.g., CIRI-full, find-circ3). `circyto`:

- runs the detector per-cell (or per-sample)
- normalizes output locations and naming
- provides standardized collectors to build a circRNA × cell matrix

### What is a “manifest”?

A manifest is a TSV file listing your samples/cells and their FASTQ paths (plus any required metadata).

Instead of describing the schema in prose (which drifts), **use the repo examples as templates**:

- `manifest.tsv`
- `manifest_2.tsv`

When you create your own manifest, start by copying one of these files and replacing the FASTQ paths.

---

## External dependencies

Circyto orchestrates detectors, but many detectors rely on standard bioinformatics executables being available on your `PATH`.

### Common executables

- `bowtie2` and `samtools` (required for find-circ3 workflows)
- `bwa` and `java` (required for CIRI-full workflows)
- `java` and `samtools` (required for CIRI3 workflows)
- `bwa` (required for CIRI3 BWA mode)
- `STAR` (required only for CIRI3 STAR mode)

Use the built-in diagnostics first:

```bash
circyto doctor
circyto detectors
```

In the default environment, these commands validate the installed Python package and report workflow-specific executables as optional warnings if they are absent.

For bundled CIRI3 direct mode, `circyto doctor` now checks both `java` presence and the detected Java major version. The bundled CIRI3 jar requires **Java >= 12**; **Java 17** is the recommended target.

If you need a manual spot check:

```bash
command -v bowtie2 samtools bwa java
```

### CIRI3 environment requirements

CIRI3 requires all of the following:

- Java 12+ for the bundled jar; Java 17 recommended
- a CIRI3 jar: either vendored under `tools/CIRI3/` or configured via environment variables
- `samtools`: required for alignment handling and inspection
- `bwa`: required for BWA alignment-first workflows
- `STAR`: required only for STAR-based workflows

Supported environment variables:

- `CIRCYTO_CIRI3_HOME`
- `CIRCYTO_CIRI3_JAR`
- `CIRCYTO_CIRI3_JAVA`
- `CIRCYTO_CIRI3_EXTRA_ARGS`

For STAR alignment-first runs, `circyto` executes STAR in a local Linux temp workspace and copies the required outputs back into the alignment cache. If your cluster expects a specific node-local scratch path, set `CIRCYTO_STAR_TMPDIR` before running `prepare-alignment-cache`.

`CIRCYTO_CIRI3_CMD_TEMPLATE` is also supported when you want explicit template execution instead of the direct `java -jar` path.

Expected vendored layout:

```text
tools/
  CIRI3/
    CIRI3_Java_*.jar
```

Minimal verification:

```bash
export CIRCYTO_CIRI3_HOME=/path/to/CIRI3
circyto doctor
circyto detectors
```

Readiness semantics:

- `READY`: circyto found a usable CIRI3 jar and Java, so the direct `java -jar` contract is available.
- For bundled direct mode, `READY` also requires a detected Java version `>=12`.
- `NOT READY`: circyto could not find a usable CIRI3 jar or execution path.
- `PARTIAL`: local CIRI3 assets were detected but the execution contract is incomplete.

---

## Workflow 1: run a detector on many cells

`run-detector` is the canonical single-detector command.

`run-batch` is still available, but it is effectively an alias of `run-detector` with a required `--detector` flag.

### Example: find-circ3 on the bundled chr21 manifest

```bash
circyto run-detector find-circ3 \
  --manifest manifest_2.tsv \
  --outdir work/find_circ3_chr21 \
  --ref-fa ref/chr21.fa \
  --threads 4 \
  --parallel 2
```

What you should see:

- an output directory under `work/find_circ3_chr21/`
- per-cell subdirectories and detector outputs (layout varies by detector)

If this fails:
- confirm `find-circ3 --help` works
- confirm `bowtie2` and `samtools` are on your PATH

### Example: `ciri-full` on your own manifest

```bash
circyto run-detector ciri-full \
  --manifest manifest.tsv \
  --outdir work/ciri_full_run \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 4 \
  --parallel 1
```

`ciri-full` layout semantics are important:

- paired-end rows (`r1` + `r2`) use the upstream bundled **CIRI-full Java Pipeline**
- single-end rows (`r1` only) use a bundled **CIRI2-based fallback path**
- both layouts are normalized to the same per-cell TSV schema for downstream collection

If you specifically need upstream CIRI-full full-length reconstruction behavior, use paired-end input.

---

## Workflow 1b: experimental SMART-Seq3 to CIRI3

For validated SMART-Seq3 pooled demux plus STAR+CIRI3 orchestration, use the experimental high-level workflow:

```bash
circyto workflow smartseq3-ciri3 \
  --read1 emtab8735/subset_100k/diySpike.R1.100k.fastq.gz \
  --read2 emtab8735/subset_100k/diySpike.R2.100k.fastq.gz \
  --index1 emtab8735/subset_100k/diySpike.I1.100k.fastq.gz \
  --index2 emtab8735/subset_100k/diySpike.I2.100k.fastq.gz \
  --annotation path/to/emtab8735_annotation.tsv \
  --cell-id-column cell_id \
  --index1-column index1 \
  --index2-column index2 \
  --ref-fa ref/chr21.fa \
  --star-index ref/star_index_chr21 \
  --outdir work/emtab8735_smartseq3_ciri3 \
  --top-n 20 \
  --export-h5ad \
  --resume
```

What it does:

- demultiplexes SMART-Seq3 transcript reads with `circyto demux smartseq3`
- writes `OUTDIR/manifests/top<N>_manifest.tsv` or `OUTDIR/manifests/all_manifest.tsv`
- runs `prepare-alignment-cache --aligner star --detector ciri3`
- runs `run-detector-from-alignments --detector ciri3`
- runs `collect-matrix`
- writes `OUTDIR/qc/cell_qc.tsv` and `OUTDIR/qc/circ_qc.tsv`
- writes `OUTDIR/anndata/circ_counts.h5ad`
- writes `workflow_summary.json`

Server-validated E-MTAB-8735 diySpike summary:

- full demux: `75,015,128` reads processed, `68,649,627` assigned, `192` cells detected
- all192 hg38 STAR+CIRI3: completed successfully
- all192 final result: `192` cells, `2503` circRNAs, `2659` nonzero entries
- top20 STAR+CIRI3: all `20` cells succeeded
- matrix: `588` circRNAs x `20` cells, `600` nonzero entries

Analysis-ready exports:

- `matrix/`
- `qc/`
- `anndata/`
- `workflow_summary.json`

`circ_counts.h5ad` is ready for Scanpy / ScGeo-style downstream analysis:

- `adata.X`: circRNA counts with shape `cells x circRNAs`
- `adata.obs`: cell QC indexed by `cell_id`
- `adata.var`: circRNA QC / feature metadata indexed by `circ_id`
- `adata.uns["circyto"]`: workflow provenance

Optional multimodal export for joint RNA + circRNA analysis:

```bash
circyto workflow smartseq3-ciri3 \
  ... \
  --export-h5ad \
  --gene-counts path/to/gene_counts.tsv \
  --gene-counts-format tsv \
  --export-mudata
```

When enabled, the workflow writes `OUTDIR/mudata/circyto_multimodal.h5mu` with:

- `mdata.mod["rna"]` = gene expression AnnData
- `mdata.mod["circ"]` = circRNA AnnData
- shared cell metadata in `mdata.obs`
- workflow provenance in `mdata.uns["circyto"]`

If `mudata` is not installed, the workflow only errors when `--export-mudata` is explicitly requested.

For lightweight downstream summaries from the circ-only AnnData:

```bash
circyto analyze summarize-h5ad \
  --input work/emtab8735_smartseq3_ciri3/anndata/circ_counts.h5ad \
  --outdir work/emtab8735_smartseq3_ciri3/anndata_summary
```

This support is still experimental:

- use local `ref/chr21.fa` and `ref/star_index_chr21` for smoke tests
- do not require full hg38 for CI or local preflight
- keep the low-level `demux`, `prepare-alignment-cache`, `run-detector-from-alignments`, and `collect-matrix` commands available for manual control

---

## Workflow 2: collect a circRNA × cell matrix

Once you have per-cell calls, build a single sparse matrix:

```bash
mkdir -p work/find_circ3_chr21_matrix

circyto collect-matrix \
  --detector find-circ3 \
  --indir work/find_circ3_chr21 \
  --matrix work/find_circ3_chr21_matrix/circ_counts.mtx \
  --circ-index work/find_circ3_chr21_matrix/circ_index.txt \
  --cell-index work/find_circ3_chr21_matrix/cell_index.txt
```

Artifacts:

- `circ_counts.mtx` – sparse MatrixMarket counts
- `circ_index.txt` – circ feature IDs (row index)
- `cell_index.txt` – cell IDs (column index)

---

## Workflow 3: host-gene annotation and multimodal export (optional)

Circyto also supports:

- `annotate-host-genes` to map circ features to host genes using a GTF
- `annotate-circs` to label known versus novel circRNAs against generic circRNA database TSVs
- `export-multimodal` to create an `.h5ad` that combines:
  - existing mRNA expression (`.X`)
  - circRNA counts (in `obsm["X_circ"]`)
  - circ feature metadata (in `uns["circ"]`)

Because these commands evolve, rely on `--help` for exact flags:

```bash
circyto annotate-host-genes --help
circyto annotate-circs --help
circyto export-multimodal --help
```

---

## Workflow 4: multi-detector runs (optional)

If you want to run multiple detectors on the same manifest:

```bash
circyto run-multidetector \
  ciri-full find-circ3 \
  --manifest manifest_2.tsv \
  --outdir work/multidetector_chr21 \
  --ref-fa ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --threads 4 \
  --parallel 2
```

Then explore:

```bash
circyto merge-detectors --help
circyto compare-detectors --help
circyto collect-multidetector --help
```

---

## CIRI3 setup

Use the alignment-first entrypoints for CIRI3:

Expected layout:

```text
tools/
  CIRI3/
    CIRI3_Java_*.jar
```

Verify the runtime first:

```bash
export CIRCYTO_CIRI3_HOME=/path/to/CIRI3
circyto doctor
circyto detectors
```

```bash
circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner bwa-mem \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/alignment_cache

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --outdir work/ciri3_run \
  --ref-fa ref/genome.fa
```

Minimal required tools by mode:

- BWA mode: `bwa`, `samtools`, `java`
- STAR mode: `STAR`, `samtools`, `java`

Validated local BWA mode:

- alignment input: unsorted SAM
- `bwa mem` parameters: `-k 15 -T 15`
- CIRI3 flags: `-S 0 -Ma 0`

STAR mode is supported in code for alignment-first workflows, but is not yet fully validated end-to-end in this release.

Reference consistency matters: the alignment cache, the alignment manifest, and the `--ref-fa` used for direct CIRI3 execution must all match the same reference build.
