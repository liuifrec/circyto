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

That environment covers the common executables used by the currently integrated workflows: `bwa`, `samtools`, `STAR`, `bowtie2`, `perl`, and Java 17. It does not include separately obtained detector artifacts such as the CIRI3 jar.

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

If you need a manual spot check:

```bash
command -v bowtie2 samtools bwa java
```

### CIRI3 environment requirements

CIRI3 requires all of the following:

- Java: mandatory
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
- `export-multimodal` to create an `.h5ad` that combines:
  - existing mRNA expression (`.X`)
  - circRNA counts (in `obsm["X_circ"]`)
  - circ feature metadata (in `uns["circ"]`)

Because these commands evolve, rely on `--help` for exact flags:

```bash
circyto annotate-host-genes --help
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
