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

If your first run fails with `command not found`, skip ahead to **External dependencies**.

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
- `STAR` (optional; relevant only for future STAR-based detectors)

Use the built-in diagnostics first:

```bash
circyto doctor
circyto detectors
```

If you need a manual spot check:

```bash
command -v bowtie2 samtools bwa java
```

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
