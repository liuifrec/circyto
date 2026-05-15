# Public Dataset Preparation

`circyto prepare-public-dataset` creates reproducible metadata templates and download plans for known benchmark datasets without downloading real FASTQ files.

It now surfaces dataset-level scientific context explicitly so organism and reference assumptions are harder to miss.

The command is intentionally planning-only:

- tests do not require network access
- dry runs only write small text artifacts
- accession expansion can be added later without changing the CLI contract

## Curated vs Placeholder Mode

When a dataset-specific curated table exists under `testdata/public_datasets/`, `circyto prepare-public-dataset` uses that TSV as the source of truth for `selected_runs.tsv`.

Current curated-table targets:

- `testdata/public_datasets/emtab8735_runs.tsv`
- `testdata/public_datasets/gse98664_runs.tsv`

If a curated table is missing, the command falls back to embedded placeholder rows in the Python module. Placeholder rows are marked with `TODO_REAL_ACCESSION` in `run_id`, `fastq_1`, `fastq_2`, or `notes` so they are easy to spot.

To replace placeholder rows with real accessions, edit or extend the curated TSV directly. Keep the header stable:

- `dataset_id`
- `run_id`
- `sample_id`
- `protocol`
- `source`
- `fastq_1`
- `fastq_2`
- `notes`

The generated `selected_runs.tsv` also includes dataset-level metadata columns:

- `organism`
- `expected_read_layout`
- `expected_reference`
- `recommended_route`

## Command

```bash
circyto prepare-public-dataset \
  --dataset-id E-MTAB-8735 \
  --protocol smartseq3 \
  --download-method ena \
  --max-runs 2 \
  --dry-run \
  --outdir work/public_emtab8735_plan
```

Dry-run output files:

- `selected_runs.tsv`
- `download_plan.sh`
- `README_next_steps.md`

Dry-run CLI output may also print a warning block when the selected dataset is scientifically incompatible with a common misuse case. For example, `GSE98664` is explicitly marked as mouse and warns against using `hg38` for biological validation.

## Supported dataset IDs

- `E-MTAB-8735`
  Protocol: `smartseq3`
  Use: existing validated Smart-seq3 / DIY spike-in / 192-cell CIRI3 benchmark
- `GSE98664`
  Protocol: `ramda`
  Use: original RamDA-seq / RNA-only full-length total RNA benchmark
- `GSE278944`
  Protocol: `scrr`
  Use: future scRepli-RamDA-seq / DNA replication plus full-length RNA benchmark
- `shin-ramda-riken`
  Protocol: `shin-ramda`
  Use: Shin-RamDA workflow reference from `rikenbit/shin-ramda-seq-paper`

## Dataset-specific behavior

### E-MTAB-8735

This target prefers `testdata/public_datasets/emtab8735_runs.tsv`. The current curated subset is intentionally small and is suitable for reproducible planning only. If confirmed ENA / ArrayExpress run accessions are not available locally, rows remain marked `TODO_REAL_ACCESSION` until someone updates the TSV manually.

Planner metadata:

- organism: `Homo sapiens`
- expected read layout: `paired-end`
- expected reference: `hg38`
- recommended route: `Smart-seq3 paired-end workflow`

### GSE98664 and GSE278944

These targets write placeholder SRA-style shell plans such as:

```bash
prefetch <RUN>
fasterq-dump <RUN> --split-files --threads 8 -O <outdir>/fastq
```

`GSE98664` prefers `testdata/public_datasets/gse98664_runs.tsv`. `GSE278944` currently uses placeholder fallback rows only. Any unresolved rows are explicitly marked with `TODO_REAL_ACCESSION` in the curated table or fallback notes.

Important for `GSE98664`:

- organism: `Mus musculus`
- expected reference: `mm10/mm39`
- expected read layout: `single-end`
- recommended route: `BWA+CIRI3 single-end`
- do not use it as a human `hg38` biological validation target

### shin-ramda-riken

This target is README-only for now. It writes the standard artifacts but does not generate download commands and instead points to the upstream GitHub workflow reference:

- https://github.com/rikenbit/shin-ramda-seq-paper

## Testing

The planner is dependency-light and designed for synthetic fixture tests only:

```bash
pytest -q tests/test_public_dataset_prepare.py
```
