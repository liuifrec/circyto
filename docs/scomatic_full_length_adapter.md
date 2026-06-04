# Full-Length RNA SComatic Adapter

## Scope

`circyto` provides a generic full-length RNA SComatic interoperability layer for:

- Smart-seq2
- Smart-seq3
- RamDA-seq
- ShinRamDA
- scRR RNA arm

This layer prepares inputs, records runnable SComatic command plans, imports external SComatic outputs, and merges normalized RNA-derived candidate variant signals into MuData.

SComatic outputs in circyto are always described as:

```text
RNA-derived candidate variant signals
```

They are not validated somatic mutations unless independent DNA validation exists.

## Environment Model

Use separate environments for circyto and SComatic:

```bash
conda activate circyto
circyto prepare-scomatic-input \
  --alignment-manifest alignment_manifest.tsv \
  --outdir scomatic_prepared \
  --protocol ramda \
  --reference-fasta hg38.fa

conda activate scomatic-smoke
bash scomatic_run/scomatic_run_plan.sh

conda activate circyto
circyto import-scomatic \
  --scomatic-output sample.step2.calling.step2.tsv \
  --cell-annotations scomatic_prepared/cell_annotations.tsv \
  --provenance-metadata scomatic_run/run_scomatic_summary.json \
  --outdir scomatic_imported

circyto merge-scomatic \
  --input rna_circ.h5mu \
  --scomatic-candidates scomatic_imported/scomatic_candidate_summary.tsv \
  --output rna_circ_candidate_snv.h5mu
```

Recommended split:

- circyto environment: preparation, documentation, import, normalization, AnnData/MuData merge
- `scomatic-smoke` or equivalent SComatic environment: external BaseCellCounter, Step1, and Step2

Local WSL/native SComatic execution has been unstable in prior testing. Prefer server, container, or HPC execution for real SComatic runs. `circyto run-scomatic` defaults to planning-only and does not invoke external SComatic unless `--execute` is explicitly supplied.

## Commands

### `circyto prepare-scomatic-input`

Prepare SComatic-ready inputs from either one-cell-per-BAM/SAM alignments or an already merged multi-cell BAM.

One-cell-per-BAM/SAM mode:

```bash
circyto prepare-scomatic-input \
  --alignment-manifest alignment_manifest.tsv \
  --outdir scomatic_prepared \
  --protocol ramda \
  --reference-fasta hg38.fa \
  --sample-id hap1_batch10 \
  --cell-type-column cell_type
```

This mode:

- reads `cell_id` plus `bam` or `sam` from the manifest
- injects or normalizes `CB:Z:<cell_id>` tags
- preserves and counts `nM` and `NH`
- writes per-cell sorted BAMs
- writes a merged SComatic BAM
- writes `cell_annotations.tsv` with `Index` and `Cell_type`

Merged BAM mode:

```bash
circyto prepare-scomatic-input \
  --merged-bam pooled.scomatic.bam \
  --cell-metadata cell_metadata.tsv \
  --outdir scomatic_prepared \
  --protocol smartseq3 \
  --cell-type-column cell_type
```

This mode assumes the merged BAM already has usable `CB` tags. circyto records the BAM path and writes SComatic cell annotations; it does not rewrite the BAM.

Outputs:

```text
adapter_rows.tsv
cell_annotations.tsv
scomatic_prepared_bams.tsv
prepare_scomatic_input_summary.json
per_cell/*.scomatic.sorted.bam      # one-cell-per-BAM/SAM mode
merged/*.scomatic.bam               # one-cell-per-BAM/SAM mode
tag_reports/*.tags.tsv              # one-cell-per-BAM/SAM mode
```

### `circyto run-scomatic`

Create an external SComatic command plan and optionally execute it.

```bash
circyto run-scomatic \
  --prepared-dir scomatic_prepared \
  --scomatic-dir /path/to/SComatic \
  --reference-fasta hg38.fa \
  --outdir scomatic_run \
  --run-step1 \
  --run-step2
```

Default behavior is dry-run planning. Add `--execute` only when the SComatic environment, reference, filters, and Panel of Normals are ready. For split-environment operation, run this command in the circyto environment to write `scomatic_run_plan.sh`, then execute that plan from the SComatic environment.

The command writes:

```text
scomatic_run_plan.sh
run_scomatic_summary.json
```

Custom SComatic arguments can be supplied with:

- `--basecellcounter-args`
- `--step1-args`
- `--step2-args`

Supported placeholders include:

- `{merged_bam}`
- `{cell_annotations}`
- `{reference_fasta}`
- `{basecellcounter_out}`
- `{step1_out}`
- `{step2_out}`
- `{threads}`

### `circyto import-scomatic`

Import SComatic Step1, Step2, or generic candidate tables into the normalized circyto schema.

```bash
circyto import-scomatic \
  --scomatic-output sample.step2.calling.step2.tsv \
  --cell-annotations scomatic_prepared/cell_annotations.tsv \
  --provenance-metadata scomatic_run/run_scomatic_summary.json \
  --outdir scomatic_imported
```

Outputs:

```text
scomatic_candidate_summary.tsv
scomatic_candidate_summary.tsv.provenance.json
normalize_scomatic_results_summary.json
import_scomatic_summary.json
```

The normalized schema is:

```text
variant_id
cell_id
chrom
pos
ref
alt
gene
filter_status
candidate_variant_class
read_support
vaf
caller
```

### `circyto merge-scomatic`

Add normalized RNA-derived candidate variant signals to MuData:

```bash
circyto merge-scomatic \
  --input rna_circ.h5mu \
  --scomatic-candidates scomatic_imported/scomatic_candidate_summary.tsv \
  --output rna_circ_candidate_snv.h5mu
```

The resulting MuData contains:

```text
MuData
|- rna
|- circ
`- candidate_snv
```

`candidate_snv.X` stores cell x candidate-variant read support.

Layers:

- `candidate_snv.layers["vaf"]`: maximum VAF per cell and candidate variant
- `candidate_snv.layers["presence"]`: binary presence of RNA-derived candidate variant signals

## Validation Status

Validated:

- one-cell-per-BAM/SAM preparation with CB tag injection or normalization
- `nM` and `NH` preservation and missing-tag counting
- merged BAM and SComatic `Index`/`Cell_type` annotation writing
- real SComatic Step1/Step2 output normalization into `scomatic_candidate_summary.tsv`
- HAP1 batch10 technical smoke through BaseCellCounter, Step1, Step2, and normalization

Partially validated:

- `candidate_snv` MuData export through `merge-scomatic`, covered by synthetic AnnData/MuData tests
- dry-run command planning for external BaseCellCounter and optional Step1/Step2
- split-environment workflow between circyto and `scomatic-smoke`

Not yet validated:

- clinically validated variants
- validated somatic mutation discovery
- genome-scale biological interpretation of SComatic outputs
- stable candidate behavior in homogeneous or one-cell-type datasets
- SComatic execution from local WSL as a supported path

All imported SComatic rows remain RNA-derived candidate variant signals unless orthogonal DNA validation is added outside circyto.

## Lessons Learned from HAP1 and IMR90

- HAP1 batch10 validated the technical path: circyto alignment outputs can be adapted for SComatic, CB tags can be injected or normalized, NH/nM tags can be preserved, and real Step1/Step2 outputs can be normalized.
- HAP1 Step2 produced no candidate calls in the single-cell-type technical smoke. This is compatible with SComatic's cell-type/group contrast design and should not be treated as an integration failure.
- IMR90 establishes that the primary scRR DNA branch is processed DNA state integration, especially CNV for IMR90, not RNA-derived variant interpretation.
- Both datasets support keeping SComatic as an optional sidecar for RNA-derived candidate variant signals while preserving CNV or RT/state as the main scRR DNA modalities.
- Future manuscripts should report SComatic as technical interoperability and hypothesis generation, not validated somatic mutation discovery.

## Validated Path

Validated in circyto:

- one-cell-per-BAM/SAM preparation with CB tag injection/normalization
- `nM` and `NH` preservation/counting
- merged BAM and SComatic annotation writing
- real SComatic Step1/Step2 output normalization
- HAP1 batch10 technical smoke through BaseCellCounter, Step1, Step2, and normalization
- synthetic `candidate_snv` MuData merge tests

## Exploratory Path

Exploratory:

- biological interpretation of RNA-derived candidate variant signals
- full genome-scale SComatic execution outside the validated smoke path
- burden comparisons for RNA-derived candidate variant signals
- candidate_snv downstream association with RNA/circ features

## Limitations

- Homogeneous datasets or one-cell-type annotations may produce no Step2 candidate calls.
- Without cell types, SComatic grouping is technical rather than biological.
- Without a Panel of Normals, output must remain exploratory.
- RNA-derived signals are affected by expression, allelic expression, RNA editing, mapping artifacts, and transcript coverage.
- `candidate_snv` is an interoperability modality, not a validated DNA mutation modality.

## Protocol Examples

Smart-seq3:

```bash
circyto prepare-scomatic-input \
  --merged-bam smartseq3.pooled.bam \
  --cell-metadata smartseq3_cells.tsv \
  --protocol smartseq3 \
  --outdir smartseq3_scomatic_prepared
```

RamDA:

```bash
circyto prepare-scomatic-input \
  --alignment-manifest ramda_alignment_manifest.tsv \
  --protocol ramda \
  --reference-fasta hg38.fa \
  --outdir ramda_scomatic_prepared
```

ShinRamDA:

```bash
circyto prepare-scomatic-input \
  --alignment-manifest shinramda_alignment_manifest.tsv \
  --protocol shinramda \
  --reference-fasta hg38.fa \
  --outdir shinramda_scomatic_prepared
```

scRR RNA arm:

```bash
circyto prepare-scomatic-input \
  --alignment-manifest scrr_alignment_manifest.tsv \
  --protocol scrr-rna \
  --reference-fasta hg38.fa \
  --outdir scrr_scomatic_prepared
```
