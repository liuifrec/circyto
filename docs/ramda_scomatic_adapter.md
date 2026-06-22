# RamDA/scRR SComatic Adapter

This note documents the local adapter that prepares one-cell-per-SRR RamDA/scRR alignments for SComatic's pooled single-cell BAM input. It does not run SComatic.

SComatic outputs from this path should be described as RNA-derived candidate variant signals unless orthogonal DNA validation exists.

## SComatic Input Contract

Local inspection of `/mnt/d/SComatic/README.md` and `/mnt/d/SComatic/scripts/SplitBam/SplitBamCellTypes.py` shows:

- Step 1 expects one coordinate-sorted BAM containing reads from all cells.
- The BAM must carry cell barcode information in the `CB` tag. `SplitBamCellTypes.py` reads `read.opt("CB")` and drops reads without it.
- `nM` and `NH` are strongly recommended by SComatic and are required when running SComatic with `--max_nM` or `--max_NH`. The adapter preserves these tags and reports missing counts, but it does not infer them.
- The cell annotation file is tab-separated and must contain at least:

```text
Index	Cell_type
cellA	cellA
cellB	cellB
```

`Index` is the cell barcode matched against the BAM `CB` tag. `Cell_type` is the SComatic grouping label. For one-cell-per-SRR RamDA/scRR inputs, the adapter defaults `Cell_type` to the `cell_id`, preserving one-cell-per-library granularity.

## Adapter Command

```bash
bash scripts/prepare_ramda_scomatic_input.sh \
  --alignment-manifest WORKDIR/align/alignment_manifest.tsv \
  --outdir WORKDIR/scomatic_adapter \
  --sample-id ramda_scrr
```

Outputs:

- `per_cell/*.scomatic.sorted.bam`
- `merged/<sample-id>.scomatic.bam`
- `merged/<sample-id>.scomatic.bam.bai`
- `cell_annotations.tsv`
- `adapter_summary.json`

Optional grouping:

```bash
bash scripts/prepare_ramda_scomatic_input.sh \
  --alignment-manifest WORKDIR/align/alignment_manifest.tsv \
  --outdir WORKDIR/scomatic_adapter \
  --sample-id ramda_scrr \
  --cell-type-column cell_type
```

If `--cell-type-column` is not supplied, `Cell_type` defaults to `cell_id`. `--default-cell-type LABEL` can be used to place all cells in one SComatic group.

## Implementation Notes

The adapter:

1. Reads circyto `alignment_manifest.tsv`.
2. Takes `cell_id` as the synthetic single-cell barcode.
3. Streams each SAM/BAM through `samtools view -h`.
4. Adds `CB:Z:<cell_id>` to alignment records that do not already have a `CB` tag, and normalizes non-matching `CB` tags to the manifest `cell_id`.
5. Converts and coordinate-sorts each per-cell alignment with `samtools sort`.
6. Merges sorted per-cell BAMs into one pooled BAM.
7. Indexes the merged BAM.
8. Writes SComatic-compatible `Index`/`Cell_type` annotations.
9. Writes `adapter_summary.json`, including missing/replaced `CB`, missing `nM`, and missing `NH` counts.

The adapter intentionally does not synthesize `nM` or `NH`. For aligners that only emit `NM`, choose SComatic parameters accordingly or add a separate, explicit tag-normalization step after validation.

## Local Safety

Do not run real SComatic locally in WSL for this project state. Real SComatic execution remains deferred to HPC, a server, or a container because native WSL dependency installation repeatedly failed with `Bus error (core dumped)`.
