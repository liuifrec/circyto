# MuData Schema

`circyto` uses MuData to keep full-length single-cell RNA, circRNA, and genome-state modalities in one scverse-compatible object.

## Expected Layouts

RNA+circRNA:

```text
MuData
|- rna
`- circ
```

HAP1 RNA+circRNA+RT:

```text
MuData
|- rna
|- circ
`- rt
```

IMR90 RNA+circRNA+CNV:

```text
MuData
|- rna
|- circ
`- cnv
```

Optional exploratory candidate signal layout:

```text
MuData
|- rna
|- circ
`- candidate_snv
```

## Cell Axis

For manuscript analyses, modality observations should use a shared biological-cell namespace whenever possible.

Expected overlap checks:

- RNA/circRNA overlap for full-length RNA+circRNA objects.
- RNA/circRNA/RT overlap for HAP1 RT analyses.
- RNA/circRNA/CNV overlap for IMR90 CNV analyses.

When an imported DNA-derived modality uses DNA sample IDs, run the scRR cell-map and remapping commands before merging.

## Provenance

Expected high-level provenance location:

```text
mdata.uns["circyto"]
```

Common fields include:

- command name
- `circyto_version`
- source workdir or input paths
- workflow summary JSON
- RNA import summary JSON
- RNA/circRNA summary JSON
- modality-specific merge summaries

## Manuscript Metrics

Scripts under `scripts/manuscript/` first look for metrics in `mdata.obs` or modality `obs`, then compute supported fallbacks from `X`.

Supported fallback metrics:

- `circRNA_count`: nonzero circRNAs per cell.
- `circRNA_total_support`: total circRNA support per cell.
- `detected_genes`: nonzero genes per cell.
- `frac_rt_pos`: fraction of RT features with positive state/value.
- `frac_rt_neg`: fraction of RT features with negative state/value, or zero state for binary 0/1 RT matrices.
- `frac_non_diploid`: fraction of CNV bins with state not equal to 2.
- `frac_loss`: fraction of CNV bins with state below 2.
- `frac_gain`: fraction of CNV bins with state above 2.
