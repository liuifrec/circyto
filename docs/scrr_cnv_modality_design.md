# scRR CNV Modality Design

## Goal

Represent processed scRR-seq DNA CNV outputs in a scverse-compatible form that can be joined to existing `rna` and `circ` modalities without rerunning DNA workflows.

## Options

## Option A: Cell x Genomic-Bin CNV Matrix

Structure:

- `obs`: cells
- `var`: genomic bins
- `X`: integer copy number parsed from upstream CNV state labels
- `layers["mappabilitynorm"]`: optional numeric mappability-normalized DNA signal

Pros:

- matches GEO processed CNV tables directly
- preserves genomic resolution
- compatible with AnnData and MuData
- supports per-bin heatmaps, chromosome summaries, and per-cell CNV burden
- keeps provenance simple because each feature is a source-table bin

Cons:

- dense by nature; most bins have a state in every cell
- multiple resolutions need separate imports or separate modalities because feature axes differ

## Option B: Cell x CNV-Segment Matrix

Structure:

- features are inferred CNV segments instead of fixed bins

Pros:

- more compact for broad arm-level CNVs
- closer to biological events when robust segmentation is available

Cons:

- GEO tables inspected here are bin matrices, not segment tables
- segment generation would require extra algorithmic choices inside circyto
- per-cell segment boundaries do not naturally share one feature axis

## Option C: Store CNV Tracks In `uns`

Structure:

- keep CNV tables as nested metadata under `mdata.uns`

Pros:

- easiest to store arbitrary external files
- avoids AnnData matrix decisions

Cons:

- not a first-class modality
- hard to analyze with scverse tools
- weak validation of cell and bin axes
- easy to lose coordinate-level provenance

## Recommendation

Use Option A as the primary representation:

```text
MuData
|- rna
|- circ
`- cnv
```

`mdata.mod["cnv"]` should be an AnnData object:

- `cnv.X`: `int16` cell x bin copy-number matrix parsed from `N-somy`
- `cnv.obs_names`: stable scRR cell IDs
- `cnv.obs["canonical_cell_id"]`: modality-prefix-free physical-cell ID
- `cnv.obs["dna_cell_id"]`: original GEO DNA sample ID
- `cnv.obs["rna_cell_id"]`: paired RNA sample ID when inferable or mapped
- `cnv.var_names`: `seqname:start-end`
- `cnv.var`: `seqname`, `start`, `end`, `bin_size`
- `cnv.layers["mappabilitynorm"]`: optional numeric DNA signal
- `cnv.uns["circyto"]`: source filenames, accession, command, bin size, and state semantics

## Cell Axis Strategy

The preferred shared cell key is the canonical scRR physical-cell ID:

```text
DNA_IMR90_A_100 -> IMR90_A_100
RNA_IMR90_A_100 -> IMR90_A_100
```

For existing circyto runs whose RNA/circ cell IDs are already `RNA_*`, future MuData export should either:

- normalize RNA/circ obs names to the canonical ID, or
- join through `canonical_cell_id` while preserving original RNA IDs in metadata

If RNA cell IDs are SRR accessions or other aliases, require an explicit mapping table.

## Multiple Resolutions

Do not mix 50 kb, 100 kb, 200 kb, and 500 kb bins in one AnnData layer because the feature axes differ.

Recommended handling:

- one selected resolution as `mdata.mod["cnv"]`
- optional additional resolutions as `mdata.mod["cnv_100kb"]`, `mdata.mod["cnv_200kb"]`, etc. only when needed
- preserve `bin_size` in `var` and provenance

## Future MuData Export

The current implementation adds a parser command that can emit `cnv.h5ad`. Full integration into `export-mudata` should be a separate step:

1. detect `WORKDIR/dna/cnv.h5ad`
2. read the CNV AnnData
3. align cells by `canonical_cell_id`
4. emit `MuData({"rna": rna, "circ": circ, "cnv": cnv})`
5. preserve original RNA, circ, and DNA identifiers in `mdata.obs`

## Non-goals

- circyto should not rerun scRepliseq CNV calling in this integration layer.
- circyto should not segment CNVs unless a future, explicit segmentation command is designed.
- circyto should not represent RNA-derived SComatic signals as validated DNA mutations.
