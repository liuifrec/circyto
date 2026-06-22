# Host-Gene Provenance

Host-gene annotation in `circyto` is intentionally provenance-aware. The final `host_gene` field is a display and analysis convenience; source-specific fields preserve how that value was chosen.

## Priority Order

1. GTF/GFF genomic overlap.
2. circAtlas table host gene if available.
3. circAtlas ID parsing fallback.

## Expected Fields

- `host_gene`: final selected host gene.
- `host_gene_source`: source of the final selected host gene.
- `host_gene_from_gtf`: host gene from GTF/GFF overlap.
- `host_gene_from_circatlas`: host gene from circAtlas table fields.
- `host_gene_from_circatlas_id`: host gene parsed from circAtlas-style IDs.
- `host_gene_id`: selected host gene ID when available.
- `host_genes_multi`: all overlapping or candidate host genes when multiple are present.
- `host_gene_ids_multi`: all corresponding host gene IDs when available.
- `host_gene_n`: number of candidate host genes.

Expected `host_gene_source` values:

- `gtf`
- `circatlas`
- `circatlas_id`
- empty or missing when no supported source is available

## Repair Existing Objects

Use `repair-host-genes` for existing `.h5ad` or `.h5mu` outputs that lack complete host-gene fields or need GTF/GFF overlap applied.

```bash
circyto repair-host-genes \
  --input /path/to/input.h5mu \
  --output /path/to/output.hostgene_fixed.h5mu \
  --circ-mod circ \
  --gtf /path/to/ref/gencode.v38.annotation.gtf
```

For AnnData circRNA objects:

```bash
circyto repair-host-genes \
  --input /path/to/circ_counts.h5ad \
  --output /path/to/circ_counts.hostgene_fixed.h5ad \
  --gtf /path/to/ref/gencode.v38.annotation.gtf
```

## Manuscript Summary

Use the inventory script to summarize recovery and source counts:

```bash
python scripts/manuscript/summarize_mudata_inventory.py \
  /path/to/full_length.hostgene_fixed.h5mu \
  --dataset-name dataset_label \
  --out /path/to/results/host_gene_inventory.tsv
```
