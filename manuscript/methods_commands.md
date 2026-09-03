# Methods Commands

These command shapes document how manuscript-scale objects are expected to be regenerated. Paths are placeholders unless explicitly described as committed repo assets.

For the current Application Note, the Smart-seq3 commands support the primary
demonstration and the IMR90 commands support a multimodal object/schema
example. The HAP1 RT analysis and IMR90 CNV association scripts below are
preserved as **deferred biological follow-up**, not as main-paper requirements.
Real HAP1 processed-file RT validation must be reconciled before those outputs
are used as evidence.

## Smart-seq3 / E-MTAB-8735

```bash
circyto workflow smartseq3-ciri3 \
  --read1 /path/to/data/E-MTAB-8735/Smartseq3.R1.fastq.gz \
  --read2 /path/to/data/E-MTAB-8735/Smartseq3.R2.fastq.gz \
  --index1 /path/to/data/E-MTAB-8735/Smartseq3.I1.fastq.gz \
  --index2 /path/to/data/E-MTAB-8735/Smartseq3.I2.fastq.gz \
  --annotation /path/to/data/E-MTAB-8735/sample_annotation.tsv \
  --cell-id-column cell_id \
  --index1-column index1 \
  --index2-column index2 \
  --ref-fa /path/to/ref/hg38.fa \
  --gtf /path/to/ref/gencode.v38.annotation.gtf \
  --star-index /path/to/ref/star_index_hg38 \
  --gene-counts /path/to/data/E-MTAB-8735/gene_counts.tsv \
  --gene-counts-format tsv \
  --outdir /path/to/work/mtab8735_smartseq3 \
  --threads 24 \
  --resume \
  --export-h5ad \
  --export-mudata
```

Repair host-gene annotations:

```bash
circyto repair-host-genes \
  --input /path/to/work/mtab8735_smartseq3/mudata/circyto_multimodal.h5mu \
  --output /path/to/work/mtab8735_smartseq3/full_length.hostgene_fixed.h5mu \
  --circ-mod circ \
  --gtf /path/to/ref/gencode.v38.annotation.gtf
```

## HAP1 RNA+circRNA+RT

**Deferred biological follow-up:** the RNA/circ route is real-data validated,
but the RT commands below currently document the integration contract. They do
not establish a completed real-file RT-circRNA analysis.

```bash
circyto workflow full-length-circrna \
  --manifest /path/to/data/HAP1/hap1_rna_manifest.tsv \
  --outdir /path/to/work/scrr_hap1/full_length \
  --protocol ramda \
  --genome-fasta /path/to/ref/hg38.fa \
  --gtf /path/to/ref/gencode.v38.annotation.gtf \
  --gene-counts /path/to/data/HAP1/gene_counts.tsv \
  --export-h5ad
```

```bash
circyto build-scrr-cell-map \
  --soft /path/to/data/GSE278952_family.soft.gz \
  --out /path/to/work/scrr_hap1/scrr_cell_map.tsv
```

```bash
circyto remap-scrr-mudata-obs \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ.h5mu \
  --cell-map /path/to/work/scrr_hap1/scrr_cell_map.tsv \
  --output /path/to/work/scrr_hap1/mudata/full_length_rna_circ.remapped.h5mu
```

```bash
circyto import-scrr-rt \
  --rt-table /path/to/data/GSE278952/HAP1_binarized_rt_state_hg38.txt.gz \
  --avg-rt-bedgraph /path/to/data/GSE278952/HAP1_midS_Avg_RT_hg38.bedGraph.gz \
  --outdir /path/to/work/scrr_hap1/rt
```

```bash
circyto merge-scrr-rt \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ.remapped.h5mu \
  --rt /path/to/work/scrr_hap1/rt/rt.h5ad \
  --output /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.h5mu
```

```bash
circyto repair-host-genes \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.h5mu \
  --output /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu \
  --circ-mod circ \
  --gtf /path/to/ref/gencode.v38.annotation.gtf
```

## IMR90 RNA+circRNA+CNV

**Supplementary interoperability only:** use these commands to demonstrate the
RNA+circ+processed-CNV object contract. CNV association or host-gene biology is
deferred biological follow-up.

```bash
circyto build-scrr-cell-map \
  --soft /path/to/data/GSE278958_family.soft.gz \
  --out /path/to/work/scrr_imr90/scrr_cell_map.tsv
```

```bash
circyto import-scrr-cnv \
  --cnv-states /path/to/data/IMR90/summary_CNV_states_bin50kb.txt.gz \
  --cnv-mappabilitynorm /path/to/data/IMR90/summary_CNV_mappabilitynorm_bin50kb.txt.gz \
  --outdir /path/to/work/scrr_imr90/cnv
```

```bash
circyto remap-scrr-mudata-obs \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ.h5mu \
  --cell-map /path/to/work/scrr_imr90/scrr_cell_map.tsv \
  --output /path/to/work/scrr_imr90/mudata/full_length_rna_circ.remapped.h5mu
```

```bash
circyto merge-scrr-cnv \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ.remapped.h5mu \
  --cnv /path/to/work/scrr_imr90/cnv/cnv.h5ad \
  --output /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.h5mu
```

```bash
circyto repair-host-genes \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.h5mu \
  --output /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu \
  --circ-mod circ \
  --gtf /path/to/ref/gencode.v38.annotation.gtf
```

## Manuscript Summary Scripts

The inventory command is current manuscript support. HAP1 RT correlations and
IMR90 CNV association outputs are deferred biological follow-up.

```bash
python scripts/manuscript/summarize_mudata_inventory.py \
  /path/to/work/mtab8735_smartseq3/full_length.hostgene_fixed.h5mu \
  /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu \
  /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu \
  --dataset-name Smart-seq3_E-MTAB-8735 \
  --dataset-name HAP1_scRR_RT \
  --dataset-name IMR90_scRR_CNV \
  --out /path/to/results/manuscript/mudata_inventory.tsv
```

```bash
python scripts/manuscript/hap1_rt_circ_analysis.py \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu \
  --correlations-out /path/to/results/manuscript/hap1_rt_circ_correlations.tsv \
  --ols-out /path/to/results/manuscript/hap1_rt_circ_ols.tsv \
  --scatter-out /path/to/results/manuscript/hap1_rt_circ_cell_metrics.tsv
```

```bash
python scripts/manuscript/imr90_cnv_circ_analysis.py \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu \
  --correlations-out /path/to/results/manuscript/imr90_cnv_circ_correlations.tsv \
  --enrichment-out /path/to/results/manuscript/imr90_cnv_high_host_genes.tsv \
  --cell-metrics-out /path/to/results/manuscript/imr90_cnv_circ_cell_metrics.tsv \
  --local-cnv-out /path/to/results/manuscript/imr90_local_cnv_at_circ_loci.tsv
```
