
# Detectors

`circyto` is designed to work with multiple circRNA detectors through small, focused adapter scripts. The goal is to provide a consistent CLI and output schema across tools.

## Detector comparison

| Detector        | Input type              | Reference needed                   | Output granularity        | Strengths                                                     | Limitations / Notes                                           | circyto support level |
|----------------|-------------------------|------------------------------------|---------------------------|----------------------------------------------------------------|---------------------------------------------------------------|------------------------|
| **CIRI-full**  | Paired-end FASTQ        | Genome FASTA + BWA index + GTF     | Full-length circRNA + AS  | Full-length reconstruction, integrates RO/AS, detailed output | Heavy, multi-stage pipeline; requires good reference indices | ✅ Fully integrated (chr21 Smart-seq2 example) |
| **CIRI-long**  | Long-read FASTQ (ONT)   | Genome FASTA + minimap2 index      | Full-length circRNA       | Designed for long-read; captures complex isoforms             | Long-read only; runtime and memory depend on read length     | 🔜 Planned integration |
| **CIRCexplorer2** | BAM (spliced align.) | Genome FASTA + annotation          | Back-splice junctions     | Widely used; works from alignment BAMs                        | No full-length reconstruction; depends on upstream aligner   | 🔜 Planned integration |
| **find_circ**  | FASTQ (or BAM)          | Genome FASTA + BWA index           | Back-splice junctions     | Simple, fast; classic circRNA detector                        | Older; higher FP rate; no full-length info                   | 🔜 Planned integration |
| **circRNA_finder** (optional) | BAM    | Genome FASTA                       | Back-splice junctions     | Works directly on STAR outputs                                | STAR-specific; less widely maintained                        | ❓ Maybe (low priority) |

Support levels:

- ✅ Fully integrated: wired through `run`/`run-manifest` + `collect` and tested end-to-end.
- 🔜 Planned: interface design ready; implementation pending.
- ❓ Maybe: candidate for future integration depending on demand.

## Adapter philosophy

Each detector is wrapped by a small adapter that is responsible for:

- Mapping manifest columns (R1/R2/BAM) to detector CLI arguments
- Handling compression/decompression
- Managing per-cell working directories
- Normalizing detector output to a common TSV schema:

```text
circ_id    chr    start    end    strand    support
```

Upstream complexity stays inside the detector; `circyto` focuses on orchestration and harmonization.

