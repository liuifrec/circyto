# Public RamDA-family dataset helpers

This page is a lightweight reference for public RamDA-family datasets that may be useful when planning future `circyto` smoke tests, benchmarks, or protocol-specific examples.

It is intentionally limited to metadata templates. It does not add download automation.

## Current recommendations

- For the first CIRI3 smoke test, prefer a small protocol-aware RamDA/Shin-RamDA example with tiny manifests and a chr21-scale reference subset.
- Do not require scRepli-RamDA-seq / scRR-seq for that first smoke test.
- scRepli-RamDA-seq is especially valuable later for DNA/RNA/circRNA correlation work because it couples genome-wide replication profiling with full-length total RNA from the same single cell.

## Referenced public sources

- Original RamDA-seq GEO series: `GSE98664`
  Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98664
- scRepli-RamDA-seq / scRR-seq related GEO accessions discussed here:
  `GSE278944`, `GSE278941`, `GSE278951`, `GSE278949`, `GSE278958`, `GSE278952`
- Shin-RamDA-seq workflow reference:
  `rikenbit/shin-ramda-seq-paper`
  Source: https://github.com/rikenbit/shin-ramda-seq-paper

## Notes by protocol

### RamDA-seq

`GSE98664` is the original RamDA-seq GEO series. GEO describes it as single-cell full-length total RNA sequencing and reports mouse ES / PrE-related material across multiple experimental settings. This makes it useful as protocol background and for understanding expected full-length total-RNA behavior, though it is not the smallest first-pass CIRI3 smoke target.

### scRepli-RamDA-seq / scRR-seq

The scRR family is the most compelling future direction for integrated analyses because the public GEO records describe genome-wide DNA replication profiling and full-length total RNA sequencing from the same single cell. That is a good fit for eventual DNA/RNA/circRNA correlation analyses, but it should remain optional for the initial CIRI3 smoke workflow.

### Shin-RamDA-seq

The RIKEN reference repository is useful even without bundling large data files because it exposes realistic analysis scopes:

- `HEK293T_benchmark`
- `Detection_of_genomic_breakpoints`
- `Lineage_estimation_using_etoposide_K562`

Those are better treated as future benchmark and workflow-shape references than as immediate smoke-test inputs.
