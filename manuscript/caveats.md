# Caveats and required wording boundaries

- circyto bridges detector output and single-cell/scverse analysis; it is not a
  new back-splice detection algorithm.
- Detector-derived features should be called circRNA candidates or detected
  circRNAs unless orthogonal validation supports stronger language.
- The Smart-seq3 MAN1A2-associated candidate is an illustrative overlay, not a
  functional or disease association.
- IMR90 CNV is imported from processed GEO summaries. The 23-cell object shows
  multimodal interoperability and does not support broad CNV biology claims.
- The validated HAP1 evidence is the real-data paired-end RNA/circRNA route.
  Real processed-file RT validation is unresolved in the current evidence
  base, so no completed RT-circRNA discovery analysis should be claimed.
- SComatic-derived layers are exploratory RNA-derived candidate variant signals,
  not orthogonally confirmed DNA variants.
- Generic ONT cDNA interoperability is alignment/QC/provenance evidence only;
  `circRNA_call=false` and non-detection is not biological absence.
- The official CIRI-long demonstration is bulk, chemistry-gated adapter
  interoperability and does not establish single-cell performance or
  biological accuracy.
- Biogenesis exports implement positive/unlabelled schemas and provenance; no
  predictive biogenesis model or true-negative semantics exist.
- MuData compatibility preserves circyto's historical synchronization policy;
  broader AnnData 0.13 string-dtype migration is deferred.
- Public HAP1 and IMR90 data are not radiation-exposure cohorts.

HAP1 RT association, IMR90 CNV-burden biology, local CNV-locus inference,
cross-dataset host-gene programs, and radiation applications are deferred
biological follow-up rather than Application Note requirements.
