# Circyto repo architecture (canonical)

## Entry points
- circyto/cli/circyto.py : main Typer app (registers subcommands)
- circyto/cli/smoke.py   : one-command smoke flows (demux -> run -> collect -> convert)

## Demux
- circyto/demux/smartseq2.py : Smart-seq2 pooled demux implementation
  - SmartSeq2DemuxParams
  - demux_smartseq2_pooled()

## Pipelines
- circyto/pipeline/run_detector.py : read_manifest() + run_detector_manifest()
- circyto/pipeline/collect*.py     : detector-specific collectors
- circyto/writers/convert.py       : convert_matrix_files() -> h5ad/loom
- circyto/pipeline/export_multimodal.py : attach circRNA matrix to AnnData

## Detectors
- circyto/detectors/ : detector engines + adapters (find-circ3, ciri2, ciri-full, ...)
