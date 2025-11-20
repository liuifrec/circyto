def test_multidetector_smoke(tmp_path):
    from circyto.pipeline.run_multidetector import run_multidetector_pipeline
    from pathlib import Path

    # minimal fake environment
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("cell_id\tr1\tr2\nc1\tx\ty\n")

    out = tmp_path / "multi"
    out.mkdir()

    meta = run_multidetector_pipeline(
        detectors=["ciri-full"],
        manifest=manifest,
        outdir=out,
        ref_fa=None,
        gtf=None,
    )

    assert "ciri-full" in meta
