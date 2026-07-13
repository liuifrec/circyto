from __future__ import annotations

import json

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def test_demo_mini_writes_self_contained_matrix(tmp_path) -> None:
    outdir = tmp_path / "demo_out"
    result = runner.invoke(app, ["demo", "mini", "--out", str(outdir)])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["external_tools_required"] is False
    assert payload["n_circs"] == 2
    assert payload["n_cells"] == 2
    assert payload["nnz"] == 3
    assert (outdir / "matrix" / "circ_counts.mtx").exists()
    assert (outdir / "matrix" / "circ_feature_table.tsv").exists()
