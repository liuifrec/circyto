import json
from pathlib import Path

import pytest

from output_audit import archive_payload, check, exact_matrix, Paragraphs
from reference import convert
from test_reference import fixture
from robustness import archive


def test_exact_matrix_large_integer_and_invalid_schema():
    data = b"%%MatrixMarket matrix coordinate integer general\n%\n1 2 2\n1 1 9007199254740993\n1 2 0\n"
    assert exact_matrix(data, ["001", "NA"], ["f"]) == {("001", "f"): 2**53+1}
    with pytest.raises(ValueError):
        exact_matrix(data.replace(b"integer", b"real"), ["001", "NA"], ["f"])
    with pytest.raises(ValueError):
        exact_matrix(data.replace(b"1 2 0", b"1 1 0"), ["001", "NA"], ["f"])


def test_summary_and_separate_matrix_corruption_are_detected(tmp_path):
    config = fixture(tmp_path)
    out = tmp_path / "out"
    convert(config, out)
    source = lambda p: Path(p).read_bytes()
    output = lambda p: (out / p).read_bytes()
    assert all(check(config, source, output, "reference").values())
    html = out / "qc/report.html"
    html.write_text(html.read_text().replace("Completed denominator 2", "Completed denominator 1"))
    result = check(config, source, output, "reference")
    assert not result["html_summary_counts_and_denominator"]
    assert result["exported_matrix_all_integer_values"]
    matrix = out / "matrix/circ_counts.mtx"
    matrix.write_text(matrix.read_text().replace("1 1 7", "1 1 8"))
    assert not check(config, source, output, "reference")["exported_matrix_all_integer_values"]
    metrics = out / "qc/metrics.json"
    m = json.loads(metrics.read_text())
    m["rna_overlap"]["n_shared"] = 1
    metrics.write_text(json.dumps(m))
    assert not check(config, source, output, "reference")["json_RNA_overlap"]


def test_html_parser_preserves_nested_paragraph_text():
    p = Paragraphs()
    p.feed('<p>Fraction: <strong>1 / 2</strong>.</p><script src="https://invalid.example/x"></script>')
    assert p.paragraphs == ["Fraction: 1 / 2."]
    assert p.remote_resources == ["https://invalid.example/x"]


def test_archived_and_live_output_audit_are_equivalent(tmp_path):
    case = tmp_path / "case"
    inputs = case / "inputs"
    inputs.mkdir(parents=True)
    config = fixture(inputs)
    (inputs / "config.json").write_text(json.dumps(config))
    out = case / "reference/output"
    convert(config, out)
    live = check(config, lambda p: Path(p).read_bytes(), lambda p: (out / p).read_bytes(), "reference")
    retained = archive(case)
    payload = archive_payload(Path(retained["archive"]))
    source = lambda p: payload["inputs/"+str(Path(p).relative_to(inputs))]
    output = lambda p: payload["reference/output/"+p]
    assert check(config, source, output, "reference") == live
    assert all(live.values())
