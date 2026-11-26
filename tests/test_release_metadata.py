"""Release metadata consistency checks, independent of PyPI availability."""
from importlib.metadata import version
from pathlib import Path
import re

import circyto


ROOT = Path(__file__).resolve().parents[1]


def test_project_cli_and_citation_versions_agree():
    project = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    project_match = re.search(r'(?m)^version\s*=\s*"([^"]+)"\s*$', project)
    citation_match = re.search(r"(?m)^version:\s*([^\s]+)\s*$", citation)
    assert project_match and citation_match
    assert project_match.group(1) == circyto.__version__
    assert citation_match.group(1) == circyto.__version__
    assert version("circyto") == circyto.__version__


def test_project_citation_uses_publication_name():
    project = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert 'name = "Yu-Chen Liu"' in project
    assert "given-names: Yu-Chen" in citation
    assert "family-names: Liu" in citation
