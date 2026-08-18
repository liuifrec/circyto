from __future__ import annotations

from typer.testing import CliRunner

from circyto.cli.doctor import doctor_app
from circyto.paths import PathResolution


runner = CliRunner()


def test_doctor_separates_core_readiness_from_optional_tools(monkeypatch) -> None:
    monkeypatch.setattr("circyto.cli.doctor.metadata.version", lambda _name: "1.2.3")
    monkeypatch.setattr("circyto.cli.doctor.shutil.which", lambda _name: None)
    monkeypatch.setattr(
        "circyto.cli.doctor.resolve_asset",
        lambda label: PathResolution(label, None, (), None),
    )
    monkeypatch.setattr(
        "circyto.cli.doctor.inspect_ciri3_runtime",
        lambda: {
            "jar": None,
            "java": None,
            "java_version": None,
            "java_version_ok": False,
            "required_java_major": 12,
            "recommended_java_major": 17,
        },
    )
    monkeypatch.setattr(
        "circyto.cli.doctor.detector_runtime_status",
        lambda _name: {"status": "NOT READY", "reason": "optional test dependency missing"},
    )
    monkeypatch.setattr(
        "circyto.cli.doctor.check_ciri_long_readiness",
        lambda: {
            "ok": False,
            "errors": ["CIRI-long executable not found on PATH"],
            "tools": {"ciri_long": None, "bwa": None},
        },
    )

    result = runner.invoke(doctor_app, [])

    assert result.exit_code == 0, result.output
    assert "Required core Python package readiness" in result.output
    assert "Bundled/configured detector assets" in result.output
    assert "Optional external executables" in result.output
    assert "Detector-specific readiness" in result.output
    assert "Protocol-specific limitations" in result.output
    assert "[OPT ] minimap2" in result.output
    assert "[OPT ] CIRI-long" in result.output
    assert "does not indicate a circyto package failure" in result.output
    assert "No network checks were performed" in result.output


def test_doctor_fails_only_for_missing_core_packages(monkeypatch) -> None:
    def fake_version(name: str) -> str:
        if name == "numpy":
            from importlib.metadata import PackageNotFoundError

            raise PackageNotFoundError(name)
        return "1.2.3"

    monkeypatch.setattr("circyto.cli.doctor.metadata.version", fake_version)
    monkeypatch.setattr("circyto.cli.doctor.shutil.which", lambda _name: None)
    monkeypatch.setattr(
        "circyto.cli.doctor.resolve_asset",
        lambda label: PathResolution(label, None, (), None),
    )
    monkeypatch.setattr(
        "circyto.cli.doctor.inspect_ciri3_runtime",
        lambda: {
            "jar": None,
            "java": None,
            "java_version": None,
            "java_version_ok": False,
            "required_java_major": 12,
            "recommended_java_major": 17,
        },
    )
    monkeypatch.setattr(
        "circyto.cli.doctor.detector_runtime_status",
        lambda _name: {"status": "NOT READY", "reason": "optional test dependency missing"},
    )
    monkeypatch.setattr(
        "circyto.cli.doctor.check_ciri_long_readiness",
        lambda: {
            "ok": False,
            "errors": ["CIRI-long executable not found on PATH"],
            "tools": {"ciri_long": None, "bwa": None},
        },
    )

    result = runner.invoke(doctor_app, [])

    assert result.exit_code == 1
    assert "[MISS] numpy" in result.output
    assert "Core Python environment is NOT ready" in result.output
