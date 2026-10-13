"""Packaging guard: fail closed if stale build assets bypass manifest exclusions."""
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist
from setuptools.errors import SetupError


def check_distribution_tree(root):
    forbidden = [p for p in Path(root).rglob("*") if p.is_file() and
                 (p.suffix.lower() in {".jar", ".class"} or
                  p.name == "CIRI-full_manual.txt")]
    if forbidden:
        raise SetupError(
            "Excluded detector assets found in build output: "
            + ", ".join(str(p) for p in forbidden)
            + ". Build from a clean source tree or remove the stale build directory; "
            "detector JARs must be installed externally."
        )


class ExternalDetectorBuild(build_py):
    def run(self):
        super().run()
        check_distribution_tree(self.build_lib)


class ExternalDetectorSdist(sdist):
    def make_release_tree(self, base_dir, files):
        super().make_release_tree(base_dir, files)
        check_distribution_tree(base_dir)


setup(cmdclass={"build_py": ExternalDetectorBuild, "sdist": ExternalDetectorSdist})
