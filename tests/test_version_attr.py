# -*- coding: utf-8 -*-
"""``abmptools.__version__`` exists and agrees with pyproject.toml and CITATION.cff.

It was missing: ``abmptools.__version__`` raised AttributeError. The version
lives in pyproject.toml; an installed wheel reads it from its metadata, a
source tree from pyproject.toml (an editable install's metadata goes stale
after a version bump).
"""
from __future__ import annotations

import pathlib
import re
from importlib import metadata

import abmptools
from abmptools import _version

ROOT = pathlib.Path(__file__).resolve().parent.parent


def _field(path, pattern):
    m = re.search(pattern, (ROOT / path).read_text(encoding="utf-8"), re.M)
    assert m, path
    return m.group(1)


def test_version_matches_pyproject_and_citation():
    want = _field("pyproject.toml", r'^version = "([^"]+)"')
    assert abmptools.__version__ == want
    assert _field("CITATION.cff", r'^version: "([^"]+)"') == want


def test_installed_copy_reads_metadata(tmp_path, monkeypatch):
    """No pyproject.toml beside the package (a wheel): the metadata's version."""
    pkg = tmp_path / "abmptools"
    pkg.mkdir()
    monkeypatch.setattr(metadata, "version", lambda name: "9.8.7")
    assert _version.get_version(pkg) == "9.8.7"


def test_another_projects_pyproject_is_not_read(tmp_path, monkeypatch):
    (tmp_path / "pyproject.toml").write_text('[project]\nname = "other"\nversion = "1.0"\n')
    pkg = tmp_path / "abmptools"
    pkg.mkdir()
    monkeypatch.setattr(metadata, "version", lambda name: "9.8.7")
    assert _version.get_version(pkg) == "9.8.7"


def test_not_installed_says_so(tmp_path, monkeypatch):
    def missing(name):
        raise metadata.PackageNotFoundError(name)
    monkeypatch.setattr(metadata, "version", missing)
    assert _version.get_version(tmp_path / "abmptools") == _version.UNKNOWN
