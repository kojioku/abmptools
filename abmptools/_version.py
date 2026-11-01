# -*- coding: utf-8 -*-
"""``abmptools.__version__``.

The version is declared once, in ``pyproject.toml``; an installed wheel carries
it in its metadata, so that is where it is read from. A source checkout is the
exception: an editable install keeps the metadata written when it was
installed, and after a version bump ``importlib.metadata`` goes on reporting
the old number (seen: 2.9.0 from a 2.19 tree). When the package is imported
from a tree that has its ``pyproject.toml`` next to it, that file is read
instead.
"""
from __future__ import annotations

import pathlib
import re
from importlib import metadata

UNKNOWN = "0+unknown"


def _from_source_tree(pkg_dir: pathlib.Path):
    pyproject = pkg_dir.parent / "pyproject.toml"
    try:
        text = pyproject.read_text(encoding="utf-8")
    except OSError:
        return None
    if not re.search(r'^name = "abmptools"', text, re.M):
        return None
    m = re.search(r'^version = "([^"]+)"', text, re.M)
    return m.group(1) if m else None


def get_version(pkg_dir: pathlib.Path = pathlib.Path(__file__).resolve().parent) -> str:
    v = _from_source_tree(pkg_dir)
    if v:
        return v
    try:
        return metadata.version("abmptools")
    except metadata.PackageNotFoundError:
        return UNKNOWN
