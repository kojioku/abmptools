# -*- coding: utf-8 -*-
"""`CITATION.cff` の版が `pyproject.toml` からずれないこと。

版は 2 か所にあり、`CITATION.cff` は release のたびに手で直す運用だった。
実際にずれる —— v2.13.0 の時点で 2.12.0 のまま放置されていたのを直したあと、
今度は 2.13.9 を出したのに 2.13.8 のまま残っていた。

これは黙って間違う類の記載である。**引用されるのは論文の中**なので、下流の
検査はどこにも無く、読者が古い版を引く形で表に出る。
"""
from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
CITATION = ROOT / "CITATION.cff"
PYPROJECT = ROOT / "pyproject.toml"


def _pyproject_version() -> str:
    m = re.search(r'^version = "([^"]+)"',
                  PYPROJECT.read_text(encoding="utf-8"), re.M)
    assert m, "pyproject.toml に version が無い"
    return m.group(1)


def _citation_field(name: str) -> str:
    m = re.search(rf'^{name}: "([^"]+)"',
                  CITATION.read_text(encoding="utf-8"), re.M)
    assert m, f"CITATION.cff に {name} が無い"
    return m.group(1)


def test_citation_version_matches_pyproject():
    assert _citation_field("version") == _pyproject_version(), (
        "CITATION.cff の version が pyproject.toml とずれている。"
        " release のたびに両方を上げること。"
    )


def test_citation_date_is_a_date():
    """`date-released` を版と一緒に上げ忘れると、日付だけ過去に残る。"""
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", _citation_field("date-released"))
