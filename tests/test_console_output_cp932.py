# -*- coding: utf-8 -*-
"""Everything that can reach a console must encode in cp932.

A Japanese Windows console is cp932 by default, and that is where J-OCTA
runs. A single character outside it — an em dash in a help string is
enough — turns `--help` into

    UnicodeEncodeError: 'cp932' codec can't encode character '\\u2014'

with a traceback and exit code 1. The Japanese text itself is fine in
cp932; it is the punctuation that is not, so the fix is ASCII punctuation
rather than forcing the stream to UTF-8, which would only move the
breakage onto the Japanese (2026-09-15, found on J-OCTA 11.1 and 12.0;
`amorphous --help` and `udf2gro --help` both died, in 2.14.1 as well).

Checked over the source rather than the output: a string only has to be
reachable to hurt, and comments cost nothing to keep clean.
"""
from __future__ import annotations

import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
PKG = ROOT / "abmptools"

#: cp932 にない記号と、置き換え先。増やすときは「画面に出す意味が
#: 変わらない」ことを確かめること。
ASCII_FOR = {"—": "--", "–": "-", "Å": "A", "·": "*", "≈": "~=",
             "≤": "<=", "≥": ">=", "²": "^2", "³": "^3", "µ": "u"}


def test_no_source_file_holds_a_character_cp932_cannot_write():
    bad = []
    for path in sorted(PKG.rglob("*.py")):
        for num, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            for ch in line:
                try:
                    ch.encode("cp932")
                except UnicodeEncodeError:
                    bad.append(f"{path.relative_to(ROOT)}:{num}: {ch!r}"
                               f" -> {ASCII_FOR.get(ch, '?')}")
    assert not bad, (
        "cp932 のコンソールで落ちる文字がある:\n  " + "\n  ".join(bad[:20])
        + (f"\n  ... 他 {len(bad) - 20} 件" if len(bad) > 20 else ""))
