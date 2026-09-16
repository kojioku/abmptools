# -*- coding: utf-8 -*-
"""Everything that can reach a console must encode in cp932.

A Japanese Windows console is cp932 by default, and that is where the downstream tooling
runs. A single character outside it — an em dash in a help string is
enough — turns `--help` into

    UnicodeEncodeError: 'cp932' codec can't encode character '\\u2014'

with a traceback and exit code 1. The Japanese text itself is fine in
cp932; it is the punctuation that is not, so the fix is ASCII punctuation
rather than forcing the stream to UTF-8, which would only move the
breakage onto the Japanese (2026-09-15, found on two releases of it;
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


def test_gro2udf_help_lists_the_options_people_come_for():
    """`--help` に `--trajectory` が出ること。

    以前は 2 行の usage を出して「全部見るなら `--from-top --help`」と
    案内していた。**ほとんどの人がここに来る理由である `--trajectory` が、
    実際に読まれる help に載っていなかった** (2026-09-16、利用者の指摘)。

    `--help` は**実際のオプション一覧**を出す。案内の中に案内を置かない。
    """
    import io
    import contextlib

    from abmptools.gro2udf.cli import _usage

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        _usage(["gro2udf"])
    out = buf.getvalue()
    for flag in ("--trajectory", "--energy", "--from-top", "--topology-only"):
        assert flag in out, f"{flag} が --help に出ていない"
    assert ".xtc" in out, "xtc を受けることが help から読み取れない"
