# -*- coding: utf-8 -*-
"""`.top` / `.itp` / `.gro` / `.mdp` は UTF-8 と明示して開くこと。

`open()` の encoding を省くとロケール既定になる。Linux は UTF-8 なので何も
起きないが、**日本語 Windows (cp932) では UTF-8 の `.top` を読んだ瞬間に落ちる**。

    UnicodeDecodeError: 'cp932' codec can't decode byte 0x87 in position 2874

`.top` にコメントで日本語を書くのは普通にあることで、受け取る側の importer 2 版
どちらでも再現した (2026-09-13、実機)。**`moldeck.hbond` の `.top` 経路は
このパーサを通る**ので、機能ごと使えなくなる。

書き側も同じ: cp932 に無い文字を書く瞬間に `UnicodeEncodeError` になる。

Linux でも `LC_ALL=C PYTHONUTF8=0` で既定を ASCII にすれば再現できる
(`PYTHONIOENCODING=utf-8` を併せないと、日本語の表示のほうが先に落ちる)。
"""
from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
#: テキストを読み書きするサブパッケージ。バイナリしか触らないものは対象外。
PKGS = ("gro2udf", "udf2gro")

_OPEN = re.compile(r'(?<![\w.])open\s*\(|\.open\s*\(')
_BINARY = re.compile(r'''["'][rwax]\+?b["']''')


def test_every_text_open_names_utf8():
    bad = []
    for pkg in PKGS:
        for path in sorted((ROOT / "abmptools" / pkg).rglob("*.py")):
            for num, line in enumerate(
                    path.read_text(encoding="utf-8").splitlines(), 1):
                code = line.split("#")[0]
                if not _OPEN.search(code) or "encoding=" in code:
                    continue
                if _BINARY.search(code):
                    continue
                bad.append(f"{path.relative_to(ROOT)}:{num}  {code.strip()}")
    assert not bad, (
        "encoding を明示していない open がある:\n  " + "\n  ".join(bad))


def test_a_top_with_japanese_comments_parses(tmp_path):
    """実際に非 ASCII を含む `.top` を食わせる。

    静的検査だけだと、新しい経路が別の関数経由で開いたときに漏れる。
    """
    from abmptools.gro2udf.top_parser import TopParser
    top = tmp_path / "jp.top"
    top.write_text(
        "; 自動生成: 日本語のコメント\n"
        "[ defaults ]\n"
        "1 2 yes 0.5 0.8333\n"
        "[ atomtypes ]\n"
        "  c3  12.011  0.0  A  0.3399  0.4577\n"
        "[ moleculetype ]\n"
        "MOL 3\n"
        "[ atoms ]\n"
        "  1  c3  1  MOL  C1  1  0.0  12.011\n"
        "[ system ]\n"
        "test\n"
        "[ molecules ]\n"
        "MOL 1\n", encoding="utf-8")
    model = TopParser().parse(str(top))
    assert model is not None
