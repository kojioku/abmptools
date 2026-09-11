# -*- coding: utf-8 -*-
"""Amber 二面角には 1-4 のスケーリングを添える。

`gro2udf` は `User_Torsion.Parameters[]` に **何も書いていなかった**。UDF は
形式として正しく、変換も GROMACS 書き出しも通る。ところが **J-OCTA で NPT を
流すと落ち、「力場を取得しなおす」と通る** (2026-09-12 に実機で確認)。
取り直した後のファイルには `SCNB` / `SCEE` が入っていた —— J-OCTA は 1-4 の
扱いをここから読む。

値は GROMACS の `[ defaults ]` の `fudgeLJ` / `fudgeQQ` をそのまま入れる
(AMBER の SCNB は本来「割る数」だが、実測で一致したのは倍率の方だった)。
多重度の 2 つめ以降 (`:1`, `:2`) には付けない。取り直し後のファイルも
先頭の項にだけ持っていた。
"""
from __future__ import annotations

from abmptools.gro2udf.top_exporter import (
    _fudge_str,
    _is_multiplicity_continuation,
)


def test_the_first_multiplicity_term_carries_the_scaling():
    assert _is_multiplicity_continuation("c3-c3-c3-c3:0") is False


def test_later_multiplicity_terms_do_not():
    assert _is_multiplicity_continuation("c3-c3-c3-c3:1") is True
    assert _is_multiplicity_continuation("c3-c3-c3-c3:2") is True


def test_a_plain_torsion_carries_the_scaling():
    assert _is_multiplicity_continuation("hc-c3-c3-hc") is False


def test_a_name_with_a_colon_but_no_number_is_not_a_continuation():
    assert _is_multiplicity_continuation("weird:name") is False


def test_fudge_values_keep_the_shape_j_octa_writes():
    """取り直し後のファイルは '0.5' と '0.8333333'。末尾の 0 は付かない。"""
    assert _fudge_str(0.5) == "0.5"
    assert _fudge_str(0.8333333) == "0.8333333"


def test_an_integral_fudge_does_not_become_empty():
    assert _fudge_str(1.0) == "1"
    assert _fudge_str(0.0) == "0"
