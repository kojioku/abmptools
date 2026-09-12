# -*- coding: utf-8 -*-
"""分子名を GROMACS 側の制約に合わせるヘルパ。

**この 2 つは 2026-09-12 まで 1 度もテストされていなかった。** リファクタで
引数名を変えたとき本体を直し忘れても、 全 1036 件が緑のままだったので
気付いた。 UDF の `Mol_Name` が素直な英数字であるかぎり例外を踏まないだけで、
空白や記号が入った瞬間に効いてくる経路。
"""
from __future__ import annotations

import pytest

from abmptools.udf2gro.udf_adapter import (
    _sanitize_gromacs_molname,
    _shorten_molname,
)


# ---------------------------------------------------------------------------
# [ moleculetype ] Name — 空白とコメント文字が入ると .top が壊れる
# ---------------------------------------------------------------------------

def test_spaces_become_underscores():
    """空白があると GROMACS が 2 つ目のフィールドとして読む。"""
    assert _sanitize_gromacs_molname("poly ethylene") == "poly_ethylene"


def test_comment_characters_are_replaced():
    """`;` 以降が黙って捨てられるので、 記号は残さない。"""
    assert _sanitize_gromacs_molname("mol;name") == "mol_name"
    assert _sanitize_gromacs_molname("a#b*c") == "a_b_c"


def test_a_leading_digit_gets_a_prefix():
    assert _sanitize_gromacs_molname("3mol") == "M_3mol"


def test_alphanumeric_and_underscore_pass_through():
    assert _sanitize_gromacs_molname("nylon6_A") == "nylon6_A"


@pytest.mark.parametrize("empty", ["", "___", "!!!"])
def test_nothing_usable_falls_back_to_mol(empty):
    """空や記号だけだと名前が消える。 落とさず MOL にする。"""
    assert _sanitize_gromacs_molname(empty) == "MOL"


# ---------------------------------------------------------------------------
# .gro の残基名欄は 5 文字
# ---------------------------------------------------------------------------

def test_short_names_are_untouched():
    assert _shorten_molname("ABC") == "ABC"
    assert _shorten_molname("ABCDE") == "ABCDE"


def test_long_names_keep_both_ends():
    """真ん中を落とす。 先頭 2 + 末尾 3 なので、 連番の区別が残る。"""
    assert _shorten_molname("nylon6") == "nyon6"
    assert _shorten_molname("polyethylene") == "poene"


def test_the_result_always_fits_the_gro_column():
    for name in ("a", "abcdef", "x"*40, "molecule_001"):
        assert len(_shorten_molname(name)) <= 5


def test_names_differing_only_at_the_end_stay_distinct():
    """末尾を残すのはこのため。 連番の系で衝突しない。"""
    a = _shorten_molname("component_01")
    b = _shorten_molname("component_02")
    assert a != b
