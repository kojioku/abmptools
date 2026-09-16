# -*- coding: utf-8 -*-
"""`[ atomtypes ]` carries the atomic number, so a reader can tell the element.

GROMACS does not need it — it works out what it needs from the parameters —
but an importer does. The names in the column are force field types ("c3",
"os", "hc"), which mean nothing outside that force field, and the atom names
in `[ atoms ]` are derived from them. Without an atomic number the only clue
left is the mass, and a reader that does not guess from mass has no way to
distinguish an all-atom system from coarse-grained beads: the importer
on the other side read a converted all-atom system as CG (2026-09-15,
11.0).

A coarse-grained bead genuinely has no element, and 0 — GROMACS' own value
for "unknown" — is the right answer there, so this is not a case of always
filling something in.
"""
from __future__ import annotations

from abmptools.udf2gro.gromacs.writers.top_writer import _atomic_number


def test_the_elements_a_force_field_actually_ships():
    """Rounded masses, as force fields write them."""
    assert _atomic_number(1.008) == 1
    assert _atomic_number(12.01) == 6      # GAFF carbon
    assert _atomic_number(12.011) == 6
    assert _atomic_number(14.01) == 7
    assert _atomic_number(16.0) == 8       # GAFF oxygen, rounded
    assert _atomic_number(15.999) == 8
    assert _atomic_number(32.06) == 16
    assert _atomic_number(35.45) == 17


def test_a_bead_gets_zero_rather_than_a_guess():
    """A CG bead has no element; 0 says so.

    Martini beads sit around 72 amu, nowhere near an element, and inventing
    one would turn "unknown" into a wrong answer.
    """
    assert _atomic_number(72.0) == 0
    assert _atomic_number(0.0) == 0
    assert _atomic_number(45.0) == 0


def test_neighbouring_elements_are_not_confused():
    """The tolerance must not let one element answer for the next."""
    assert _atomic_number(14.007) == 7
    assert _atomic_number(15.999) == 8
    assert _atomic_number(18.998) == 9
    assert _atomic_number(19.5) == 0       # 間はどちらでもない


def test_atoms_carry_the_mass_column():
    """`[ atoms ]` の 8 列目に質量を出すこと。

    GROMACS では省略可 —— 無ければ `[ atomtypes ]` から引く、という規約。
    **受け取る側の importer はその規約を実装していない。**
    `convert_gromacs_udf.py:238-246` は

        if len(stmp) > 7:
            mass = float(stmp[7])
            ... mass_map[round(mass)] で元素名を決める

    と**この列だけ**を見ており、無ければ力場の型名 (`hc1`) がそのまま
    原子名になる。結果、**全原子系が粗視化として読まれる** (2026-09-15、
    その importer の 2 つの版で両方向を確認)。

    `[ atomtypes ]` の `at.num` は同じ importer の中で**コメントアウト
    されている**ので、そちらを足しても解決しない —— 最初にそう直して
    外した。
    """
    import pathlib
    ref = (pathlib.Path(__file__).parent / "regression" / "reference"
           / "prerefactor" / "udf2gro" / "test.top")
    lines = ref.read_text(encoding="utf-8").splitlines()
    i = next(n for n, s in enumerate(lines) if s.startswith("[ atoms ]"))
    assert "mass" in lines[i + 1], "[ atoms ] のヘッダに mass が無い"
    row = lines[i + 2].split()
    assert len(row) == 8, f"[ atoms ] が 8 列でない: {row}"
    assert float(row[7]) > 0, f"質量が入っていない: {row}"
