#!/usr/bin/env python
"""abmptools.udfcharge の「形式電荷の復元 (restore)」デモ。

MD では系を中性 (Σq=0) にして走らせるため、 元の形式電荷 (整数) を持つ分子の
電荷が中和されて UDF に入っていることがある。 ``restore`` はその中和電荷と
目標形式電荷から、 中和前の per-atom 電荷 (Σ=形式電荷) を復元する。

ここでは methylammonium 様の +1 カチオン (8 atom) を題材に:
1. 元電荷 A (Σ=+1) を |A| 比例で中和して B (Σ=0) にした UDF を作る
2. restore で B + 形式電荷 +1 → A を復元する

    python restore_example.py
    python -m abmptools.udfcharge restore \
        --udf methylammonium_neutral.udf --formal-charge 1 \
        --out methylammonium_restored.udf -v
"""
from __future__ import annotations

import shutil
from pathlib import Path

from UDFManager import UDFManager

import abmptools.gro2udf as _g
from abmptools.udfcharge import CHARGE_UNIT, read_molecule_charges, restore_formal_charge

TEMPLATE = Path(_g.__file__).parent / "default_template.udf"

MOL_NAME = "MAM"  # methylammonium CH3-NH3+
#           N      HN     HN     HN     C      HC     HC     HC
ELEMS = ["N", "H", "H", "H", "C", "H", "H", "H"]
TYPES = ["n4", "hn", "hn", "hn", "c3", "hc", "hc", "hc"]
# 元の per-atom 電荷 A (形式電荷 +1 に和が一致)
A = [-0.80, 0.45, 0.45, 0.45, -0.20, 0.20, 0.20, 0.25]
FORMAL_CHARGE = 1


def _neutralize(charges, S):
    """A (Σ=S) を |A| 比例で中和して B (Σ≈0) を返す (forward)。"""
    sa = sum(abs(a) for a in charges)
    return [a - S * abs(a) / sa for a in charges]


def build_neutral(out_path: Path) -> None:
    B = _neutralize(A, FORMAL_CHARGE)
    shutil.copy(TEMPLATE, out_path)
    u = UDFManager(str(out_path))
    u.jump(-1)
    u.put(MOL_NAME, "Set_of_Molecules.molecule[].Mol_Name", [0])
    for i, (el, ty, b) in enumerate(zip(ELEMS, TYPES, B)):
        u.put(i, "Set_of_Molecules.molecule[].atom[].Atom_ID", [0, i])
        u.put(el, "Set_of_Molecules.molecule[].atom[].Atom_Name", [0, i])
        u.put(ty, "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [0, i])
        u.put("POINT_CHARGE",
              "Set_of_Molecules.molecule[].electrostatic_Site[].Type_Name", [0, i])
        u.put(float(b) * CHARGE_UNIT,
              "Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element", [0, i])
        u.put(i, "Set_of_Molecules.molecule[].electrostatic_Site[].atom[]", [0, i, 0])
    u.write(str(out_path))
    print(f"wrote {out_path.name}: {MOL_NAME} 8 atoms, Σq={sum(B):+.6f} (中和済み)")


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    neutral = here / "methylammonium_neutral.udf"
    build_neutral(neutral)

    res = restore_formal_charge(neutral, FORMAL_CHARGE,
                                here / "methylammonium_restored.udf")
    print(f"restore: Σq {res.input_total:+.6f} → {res.output_total:+.6f} "
          f"(formal {res.formal_charge:+d}, λ={res.lam:.8f})")
    restored = read_molecule_charges(res.out_path).charges
    print("  atom |   A(元)   | restored")
    for i, (a, r) in enumerate(zip(A, restored)):
        print(f"  {ELEMS[i]:>4} | {a:+8.4f} | {r:+8.4f}")
    print(f"  max|restored - A| = {max(abs(a-r) for a,r in zip(A,restored)):.2e}")
