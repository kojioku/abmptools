#!/usr/bin/env python
"""abmptools.udfcharge のデモ用 UDF を生成する。

methanol (CH3OH) を題材に:
- ``methanol_single_charged.udf`` — 1 分子、 per-atom 電荷あり (template)
- ``methanol_bulk_uncharged.udf`` — 同じ methanol × 8、 電荷なし (割り当て先)

その後:
    python -m abmptools.udfcharge \
        --template methanol_single_charged.udf \
        --bulk methanol_bulk_uncharged.udf \
        --out methanol_bulk_charged.udf -v
で全 8 分子へ電荷を転写できる。
"""
from __future__ import annotations

import shutil
from pathlib import Path

from UDFManager import UDFManager

import abmptools.gro2udf as _g
from abmptools.udfcharge import CHARGE_UNIT

TEMPLATE = Path(_g.__file__).parent / "default_template.udf"

MOL_NAME = "MeOH"
#            C        O        Ho       H1       H2       H3
ATOM_TYPES = ["c3", "oh", "ho", "h1", "h1", "h1"]
ATOM_ELEMS = ["C", "O", "H", "H", "H", "H"]
CHARGES = [0.145, -0.683, 0.418, 0.040, 0.040, 0.040]  # 和 = 0


def build(out_path: Path, n_copies: int, with_charges: bool) -> None:
    shutil.copy(TEMPLATE, out_path)
    u = UDFManager(str(out_path))
    u.jump(-1)
    for imol in range(n_copies):
        u.put(MOL_NAME, "Set_of_Molecules.molecule[].Mol_Name", [imol])
        for i, (el, ty) in enumerate(zip(ATOM_ELEMS, ATOM_TYPES)):
            u.put(i, "Set_of_Molecules.molecule[].atom[].Atom_ID", [imol, i])
            u.put(el, "Set_of_Molecules.molecule[].atom[].Atom_Name", [imol, i])
            u.put(ty, "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [imol, i])
            if with_charges:
                u.put("POINT_CHARGE",
                      "Set_of_Molecules.molecule[].electrostatic_Site[].Type_Name", [imol, i])
                u.put(float(CHARGES[i] * CHARGE_UNIT),
                      "Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element", [imol, i])
                u.put(i,
                      "Set_of_Molecules.molecule[].electrostatic_Site[].atom[]", [imol, i, 0])
    u.write(str(out_path))
    print(f"wrote {out_path.name}: {n_copies} x {MOL_NAME}, charges={'yes' if with_charges else 'no'}")


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    build(here / "methanol_single_charged.udf", n_copies=1, with_charges=True)
    build(here / "methanol_bulk_uncharged.udf", n_copies=8, with_charges=False)
