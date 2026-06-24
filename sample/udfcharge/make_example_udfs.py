#!/usr/bin/env python
"""abmptools.udfcharge のデモ用 UDF を生成する (topology + 電荷 + **座標**)。

methanol (CH3OH) を題材に:
- ``methanol_single_charged.udf`` — 1 分子、 per-atom 電荷あり (template)
- ``methanol_bulk_uncharged.udf`` — 同じ methanol × 8 (2×2×2 格子)、 電荷なし

両ファイルとも **Structure.Position に実座標 + Unit_Cell** を持つので OCTA
viewer 等でそのまま開ける。 電荷は ``Set_of_Molecules`` (static record)、
座標は dynamic record に入る (gro2udf と同じ構成)。

その後:
    python -m abmptools.udfcharge \
        --template methanol_single_charged.udf \
        --bulk methanol_bulk_uncharged.udf \
        --out methanol_bulk_charged.udf -v
で全 8 分子へ電荷を転写できる (座標は保持される)。
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

# methanol の分子内座標 [nm] (C を原点付近に置いた近似配座)
#   atom index: 0=C 1=O 2=Ho 3=H1 4=H2 5=H3
MOL_XYZ = [
    (0.000, 0.000, 0.000),   # 0 C
    (0.140, 0.000, 0.000),   # 1 O  (C-O ~1.4 Å)
    (0.176, 0.091, 0.000),   # 2 Ho (O-H)
    (-0.036, 0.103, 0.000),  # 3 H1 (methyl)
    (-0.036, -0.051, 0.089),  # 4 H2
    (-0.036, -0.051, -0.089),  # 5 H3
]

# 結合トポロジー (atom index は 0-based、 分子の「形」= viewer の骨格描画に必須)
BONDS = [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5)]                 # C-O, O-Ho, C-H×3
ANGLES = [(1, 0, 3), (1, 0, 4), (1, 0, 5),                       # O-C-H
          (3, 0, 4), (3, 0, 5), (4, 0, 5),                       # H-C-H
          (0, 1, 2)]                                             # C-O-Ho
DIHEDRALS = [(3, 0, 1, 2), (4, 0, 1, 2), (5, 0, 1, 2)]           # H-C-O-Ho (methyl 回転)


def _placements(n_copies: int):
    """分子配置のオフセット [nm] を返す。 1 個なら箱中心、 8 個なら 2×2×2 格子。"""
    if n_copies == 1:
        return [(1.0, 1.0, 1.0)], 2.0          # box 2.0 nm
    # 2×2×2 = 8 個、 格子点 {0.6, 1.6}、 box 2.4 nm (methanol ~0.3 nm で非重複)
    grid = [0.6, 1.6]
    offs = [(x, y, z) for x in grid for y in grid for z in grid]
    return offs[:n_copies], 2.4


def build(out_path: Path, n_copies: int, with_charges: bool) -> None:
    offs, box = _placements(n_copies)
    shutil.copy(TEMPLATE, out_path)
    u = UDFManager(str(out_path))

    # --- topology + charge は static record (-1) ---
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

        # bond / angle / dihedral (分子の形 = viewer 骨格)。 atom は 0-based local index
        for k, (a1, a2) in enumerate(BONDS):
            u.put("", "Set_of_Molecules.molecule[].bond[].Potential_Name", [imol, k])
            u.put(a1, "Set_of_Molecules.molecule[].bond[].atom1", [imol, k])
            u.put(a2, "Set_of_Molecules.molecule[].bond[].atom2", [imol, k])
            u.put(1.0, "Set_of_Molecules.molecule[].bond[].Order", [imol, k])
        for k, (a1, a2, a3) in enumerate(ANGLES):
            u.put("", "Set_of_Molecules.molecule[].angle[].Potential_Name", [imol, k])
            u.put(a1, "Set_of_Molecules.molecule[].angle[].atom1", [imol, k])
            u.put(a2, "Set_of_Molecules.molecule[].angle[].atom2", [imol, k])
            u.put(a3, "Set_of_Molecules.molecule[].angle[].atom3", [imol, k])
        for k, (a1, a2, a3, a4) in enumerate(DIHEDRALS):
            u.put("", "Set_of_Molecules.molecule[].torsion[].Potential_Name", [imol, k])
            u.put(a1, "Set_of_Molecules.molecule[].torsion[].atom1", [imol, k])
            u.put(a2, "Set_of_Molecules.molecule[].torsion[].atom2", [imol, k])
            u.put(a3, "Set_of_Molecules.molecule[].torsion[].atom3", [imol, k])
            u.put(a4, "Set_of_Molecules.molecule[].torsion[].atom4", [imol, k])

    # --- 座標 + cell は dynamic record (template の placeholder を erase して 1 frame 追加) ---
    u.eraseRecord(0, u.totalRecord())
    u.newRecord()
    for imol in range(n_copies):
        ox, oy, oz = offs[imol]
        for i, (x, y, z) in enumerate(MOL_XYZ):
            u.put(ox + x, "Structure.Position.mol[].atom[].x", [imol, i], "[nm]")
            u.put(oy + y, "Structure.Position.mol[].atom[].y", [imol, i], "[nm]")
            u.put(oz + z, "Structure.Position.mol[].atom[].z", [imol, i], "[nm]")
    for ax in ("a", "b", "c"):
        u.put(float(box), f"Structure.Unit_Cell.Cell_Size.{ax}", "[nm]")
    for ax in ("alpha", "beta", "gamma"):
        u.put(90.0, f"Structure.Unit_Cell.Cell_Size.{ax}")
    u.put(0, "Steps")
    u.put(0.0, "Time")

    u.write(str(out_path))
    print(f"wrote {out_path.name}: {n_copies} x {MOL_NAME}, "
          f"charges={'yes' if with_charges else 'no'}, coords=yes, "
          f"bonds={len(BONDS)}/angles={len(ANGLES)}/dih={len(DIHEDRALS)}, box={box} nm")


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    build(here / "methanol_single_charged.udf", n_copies=1, with_charges=True)
    build(here / "methanol_bulk_uncharged.udf", n_copies=8, with_charges=False)
