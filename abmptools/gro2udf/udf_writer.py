# -*- coding: utf-8 -*-
"""
udf_writer.py
-------------
Writes one GRO frame's worth of structural data into the current UDF record.

Contract (gro2udf context):
  AtomPosition.mol_id   = 0-based UDF mol index  (index into mol[])
  AtomPosition.atom_id  = 0-based atom index within mol (index into atom[])
  AtomPosition.x/y/z    = position  [nm]
  AtomPosition.vx/vy/vz = velocity  [nm/ps] -> converted to m/s here

UDF paths written per frame:
  Steps
  Time                                          (value in [ps])
  Structure.Position.mol[m].atom[a].{x,y,z}    [nm]
  Structure.Velocity.mol[m].atom[a].{x,y,z}    [m/s]
  Structure.Unit_Cell.Cell_Size.{a,b,c}         [nm]
  Structure.Unit_Cell.Cell_Size.{alpha,beta,gamma}
  Structure.Unit_Cell.Shear_Strain

Each UDF put first tries with an explicit unit string; on failure it falls
back to the no-unit form.  This replicates the try/except pattern from the
original gro2udf.py.
"""
from __future__ import annotations
import logging

logger = logging.getLogger(__name__)


from typing import List

from ..core.system_model import AtomPosition, CellGeometry

_VELOCITY_UNIT = 1000   # nm/ps → m/s  (1 nm/ps = 1000 m/s)
_SHEAR_STRAIN  = 0.0    # constant, same as original


#: `Unit_Parameter.Comment` に書く力場 ID。 番号は
#: 力場定義ファイルの並び (0 始まり) に対応する。 実データ 3 件で
#: 検算済み: GAFF 型の UDF が FF=2、 DREIDING 型 (C_3 / C_33) の参照
#: サンプルが FF=4。
FF_IDS = {
    "amber": 0, "amber20": 1, "gaff": 2, "gaff2": 3, "dreiding": 4,
    "uff": 5, "oplsaa": 6, "loplsaa": 7, "loplsaa2023": 8, "pcff": 9,
}


def ff_comment(force_field) -> str:
    """力場の名前か番号を ``FF=n`` にする。"""
    if force_field is None:
        return ""
    s = str(force_field).strip()
    if not s:                       # "" = 名乗らない (既存の値をそのまま使う)
        return ""
    if s.upper().startswith("FF="):
        return s
    if s.isdigit():
        return "FF=%s" % s
    key = s.lower().replace("-", "").replace("_", "")
    if key not in FF_IDS:
        raise ValueError(
            "unknown force field %r; use one of %s, or a number"
            % (force_field, ", ".join(sorted(FF_IDS))))
    return "FF=%d" % FF_IDS[key]


def set_force_field_comment(udf, force_field, overwrite: bool = False) -> None:
    """``Unit_Parameter.Comment`` に力場 ID を書く。

    下流はここを見て力場を決める。 **空だと「力場が分からない」扱いになり、
    GROMACS への書き出しが通らない** (2026-09-12 に実機で確認。 下流で
    力場を取得しなおすと ``FF=2`` が入り、 それで通るようになった)。
    ``.top`` には力場の種類が書かれていないので、 こちらで名乗るしかない。

    既存の値は既定では残す。 テンプレート由来の ``FF=n`` を上書きしないため。
    """
    comment = ff_comment(force_field)
    if not comment:
        return
    try:
        udf.jump(-1)
        current = udf.get("Unit_Parameter.Comment")
    except Exception:                                    # noqa: BLE001
        current = None
    if current and not overwrite:
        if str(current).strip() != comment:
            logger.info("Unit_Parameter.Comment: keeping %r (not %r)",
                        current, comment)
        return
    udf.put(comment, "Unit_Parameter.Comment")
    logger.info("Unit_Parameter.Comment = %r (force field id)", comment)


def write_static_cell_abc(udf, a_nm: float, b_nm: float, c_nm: float) -> None:
    """Put a box into the static cell and Initial_Unit_Cell (both in nm).

    ``write_frame`` only writes the cell of the record it is filling, so
    without this the two static fields keep whatever the template had. For
    the bundled default that is a 20 A cube and a 100 A cube, neither of
    which has anything to do with the system. Everything downstream that
    reads a record is fine, which is why it went unnoticed: a reader that
    takes the static section instead gets coordinates from one box and a
    cell from another. A downstream GROMACS converter accepted such a file and
    produced an infinite cell at run time.

    NVT hides it, because then every box is the same one.
    """
    a, b, c = float(a_nm), float(b_nm), float(c_nm)
    udf.jump(-1)
    for field, value in zip("abc", (a, b, c)):
        for root in ("Structure.Unit_Cell",
                     "Initial_Structure.Initial_Unit_Cell"):
            try:
                udf.put(value, "%s.Cell_Size.%s" % (root, field), "[nm]")
            except TypeError:
                udf.put(value, "%s.Cell_Size.%s" % (root, field))
    for field in ("alpha", "beta", "gamma"):
        for root in ("Structure.Unit_Cell",
                     "Initial_Structure.Initial_Unit_Cell"):
            udf.put(90.0, "%s.Cell_Size.%s" % (root, field))


def warn_if_template_box_differs(udf, template_path,
                                 a_nm: float, b_nm: float, c_nm: float) -> None:
    """Say so when the template's static box is not the box of the data.

    The per-record cells are always right, so a mismatch is invisible until
    a reader takes the static section. Worth a line either way, because it
    also catches a template picked up by accident from the working
    directory.
    """
    try:
        udf.jump(-1)
        old = udf.get("Structure.Unit_Cell.Cell_Size")
    except Exception:                                    # noqa: BLE001
        return
    if not old:
        return
    new = [float(a_nm) * 10.0, float(b_nm) * 10.0, float(c_nm) * 10.0]
    if all(abs(float(o) - n) <= 1e-6 * max(1.0, n)
           for o, n in zip(old[:3], new)):
        return
    logger.warning(
        "template %s declares a static cell of %.4f x %.4f x %.4f A, but the "
        "data is %.4f x %.4f x %.4f A. The static cell and Initial_Unit_Cell "
        "are being set from the first frame so the two agree. Check that this "
        "template is the one you meant: a template supplies the static "
        "structure, so a pre-MD UDF brings the box the system had before it "
        "ran.",
        template_path, float(old[0]), float(old[1]), float(old[2]),
        new[0], new[1], new[2],
    )


class UDFWriter:
    """Writes positional / cell data for one GRO frame into the UDF record.

    Caller is responsible for calling ``udf.newRecord()`` before invoking
    ``write_frame()``.
    """

    def write_frame(
        self,
        udf,
        positions: List[AtomPosition],
        cell: CellGeometry,
        step: int,
        time: float,    # [ps]
    ) -> None:
        """Write Steps, Time, Position, Velocity and Cell to the current record."""
        # --- Steps and Time ---
        udf.put(step, "Steps")
        try:
            udf.put(time, "Time", "[ps]")
        except TypeError:
            udf.put(time, "Time")

        # --- Atom positions and velocities ---
        for ap in positions:
            m = ap.mol_id       # 0-based mol index
            a = ap.atom_id      # 0-based atom index within mol
            vx = ap.vx * _VELOCITY_UNIT
            vy = ap.vy * _VELOCITY_UNIT
            vz = ap.vz * _VELOCITY_UNIT
            try:
                udf.put(ap.x, "Structure.Position.mol[].atom[].x", [m, a], "[nm]")
                udf.put(ap.y, "Structure.Position.mol[].atom[].y", [m, a], "[nm]")
                udf.put(ap.z, "Structure.Position.mol[].atom[].z", [m, a], "[nm]")
                udf.put(vx,   "Structure.Velocity.mol[].atom[].x", [m, a], "[m/s]")
                udf.put(vy,   "Structure.Velocity.mol[].atom[].y", [m, a], "[m/s]")
                udf.put(vz,   "Structure.Velocity.mol[].atom[].z", [m, a], "[m/s]")
            except TypeError:
                udf.put(ap.x, "Structure.Position.mol[].atom[].x", [m, a])
                udf.put(ap.y, "Structure.Position.mol[].atom[].y", [m, a])
                udf.put(ap.z, "Structure.Position.mol[].atom[].z", [m, a])
                udf.put(vx,   "Structure.Velocity.mol[].atom[].x", [m, a])
                udf.put(vy,   "Structure.Velocity.mol[].atom[].y", [m, a])
                udf.put(vz,   "Structure.Velocity.mol[].atom[].z", [m, a])

        # --- Unit cell ---
        # Cast to Python float — UDFManager.put silently writes 0 for
        # numpy float32 / float64 values (see trajectory_ingest fix).
        ca, cb, cc = float(cell.a), float(cell.b), float(cell.c)
        try:
            udf.put(ca, "Structure.Unit_Cell.Cell_Size.a", "[nm]")
            udf.put(cb, "Structure.Unit_Cell.Cell_Size.b", "[nm]")
            udf.put(cc, "Structure.Unit_Cell.Cell_Size.c", "[nm]")
        except TypeError:
            udf.put(ca, "Structure.Unit_Cell.Cell_Size.a")
            udf.put(cb, "Structure.Unit_Cell.Cell_Size.b")
            udf.put(cc, "Structure.Unit_Cell.Cell_Size.c")

        udf.put(float(cell.alpha), "Structure.Unit_Cell.Cell_Size.alpha")
        udf.put(float(cell.beta),  "Structure.Unit_Cell.Cell_Size.beta")
        udf.put(float(cell.gamma), "Structure.Unit_Cell.Cell_Size.gamma")
        udf.put(_SHEAR_STRAIN, "Structure.Unit_Cell.Shear_Strain")
