"""OCTA/COGNAC UDF への per-atom 電荷割り当て (charge transfer)。

単分子 UDF (電荷あり) から per-atom partial charge を抽出し、 バルク系 UDF
(同名分子、 電荷なし) の **同名分子すべて**へ同じ電荷を割り当てる。

COGNAC UDF の電荷規約 (gro2udf / udf2gro と共通):

- 電荷は ``Set_of_Molecules.molecule[i].electrostatic_Site[j]`` に格納:
  - ``.Type_Name``  = ``"POINT_CHARGE"``
  - ``.ES_Element`` = **partial charge [e] × CHARGE_UNIT (18.224159264)**
    (COGNAC 内部の電荷スケール単位)
  - ``.atom[0]``    = この site が属する atom の (分子内 0-origin) index
- 読み戻しは ``charge[e] = ES_Element / CHARGE_UNIT``。

注意点:
- ``UDFManager.put`` は numpy float32 をサイレントに 0 化するため、 ES_Element は
  必ず Python ``float()`` で渡す ([[reference_udfmanager_put_numpy_silent_zero]])。
- Set_of_Molecules は static record (``jump(-1)``) に格納される。
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

__all__ = [
    "CHARGE_UNIT",
    "POINT_CHARGE",
    "MoleculeChargeTemplate",
    "AssignResult",
    "read_molecule_charges",
    "assign_charges_to_bulk",
]

# COGNAC: ES_Element = charge[e] * CHARGE_UNIT  (gro2udf/udf2gro と同一定数)
CHARGE_UNIT = 18.224159264
POINT_CHARGE = "POINT_CHARGE"

_MOL = "Set_of_Molecules.molecule[]"
_ATOM = "Set_of_Molecules.molecule[].atom[]"
_ES = "Set_of_Molecules.molecule[].electrostatic_Site[]"


def _open(udf_path) -> "object":
    """UDFManager を開いて static record (jump(-1)) に移動して返す。"""
    from UDFManager import UDFManager  # OCTA 同梱 (PyPI 非配布)

    u = UDFManager(str(udf_path))
    u.jump(-1)
    return u


@dataclass
class MoleculeChargeTemplate:
    """単分子 UDF から抽出した per-atom 電荷テンプレート。

    Attributes
    ----------
    mol_name
        抽出元分子の ``Mol_Name``。 割り当て対象の選別キー。
    n_atoms
        分子内 atom 数。
    charges
        per-atom partial charge [e] (atom 0..n-1 順)。
    atom_type_names
        per-atom ``Atom_Type_Name`` (割り当て時の整合検証に使う)。
    atom_names
        per-atom ``Atom_Name`` (element symbol、 参考)。
    net_charge
        ``sum(charges)`` (中性確認の目安)。
    source_index
        抽出元 UDF 内での molecule index。
    """

    mol_name: str
    n_atoms: int
    charges: List[float]
    atom_type_names: List[str] = field(default_factory=list)
    atom_names: List[str] = field(default_factory=list)
    source_index: int = 0

    @property
    def net_charge(self) -> float:
        return float(sum(self.charges))


@dataclass
class AssignResult:
    """:func:`assign_charges_to_bulk` の結果サマリ。"""

    out_path: str
    mol_name: str
    n_atoms_per_mol: int
    n_molecules_assigned: int
    n_molecules_total: int
    skipped_indices: List[int] = field(default_factory=list)


def _find_mol_index(u, *, mol_name: Optional[str], mol_index: Optional[int]) -> int:
    """読み出し対象の molecule index を決める。"""
    nmol = u.size(_MOL)
    if nmol is None or nmol == 0:
        raise ValueError("UDF に Set_of_Molecules.molecule[] が存在しません")
    if mol_index is not None:
        if not (0 <= mol_index < nmol):
            raise IndexError(f"mol_index={mol_index} は範囲外 (nmol={nmol})")
        return mol_index
    if mol_name is not None:
        for i in range(nmol):
            if u.get(f"{_MOL}.Mol_Name", [i]) == mol_name:
                return i
        raise ValueError(f"Mol_Name={mol_name!r} の分子が見つかりません")
    return 0  # default: 先頭分子 (単分子 UDF を想定)


def read_molecule_charges(
    udf_path,
    *,
    mol_name: Optional[str] = None,
    mol_index: Optional[int] = None,
) -> MoleculeChargeTemplate:
    """単分子 (template) UDF から per-atom 電荷を抽出する。

    Parameters
    ----------
    udf_path
        電荷を持つ UDF へのパス (通常は単分子)。
    mol_name
        抽出する分子の ``Mol_Name``。 指定時は最初に一致した分子。
    mol_index
        抽出する分子の index。 ``mol_name`` より優先。 どちらも未指定なら 0。

    Returns
    -------
    MoleculeChargeTemplate
        ``electrostatic_Site[].atom[0]`` で原子対応を取り、 ES site の無い原子は
        電荷 0.0 とする (atom 数は ``atom[]`` のサイズが基準)。
    """
    u = _open(udf_path)
    idx = _find_mol_index(u, mol_name=mol_name, mol_index=mol_index)

    name = u.get(f"{_MOL}.Mol_Name", [idx])
    n_atoms = int(u.size(_ATOM, [idx]) or 0)
    if n_atoms == 0:
        raise ValueError(f"molecule[{idx}] ({name!r}) に atom がありません")

    atom_types = [u.get(f"{_ATOM}.Atom_Type_Name", [idx, i]) for i in range(n_atoms)]
    atom_names = [u.get(f"{_ATOM}.Atom_Name", [idx, i]) for i in range(n_atoms)]

    charges = [0.0] * n_atoms
    n_es = int(u.size(_ES, [idx]) or 0)
    if n_es == 0:
        logger.warning(
            "molecule[%d] (%r) に electrostatic_Site が無く、 全電荷 0 として読み出します",
            idx, name,
        )
    for j in range(n_es):
        at = u.get(f"{_ES}.atom[]", [idx, j, 0])
        es = u.get(f"{_ES}.ES_Element", [idx, j])
        ai = int(at) if at is not None else j
        if not (0 <= ai < n_atoms):
            logger.warning("ES site %d -> atom %d が範囲外、 skip", j, ai)
            continue
        charges[ai] = float(es) / CHARGE_UNIT

    tmpl = MoleculeChargeTemplate(
        mol_name=str(name),
        n_atoms=n_atoms,
        charges=charges,
        atom_type_names=atom_types,
        atom_names=atom_names,
        source_index=idx,
    )
    logger.info(
        "read template: mol=%r n_atoms=%d net_charge=%.4f (from %s)",
        tmpl.mol_name, tmpl.n_atoms, tmpl.net_charge, Path(udf_path).name,
    )
    return tmpl


def assign_charges_to_bulk(
    bulk_udf_path,
    template: MoleculeChargeTemplate,
    out_path=None,
    *,
    verify_atom_types: bool = True,
    strict: bool = True,
) -> AssignResult:
    """バルク UDF の ``template.mol_name`` 分子すべてに per-atom 電荷を割り当てる。

    電荷は **atom index 対応** (template の atom i → 対象分子の atom i) で書き込む。
    同一分子であることの担保として atom 数 (必須) と ``Atom_Type_Name`` 列
    (``verify_atom_types``) を検証する。

    Parameters
    ----------
    bulk_udf_path
        割り当て先のバルク UDF (同名分子が複数、 通常 electrostatic_Site なし)。
    template
        :func:`read_molecule_charges` で得たテンプレート。
    out_path
        出力 UDF パス。 ``None`` の場合 ``<bulk>_charged.udf`` を同ディレクトリに作る
        (入力を上書きしない)。
    verify_atom_types
        対象分子の ``Atom_Type_Name`` 列が template と一致するか検証する。
    strict
        検証不一致のとき ``True`` なら例外、 ``False`` なら warning で skip。

    Returns
    -------
    AssignResult
    """
    bulk_udf_path = Path(bulk_udf_path)
    if out_path is None:
        out_path = bulk_udf_path.with_name(bulk_udf_path.stem + "_charged" + bulk_udf_path.suffix)
    out_path = Path(out_path)

    u = _open(bulk_udf_path)
    nmol = int(u.size(_MOL) or 0)
    if nmol == 0:
        raise ValueError("bulk UDF に分子がありません")

    assigned: List[int] = []
    skipped: List[int] = []
    for i in range(nmol):
        if u.get(f"{_MOL}.Mol_Name", [i]) != template.mol_name:
            continue

        n_atoms = int(u.size(_ATOM, [i]) or 0)
        if n_atoms != template.n_atoms:
            msg = (f"molecule[{i}] ({template.mol_name!r}) の atom 数 {n_atoms} が "
                   f"template の {template.n_atoms} と不一致")
            if strict:
                raise ValueError(msg)
            logger.warning("%s → skip", msg)
            skipped.append(i)
            continue

        if verify_atom_types:
            tgt_types = [u.get(f"{_ATOM}.Atom_Type_Name", [i, k]) for k in range(n_atoms)]
            if tgt_types != template.atom_type_names:
                msg = (f"molecule[{i}] の Atom_Type_Name 列が template と不一致 "
                       f"(例: {tgt_types[:3]} vs {template.atom_type_names[:3]})")
                if strict:
                    raise ValueError(msg)
                logger.warning("%s → skip", msg)
                skipped.append(i)
                continue

        for k in range(n_atoms):
            u.put(POINT_CHARGE, f"{_ES}.Type_Name", [i, k])
            # float() cast 必須: numpy 値だと UDFManager が silent 0 化する
            u.put(float(template.charges[k] * CHARGE_UNIT), f"{_ES}.ES_Element", [i, k])
            u.put(int(k), f"{_ES}.atom[]", [i, k, 0])
        assigned.append(i)

    if not assigned and strict:
        raise ValueError(
            f"Mol_Name={template.mol_name!r} に一致する分子が bulk UDF にありません"
        )

    u.write(str(out_path))
    res = AssignResult(
        out_path=str(out_path),
        mol_name=template.mol_name,
        n_atoms_per_mol=template.n_atoms,
        n_molecules_assigned=len(assigned),
        n_molecules_total=nmol,
        skipped_indices=skipped,
    )
    logger.info(
        "assigned charges to %d/%d molecules (%r) → %s",
        res.n_molecules_assigned, res.n_molecules_total, res.mol_name, out_path.name,
    )
    return res
