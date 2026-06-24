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
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

__all__ = [
    "CHARGE_UNIT",
    "POINT_CHARGE",
    "MoleculeChargeTemplate",
    "AssignResult",
    "RestoreResult",
    "read_molecule_charges",
    "assign_charges_to_bulk",
    "restore_formal_charge",
]

# COGNAC: ES_Element = charge[e] * CHARGE_UNIT  (gro2udf/udf2gro と同一定数)
CHARGE_UNIT = 18.224159264
POINT_CHARGE = "POINT_CHARGE"

_MOL = "Set_of_Molecules.molecule[]"
_ATOM = "Set_of_Molecules.molecule[].atom[]"
_ES = "Set_of_Molecules.molecule[].electrostatic_Site[]"


def _put_molecule_charges(u, imol, charges) -> None:
    """molecule[imol] の各 atom に electrostatic_Site (POINT_CHARGE) を書く。

    float() cast 必須: numpy 値だと UDFManager が silent 0 化する。
    """
    for k, q in enumerate(charges):
        u.put(POINT_CHARGE, f"{_ES}.Type_Name", [imol, k])
        u.put(float(q) * CHARGE_UNIT, f"{_ES}.ES_Element", [imol, k])
        u.put(int(k), f"{_ES}.atom[]", [imol, k, 0])


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

        _put_molecule_charges(u, i, template.charges)
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


# ---------------------------------------------------------------------------
# 形式電荷の復元 (neutralize の逆変換)
# ---------------------------------------------------------------------------
#
# 前提となる中和ルール (forward): 元電荷 A_i (Σ A = S = 形式電荷) を、 過剰分 S を
# ``|A_i|`` 比例で各原子に分散して中和し、 B_i (Σ B ≈ 0) を得る:
#
#     B_i = A_i - S·|A_i| / Σ|A|          (λ = S/Σ|A| とおくと)
#     B_i = (1-λ)A_i   (A_i>0)
#     B_i = (1+λ)A_i   (A_i<0)
#
# reverse (本機能): B_i と目標形式電荷 S だけから A_i を復元する。 符号は保存される
# (λ<1 のとき) ので:
#
#     A_i = B_i/(1-λ)   (B_i>0)
#     A_i = B_i/(1+λ)   (B_i<0)
#
# λ は Σ A = S の制約から、 P=Σ_{B>0}B / N=Σ_{B<0}B を用いて
#
#     S·λ² + (P-N)·λ + (P+N-S) = 0
#
# の |λ|<1 の根として求まる (詳細は docs/udfcharge.md / SI/A列再現方法.md)。
#
# もう一つの中和ルール ``mode="uniform"`` (均等分配): 過剰分 S を全原子に均等に
# 分散して中和する。
#
#     forward:  B_i = A_i - S/N        (N = 原子数、 Σ B = 0)
#     reverse:  A_i = B_i + S/N
#
# こちらは二次方程式も符号反転問題も無く、 常に厳密・一意に復元できる。


@dataclass
class RestoreResult:
    """:func:`restore_formal_charge` の結果サマリ。"""

    out_path: str
    mol_name: str
    n_atoms: int
    formal_charge: int
    input_total: float       # 入力 UDF の総電荷 (≈0 を想定)
    output_total: float      # 出力 UDF の総電荷 (= formal_charge)
    mode: str = "proportional"          # "proportional" | "uniform"
    lam: Optional[float] = None         # proportional 時の λ (uniform では None)
    shift: Optional[float] = None       # uniform 時の S/N (proportional では None)


def _solve_neutralization_lambda(charges: List[float], formal_charge: float) -> float:
    """中和電荷 (B) と目標形式電荷 S から逆変換用 λ を解く。

    ``S·λ² + (P-N)·λ + (P+N-S) = 0`` (P=Σ_{B>0}B, N=Σ_{B<0}B) の |λ|<1 の根。
    S≈0 のときは線形に縮退する。
    """
    P = sum(b for b in charges if b > 0.0)
    N = sum(b for b in charges if b < 0.0)
    S = float(formal_charge)
    a, b, c = S, (P - N), (P + N - S)

    if abs(a) < 1e-12:                       # S ≈ 0 → 線形 (P-N)λ + (P+N) = 0
        if abs(b) < 1e-12:
            return 0.0
        lam = -(P + N) / b
    else:
        disc = b * b - 4.0 * a * c
        if disc < 0:
            raise ValueError(
                f"形式電荷 {S:g} は |q| 比例分配ルールで到達不能です (判別式 < 0)"
            )
        sq = math.sqrt(disc)
        roots = [(-b + sq) / (2.0 * a), (-b - sq) / (2.0 * a)]
        valid = [r for r in roots if abs(r) < 1.0]
        if not valid:
            raise ValueError(
                f"|λ|<1 の根が無く電荷符号が反転します (formal_charge={S:g} が大きすぎ、 "
                "B 列と総和だけからは一意に復元できないケース)"
            )
        lam = min(valid, key=abs)
    return lam


def restore_formal_charge(
    udf_path,
    formal_charge: int,
    out_path=None,
    *,
    mol_index: int = 0,
    mol_name: Optional[str] = None,
    mode: str = "proportional",
) -> RestoreResult:
    """中和 (Σq≈0) された 1 分子 UDF の電荷を、 指定形式電荷になるよう逆変換して出力。

    過剰電荷を分散して中和した UDF (Σ電荷≈0) を入力に、 目標の **形式電荷 (整数) S**
    を与えると、 中和前の per-atom 電荷 (Σ=S) を復元して別 UDF に書き出す。
    ``electrostatic_Site`` のみ更新し、 座標・結合等は無改変。

    **中和ルール (``mode``) を正しく選ぶこと** — UDF を中和した方法に一致させる:

    - ``"proportional"`` (既定): 過剰分 S を ``|q|`` 比例で分散
      (`B_i = A_i − S·|A_i|/Σ|A|`)。 逆変換は二次方程式で λ を解く。
      ``|S| ≥ Σ|q|`` (符号反転) のケースは ``ValueError``。
    - ``"uniform"``: 過剰分 S を全原子に**均等**に分散 (`B_i = A_i − S/N`)。
      逆変換は `A_i = B_i + S/N` で常に厳密・一意 (符号問題なし)。

    Parameters
    ----------
    udf_path
        中和済み電荷を持つ 1 分子 UDF。
    formal_charge
        目標の形式電荷 (整数)。 出力 UDF の総電荷がこの値になる。
    out_path
        出力 UDF。 省略時 ``<udf>_q<±S>.udf``。
    mol_index / mol_name
        対象分子 (既定は先頭=0、 単分子 UDF を想定)。
    mode
        中和ルール: ``"proportional"`` (既定) または ``"uniform"``。

    Returns
    -------
    RestoreResult
    """
    udf_path = Path(udf_path)
    S = int(formal_charge)
    if out_path is None:
        out_path = udf_path.with_name(f"{udf_path.stem}_q{S:+d}{udf_path.suffix}")
    out_path = Path(out_path)

    tmpl = read_molecule_charges(udf_path, mol_index=mol_index, mol_name=mol_name)
    B = tmpl.charges
    input_total = float(sum(B))

    lam: Optional[float] = None
    shift: Optional[float] = None
    if mode == "proportional":
        lam = _solve_neutralization_lambda(B, S)
        A = [
            (b / (1.0 - lam) if b > 0.0 else (b / (1.0 + lam) if b < 0.0 else 0.0))
            for b in B
        ]
    elif mode == "uniform":
        shift = float(S) / tmpl.n_atoms
        A = [b + shift for b in B]
    else:
        raise ValueError(
            f"mode={mode!r} は未対応です ('proportional' / 'uniform' から選択)"
        )
    output_total = float(sum(A))

    u = _open(udf_path)
    _put_molecule_charges(u, tmpl.source_index, A)
    u.write(str(out_path))

    res = RestoreResult(
        out_path=str(out_path),
        mol_name=tmpl.mol_name,
        n_atoms=tmpl.n_atoms,
        formal_charge=S,
        input_total=input_total,
        output_total=output_total,
        mode=mode,
        lam=lam,
        shift=shift,
    )
    detail = f"λ={lam:.8f}" if mode == "proportional" else f"shift=S/N={shift:.8f}"
    logger.info(
        "restore[%s]: %r %d atoms, Σq %.6f → %.6f (target %d, %s) → %s",
        mode, res.mol_name, res.n_atoms, input_total, output_total, S, detail,
        out_path.name,
    )
    return res
