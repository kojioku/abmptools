# -*- coding: utf-8 -*-
"""
udf_adapter.py
--------------
Adapter layer: reads a COGNAC UDF file via UDFManager and builds a SystemModel.

The public entry point is::

    model = UdfAdapter(udf_manager).build()

All UDFManager calls are contained here; no UDF access occurs in the Writers.
"""
from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass
from typing import List, Optional, Tuple

from ..core.system_model import (
    AtomType, AtomRecord, BondRecord, PairRecord, AngleRecord, DihedralRecord,
    MoleculeTopology, AtomPosition, CellGeometry, NdxData,
    SimulationParams, SystemModel,
)

logger = logging.getLogger(__name__)

# Charge conversion constant (same as original script)
_CHARGE_UNIT = 18.224159264


# ---------------------------------------------------------------------------
# Helpers (ported from udf2gro.py, but now pure functions / instance methods)
# ---------------------------------------------------------------------------

def _sanitize_gromacs_molname(name: str) -> str:
    """Sanitize a UDF Mol_Name into a string usable as GROMACS
    ``[ moleculetype ] Name``.

    GROMACS doesn't enforce a strict length on mol type names, but the
    name must not contain whitespace or commenting characters. We
    replace anything outside ``[A-Za-z0-9_]`` with ``_`` and avoid
    leading-digit names by prepending ``M_`` if needed. Empty or all-
    underscore inputs fall back to ``MOL``.
    """
    if not name:
        return "MOL"
    out = "".join(c if (c.isalnum() or c == "_") else "_" for c in name)
    if out and out[0].isdigit():
        out = "M_" + out
    if not out.replace("_", ""):
        return "MOL"
    return out


def _shorten_molname(name: str) -> str:
    """.gro の残基名欄は 5 文字。 溢れる分は真ん中を落として端を残す。"""
    if len(name) > 5:
        return name[0] + name[1] + name[-3] + name[-2] + name[-1]
    return name


def _float2str(value: float, ndigit: int) -> str:
    return "{:.15f}".format(value)[:ndigit]


def _is_rectangular(cell_raw, thres: float = 1e-5) -> bool:
    for check in [cell_raw[3], cell_raw[4], cell_raw[5]]:
        if abs(check - 90.0) > thres:
            return False
    return True


def _cosine_polynomial_to_proper_dihedral(udf, torsion_index: int,
                                          n_terms: int, energy_scale: float):
    """COGNAC の余弦多項式を GROMACS の proper dihedral (funct 9) に直す。

    COGNAC は ``Cosine_Polynomial.p[]`` の係数 ``A_0..A_{n-1}`` で
    ``sum A_i cos^i(phi)`` と書く。 GROMACS は ``k(1 + cos(n*phi - phi_s))``
    なので、 次数ごとに ``(phi_s, k, n)`` へ読み替える。 位相は係数の符号で
    0 度か 180 度のどちらかになり、 係数の分母は倍角公式から出る。

    Parameters
    ----------
    torsion_index : ``Torsion_Potential[]`` の何番目か
    n_terms       : 係数の数。 多重度は ``n_terms - 1``
    energy_scale  : UDF のエネルギー単位 -> kJ/mol の換算係数

    Returns
    -------
    ``(phi_s [deg], k [kJ/mol], multiplicity)``
    """
    j, n, kk = torsion_index, n_terms, energy_scale
    loc = "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.p[]"
    gro_k = 0.0
    gro_phi = 0.0
    gro_mult = 0
    if n == 1:
        gro_mult = 0
        A0 = udf.get(loc, [j, 0])
        gro_phi = 0.0 if A0 != 0 else 180.0
        gro_k = A0 / 2 * kk
    elif n == 2:
        gro_mult = 1
        A0 = udf.get(loc, [j, 0])
        A1 = udf.get(loc, [j, 1])
        A01 = A0 * A1
        gro_phi = 0.0 if A01 < 0 else 180.0
        gro_k = A0 * kk
    elif n == 3:
        gro_mult = 2
        A0 = udf.get(loc, [j, 0])
        A2 = udf.get(loc, [j, 2])
        if A0 == 0:
            gro_phi = 0.0
            gro_k = A2 / 2 * kk
        else:
            gro_phi = 180.0
            gro_k = A0 / 2 * kk
    elif n == 4:
        gro_mult = 3
        A0 = udf.get(loc, [j, 0])
        A1 = udf.get(loc, [j, 1])
        A01 = A0 * A1
        gro_phi = 0.0 if A01 >= 0 else 180.0
        gro_k = A0 * kk
    elif n == 5:
        gro_mult = 4
        A0 = udf.get(loc, [j, 0])
        A2 = udf.get(loc, [j, 2])
        if A0 != 0:
            gro_phi = 0.0
            gro_k = A0 / 2 * kk
        else:
            gro_phi = 180.0
            gro_k = A2 / 8 * kk
    elif n == 6:
        gro_mult = 5
        A0 = udf.get(loc, [j, 0])
        A1 = udf.get(loc, [j, 1])
        A01 = A0 * A1
        gro_phi = 0.0 if A01 < 0 else 180.0
        gro_k = A0 * kk
    elif n == 7:
        gro_mult = 6
        A0 = udf.get(loc, [j, 0])
        A2 = udf.get(loc, [j, 2])
        if A0 == 0:
            gro_phi = 0.0
            gro_k = A2 / 18 * kk
        else:
            gro_phi = 180.0
            gro_k = A0 / 2 * kk
    return gro_phi, gro_k, gro_mult


# ---------------------------------------------------------------------------
# Main adapter class
# ---------------------------------------------------------------------------

#: Q を持つ Nose-Hoover 系のアルゴリズム。 いずれも同じ式で tau_t を出す。
_NOSE_HOOVER_ALGORITHMS = (
    "NVT_Nose_Hoover",
    "NPT_Parrinello_Rahman_Nose_Hoover",
    "NPT_Andersen_Nose_Hoover",
)

#: Boltzmann 定数 [amu nm^2 / (ps^2 K)]。 amu/nm/ps 系では energy = kJ/mol。
_KB_AMU_NM2_PS2_K = 0.0083144626

#: 水の等温圧縮率 [bar^-1]。 GROMACS の compressibility の慣用既定値。
_COMPRESSIBILITY_BAR = 4.5e-5

#: 1 bar を kJ/mol/nm^3 に直す係数。
_BAR_TO_KJ_MOL_NM3 = 0.06022140762

#: tau_p の既定 [ps]。 Parrinello-Rahman の実用域 (2-5 ps) の下端。
DEFAULT_TAU_P_PS = 2.0

#: ``--barostat`` で受ける名前 -> GROMACS の ``pcoupl`` に書く綴り。
#: 綴りを間違えると grompp が落ちるだけなので、 ここで正規化して弾く。
#: C-rescale は GROMACS 2021+。 COGNAC 側に対応概念が無いので指定でしか選べない。
BAROSTAT_NAMES = {
    "c-rescale":          "C-rescale",
    "crescale":           "C-rescale",
    "parrinello-rahman":  "Parrinello-Rahman",
    "parrinello_rahman":  "Parrinello-Rahman",
    "berendsen":          "Berendsen",
    "mttk":               "MTTK",
    "no":                 "no",
}


#: GROMACS で実際に使われる tau_p の範囲 [ps]。 目安であって制約ではない。
PRACTICAL_TAU_P_PS = (2.0, 5.0)


def _check_tau_pair(t_coupl, tau_t, p_coupl, tau_p) -> None:
    """``tau_t`` と ``tau_p`` の組合せを見て、 実行前に言えることを言う。

    UDF から出した値は「その UDF が記述している計算」であって、 GROMACS で
    普通に選ぶ値とは限らない。 ここで気付けないと、 短い NPT で箱がまだ緩和
    していないのに平衡とみなす、 といった形で静かに効く。
    """
    if str(p_coupl).lower() in ("", "no"):
        return

    lo, hi = PRACTICAL_TAU_P_PS
    if tau_p > hi:
        # Cell_Mass = 系の全質量 という COGNAC の慣用値のせいで、 tau_p は
        # 箱の一辺に比例して伸びる (5 nm 立方で PR 11 ps / Andersen 19 ps)。
        # 値としては UDF に忠実だが、 密度の緩和には数 x tau_p かかる。
        logger.warning(
            "tau_p = %.2f ps is longer than the %g-%g ps usually used in "
            "GROMACS. It is what the UDF's Cell_Mass asks for (COGNAC sets "
            "it to the system's total mass, so tau_p grows with the box). "
            "The box will relax slowly -- allow several times tau_p, or "
            "pass --tau-p %g.", tau_p, lo, hi, hi)

    # grompp が出すのと同じ条件。 こちらは両方の値を持っているので先に言える。
    if str(t_coupl).lower() == "nose-hoover" and tau_p < 2.0 * tau_t:
        logger.warning(
            "tau_p = %.3f ps is less than twice tau_t = %.3f ps; with "
            "nose-hoover this can resonate and grompp will say so. Raise "
            "tau_p with --tau-p, or lower tau_t with --tau-t.", tau_p, tau_t)


@dataclass
class ThermostatSettings:
    """``.mdp`` の熱浴まわり。 ``_build_thermostat`` が返す。"""
    t_coupl: str                      #: GROMACS の ``tcoupl``
    tau_t: float                      #: [ps]。 nose-hoover では振動の周期
    ref_t: float                      #: [K]


@dataclass
class BarostatSettings:
    """``.mdp`` の圧力浴まわり。 ``_build_barostat`` が返す。

    以前は 7 要素のタプルで返していたので、 呼び出し側が位置で数えていた。
    要素はどれも float か str なので、 順序を取り違えても型エラーにならず、
    ``.mdp`` に違う値が書かれるだけになる。
    """
    p_coupl: str                      #: GROMACS の ``pcoupl``
    pcoupltype: str                   #: isotropic / anisotropic
    tau_p: float                      #: [ps]
    ref_p: float                      #: [bar]
    ref_p_tensor: Optional[list]      #: anisotropic のときの 6 成分
    compressibility: float            #: [bar^-1]
    compressibility_tensor: Optional[list]

    @classmethod
    def none(cls, compressibility: float) -> "BarostatSettings":
        """圧力浴なし (NVE / NVT)。"""
        return cls("no", "isotropic", DEFAULT_TAU_P_PS, 1.0, None,
                   compressibility, None)


def canonical_barostat(name: str) -> str:
    """``--barostat`` の値を GROMACS の綴りに直す。 未知なら ValueError。"""
    try:
        return BAROSTAT_NAMES[str(name).strip().lower()]
    except KeyError:
        raise ValueError(
            "unknown barostat %r; choose one of %s"
            % (name, ", ".join(sorted(set(BAROSTAT_NAMES.values())))))


class UdfAdapter:
    """Reads a UDFManager object and produces a SystemModel."""

    #: Unit_Parameter が無い UDF に当てる既定のスケール。
    #: 全原子 UDF (AMBER / GAFF 系) の慣用単位 = 長さ Å, エネルギー kcal/mol。
    #: (Mass [amu], Energy [kJ/mol], Length [nm])
    ALL_ATOM_UNIT = (1.0, 4.184, 0.1)

    def __init__(self, udf, unit_parameter=None, tau_t=None, tau_p=None,
                 barostat=None):
        """
        Parameters
        ----------
        udf : UDFManager
        unit_parameter : tuple | str | None
            ``Unit_Parameter`` が UDF に無いときに使うスケール
            ``(Mass[amu], Energy[kJ/mol], Length[nm])``。
            ``"all_atom"`` で :data:`ALL_ATOM_UNIT` (Å / kcal/mol)。
            ``None`` かつ UDF にも無ければ **エラーにする** (黙って
            無次元値を GROMACS 単位として書き出さないため)。
        tau_t, tau_p : float | None
            熱浴 / 圧力浴の時定数 [ps] を直接指定する。 ``None`` なら
            ``tau_t`` は ``Q`` から算出し、 ``tau_p`` は既定値を使う。
        """
        #: --tau-t / --tau-p による上書き
        self._tau_t_override = tau_t
        self._tau_p_override = tau_p
        #: 圧力浴の上書き。 UDF に対応する概念が無いので指定でのみ効く。
        #: C-rescale (GROMACS 2021+) は Berendsen 並に安定で、
        #: かつ正しい NPT アンサンブルを与える。 拘束とも併用できる。
        self._barostat_override = barostat
        self._udf = udf
        if isinstance(unit_parameter, str):
            if unit_parameter != "all_atom":
                raise ValueError(
                    f"unit_parameter={unit_parameter!r} は不明です "
                    "('all_atom' か (mass, energy, length) のタプル)"
                )
            unit_parameter = self.ALL_ATOM_UNIT
        self._unit_parameter = unit_parameter

    # ------------------------------------------------------------------
    # Unit system
    # ------------------------------------------------------------------

    def _ensure_unit_parameter(self):
        """``Unit_Parameter`` が宣言されているか確かめる。

        本アダプタは値を ``udf.get(..., "[nm]")`` のように **単位を指定して**
        読んでおり、換算は UDFManager に任せている。その換算は UDF の
        ``Unit_Parameter`` (1 sigma が何 nm か、1 epsilon が何 kJ/mol か) に
        基づく。

        **``Unit_Parameter`` が無いと換算は黙って素通りする。** その場合
        Å の座標が nm、kcal/mol の epsilon が kJ/mol として書き出され、
        ``gmx grompp`` は形式が正しいので通してしまう。箱が 10 倍
        (= 密度 1/1000) の系が警告なしに走るので、極めて気付きにくい。

        そこで:

        * UDF が持っていればそれを使う (何もしない)
        * ``unit_parameter`` が渡されていればメモリ上に注入する
        * どちらも無ければ **エラーにする**
        """
        udf = self._udf
        declared = None
        try:
            declared = udf.get("Unit_Parameter.Length")
        except Exception:                       # noqa: BLE001 - 古い定義には無い
            declared = None

        if declared:
            logger.info(
                "Unit_Parameter: Mass=%s amu, Energy=%s kJ/mol, Length=%s nm (UDF 宣言値)",
                udf.get("Unit_Parameter.Mass"),
                udf.get("Unit_Parameter.Energy"),
                declared,
            )
            return

        if self._unit_parameter is None:
            raise RuntimeError(
                "この UDF は Unit_Parameter を宣言していません。\n"
                "UDFManager の単位換算はこれを基準に行うため、このまま変換すると\n"
                "Å の値が nm、kcal/mol の値が kJ/mol として書き出されます\n"
                "(grompp は通ってしまい、箱が 10 倍 = 密度 1/1000 の系が走ります)。\n"
                "\n"
                "全原子 UDF (GAFF 系、長さ Å・エネルギー kcal/mol) なら:\n"
                "    Exporter().export(udf, prefix, unit_parameter='all_atom')\n"
                "    python -m abmptools.udf2gro in.udf out --unit all_atom\n"
                "別のスケールなら (Mass[amu], Energy[kJ/mol], Length[nm]) を渡すか、\n"
                "UDF 側に Unit_Parameter を書いてください。"
            )

        mass, energy, length = self._unit_parameter
        udf.put(float(mass), "Unit_Parameter.Mass")
        udf.put(float(energy), "Unit_Parameter.Energy")
        udf.put(float(length), "Unit_Parameter.Length")
        logger.warning(
            "UDF に Unit_Parameter が無いので、指定された値を当てました "
            "(Mass=%s amu, Energy=%s kJ/mol, Length=%s nm)。"
            "UDF 側に書いておくと以後この指定は不要です",
            mass, energy, length,
        )

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def build(self) -> SystemModel:
        """UDFManagerからSystemModelを構築して返す。

        UDFの全データ (原子型、トポロジ、座標、シミュレーション条件) を読み取り、
        GROMACS変換用の中間表現にまとめる。
        """
        udf = self._udf
        self._ensure_unit_parameter()
        udf.jump(udf.totalRecord() - 1)
        logger.info("Data output: Record number = %s", udf.totalRecord() - 1)

        udf_path = str(udf.udfFile()).strip()
        title = udf_path.split("/")[-1] if "/" in udf_path else udf_path

        calcQQ = udf.get("Simulation_Conditions.Calc_Potential_Flags.Electrostatic")
        logger.info("Electrostatic is %s.", "ON" if calcQQ == 1 else "OFF")

        # --- counts ---
        atm_type_num = udf.size(udf.rlocation("Molecular_Attributes.Atom_Type[]", []))
        mol_num      = udf.size(udf.rlocation("Set_of_Molecules.molecule[]",       []))
        int_num      = udf.size(udf.rlocation("Interactions.Pair_Interaction[]",   []))
        bnd_type_num = udf.size(udf.rlocation("Molecular_Attributes.Bond_Potential[]",     []))
        agl_type_num = udf.size(udf.rlocation("Molecular_Attributes.Angle_Potential[]",    []))
        tor_type_num = udf.size(udf.rlocation("Molecular_Attributes.Torsion_Potential[]",  []))

        all_atm_num = sum(
            udf.size(udf.rlocation("Set_of_Molecules.molecule[].atom[]", [i]))
            for i in range(mol_num)
        )

        # --- OPLS detection ---
        totallyOPLS, partiallyOPLS = self._detect_opls()
        comb_rule = 3 if totallyOPLS else 2
        if totallyOPLS:
            logger.info("OPLS FF is detected. LJ mixing-rule -> geometric")
        elif partiallyOPLS:
            logger.warning("LJ mixing-rule cannot be identified!")

        # --- fudge factors ---
        flags14 = udf.get("Simulation_Conditions.Calc_Potential_Flags.Non_Bonding_1_4")
        fudgeLJ = udf.getArray("Interactions.Pair_Interaction[].Scale_1_4_Pair", [0])
        fudgeQQ = 0.0
        if calcQQ == 1:
            fudgeQQ = udf.getArray("Interactions.Electrostatic_Interaction[].Scale_1_4_Pair", [0])

        # --- atomname_in_gro map ---
        atomname_in_gro = self._build_atomname_map()

        # --- mol name mappings ---
        mol_name_list, molname_map, mol_name = self._build_mol_name_mappings(mol_num)

        mol_type_num = len(mol_name)
        logger.info("Molecular type's num = %s", mol_type_num)
        logger.info("Molecular's num      = %s", mol_num)
        logger.info("Atom's num   = %s", all_atm_num)
        logger.info("Atomtype num = %s", atm_type_num)
        logger.info("molecule name map (udf -> top)")
        for name in molname_map.keys():
            logger.info("%-15s -> %-5s", name, molname_map[name])

        # --- bond name -> potential type map ---
        bond_name_potmap = self._build_bond_name_potmap()

        # --- atom types with LJ params ---
        atom_types, lj_cutoff = self._build_atom_types(
            atm_type_num, int_num, totallyOPLS
        )

        # --- mol topologies ---
        mol_topologies = self._build_mol_topologies(
            mol_type_num, mol_name, mol_name_list, mol_num, molname_map,
            atomname_in_gro, calcQQ, bond_name_potmap,
            bnd_type_num, agl_type_num, tor_type_num,
            totallyOPLS, int_num
        )

        # --- mol sequence for [ molecules ] ---
        mol_sequence = self._build_mol_sequence(mol_num, mol_name_list, molname_map)

        # --- atom positions (gro structure) ---
        atom_positions, cell, vel_gen, udf_gro = self._build_atom_positions(
            mol_num, atomname_in_gro
        )

        # --- NDX data (optional) ---
        ndx_data = self._build_ndx_data(mol_num, mol_name_list, molname_map)

        # --- simulation params ---
        sim_params = self._build_sim_params(
            title, calcQQ, lj_cutoff, all_atm_num,
            vel_gen, cell, udf_gro
        )

        return SystemModel(
            title=title,
            udf_path=udf_path,
            comb_rule=comb_rule,
            flags14=flags14,
            fudgeLJ=fudgeLJ,
            fudgeQQ=fudgeQQ,
            calcQQ=calcQQ,
            atom_types=atom_types,
            mol_topologies=mol_topologies,
            mol_sequence=mol_sequence,
            atom_positions=atom_positions,
            cell=cell,
            sim_params=sim_params,
            ndx_data=ndx_data,
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _detect_opls(self):
        atypes = self._udf.get("Molecular_Attributes.Atom_Type[].Name")
        totallyOPLS = True
        partiallyOPLS = False
        for atype in atypes:
            if '&' in atype:
                partiallyOPLS = True
            else:
                totallyOPLS = False
        return totallyOPLS, partiallyOPLS

    def _build_atomname_map(self):
        udf = self._udf
        atomtypedata = udf.get("Molecular_Attributes.Atom_Type[].Name")
        atomname_in_gro = {}
        for i, name in enumerate(atomtypedata):
            atomname_in_gro[name] = name[:2] + hex(i)[2:]
        return atomname_in_gro

    def _build_mol_name_mappings(self, mol_num: int):
        udf = self._udf
        mol_name_list = list(udf.get("Set_of_Molecules.molecule[].Mol_Name"))
        map_molname_atypes = {}
        ncount_diff = 1

        for j, molname in enumerate(mol_name_list):
            atypelist = udf.get(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [j]
            )
            if molname in map_molname_atypes:
                if atypelist == map_molname_atypes[molname]:
                    pass
                else:
                    molname_new = self._molname_with_same_topology(
                        atypelist, map_molname_atypes
                    )
                    if molname_new is None:
                        molname_new = "".join(a[0] for a in atypelist)
                        map_molname_atypes[molname_new] = atypelist
                        ncount_diff += 1
                    logger.info("mol[%s], name=%s is renamed to %s",
                        j, molname, molname_new
                    )
                    mol_name_list[j] = molname_new
            else:
                map_molname_atypes[molname] = atypelist

        # GROMACS [ moleculetype ] Name は記号や空白を避ければ長さ制限は
        # 実質ない。元の UDF Mol_Name をそのまま使い、不正文字だけ '_' に
        # 置換する (これにより fcewsmb 等の下流が UDF Mol_Name と top の
        # mol type 名を直接一致させられる)。
        # 旧挙動 ("M0000", "M0001", ...) は同名ぶつかり時の自動 dedup と、
        # 5 chars 想定の古いツール用に残す形でフォールバック。
        molname_map = {}
        mol_name = []
        used = set()
        for j, molname in enumerate(mol_name_list):
            if molname in mol_name:
                continue
            sanitized = _sanitize_gromacs_molname(molname)
            if sanitized in used:
                # collision (元の名前が衝突する稀ケース) → "M{:04d}" にフォールバック
                sanitized = "M{:04d}".format(j)
            used.add(sanitized)
            molname_map[molname] = sanitized
            mol_name.append(molname)

        return mol_name_list, molname_map, mol_name

    @staticmethod
    def _molname_with_same_topology(atypes, map_name_atypes):
        for name, atypes_check in map_name_atypes.items():
            if atypes == atypes_check:
                return name
        return None

    def _build_bond_name_potmap(self):
        udf = self._udf
        location = "Molecular_Attributes.Bond_Potential[]"
        ntypes = udf.size(location)
        result = {}
        for itype in range(ntypes):
            name    = udf.get(location + ".Name",           [itype])
            pottype = udf.get(location + ".Potential_Type", [itype])
            result[name] = pottype
        return result

    def _build_atom_types(self, atm_type_num: int, int_num: int,
                          totallyOPLS: bool):
        udf = self._udf
        atom_types = []
        lj_cutoff = 0.0

        for i in range(atm_type_num):
            name = udf.getArray("Molecular_Attributes.Atom_Type[].Name", [i])
            mass = udf.getArray("Molecular_Attributes.Atom_Type[].Mass", [i])

            sigma   = 0.0
            epsilon = 0.0
            found   = False
            for j in range(int_num):
                site1   = udf.getArray("Interactions.Pair_Interaction[].Site1_Name", [j])
                site2   = udf.getArray("Interactions.Pair_Interaction[].Site2_Name", [j])
                pottype = udf.get("Interactions.Pair_Interaction[].Potential_Type", [j])
                if pottype != "Lennard_Jones":
                    logger.error("non LJ(12,6) type potential is detected")
                    raise RuntimeError("non-LJ potential")
                if site1 == name and site1 == site2:
                    sigma   = udf.get("Interactions.Pair_Interaction[].Lennard_Jones.sigma",   [j], "[nm]")
                    epsilon = udf.get("Interactions.Pair_Interaction[].Lennard_Jones.epsilon",  [j], "[kJ/mol]")
                    cc = udf.get("Interactions.Pair_Interaction[].Cutoff", [j], "[nm]")
                    if cc > lj_cutoff:
                        lj_cutoff = cc
                    found = True
                    break

            atom_types.append(AtomType(
                name=name, mass=mass, sigma=sigma, epsilon=epsilon
            ))

        return atom_types, lj_cutoff

    # ------------------------------------------------------------------
    # Molecule topology building
    # ------------------------------------------------------------------

    def _build_mol_topologies(
        self, mol_type_num, mol_name, mol_name_list, mol_num, molname_map,
        atomname_in_gro, calcQQ, bond_name_potmap,
        bnd_type_num, agl_type_num, tor_type_num,
        totallyOPLS, int_num
    ):
        udf = self._udf
        topologies = []

        for m in range(mol_type_num):
            udf_name = mol_name[m]
            gro_raw  = molname_map[udf_name]
            # [ moleculetype ] Name は長さ制限が実質ない (whitespace と
            # コメント文字だけ避ければ良い)。`_sanitize_gromacs_molname`
            # で既に GROMACS-safe にしてあるのでそのまま使う。
            # ``_shorten_molname`` (5 chars) は .gro 残基列向けで、
            # AtomPosition 構築時に別途適用される (lines 678-688 参照)。
            gro_name = gro_raw
            logger.info("moltype %d %s", m, gro_name)

            # find representative molecule index
            mol_no = next(
                i for i in range(mol_num)
                if mol_name_list[i] == udf_name
            )

            topo = MoleculeTopology(
                udf_name=udf_name,
                gro_name=gro_name,
                nrexcl=3,
            )

            # --- atoms ---
            topo.atoms = self._build_atoms(
                mol_no, gro_name, atomname_in_gro, calcQQ
            )

            # --- bonds ---
            topo.bonds = self._build_bonds(
                mol_no, bond_name_potmap, bnd_type_num
            )

            # --- pairs (1-4) ---
            topo.pairs = self._build_pairs(mol_no, int_num)

            # --- angles ---
            topo.angles = self._build_angles(mol_no, agl_type_num)

            # --- dihedrals ---
            topo.dihedrals = self._build_dihedrals(mol_no, tor_type_num)

            topologies.append(topo)

        return topologies

    def _build_atoms(self, mol_no, gro_name, atomname_in_gro, calcQQ):
        udf = self._udf
        atm_num = udf.size("Set_of_Molecules.molecule[].atom[]", [mol_no])
        ele_num = udf.size("Set_of_Molecules.molecule[].electrostatic_Site[]", [mol_no])
        logger.info("atom_num     = %d", atm_num)

        atoms = []
        for i in range(atm_num):
            atype = udf.get(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [mol_no, i]
            )
            gro_aname = atomname_in_gro[atype]

            if ele_num == 0 or calcQQ == 0:
                charge = 0.0
            else:
                q_raw = udf.get(
                    "Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element",
                    [mol_no, i]
                )
                charge = round(q_raw / _CHARGE_UNIT, 8)

            atoms.append(AtomRecord(
                index=i + 1,
                type_name=atype,
                gro_name=gro_aname,
                charge=charge,
            ))
        return atoms

    def _build_bonds(self, mol_no, bond_name_potmap, bnd_type_num):
        udf = self._udf
        bnd_num = udf.size(udf.rlocation("Set_of_Molecules.molecule[].bond[]", [mol_no]))
        logger.info("bond_num     = %d", bnd_num)

        bonds = []
        for i in range(bnd_num):
            a1 = udf.getArray("Set_of_Molecules.molecule[].bond[].atom1", [mol_no, i]) + 1
            a2 = udf.getArray("Set_of_Molecules.molecule[].bond[].atom2", [mol_no, i]) + 1
            pot_name = udf.getArray(
                "Set_of_Molecules.molecule[].bond[].Potential_Name", [mol_no, i]
            )
            pottype = bond_name_potmap[pot_name]

            funct = "1"
            r0 = 0.0
            kb = 0.0
            for j in range(bnd_type_num):
                name = str(udf.getArray("Molecular_Attributes.Bond_Potential[].Name", [j]))
                if name == pot_name:
                    if pottype == "FENE_LJ":
                        funct = "7"
                        r0 = udf.get(
                            "Molecular_Attributes.Bond_Potential[].FENE_LJ.R_max",
                            [j], "[sigma]"
                        )
                        kb = udf.get(
                            "Molecular_Attributes.Bond_Potential[].FENE_LJ.K",
                            [j], "[epsilon/(sigma^2)]"
                        )
                    elif pottype == "Harmonic":
                        funct = "1"
                        r0 = udf.get("Molecular_Attributes.Bond_Potential[].R0",
                                     [j], "[nm]")
                        kb = udf.get("Molecular_Attributes.Bond_Potential[].Harmonic.K",
                                     [j], "[kJ/mol/nm^2]")
                    break

            bonds.append(BondRecord(atom1=a1, atom2=a2, funct=funct, r0=r0, kb=kb))
        return bonds

    def _build_pairs(self, mol_no: int, int_num: int) -> List[PairRecord]:
        """Build 1-4 pair list using the same algorithm as get14() in original."""
        udf = self._udf
        atom_num = udf.size("Set_of_Molecules.molecule[].atom[]", [mol_no])

        bond_map  = [[] for _ in range(atom_num)]
        angle_map = [[] for _ in range(atom_num)]

        bond_num = udf.size("Set_of_Molecules.molecule[].bond[]", [mol_no])
        for j in range(bond_num):
            a1 = udf.get("Set_of_Molecules.molecule[].bond[].atom1", [mol_no, j])
            a2 = udf.get("Set_of_Molecules.molecule[].bond[].atom2", [mol_no, j])
            bond_map[a1].append(a2)
            bond_map[a2].append(a1)

        angle_num = udf.size("Set_of_Molecules.molecule[].angle[]", [mol_no])
        for j in range(angle_num):
            a1 = udf.get("Set_of_Molecules.molecule[].angle[].atom1", [mol_no, j])
            a3 = udf.get("Set_of_Molecules.molecule[].angle[].atom3", [mol_no, j])
            angle_map[a1].append(a3)
            angle_map[a3].append(a1)

        pair1: List[int] = []
        pair2: List[int] = []

        tors_num = udf.size("Set_of_Molecules.molecule[].torsion[]", [mol_no])
        for j in range(tors_num):
            a1 = udf.get("Set_of_Molecules.molecule[].torsion[].atom1", [mol_no, j])
            a3 = udf.get("Set_of_Molecules.molecule[].torsion[].atom3", [mol_no, j])
            a4 = udf.get("Set_of_Molecules.molecule[].torsion[].atom4", [mol_no, j])
            if a3 not in bond_map[a1] and a4 not in angle_map[a1]:
                chk = True
                for m in range(len(pair1)):
                    if ((pair1[m] == a1 and pair2[m] == a4) or
                            (pair2[m] == a1 and pair1[m] == a4)):
                        chk = False
                        break
                if chk:
                    pair1.append(a1)
                    pair2.append(a4)

        del bond_map
        logger.info("1-4 pair_num = %s", len(pair1))

        site1names = udf.get("Interactions.Pair_Interaction[].Site1_Name")
        site2names = udf.get("Interactions.Pair_Interaction[].Site2_Name")

        pairs = []
        for i in range(len(pair1)):
            name1 = udf.getArray(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name",
                [mol_no, pair1[i]]
            )
            name2 = udf.getArray(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name",
                [mol_no, pair2[i]]
            )
            found = False
            for j in range(int_num):
                s1 = site1names[j]
                s2 = site2names[j]
                if (s1 == name1 and s2 == name2) or (s1 == name2 and s2 == name1):
                    found = True
                    break
            # Always add pair (found or not) — matches original behaviour
            pairs.append(PairRecord(atom1=pair1[i] + 1, atom2=pair2[i] + 1))

        return pairs

    def _build_angles(self, mol_no: int, agl_type_num: int) -> List[AngleRecord]:
        udf = self._udf
        agl_num = udf.size(udf.rlocation("Set_of_Molecules.molecule[].angle[]", [mol_no]))
        logger.info("angle_num    = %d", agl_num)

        angles = []
        for i in range(agl_num):
            a1 = udf.getArray("Set_of_Molecules.molecule[].angle[].atom1", [mol_no, i]) + 1
            a2 = udf.getArray("Set_of_Molecules.molecule[].angle[].atom2", [mol_no, i]) + 1
            a3 = udf.getArray("Set_of_Molecules.molecule[].angle[].atom3", [mol_no, i]) + 1
            pot_name = udf.getArray(
                "Set_of_Molecules.molecule[].angle[].Potential_Name", [mol_no, i]
            )
            theta0 = 0.0
            k = 0.0
            for j in range(agl_type_num):
                name = str(udf.getArray("Molecular_Attributes.Angle_Potential[].Name", [j]))
                if name == pot_name:
                    theta0 = 180.0 - float(udf.get(
                        "Molecular_Attributes.Angle_Potential[].theta0", [j], "[degree]"
                    ))
                    k = udf.get(
                        "Molecular_Attributes.Angle_Potential[].Theta.K",
                        [j], "[kJ/mol/rad^2]"
                    )
                    break
            angles.append(AngleRecord(atom1=a1, atom2=a2, atom3=a3,
                                       theta0=theta0, k=k))
        return angles

    def _build_dihedrals(self, mol_no: int, tor_type_num: int) -> List[DihedralRecord]:
        udf = self._udf
        tor_num = udf.size(udf.rlocation("Set_of_Molecules.molecule[].torsion[]", [mol_no]))
        logger.info("torsion_num  = %d", tor_num)

        tors_pot_type_list = udf.get("Molecular_Attributes.Torsion_Potential[].Potential_Type")
        use_amber_header = "Amber" in tors_pot_type_list

        dihedrals = []
        for i in range(tor_num):
            a1 = udf.getArray("Set_of_Molecules.molecule[].torsion[].atom1", [mol_no, i]) + 1
            a2 = udf.getArray("Set_of_Molecules.molecule[].torsion[].atom2", [mol_no, i]) + 1
            a3 = udf.getArray("Set_of_Molecules.molecule[].torsion[].atom3", [mol_no, i]) + 1
            a4 = udf.getArray("Set_of_Molecules.molecule[].torsion[].atom4", [mol_no, i]) + 1
            pot_name = udf.getArray(
                "Set_of_Molecules.molecule[].torsion[].Potential_Name", [mol_no, i]
            )

            funct = ""
            params: List[float] = []

            for j in range(tor_type_num):
                name    = str(udf.getArray("Molecular_Attributes.Torsion_Potential[].Name",          [j]))
                pottype = str(udf.getArray("Molecular_Attributes.Torsion_Potential[].Potential_Type", [j]))
                if name != pot_name:
                    continue

                if pottype == "Cosine_Polynomial":
                    kk = float(udf.get(
                        "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.K",
                        [j], "[kJ/mol]"
                    ))
                    n = udf.getArray(
                        "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.N", [j]
                    )
                    if "oopa" in name:
                        funct = "4"
                        gro_phi, gro_k, gro_mult = _cosine_polynomial_to_proper_dihedral(udf, j, n, kk)
                        params = [gro_phi, gro_k, gro_mult]
                    elif n <= 6:
                        funct = "3"
                        raw_params = udf.get(
                            "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.p[]",
                            [j]
                        )
                        for k_idx in range(6):
                            if k_idx < len(raw_params):
                                params.append(kk * float(raw_params[k_idx]))
                            else:
                                params.append(0.0)
                    else:
                        funct = "1"
                        gro_phi, gro_k, gro_mult = _cosine_polynomial_to_proper_dihedral(udf, j, n, kk)
                        params = [gro_phi, gro_k, gro_mult]
                    break

                elif pottype == "Amber":
                    if ":" in name:
                        funct = "9"
                    else:
                        funct = "1"
                    pk    = udf.get("Molecular_Attributes.Torsion_Potential[].Amber.PK",    [j], "[kJ/mol]")
                    idivf = udf.get("Molecular_Attributes.Torsion_Potential[].Amber.IDIVF", [j])
                    pn    = udf.get("Molecular_Attributes.Torsion_Potential[].Amber.PN",    [j])
                    phase = udf.get("Molecular_Attributes.Torsion_Potential[].Amber.PHASE", [j])
                    k_val = pk / idivf
                    params = [phase, k_val, pn]
                    break

            dihedrals.append(DihedralRecord(
                atom1=a1, atom2=a2, atom3=a3, atom4=a4,
                funct=funct, params=params
            ))
        return dihedrals

    # ------------------------------------------------------------------
    # Mol sequence ([ molecules ] section)
    # ------------------------------------------------------------------

    def _build_mol_sequence(self, mol_num, mol_name_list, molname_map):
        # [ molecules ] section の mol 名は [ moleculetype ] Name と一致させる
        # 必要があるため、ここでは shorten しない (moleculetype 側も
        # _build_mol_topologies で full name を使うようにした)。
        # ``_shorten_molname`` は .gro 残基列 (5 chars) 専用に残してある。
        sequence = []
        if mol_num > 1:
            ncount = 1
            for i in range(1, mol_num):
                if mol_name_list[i] == mol_name_list[i - 1]:
                    ncount += 1
                    if i == mol_num - 1:
                        gn = molname_map[mol_name_list[i]]
                        sequence.append((gn, ncount))
                else:
                    gn = molname_map[mol_name_list[i - 1]]
                    sequence.append((gn, ncount))
                    ncount = 1
                    if i == mol_num - 1:
                        gn = molname_map[mol_name_list[i]]
                        sequence.append((gn, ncount))
        else:
            gn = molname_map[mol_name_list[0]]
            sequence.append((gn, 1))
        return sequence

    # ------------------------------------------------------------------
    # Atom positions (GRO structure)
    # ------------------------------------------------------------------

    def _build_atom_positions(self, mol_num, atomname_in_gro):
        udf = self._udf
        vel_gen = False
        udf_gro = udf

        if udf.get("Initial_Structure.Generate_Method.Method") == "Restart":
            rest_udfname = udf.get("Initial_Structure.Generate_Method.Restart.UDF_Name")
            rest_record  = udf.get("Initial_Structure.Generate_Method.Restart.Record")
            if len(rest_udfname) > 0:
                if len(os.path.dirname(rest_udfname)) == 0:
                    udfdir = udf.udfDirectory()
                    rest_udfpath = os.path.join(udfdir, rest_udfname)
                else:
                    rest_udfpath = rest_udfname
                logger.info("Restart geometry is read from record %s in '%s'",
                    rest_record, rest_udfname
                )
                if not os.path.exists(rest_udfpath):
                    raise RuntimeError(
                        "Error! UDF file for restart does not exist!\n"
                        "       Check UDF path 'Initial_Structure.Generate_Method.Restart.UDF_Name'"
                    )
                from UDFManager import UDFManager
                udf_gro = UDFManager(rest_udfpath)
                nrecs = udf_gro.totalRecord()
                if rest_record == -1:
                    udf_gro.jump(nrecs - 1)
                else:
                    udf_gro.jump(rest_record)

        positions = []
        for i in range(mol_num):
            atm_num = udf_gro.size(udf_gro.rlocation("Structure.Position.mol[].atom[]", [i]))
            mol_id = min(i + 1, 99999)
            ss = udf_gro.get("Set_of_Molecules.molecule[].Mol_Name", [i])
            if len(ss) > 5:
                ss = ss[0] + ss[1] + ss[2] + ss[len(ss)-2] + ss[len(ss)-1]

            for j in range(atm_num):
                x  = udf_gro.get("Structure.Position.mol[].atom[].x",  [i, j], "[nm]")
                y  = udf_gro.get("Structure.Position.mol[].atom[].y",  [i, j], "[nm]")
                z  = udf_gro.get("Structure.Position.mol[].atom[].z",  [i, j], "[nm]")
                vx = udf_gro.get("Structure.Velocity.mol[].atom[].x",  [i, j], "[nm/ps]")
                vy = udf_gro.get("Structure.Velocity.mol[].atom[].y",  [i, j], "[nm/ps]")
                vz = udf_gro.get("Structure.Velocity.mol[].atom[].z",  [i, j], "[nm/ps]")
                if vx is None or vy is None or vz is None:
                    vel_gen = True
                    vx = vy = vz = 0.0

                atom_type_name = udf_gro.get(
                    "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [i, j]
                )
                atom_gro_name = atomname_in_gro[atom_type_name]
                atom_id = min(
                    udf_gro.get("Set_of_Molecules.molecule[].atom[].Atom_ID", [i, j]) + 1,
                    99999
                )

                positions.append(AtomPosition(
                    mol_id=mol_id,
                    mol_name_short=ss,
                    atom_gro_name=atom_gro_name,
                    atom_id=atom_id,
                    x=x, y=y, z=z,
                    vx=vx, vy=vy, vz=vz,
                ))

        # Cell geometry
        cell_raw = udf_gro.get("Structure.Unit_Cell.Cell_Size")
        if _is_rectangular(cell_raw):
            a = udf_gro.get("Structure.Unit_Cell.Cell_Size.a", "[nm]")
            b = udf_gro.get("Structure.Unit_Cell.Cell_Size.b", "[nm]")
            c = udf_gro.get("Structure.Unit_Cell.Cell_Size.c", "[nm]")
            cell = CellGeometry(a=a, b=b, c=c)
        else:
            cell = CellGeometry(
                a=cell_raw[0] * 0.1,
                b=cell_raw[1] * 0.1,
                c=cell_raw[2] * 0.1,
                alpha=cell_raw[3],
                beta=cell_raw[4],
                gamma=cell_raw[5],
            )

        return positions, cell, vel_gen, udf_gro

    # ------------------------------------------------------------------
    # NDX data
    # ------------------------------------------------------------------

    def _build_ndx_data(self, mol_num, mol_name_list, molname_map) -> Optional[NdxData]:
        udf = self._udf
        loc_constr = "Simulation_Conditions.Constraint_Conditions.Constraint_Atom[]"
        ndata_constr = udf.size(loc_constr)
        if ndata_constr == 0:
            return None

        constr_axis_old = udf.get(loc_constr + ".Constraint_Axis", [0])
        constr_aid_list = []
        constr_axis = None

        for idat in range(ndata_constr):
            idx_mol_atom  = udf.get(loc_constr + ".Index",           [idat])
            constr_axis   = udf.get(loc_constr + ".Constraint_Axis", [idat])
            constr_method = udf.get(loc_constr + ".Method",          [idat])
            if constr_method != "Steady":
                raise RuntimeError(
                    "Error! constraint method {} is not supported!".format(constr_method)
                )
            constr_veloc = udf.get(loc_constr + ".Steady.Velocity", [idat])
            molidx  = idx_mol_atom[0]
            atomidx = idx_mol_atom[1]
            atom_id = udf.get(
                "Set_of_Molecules.molecule[].atom[].Atom_ID", [molidx, atomidx]
            ) + 1
            if constr_veloc[0] == 0.0 and constr_veloc[1] == 0.0 and constr_veloc[2] == 0.0:
                if constr_axis == constr_axis_old:
                    constr_axis_old = constr_axis
                    constr_aid_list.append(atom_id)
                else:
                    raise RuntimeError("Error! Not supported constraint conditions!")
            else:
                raise RuntimeError("Error! constraint velocity other than 0 is not supported!")

        atm_num_last = udf.size("Structure.Position.mol[].atom[]", [mol_num - 1])
        atom_id_max = udf.get(
            "Set_of_Molecules.molecule[].atom[].Atom_ID",
            [mol_num - 1, atm_num_last - 1]
        ) + 1

        molnames = set(mol_name_list)
        mol_groups = {}
        for molname in molnames:
            gn = _shorten_molname(molname_map[molname])
            aid_list = []
            for imol in range(mol_num):
                if udf.get("Set_of_Molecules.molecule[].Mol_Name", [imol]) == molname:
                    natoms = udf.size("Structure.Position.mol[].atom[]", [imol])
                    for iatom in range(natoms):
                        aid = udf.get(
                            "Set_of_Molecules.molecule[].atom[].Atom_ID", [imol, iatom]
                        ) + 1
                        aid_list.append(aid)
            mol_groups[gn] = aid_list

        return NdxData(
            atom_id_max=atom_id_max,
            mol_groups=mol_groups,
            constraint_atom_ids=constr_aid_list,
            constr_axis=constr_axis,
        )

    # ------------------------------------------------------------------
    # Simulation parameters
    # ------------------------------------------------------------------

    def _build_sim_params(
        self, title, calcQQ, lj_cutoff, all_atm_num,
        vel_gen, cell, udf_gro
    ) -> SimulationParams:
        udf = self._udf
        UDFVer = udf.getEngineVersion()

        algorithm = udf.get("Simulation_Conditions.Solver.Dynamics.Dynamics_Algorithm")
        if len(algorithm) == 0:
            raise RuntimeError("Error! Specify MD algorithm!")

        tail_correction = udf.get("Simulation_Conditions.Calc_Potential_Flags.Tail_Correction")
        pbc_a = udf.get("Simulation_Conditions.Boundary_Conditions.a_axis")
        pbc_b = udf.get("Simulation_Conditions.Boundary_Conditions.b_axis")
        pbc_c = udf.get("Simulation_Conditions.Boundary_Conditions.c_axis")

        fix_cell  = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Parrinello_Rahman_Nose_Hoover.Fix_Cell_Length")
        fix_angle = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Parrinello_Rahman_Nose_Hoover.Fix_Angle")

        deform, deform_npt, deform_vel = self._build_deformation(algorithm, UDFVer, cell)

        # --- vel_gen override from restart ---
        gen_method = udf.get("Initial_Structure.Generate_Method.Method")
        logger.debug("%s", gen_method)
        if gen_method == "Restart":
            restore_vel = udf.get("Initial_Structure.Generate_Method.Restart.Restore_Velocity")
            logger.debug("restore_vel %s", restore_vel)
            vel_gen = (restore_vel == 0)

        # --- integrator ---
        integrator, ld_seed = self._build_integrator(algorithm, deform_npt)

        # --- output intervals ---
        outputinterval = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps"
        )
        outputinterval2 = outputinterval
        if outputinterval >= 10000:
            outputinterval2 = int(outputinterval / 10)

        nsteps = udf.get("Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps")
        dt     = udf.get("Simulation_Conditions.Dynamics_Conditions.Time.delta_T", "[ps]")

        # --- constraints ---
        rattle_bond  = bool(udf.get("Simulation_Conditions.Dynamics_Conditions.RATTLE.Bond"))
        rattle_angle = bool(udf.get("Simulation_Conditions.Dynamics_Conditions.RATTLE.Angle"))

        # --- electrostatics ---
        qq_algorithm = ""
        if calcQQ == 1:
            qq_algorithm = udf.getArray("Interactions.Electrostatic_Interaction[].Algorithm", [0])

        # --- cutoffs ---
        cutoff_l = lj_cutoff
        cutoff_cl = cutoff_l
        if calcQQ == 1:
            cutoff_c = udf.get(
                "Interactions.Electrostatic_Interaction[].Cutoff_Coulomb.cutoff", [0], "[nm]"
            )
            if cutoff_c <= 0.0:
                cutoff_c = udf.get(
                    "Interactions.Electrostatic_Interaction[].Ewald.R_cutoff", [0], "[nm]"
                )
                if cutoff_c > 0.0:
                    logger.info("Ewald.R_cutoff was substituted for rcoulomb.")
            if cutoff_c <= cutoff_l:
                cutoff_c = cutoff_l
            cutoff_cl = max(cutoff_c, cutoff_l)
        else:
            cutoff_cl = cutoff_l

        # --- 熱浴・圧力浴 ---
        thermostat = self._build_thermostat(algorithm, all_atm_num)
        barostat = self._build_barostat(algorithm, fix_cell, fix_angle,
                                        deform_npt, deform_vel, cell)
        _check_tau_pair(thermostat.t_coupl, thermostat.tau_t,
                        barostat.p_coupl, barostat.tau_p)

        # --- pbc ---
        if pbc_a == "NONE" and pbc_b == "NONE" and pbc_c == "NONE":
            pbc = "no"
        elif pbc_a == "PERIODIC" and pbc_b == "PERIODIC" and pbc_c == "PERIODIC":
            pbc = "xyz"
        else:
            logger.error("Please Check Your Boundary Conditions")
            pbc = "xyz"

        periodic_mol = bool(udf.get("Simulation_Conditions.Boundary_Conditions.Periodic_Bond"))

        # --- constraint freeze ---
        freeze_grps = None
        freeze_dim  = None
        loc_constr  = "Simulation_Conditions.Constraint_Conditions.Constraint_Atom[]"
        ndata_constr = udf.size(loc_constr)
        if ndata_constr > 0:
            constr_axis = udf.get(loc_constr + ".Constraint_Axis", [0])
            freeze_grps = "Constraint"
            parts = []
            for each_axis in constr_axis:
                parts.append("Y" if each_axis == "YES" else "N")
            freeze_dim = " ".join(parts)

        # --- vel_gen temp ---
        gen_temp = None
        if vel_gen:
            gen_temp = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Temperature.Temperature", "[K]"
            )

        return SimulationParams(
            title=title,
            algorithm=algorithm,
            nsteps=nsteps,
            dt=dt,
            outputinterval=outputinterval,
            outputinterval2=outputinterval2,
            integrator=integrator,
            ld_seed=ld_seed,
            vel_gen=vel_gen,
            gen_temp=gen_temp,
            rattle_bond=rattle_bond,
            rattle_angle=rattle_angle,
            calcQQ=calcQQ,
            qq_algorithm=qq_algorithm,
            lj_cutoff=cutoff_cl,
            coulomb_cutoff=cutoff_cl,
            t_coupl=thermostat.t_coupl,
            tau_t=thermostat.tau_t,
            ref_t=thermostat.ref_t,
            p_coupl=barostat.p_coupl,
            pcoupltype=barostat.pcoupltype,
            tau_p=barostat.tau_p,
            ref_p=barostat.ref_p,
            ref_p_tensor=barostat.ref_p_tensor,
            compressibility=barostat.compressibility,
            compressibility_tensor=barostat.compressibility_tensor,
            tail_correction=tail_correction,
            pbc=pbc,
            periodic_mol=periodic_mol,
            deform_vel=deform_vel,
            freeze_grps=freeze_grps,
            freeze_dim=freeze_dim,
        )

    def _build_integrator(self, algorithm: str, deform_npt: bool):
        ld_seed = None
        if "NPT" in algorithm:
            if deform_npt:
                integrator = "md"
            elif "NPT_Andersen" in algorithm:
                integrator = "md-vv"
            else:
                integrator = "md"
        elif "Kremer_Grest" in algorithm:
            integrator = "sd"
            ld_seed = 1993
        else:
            integrator = "md-vv"
        return integrator, ld_seed

    def _build_deformation(self, algorithm, UDFVer, cell):
        udf = self._udf
        deform = udf.get("Simulation_Conditions.Dynamics_Conditions.Deformation.Method")
        list_supported_deform = ["Cell_Deformation"]
        list_supported_deform_method = ["Simple_Elongation", "Deformation_Rate"]
        deform_npt = False
        deform_vel = None

        # COGNAC writes "None" (not an empty string) when no cell deformation
        # is applied; treat both as "no deformation".
        if len(deform) == 0 or deform == "None":
            return deform, deform_npt, deform_vel

        cellsize = udf.get("Structure.Unit_Cell.Cell_Size")

        if deform in list_supported_deform:
            if "NPT" in algorithm:
                deform_npt = True
        else:
            raise RuntimeError(
                "Error!! deformation type '{}' is not supported.".format(deform)
            )

        deform_vel = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        deform_method = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Deformation.Cell_Deformation.Method"
        )
        if deform_method not in list_supported_deform_method:
            raise RuntimeError(
                "Error!! deformation method {} is not supported.".format(deform_method)
            )

        if deform_method == "Simple_Elongation":
            deform_axis = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Deformation"
                ".Cell_Deformation.Simple_Elongation.Axis"
            )
            poisson_ratio = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Deformation"
                ".Cell_Deformation.Simple_Elongation.Poisson_Ratio"
            )
            if UDFVer == "Ver101":
                input_method = udf.get(
                    "Simulation_Conditions.Dynamics_Conditions.Deformation"
                    ".Cell_Deformation.Simple_Elongation.Input_Method"
                )
                if input_method == "Deformation_Speed":
                    deform_rate = udf.get(
                        "Simulation_Conditions.Dynamics_Conditions.Deformation"
                        ".Cell_Deformation.Simple_Elongation.Deformation_Speed.Speed",
                        "[m/s]"
                    )
                elif input_method == "Initial_Strain_Rate":
                    if deform_axis == "z":
                        lz = cellsize[2] * 1e-10
                        rate = udf.get(
                            "Simulation_Conditions.Dynamics_Conditions.Deformation"
                            ".Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate",
                            "[1/s]"
                        )
                        deform_rate = lz * rate
            else:
                deform_rate = udf.get(
                    "Simulation_Conditions.Dynamics_Conditions.Deformation"
                    ".Cell_Deformation.Simple_Elongation.Elongation_Rate",
                    "[m/s]"
                )

            if poisson_ratio == 0.0:
                if deform_axis == "z":
                    deform_vel[2] = deform_rate * 0.001
                elif deform_axis == "xy":
                    deform_vel[0] = deform_rate * 0.001
                    deform_vel[1] = deform_rate * 0.001
                else:
                    raise RuntimeError(
                        'Error!! deform axis "{}" is unknown!'.format(deform_axis)
                    )
            else:
                raise RuntimeError("Error!! Poisson ratio other than 0.0 is not supported!")

        if deform_method == "Deformation_Rate":
            loc = "Simulation_Conditions.Dynamics_Conditions.Deformation.Cell_Deformation.Deformation_Rate"
            La = cellsize[0] * 1e-10
            Lb = cellsize[1] * 1e-10
            Lc = cellsize[2] * 1e-10
            deform_vel[0] = La * udf.get(loc + ".xx", "[Hz]") * 1e-3
            deform_vel[1] = Lb * udf.get(loc + ".yy", "[Hz]") * 1e-3
            deform_vel[2] = Lc * udf.get(loc + ".zz", "[Hz]") * 1e-3
            deform_vel[3] = Lb * udf.get(loc + ".xy", "[Hz]") * 1e-3
            deform_vel[4] = Lc * udf.get(loc + ".zx", "[Hz]") * 1e-3
            deform_vel[5] = Lc * udf.get(loc + ".yz", "[Hz]") * 1e-3

            if deform_npt:
                if deform_vel[3] != 0 or deform_vel[4] != 0 or deform_vel[5] != 0:
                    logger.error("shear deformation with NPT ensemble is not supported!")

        return deform, deform_npt, deform_vel

    def _nose_hoover_tau_t(self, Q, T_d, unit_Mass, unit_L, n_atoms):
        """COGNAC の熱浴質量 Q から GROMACS の ``tau_t`` [ps] を出す。

        **式**::

            tau_t = 2*pi * sqrt( Q_d / (g * k_B * T) )

            Q_d = Q * Unit_Parameter.Mass * Unit_Parameter.Length^2   [amu nm^2]
            g   = 3N   (comm-mode = None)  /  3N-3  (comm-mode = Linear)
            k_B = 0.0083144626 amu nm^2 / (ps^2 K)

        **根拠**。 COGNAC の Q は Nose-Hoover の熱浴質量 ``Q = g k_B T tau^2``
        で、 ``tau`` が応答時間。 GROMACS の ``tau_t`` は Nose-Hoover では
        緩和時間ではなく**運動エネルギー振動の周期**なので、 ``2*pi*tau`` に
        なる。 したがって ``tau_t = 2*pi*sqrt(Q_d/(g k_B T))``。

        **``g * k_B`` を落とすと ``tau_t`` が ``sqrt(3N)`` に比例して増大する**
        (3050 原子で 5.48 ps、 本式なら 0.63 ps)。 系を大きくするほど熱浴が
        鈍くなるので、 変換の検算はサイズを変えた 2 系で行うとよい。
        詳細は ``docs/udf2gro.md``。

        ``g`` は UDF の ``Dynamics_Conditions.Moment`` から決める。 重心運動を
        止める設定なら GROMACS 側も ``comm-mode = Linear`` になるので 3 を引く。
        """
        udf = self._udf
        Q_d = float(Q) * float(unit_Mass) * float(unit_L) ** 2
        g = max(1, 3 * int(n_atoms) - (3 if self._removes_com_motion() else 0))
        denom = g * _KB_AMU_NM2_PS2_K * float(T_d)
        if denom <= 0.0 or Q_d <= 0.0:
            return 0.1
        return 2.0 * math.pi * math.sqrt(Q_d / denom)

    def _removes_com_motion(self) -> bool:
        """GROMACS 側が ``comm-mode = Linear`` になる設定か。"""
        udf = self._udf
        base = "Simulation_Conditions.Dynamics_Conditions.Moment."
        try:
            return bool(udf.get(base + "Calc_Moment")) and \
                bool(udf.get(base + "Stop_Translation"))
        except Exception:                                # noqa: BLE001
            return False

    def _build_thermostat(self, algorithm, all_atm_num) -> ThermostatSettings:
        udf = self._udf
        t_coupl_map = {
            "NVE":                              "no",
            "NVT_Nose_Hoover":                  "nose-hoover",
            "NPT_Parrinello_Rahman_Nose_Hoover": "nose-hoover",
            "NPT_Andersen_Nose_Hoover":          "nose-hoover",
            "NVT_Berendsen":                    "berendsen",
            "NPT_Berendsen":                    "berendsen",
            "NVT_Kremer_Grest":                 "no",
        }
        t_coupl_str = t_coupl_map.get(algorithm, "no")
        if t_coupl_str == "no" and algorithm not in t_coupl_map:
            logger.error("algorithm %s is not supported !!", algorithm)

        if t_coupl_str == "no":
            return ThermostatSettings(t_coupl_str, 0.1, 300.0)

        T_d = udf.get("Simulation_Conditions.Dynamics_Conditions.Temperature.Temperature", "[K]")
        unit_L    = udf.get("Unit_Parameter.Length", "[nm]")
        unit_Mass = udf.get("Unit_Parameter.Mass",   "[amu]")

        tau_t = 0.1
        if algorithm in _NOSE_HOOVER_ALGORITHMS:
            Q = udf.get("Simulation_Conditions.Solver.Dynamics.%s.Q" % algorithm)
            tau_t = self._nose_hoover_tau_t(Q, T_d, unit_Mass, unit_L,
                                            all_atm_num)
        elif algorithm == "NVT_Berendsen":
            tau_t = udf.get("Simulation_Conditions.Solver.Dynamics.NVT_Berendsen.tau_T", "[ps]")
        elif algorithm == "NPT_Berendsen":
            tau_t = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Berendsen.tau_T", "[ps]")

        if self._tau_t_override is not None:
            logger.info("tau_t = %.4f ps (--tau-t で指定。 算出値 "
                        "%.4f ps を上書き)", self._tau_t_override, tau_t)
            tau_t = self._tau_t_override
        return ThermostatSettings(t_coupl_str, tau_t, T_d)

    #: COGNAC のアルゴリズム -> GROMACS の ``pcoupl``。
    _PCOUPL_BY_ALGORITHM = {
        "NVE":                               "no",
        "NVT_Nose_Hoover":                   "no",
        "NVT_Berendsen":                     "no",
        "NVT_Kremer_Grest":                  "no",
        "NPT_Parrinello_Rahman_Nose_Hoover": "Parrinello-Rahman",
        # Andersen は体積だけを動かすので、 セル行列全体を動かす PR ではなく
        # MTTK (isotropic) に対応する。 ただし MTTK は LINCS / SETTLE と
        # 併用できないので、 拘束のある系では --barostat で替える。
        "NPT_Andersen_Nose_Hoover":          "MTTK",
        "NPT_Berendsen":                     "berendsen",
    }

    def _barostat_name(self, algorithm) -> str:
        """``pcoupl`` に書く名前を決める。 ``--barostat`` があればそれが勝つ。"""
        name = self._PCOUPL_BY_ALGORITHM.get(algorithm, "no")

        if not self._barostat_override:
            return name

        # UDF (COGNAC) には C-rescale に対応する概念が無いので、 指定でしか
        # 選べない。 NVE/NVT の UDF に付けるとアンサンブルが変わる。
        chosen = canonical_barostat(self._barostat_override)
        if chosen != "no" and name == "no":
            logger.warning(
                "the UDF asks for no pressure coupling, but --barostat %s "
                "was given: the run will be NPT, not %s.",
                chosen, algorithm or "NVE")
        logger.info("pcoupl = %s (--barostat で指定)", chosen)
        return chosen

    def _barostat_tau_p(self, algorithm, unit_Mass, cell, beta_mdp):
        """COGNAC の ``Cell_Mass`` から GROMACS の ``tau_p`` [ps] を出す。

        **式**::

            PR       : tau_p = 2*pi * sqrt( C * beta / (3 * L) )
            Andersen : tau_p = 2*pi * sqrt( C * beta / L )        = sqrt(3) * PR

            C    = Cell_Mass * Unit_Parameter.Mass    [amu]
            L    = max(a, b, c)                       [nm]
            beta = .mdp に書く compressibility         [nm^3 mol / kJ]

        ``Cell_Mass`` はスキーマ上 ``[mass]`` (``def_udf/cognac*.udf``)。
        ``Q`` と違って ``sigma^2`` は掛けない。

        **根拠**。 どちらの側も ``tau_p`` は箱の振動の周期なので、 周期どうしを
        等置すれば出る。

        *GROMACS Parrinello-Rahman*::

            W^-1 = 4 pi^2 beta / (3 tau_p^2 L)      L = 最長のセル辺
            b_ddot = V W^-1 b'^-1 (P - P_ref)

        直方体 ``b = diag(a,b,c)`` で ``dP = -dV/(V beta)`` と線形化すると
        ``a_ddot = -4 pi^2 da / tau_p^2``、 つまり **``tau_p`` そのものが周期**。
        ``L`` が最長辺なのは GROMACS の定義どおりなので、 立方体でなくてよい。

        *COGNAC Parrinello-Rahman* (``COGNAC1124/src/PRsystem.cpp``)::

            cellMass = cellMassFactor;                                   (l.7)
            h2 = ((-currentStress - press0Tensor)*hinv*volume)/cellMass; (l.58)

        同じ線形化で係数を突き合わせると ``4 pi^2 beta/(3 tau_p^2 L) = 1/C``。

        *COGNAC Andersen* (``COGNAC1124/src/Anphsystem.cpp``)::

            cellMass = cellMassFactor * pow(volume, -4./3.);   (l.17)
            aVolume  = (currentPress - pressSum)/cellMass;     (l.122)

        体積を座標にした ``W V_ddot = dP``。 ``V = L^3`` として ``dL`` で
        書き直すと ``omega^2 = L/(C beta)`` で、 **PR より sqrt(3) 倍遅い**。
        言い換えると、 同じ周期を PR で出すには ``Cell_Mass`` を **3 倍**する。

        > 旧実装はここを ``1/3`` 倍にしていた (向きが逆)。 ``unit_L^2`` の
        > 取り違えと重なって下限 2.0 に丸められるため、 値としては表面化して
        > いなかった。 ``docs/udf2gro.md`` 参照。

        **実際の圧縮率は式から消える。** GROMACS 側の周期は
        ``tau_p * sqrt(beta_true/beta_mdp)`` なので、 等置すると ``beta_true``
        が両辺で相殺し、 ``.mdp`` に書く ``beta_mdp`` だけが残る。 系の本当の
        圧縮率を知らなくてよい。

        **近似**。 ``pcoupl = MTTK`` (Andersen の既定の行き先) はバロスタット
        質量の定義が違うので目安にとどまる。 ``Cell_Mass`` が無い / 0 の UDF
        では ``None`` を返す。
        """
        try:
            W = self._udf.get(
                "Simulation_Conditions.Solver.Dynamics.%s.Cell_Mass" % algorithm)
        except Exception:                                # noqa: BLE001
            return None
        if not W:
            return None

        C = float(W) * float(unit_Mass)                       # amu
        beta = float(beta_mdp) / _BAR_TO_KJ_MOL_NM3           # nm^3 mol / kJ
        max_L = max(float(cell.a), float(cell.b), float(cell.c))
        if C <= 0.0 or beta <= 0.0 or max_L <= 0.0:
            return None

        tau_p = 2.0 * math.pi * math.sqrt(C * beta / (3.0 * max_L))
        if algorithm == "NPT_Andersen_Nose_Hoover":
            tau_p *= math.sqrt(3.0)
        return tau_p

    def _build_barostat(self, algorithm, fix_cell, fix_angle,
                        deform_npt, deform_vel, cell) -> BarostatSettings:
        udf = self._udf
        commp = 0.000045  # bar^-1

        p_coupl_str = self._barostat_name(algorithm)

        if p_coupl_str == "no":
            return BarostatSettings.none(commp)

        # pcoupltype
        if self._barostat_override:
            pcoupltype = "isotropic"
        elif algorithm == "NPT_Parrinello_Rahman_Nose_Hoover":
            if fix_cell in ("", "xy", "yz", "zx", "x", "y", "z") or deform_npt:
                pcoupltype = "anisotropic"
            else:
                logger.error("Fix_Cell_Length %s is not supported!", fix_cell)
                pcoupltype = "isotropic"
        else:
            pcoupltype = "isotropic"

        unit_L    = udf.get("Unit_Parameter.Length", "[nm]")
        unit_Mass = udf.get("Unit_Parameter.Mass",   "[amu]")

        # cell dims for tau_p
        cell_x = cell.a
        cell_y = cell.b
        cell_z = cell.c
        max_L = max(cell_x, cell_y, cell_z)

        # tau_p は **どのバロスタットでも UDF から変換する**。 経路は 2 つ:
        #
        #   NPT_Berendsen  : UDF が tau_P を時間で持つ -> 単位換算だけ
        #   Andersen / PR  : UDF が Cell_Mass を質量で持つ -> 運動方程式経由
        #                    (_barostat_tau_p。 導出はそちらの docstring)
        #
        # 旧実装は後者を「Cell_Mass に Q 用の [mass*sigma^2] の換算係数を
        # 掛ける」という取り違えで出していた。 Cell_Mass はスキーマ上 [mass]
        # なので unit_L^2 (all-atom で 0.01) が余計で、 現実的な系では下限
        # 2.0 に丸められて値が表面化していなかった。 運動方程式を経由する
        # 形に直したので、 実機 (3050 原子) で Andersen 11.74 ps /
        # PR 6.78 ps が出る。
        #
        # 変換できないとき (Cell_Mass が無い) だけ既定 2.0 ps に落とす。
        # --tau-p はこのブロックの後で一度だけ適用する。
        tau_p = DEFAULT_TAU_P_PS

        if algorithm in ("NPT_Parrinello_Rahman_Nose_Hoover",
                         "NPT_Andersen_Nose_Hoover"):
            _converted = self._barostat_tau_p(algorithm, unit_Mass, cell, commp)
            if _converted is None:
                logger.warning(
                    "%s has no Cell_Mass; tau_p falls back to %.1f ps. "
                    "Pass --tau-p if the run needs a particular value.",
                    algorithm, DEFAULT_TAU_P_PS)
            else:
                tau_p = _converted
                logger.info("tau_p = %.2f ps (Cell_Mass から換算)", tau_p)
                if not 0.5 <= tau_p <= 20.0:
                    # 箱の応答が遅すぎる / 速すぎる。 Cell_Mass が既定
                    # (系の全質量) のままだと大きな系で数十 ps に伸びる。
                    logger.warning(
                        "tau_p = %.2f ps is outside the usual 0.5-20 ps. "
                        "This is what Cell_Mass = %s asks for; override with "
                        "--tau-p if it is not what you want.",
                        tau_p,
                        self._udf.get(
                            "Simulation_Conditions.Solver.Dynamics.%s.Cell_Mass"
                            % algorithm))
        elif algorithm == "NPT_Berendsen":
            # COGNAC が時間で持っているので、 単位換算だけして採用する。
            unit_P = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure", "[bar]"
            ) / udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure", "[P]"
            )
            tau_p_raw = udf.get(
                "Simulation_Conditions.Solver.Dynamics.NPT_Berendsen.tau_P", "[P*ps]"
            )
            tau_p = tau_p_raw * unit_P * commp

        if self._tau_p_override is not None:
            logger.info("tau_p = %.4f ps (--tau-p で指定。 算出値 "
                        "%.4f ps を上書き)", self._tau_p_override, tau_p)
            tau_p = self._tau_p_override

        ref_p, ref_p_tensor = self._reference_pressure(algorithm)
        comp_tensor = self._compressibility_tensor(
            algorithm, fix_cell, fix_angle, deform_npt, deform_vel, commp)

        return BarostatSettings(p_coupl_str, pcoupltype, tau_p, ref_p,
                                ref_p_tensor, commp, comp_tensor)

    def _reference_pressure(self, algorithm):
        """``ref_p`` と、 anisotropic 用の 6 成分を返す。

        UDF は目標圧力と応力を別に持つ。 GROMACS の ``ref-p`` は
        anisotropic では 6 成分 (xx, yy, zz, xy, zx, yz) で、 対角は
        ``P - sigma_ii``、 非対角は ``-sigma_ij``。
        """
        base = "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress."
        pressure = self._udf.get(base + "Pressure", "[bar]")
        if algorithm != "NPT_Parrinello_Rahman_Nose_Hoover":
            return pressure, None

        stress = [self._udf.get(base + "Stress." + e, "[bar]")
                  for e in ("xx", "yy", "zz", "yz", "zx", "xy")]
        # UDF は yz, zx, xy の順、 GROMACS は xy, zx, yz の順
        order = [0, 1, 2, 5, 4, 3]
        tensor = [pressure - stress[i] if i < 3 else -stress[i] for i in order]
        return pressure, tensor

    @staticmethod
    def _compressibility_tensor(algorithm, fix_cell, fix_angle,
                                deform_npt, deform_vel, commp):
        """anisotropic 用の圧縮率 6 成分。 固定した辺は 0 にする。

        0 を入れた方向は箱が動かない。 ``Fix_Cell_Length`` / ``Fix_Angle``
        をここに写す。 isotropic では ``None`` (スカラーだけ使う)。
        """
        if algorithm != "NPT_Parrinello_Rahman_Nose_Hoover":
            return None

        if fix_angle == 0:
            offdiag = [commp, commp, commp]
        elif fix_angle == 1:
            offdiag = [0.0, 0.0, 0.0]
        else:
            raise RuntimeError("Error: Please Check Your Fix_Cell_Angle")

        #: Fix_Cell_Length -> 動かせる方向 (1 = 動く)
        free_by_fix = {
            "":   (1, 1, 1),
            "xy": (0, 0, 1),
            "yz": (1, 0, 0),
            "zx": (0, 1, 0),
            "x":  (0, 1, 1),
            "y":  (1, 0, 1),
            "z":  (1, 1, 0),
        }
        if fix_cell in free_by_fix:
            diag = [commp if f else 0.0 for f in free_by_fix[fix_cell]]
            return diag + offdiag

        if deform_npt and deform_vel is not None:
            if fix_angle != 1:
                raise RuntimeError("Error: Please Check Your Fix_Cell_Angle")
            # 変形させる方向は圧力浴に触らせない
            diag = [0.0 if dv != 0.0 else commp for dv in deform_vel[:3]]
            return diag + [0.0, 0.0, 0.0]

        raise RuntimeError("Error: Please Check Your Fix_Cell_Length")
