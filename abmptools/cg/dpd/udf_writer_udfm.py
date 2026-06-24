# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.udf_writer_udfm
--------------------------------
**UDFManager ベース** の Cognac DPD 入力 UDF writer。

旧 ``udf_writer.py`` は positional plain-text で UDF を手書きしていたが、
Pair_Interaction (14 union サブフィールドを持つ class) を positional で書く方式は
脆く、生成 UDF が **cognac112 スキーマに非準拠で UDFManager ロード不可** だった
(Pair_Interaction を `Molecular_Attributes` 配下に置く / `Interaction_Site_Type[]`
未定義、等)。

本モジュールは dpdgen `udfdpd_io` と同じく **UDFManager の named-path put** で
組み立てるため、cognac112 スキーマ準拠が保証される (= 実際に cognac でロード可能)。

正しい配置 (cognac112.udf スキーマ準拠):
    Molecular_Attributes.Atom_Type[].Name / Mass
    Molecular_Attributes.Interaction_Site_Type[].Name / Num_of_Atoms / Range
    Molecular_Attributes.Bond_Potential[].Name / Potential_Type / R0 / Harmonic.K
    Molecular_Attributes.Angle_Potential[].Name / Potential_Type / theta0 / Theta.K
    Interactions.Pair_Interaction[].Name / Potential_Type / Site1_Name / Site2_Name
                                   / Cutoff / Scale_1_4_Pair / DPD.a / DPD.gamma
    Set_of_Molecules.molecule[].Mol_Name / atom[] / bond[] / angle[] / interaction_Site[]

``UDFManager`` は OCTA 同梱 (PyPI 非配布) で **lazy import**。`put` は numpy 値を
サイレント 0 化するため全数値 `float()`/`int()` cast 済 (`reference_udfmanager_
put_numpy_silent_zero`)。
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Union

from .models import DpdSystemSpec

logger = logging.getLogger(__name__)

DEFAULT_INCLUDE = "cognac112.udf"
_MA = "Molecular_Attributes"
_PI = "Interactions.Pair_Interaction"
_SM = "Set_of_Molecules.molecule"
# cognac90 互換 calc flag [Bond, Angle, Torsion, Non_Bonding, ?, Non_Bonding_1_3, ...]
_CALC_FLAGS = [1, 1, 0, 1, 1, 0, 1, 0, 0, 0]
_DPD_GAMMA = 4.5


def _skeleton(include_file: str, project_name: str, comment: str) -> str:
    """UDFManager が開ける最小 UDF (──include + header + 空 data)。"""
    return (
        f'\\include{{"{include_file}"}}\n'
        "\\begin{header}\n\\begin{data}\n"
        'EngineType:"COGNAC"\n'
        'EngineVersion:"Ver112"\n'
        'IOType:"IN"\n'
        f'ProjectName:"{project_name}"\n'
        f'Comment:"{comment}"\n'
        'Action:""\n'
        "\\end{data}\n\\end{header}\n"
        "\\begin{data}\n\\end{data}\n"
    )


def write_dpd_udf_udfm(
    spec: DpdSystemSpec,
    output_path: Union[str, Path],
    *,
    include_file: str = DEFAULT_INCLUDE,
) -> Path:
    """DpdSystemSpec から cognac DPD 入力 UDF を UDFManager で書き出す。

    Returns
    -------
    Path
        出力 UDF path。
    """
    from UDFManager import UDFManager  # OCTA 同梱 (PyPI 非配布)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        _skeleton(include_file, spec.project_name, spec.comment), encoding="utf-8"
    )

    u = UDFManager(str(output_path))

    _put_simulation_conditions(u, spec)
    _put_structure(u, spec)
    _put_molecular_attributes(u, spec)
    _put_interactions(u, spec)
    _put_set_of_molecules(u, spec)

    u.write(str(output_path))
    logger.info(
        "write_dpd_udf_udfm: %s (%d segment(s), %d monomer(s), %d pair(s))",
        output_path, len(spec.segment_names()), len(spec.monomers),
        len(spec.aij.to_a_values()),
    )
    return output_path


def _box(spec: DpdSystemSpec) -> float:
    calc = spec.calc_sett
    if calc and calc.box_size:
        return float(calc.box_size[0])
    return 10.0


def _put_simulation_conditions(u, spec: DpdSystemSpec) -> None:
    calc = spec.calc_sett
    dt = float(calc.dt) if calc else 0.05
    steps = int(calc.step_list[0]) if calc and calc.step_list else 10000
    out_interval = int(calc.output_interval) if calc else 100
    lam = float(calc.phys_param.get("lambda", 0.65)) if calc and calc.phys_param else 0.65

    u.put("DPD", "Simulation_Conditions.Solver.Dynamics.Dynamics_Algorithm")
    u.put(0.65 if lam is None else lam, "Simulation_Conditions.Solver.Dynamics.DPD.lambda")
    u.put(1.0e12, "Simulation_Conditions.Dynamics_Conditions.Max_Force")
    u.put(dt, "Simulation_Conditions.Dynamics_Conditions.Time.delta_T")
    u.put(steps, "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps")
    u.put(out_interval, "Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps")
    u.put(1.0, "Simulation_Conditions.Dynamics_Conditions.Temperature.Temperature")
    u.put(1.0, "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure")

    # calc flags (cognac112: Non_Bonding は Interchain / Intrachain に分割)
    cf = "Simulation_Conditions.Calc_Potential_Flags"
    u.put(1, f"{cf}.Bond")
    u.put(1, f"{cf}.Angle")
    u.put(0, f"{cf}.Torsion")
    u.put(1, f"{cf}.Non_Bonding_Interchain")
    u.put(1, f"{cf}.Non_Bonding_Intrachain")
    u.put(0, f"{cf}.Non_Bonding_1_3")


def _put_structure(u, spec: DpdSystemSpec) -> None:
    box = _box(spec)
    # initial structure: Random 生成 + density
    u.put("Random", "Initial_Structure.Generate_Method.Method")
    u.put(3.0, "Initial_Structure.Initial_Unit_Cell.Density")
    u.put(box, "Initial_Structure.Initial_Unit_Cell.Cell_Size.a")
    u.put(box, "Initial_Structure.Initial_Unit_Cell.Cell_Size.b")
    u.put(box, "Initial_Structure.Initial_Unit_Cell.Cell_Size.c")
    u.put(0, "Initial_Structure.Relaxation.Relaxation")
    # unit cell
    u.put(3.0, "Structure.Unit_Cell.Density")
    u.put(box, "Structure.Unit_Cell.Cell_Size.a")
    u.put(box, "Structure.Unit_Cell.Cell_Size.b")
    u.put(box, "Structure.Unit_Cell.Cell_Size.c")
    # unit parameter (DPD reduced units)
    u.put(1.0, "Unit_Parameter.Energy")
    u.put(0.71137866, "Unit_Parameter.Length")


def _bond_table(spec: DpdSystemSpec) -> Dict[Tuple[str, str], Tuple[float, float]]:
    """(name_i, name_j) -> (R0, K) を全 monomer から集約 (重複は先勝ち)。"""
    seen: Dict[Tuple[str, str], Tuple[float, float]] = {}
    for mono in spec.monomers:
        for entry in mono.bond12h:
            if len(entry) < 5:
                continue
            _, i, j, dist, stiff = entry[:5]
            ni = mono.particle_names[i] if i < mono.n_particles else f"P{i}"
            nj = mono.particle_names[j] if j < mono.n_particles else f"P{j}"
            key = tuple(sorted([ni, nj]))
            seen.setdefault(key, (float(dist), float(stiff)))
    return seen


def _angle_table(spec: DpdSystemSpec) -> Dict[Tuple[str, str, str], Tuple[float, float]]:
    seen: Dict[Tuple[str, str, str], Tuple[float, float]] = {}
    for mono in spec.monomers:
        for entry in mono.angle13data:
            if len(entry) < 5:
                continue
            a, b, c, eq, stiff = entry[:5]
            na = mono.particle_names[a] if a < mono.n_particles else f"P{a}"
            nb = mono.particle_names[b] if b < mono.n_particles else f"P{b}"
            nc = mono.particle_names[c] if c < mono.n_particles else f"P{c}"
            seen.setdefault((na, nb, nc), (float(eq), float(stiff)))
    return seen


def _put_molecular_attributes(u, spec: DpdSystemSpec) -> None:
    segs = spec.segment_names()
    for n, s in enumerate(segs):
        u.put(str(s), f"{_MA}.Atom_Type[{n}].Name")
        u.put(1.0, f"{_MA}.Atom_Type[{n}].Mass")
        u.put(str(s), f"{_MA}.Interaction_Site_Type[{n}].Name")
        u.put(1, f"{_MA}.Interaction_Site_Type[{n}].Num_of_Atoms")
        u.put(0.6, f"{_MA}.Interaction_Site_Type[{n}].Range")

    for n, ((ni, nj), (R0, K)) in enumerate(_bond_table(spec).items()):
        u.put(f"{ni}-{nj}", f"{_MA}.Bond_Potential[{n}].Name")
        u.put("Harmonic", f"{_MA}.Bond_Potential[{n}].Potential_Type")
        u.put(float(R0), f"{_MA}.Bond_Potential[{n}].R0")
        u.put(float(K), f"{_MA}.Bond_Potential[{n}].Harmonic.K")

    for n, ((na, nb, nc), (eq, K)) in enumerate(_angle_table(spec).items()):
        u.put(f"{na}-{nb}-{nc}", f"{_MA}.Angle_Potential[{n}].Name")
        u.put("Theta", f"{_MA}.Angle_Potential[{n}].Potential_Type")
        u.put(float(eq), f"{_MA}.Angle_Potential[{n}].theta0")
        u.put(float(K), f"{_MA}.Angle_Potential[{n}].Theta.K")


def _put_interactions(u, spec: DpdSystemSpec) -> None:
    for n, (i, j, a) in enumerate(spec.aij.to_a_values()):
        u.put(f"{i}-{j}", f"{_PI}[{n}].Name")
        u.put("DPD", f"{_PI}[{n}].Potential_Type")
        u.put(str(i), f"{_PI}[{n}].Site1_Name")
        u.put(str(j), f"{_PI}[{n}].Site2_Name")
        u.put(1.0, f"{_PI}[{n}].Cutoff")
        u.put(1.0, f"{_PI}[{n}].Scale_1_4_Pair")
        u.put(float(a), f"{_PI}[{n}].DPD.a")
        u.put(float(_DPD_GAMMA), f"{_PI}[{n}].DPD.gamma")


def _put_set_of_molecules(u, spec: DpdSystemSpec) -> None:
    """各 monomer 種を 1 molecule テンプレートとして登録 (atoms/bonds/angles)。"""
    bond_names = _bond_table(spec)
    angle_names = _angle_table(spec)
    for mid, mono in enumerate(spec.monomers):
        u.put(str(mono.name), f"{_SM}[{mid}].Mol_Name")
        names = list(mono.particle_names) or [f"P{k}" for k in range(mono.n_particles)]
        for k, pn in enumerate(names):
            u.put(int(k), f"{_SM}[{mid}].atom[{k}].Atom_ID")
            u.put(str(pn), f"{_SM}[{mid}].atom[{k}].Atom_Name")
            u.put(str(pn), f"{_SM}[{mid}].atom[{k}].Atom_Type_Name")
            u.put(0, f"{_SM}[{mid}].atom[{k}].Chirality")
            u.put(1, f"{_SM}[{mid}].atom[{k}].Main_Chain")
            u.put(str(pn), f"{_SM}[{mid}].interaction_Site[{k}].Type_Name")
            u.put(int(k), f"{_SM}[{mid}].interaction_Site[{k}].atom[0]")
        for bk, (i, j) in enumerate(mono.bond12):
            ni = names[i] if i < len(names) else f"P{i}"
            nj = names[j] if j < len(names) else f"P{j}"
            key = tuple(sorted([ni, nj]))
            pname = f"{key[0]}-{key[1]}" if key in bond_names else f"{ni}-{nj}"
            u.put(str(pname), f"{_SM}[{mid}].bond[{bk}].Potential_Name")
            u.put(int(i), f"{_SM}[{mid}].bond[{bk}].atom1")
            u.put(int(j), f"{_SM}[{mid}].bond[{bk}].atom2")
            u.put(1, f"{_SM}[{mid}].bond[{bk}].Order")
        for ak, (a, b, c) in enumerate(mono.angle13):
            na = names[a] if a < len(names) else f"P{a}"
            nb = names[b] if b < len(names) else f"P{b}"
            nc = names[c] if c < len(names) else f"P{c}"
            pname = f"{na}-{nb}-{nc}"
            u.put(str(pname), f"{_SM}[{mid}].angle[{ak}].Potential_Name")
            u.put(int(a), f"{_SM}[{mid}].angle[{ak}].atom1")
            u.put(int(b), f"{_SM}[{mid}].angle[{ak}].atom2")
            u.put(int(c), f"{_SM}[{mid}].angle[{ak}].atom3")
