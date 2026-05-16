# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.udf_writer
---------------------------
**Route R1**: Cognac DPD 入力 UDF (``*_uin.udf``) を **plain text で直接書き出す**。

設計指針
~~~~~~~~
- ``\\include{"cognac112.udf"}`` 1 行で class 定義を J-OCTA install 環境に委譲
  (abmptools 側に OCTA spec を持たない、 権利配慮)
- UDFManager (OCTA) **に依存しない** plain text writer
- 粒子位置 / 速度は **空 (Position_Generation_From_Cognac で COGNAC が随時発生)** で
  skeleton 出力、 user は COGNAC 内部で random 配置するか、 別途 manual で埋める
- Boilerplate な数値フィールドは dpm-sample_uin.udf を参考にした default 値で埋める

dpdgen 既存 logic (`udfdpd_io.py:setupudfdpd` 経由 ``UDFManager.put*``) と機能
等価だが、 plain text 書き出しなので **UDFManager 非依存**、 abmptoolsenv
(non-OCTA 環境) でも動く。
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, List, Tuple, Union

from .models import DpdSystemSpec, MonomerSpec

logger = logging.getLogger(__name__)

# `\include` 1 行で class 定義は J-OCTA install dir 経由で resolve
DEFAULT_INCLUDE = 'cognac112.udf'


def write_dpd_udf(
    spec: DpdSystemSpec,
    output_path: Union[str, Path],
    *,
    include_file: str = DEFAULT_INCLUDE,
) -> Path:
    """Cognac DPD 入力 UDF skeleton を ``output_path`` に書き出す。

    Parameters
    ----------
    spec : DpdSystemSpec
        cg_segmenter monomer + fcews aij + calc_sett から構築済の system spec。
        ``spec.calc_sett`` が None なら default 値で fallback。
    output_path : str | Path
        出力 UDF path (例: ``chol_uin.udf``)。
    include_file : str
        冒頭 ``\\include`` で参照する Cognac class 定義 file (default
        ``"cognac112.udf"``、 J-OCTA 11.x 同梱)。

    Returns
    -------
    Path
        出力 path。
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    blocks: List[str] = [
        f'\\include{{"{include_file}"}}',
        _header(spec),
        '\\begin{data}',
        _simulation_conditions(spec),
        _initial_structure(spec),
        _molecular_attributes(spec),
        _set_of_molecules(spec),
        '\\end{data}',
        '',
    ]
    output_path.write_text("\n".join(blocks) + "\n", encoding="utf-8")
    logger.info(
        "write_dpd_udf: %s (%d segments, %d monomers, %d pair, %d bond, %d angle)",
        output_path,
        len(spec.segment_names()),
        len(spec.monomers),
        len(spec.aij.pairs),
        sum(len(m.bond12) for m in spec.monomers),
        sum(len(m.angle13) for m in spec.monomers),
    )
    return output_path


def _header(spec: DpdSystemSpec) -> str:
    """\\begin{header}...\\end{header} (EngineType: COGNAC, Comment)。"""
    return (
        '\\begin{header}\n'
        '\\begin{data}\n'
        'EngineType:"COGNAC"\n'
        'EngineVersion:"Ver112"\n'
        'IOType:"IN"\n'
        f'ProjectName:"{spec.project_name}"\n'
        f'Comment:"{spec.comment}"\n'
        'Action:""\n'
        '\\end{data}\n'
        '\\end{header}'
    )


def _simulation_conditions(spec: DpdSystemSpec) -> str:
    """Simulation_Conditions:{ ... } (Dynamics=DPD、 default values)。"""
    calc = spec.calc_sett
    total_num = calc.total_num_list[0] if calc and calc.total_num_list else 100000
    total_steps = calc.step_list[0] if calc and calc.step_list else 10000
    dt = calc.dt if calc else 0.05
    out_interval = calc.output_interval if calc else 100
    lambda_val = (
        calc.phys_param.get("lambda", 0.65) if calc and calc.phys_param else 0.65
    )

    lines = [
        "Simulation_Conditions:{",
        " {",
        f"  {total_num}.000000000000,",
        f"  {{{dt:.14e},{total_steps},{out_interval}}}",
        "  {1.00000000000000,0}",
        "  {1.00000000000000,{0.0,0.0,0.0,0.0,0.0,0.0}}",
        '  {"",{"",{0.0}{0.0,0.0}}{"",{0.0,0.0,0.0,0.0,0.0,0.0}{"",{0.0}{0.0}0.50000000000000,""}{0.0,0.0,0.0}100,0}}',
        f"  {{{out_interval},0,0,0,0}}",
        "  {",
        "   0,",
        "   []",
        "   0,",
        "   []",
        "   1.00000000000000e-05",
        "  }",
        " }",
        " {",
        '  "Dynamics",',
        # DPD: gamma=20.0, cutoff=1.0, etc. (defaults from dpm-sample_uin.udf)
        f'  {{"DPD",{{20.0000000000000}}{{10.0000000000000}}{{0.50000000000000}}{{25.0000000000000}}{{25.0000000000000,0,""}}{{80.0000000000000,0,""}}{{200.000000000000,20.0000000000000}}{{20.0000000000000,0.50000000000000}}{{20.0000000000000,0,"",25.0000000000000}}{{20.0000000000000,0,"",0.50000000000000}}{{20.0000000000000,1.00000000000000}}{{20.0000000000000,0,"",1.00000000000000}}{{0.0}}{{0.0,0.0}}{{{lambda_val:.14f}}}}}',
        "  {",
        '   "",',
        "   0,",
        "   0.0,",
        "   0,",
        "   {0,0.0}",
        "   []",
        "  }",
        " }",
        ' {"PERIODIC","PERIODIC","PERIODIC",0}',
        " {1,1,0,1,1,1,0,0,0,0}",
        " {{1,1,1,1,1,1,1,1,1,0,0,0}{1,1,1}{0,0,0}{0,0,0}}",
        " {0,0,0,{0.0,0.0,0.0}0.0,0}",
        " {",
        '  "NO",',
        "  []",
        " }",
        "}",
    ]
    return "\n".join(lines)


def _initial_structure(spec: DpdSystemSpec) -> str:
    """Initial_Structure (cell + Restart で粒子位置は空 = COGNAC が後で生成)。"""
    calc = spec.calc_sett
    if calc and calc.box_size:
        bx, by, bz = calc.box_size
    else:
        bx, by, bz = 10.0, 10.0, 10.0

    lines = [
        "Initial_Structure:{",
        f" {{{bx:.14f},{{0.0,0.0,0.0,90.0000000000000,90.0000000000000,90.0000000000000}}0.0}}",
        ' {"",-1}',
        " {",
        '  "Restart",',
        '  {"",-1,1,1}',
        "  {",
        "   {",
        '    "",',
        "    []",
        "   }",
        "   {0,0.0,0.0}",
        "   []",
        "   []",
        "   0.0,",
        "   []",
        "  }",
        "  {",
        '   {"",{0.0,0.0}}',
        "   0,",
        "   0,",
        "   0,",
        "   []",
        "   0.0,",
        "   0.0,",
        "   []",
        "   0.0,",
        '   ""',
        "  }",
        "  {",
        "   0.0,",
        "   0,",
        "   0.0,",
        "   {0.0,0,0.0,0.0,0.0}",
        "   {",
        "    0.0,",
        "    []",
        "    0.0",
        "   }",
        "  }",
        '  {"",0,0,0,0.0}',
        " }",
        ' {1,"DYNAMICS",300.000000000000,5000,0}',
        "}",
    ]
    _ = (by, bz)  # 現状の Cognac UDF は cell.a に総合値、 b/c は ratio
    return "\n".join(lines)


def _molecular_attributes(spec: DpdSystemSpec) -> str:
    """Molecular_Attributes:{Atom_Type[] Bond_Potential[] Angle_Potential[] Pair_Interaction[]}.

    cg_segmenter monomer の bond12 / angle13 を平坦化し、 fcews aij の各 pair を
    Pair_Interaction (DPD type) として出力する。
    """
    atom_types = _atom_types_block(spec)
    bond_pots = _bond_potentials_block(spec)
    angle_pots = _angle_potentials_block(spec)
    pair_inters = _pair_interactions_block(spec)
    torsion_pots = " []"  # angle ポテンシャル運用なので torsion なし

    return (
        "Molecular_Attributes:{\n"
        + atom_types + "\n"
        + bond_pots + "\n"
        + angle_pots + "\n"
        + torsion_pots + "\n"
        + pair_inters + "\n"
        + "}"
    )


def _atom_types_block(spec: DpdSystemSpec) -> str:
    """[{"segA",1.0}{"segB",1.0}...] 形式 (各 segment の mass=1.0)。"""
    segs = spec.segment_names()
    body = "".join(f'{{"{s}",1.00000000000000}}' for s in segs)
    return f" [{body}]"


def _bond_potentials_block(spec: DpdSystemSpec) -> str:
    """Bond_Potential[]: 各 bond12 を ``"seg_i-seg_j"`` 名で登録。

    cg_segmenter の各 bond12h (= [type, i, j, dist, stiff, 0]) を Harmonic に変換。
    重複する seg_i-seg_j ペアは集約。
    """
    seen: dict = {}  # key=(seg_i, seg_j) value=(R0, K)
    for mono in spec.monomers:
        for entry in mono.bond12h:
            # entry = [type, i, j, dist, stiff, 0]
            if len(entry) < 5:
                continue
            _, i, j, dist, stiff = entry[:5]
            i_name = mono.particle_names[i] if i < mono.n_particles else f"P{i}"
            j_name = mono.particle_names[j] if j < mono.n_particles else f"P{j}"
            key = tuple(sorted([i_name, j_name]))
            seen.setdefault(key, (float(dist), float(stiff)))

    lines = [" ["]
    for (i_name, j_name), (R0, K) in seen.items():
        name = f"{i_name}-{j_name}"
        lines.append(f'  {{')
        lines.append(f'   "{name}",')
        lines.append(f'   "Harmonic",')
        lines.append(f'   {R0:.14f},')
        lines.append(f'   {{{K:.14f}}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0}}')
        lines.append(f'   {{0.0}}')
        lines.append(f'   {{0.0,0.0,""}}')
        lines.append(f'   {{')
        lines.append(f'    0,')
        lines.append(f'    []')
        lines.append(f'    ""')
        lines.append(f'   }}')
        lines.append(f'   {{0.0}}')
        lines.append(f'   {{"","",{{""}}{{0}}}}')
        lines.append(f'   {{')
        lines.append(f'    0,')
        lines.append(f'    []')
        lines.append(f'   }}')
        lines.append(f'  }}')
    lines.append(" ]")
    return "\n".join(lines)


def _angle_potentials_block(spec: DpdSystemSpec) -> str:
    """Angle_Potential[]: cg_segmenter angle13 を Theta_Polynomial で登録。

    cognac 余角 convention (eq=0 が 180°、 eq=30 が 150° 等) を ``theta0`` フィールド
    に直接書き出す。 stiffness は ``{K}`` (Theta フォーマット)。
    """
    seen: dict = {}  # key=(a, b, c) value=(theta0_余角, K)
    for mono in spec.monomers:
        for entry in mono.angle13data:
            # entry = [a, b, c, eq余角, stiff]
            if len(entry) < 5:
                continue
            a, b, c, eq, stiff = entry[:5]
            a_name = mono.particle_names[a] if a < mono.n_particles else f"P{a}"
            b_name = mono.particle_names[b] if b < mono.n_particles else f"P{b}"
            c_name = mono.particle_names[c] if c < mono.n_particles else f"P{c}"
            key = (a_name, b_name, c_name)
            seen.setdefault(key, (float(eq), float(stiff)))

    lines = [" ["]
    for (a_name, b_name, c_name), (theta0, K) in seen.items():
        name = f"{a_name}-{b_name}-{c_name}"
        lines.append(f'  {{')
        lines.append(f'   "{name}",')
        lines.append(f'   "Theta",')
        lines.append(f'   {theta0:.14f},')
        lines.append(f'   {{{K:.14f}}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0,0.0,0.0}}')
        lines.append(f'   {{0.0,0,0,0}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0}}')
        lines.append(f'  }}')
    lines.append(" ]")
    return "\n".join(lines)


def _pair_interactions_block(spec: DpdSystemSpec) -> str:
    """Pair_Interaction[]: fcews aij の各 pair を ``"DPD"`` type で登録。

    a 値は ``aij.to_a_values()`` で chi モードからの自動変換含む。
    cutoff (range) = 1.0 (DPD 標準) を仮定。
    """
    a_values = spec.aij.to_a_values()
    lines = [" ["]
    for (i, j, a) in a_values:
        name = f"{i}-{j}"
        lines.append(f'  {{')
        lines.append(f'   "{name}",')
        lines.append(f'   "DPD",')
        lines.append(f'   "{i}",')
        lines.append(f'   "{j}",')
        lines.append(f'   1.00000000000000,')
        lines.append(f'   1.00000000000000,')
        lines.append(f'   {{0.0,0.0}}')
        lines.append(f'   {{0.0,0.0,0.0}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0,0,0}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0,0.0,0.0,0.0}}')
        # DPD-specific: {a, range}
        lines.append(f'   {{{a:.14f},4.50000000000000}}')
        lines.append(f'   {{0.0,0.0,0.0}}')
        lines.append(f'   {{0.0,0.0,0.0}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}')
        lines.append(f'   {{0.0,0.0,0.0,0.0}}')
        lines.append(f'   {{0.0,0.0,0.0}}')
        lines.append(f'   {{"","",{{""}}{{0}}}}')
        lines.append(f'   {{')
        lines.append(f'    0,')
        lines.append(f'    []')
        lines.append(f'   }}')
        lines.append(f'  }}')
    lines.append(" ]")
    return "\n".join(lines)


def _set_of_molecules_block_inner(mono: MonomerSpec, mol_idx: int) -> str:
    """1 molecule の Set_of_Molecules entry (atom_list + bond_list + angle_list)."""
    n = mono.n_particles
    lines = [" {"]
    lines.append(f'  "{mono.name}",')
    lines.append(f'  {n},')
    lines.append("  [")
    # atom_list: [{atom_id, Atom_Type, charge, position}]
    for k in range(n):
        a_type = mono.particle_names[k]
        lines.append(f'   {{{k},"{a_type}",0.0,{{0.0,0.0,0.0}}}}')
    lines.append("  ]")
    # bond_list: [{atom1, atom2, Bond_Potential}]
    lines.append("  [")
    for (i, j) in mono.bond12:
        i_name = mono.particle_names[i]
        j_name = mono.particle_names[j]
        bp = f"{min(i_name, j_name)}-{max(i_name, j_name)}"
        lines.append(f'   {{{i},{j},"{bp}"}}')
    lines.append("  ]")
    # angle_list: [{atom1, atom2, atom3, Angle_Potential}]
    lines.append("  [")
    for (a, b, c) in mono.angle13:
        a_name = mono.particle_names[a]
        b_name = mono.particle_names[b]
        c_name = mono.particle_names[c]
        ap = f"{a_name}-{b_name}-{c_name}"
        lines.append(f'   {{{a},{b},{c},"{ap}"}}')
    lines.append("  ]")
    lines.append("  []")  # torsion_list (空)
    lines.append("  []")  # exclude_pair (空)
    lines.append(" }")
    return "\n".join(lines)


def _set_of_molecules(spec: DpdSystemSpec) -> str:
    """Set_of_Molecules:{ [{...mol1...},{...mol2...}] }"""
    lines = ["Set_of_Molecules:{", " ["]
    for k, mono in enumerate(spec.monomers):
        lines.append(_set_of_molecules_block_inner(mono, k))
    lines.append(" ]")
    lines.append("}")
    return "\n".join(lines)
