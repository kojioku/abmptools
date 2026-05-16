# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.dpm_writer
---------------------------
**Route R2**: OCTA / J-OCTA DPD Modeler 用 ``*.dpm`` (GUI 用 modeling file) を
**B 案 (user template + patch)** で生成する。

権利配慮
~~~~~~~~
``*.dpm`` の class 定義 (``\\begin{def}`` 800+ 行、 J-OCTA 商用 spec) は
**abmptools 側で複製しない**。 代わりに以下のフロー:

1. ユーザーが J-OCTA で **空 dpm template** を 1 回だけ作成 (ModelType:4 / SegmentModel
   が空の状態)。
2. abmptools が :func:`patch_dpm` で **``\\begin{data}`` セクション内の特定ブロック**
   (SegmentModel / SegmentPairModel / PolymerModel / DpdInput / FcewsParam) を
   生成内容で差し替える。
3. ``Virtual.mom`` も user 提供 (J-OCTA Monomer Modeler で 1 回作成) を
   :func:`propagate_virtual_mom` で全 segment dir に copy する。

これにより abmptools のリポジトリには J-OCTA 商用 spec の内容が含まれない。

参考ファイル
~~~~~~~~~~~~
``man/octa/dpdfile-test/dpm-sample.dpm`` は **構造把握用 reference のみ**、
本パッケージに複製しない。
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Union

from .models import DpdSystemSpec, MonomerSpec

logger = logging.getLogger(__name__)

# `#Message.txt` の標準内容 (J-OCTA Monomer Modeler 流儀: 1 行 "RUN")
MESSAGE_TXT_CONTENT = "RUN\n"

# patch_dpm が差し替える default フィールド
DEFAULT_PATCH_FIELDS = (
    "SegmentModel[]",
    "SegmentPairModel[]",
    "PolymerModel[]",
    "DpdInput",
    "FcewsParam",
)


# --- B 案: Virtual.mom propagator ---------------------------------------

def propagate_virtual_mom(
    template_path: Union[str, Path],
    segments: List[str],
    output_dir: Union[str, Path],
) -> List[Path]:
    """User-provided ``Virtual.mom`` を全 segment dir に copy する。

    Parameters
    ----------
    template_path : str | Path
        ユーザーが J-OCTA Monomer Modeler で 1 回だけ作成した ``Virtual.mom``。
    segments : List[str]
        対象 segment 名のリスト (cg_segmenter で生成された segment 名と一致)。
    output_dir : str | Path
        ``monomer-lib`` 配下に相当するディレクトリ。 各 segment 名のサブディレクトリ
        を作成し、 そこに ``Virtual.mom`` を配置する。

    Returns
    -------
    List[Path]
        生成された ``Virtual.mom`` の path リスト (segment 順)。

    Notes
    -----
    全 segment dir の ``Virtual.mom`` は **完全に同一の内容** で配置される
    (J-OCTA 流儀、 ``dpm-sample/monomer-lib/*/Virtual.mom`` が ``diff`` で
    差分なしであることを確認済)。
    """
    template_path = Path(template_path)
    if not template_path.exists():
        raise FileNotFoundError(
            f"Virtual.mom template not found: {template_path}\n"
            f"User must provide a Virtual.mom by 1 round-trip in J-OCTA Monomer Modeler."
        )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    content = template_path.read_text(encoding="utf-8")
    paths: List[Path] = []
    for seg in segments:
        seg_dir = output_dir / seg
        seg_dir.mkdir(parents=True, exist_ok=True)
        out = seg_dir / "Virtual.mom"
        out.write_text(content, encoding="utf-8")
        paths.append(out)
    logger.info(
        "propagate_virtual_mom: copied %s to %d segment dir(s) under %s",
        template_path.name, len(paths), output_dir,
    )
    return paths


def write_message_txt(output_dir: Union[str, Path], content: str = MESSAGE_TXT_CONTENT) -> Path:
    """``#Message.txt`` (J-OCTA 流儀の sentinel file) を書き出す。

    Parameters
    ----------
    output_dir : str | Path
        dpm 連携ディレクトリ (例: ``./chol_project/``)。
    content : str
        ファイル内容 (default = ``"RUN\\n"``、 J-OCTA dpm-sample に倣う)。
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / "#Message.txt"
    out.write_text(content, encoding="utf-8")
    logger.info("write_message_txt: %s (%d bytes)", out, len(content))
    return out


# --- B 案: dpm template の data section patcher --------------------------

def patch_dpm(
    template_path: Union[str, Path],
    output_path: Union[str, Path],
    spec: DpdSystemSpec,
    *,
    replace_fields: Optional[List[str]] = None,
) -> Path:
    """User-provided 空 dpm template の data section を生成内容で差し替える。

    Parameters
    ----------
    template_path : str | Path
        J-OCTA で作成した空 dpm template (ModelType:4 / SegmentModel[] が空でも可)。
    output_path : str | Path
        出力 dpm path。
    spec : DpdSystemSpec
        生成元データ (monomers + aij + project_name など)。
    replace_fields : List[str], optional
        差し替える top-level field のリスト。 default は :data:`DEFAULT_PATCH_FIELDS`
        (``["SegmentModel[]", "SegmentPairModel[]", "PolymerModel[]", "DpdInput",
        "FcewsParam"]``)。

    Returns
    -------
    Path
        出力 path。

    Notes
    -----
    実装は ``\\begin{data} ... \\end{data}`` セクション内の **top-level field**
    (`{` を伴う / `[]:` を伴う) を brace 数えで識別し、 該当範囲を生成内容で
    置換する。 dpm の class 定義 (``\\begin{def}`` 800+ 行) や header は
    template のまま温存される。
    """
    template_path = Path(template_path)
    output_path = Path(output_path)
    if not template_path.exists():
        raise FileNotFoundError(
            f"dpm template not found: {template_path}\n"
            f"User must provide an empty dpm template from J-OCTA DPD Modeler."
        )
    if replace_fields is None:
        replace_fields = list(DEFAULT_PATCH_FIELDS)

    text = template_path.read_text(encoding="utf-8")

    blocks: Dict[str, str] = {
        "SegmentModel[]": _build_segment_models(spec),
        "SegmentPairModel[]": _build_segment_pair_models(spec),
        "PolymerModel[]": _build_polymer_models(spec),
        "DpdInput": _build_dpd_input(spec),
        "FcewsParam": _build_fcews_param(spec),
    }
    for field in replace_fields:
        if field not in blocks:
            raise ValueError(
                f"unknown field {field!r}; valid: {list(blocks.keys())}"
            )
        text = _replace_dpm_field(text, field, blocks[field])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")
    logger.info(
        "patch_dpm: %s -> %s (patched %d fields: %s)",
        template_path.name, output_path, len(replace_fields), replace_fields,
    )
    return output_path


# --- low-level: brace-aware field replacement -----------------------------

def _replace_dpm_field(text: str, field_name: str, new_value: str) -> str:
    """dpm の data section 内 top-level field を新値で置換。

    field の存在パターン:
    - ``Field[]:[ ... ]`` (リスト型、 末尾 ``]``)
    - ``Field:{ ... }`` (struct 型、 末尾 ``}``)

    ``\\begin{def}...\\end{def}`` 内の **型宣言** (例 ``DpdInput:DPDInput``) を
    誤って matching しないよう、 ``\\end{def}`` 以降の **main data section** から
    検索する。 field が見つからない場合は ``\\end{data}`` の直前に追加する。
    """
    lines = text.split("\n")

    # main data section の開始位置を確定。 dpm は入れ子の \begin{def}
    # (header 内の小さい def + main の class 定義) を持つので、 **最後の**
    # \end{def} を探し、 その後の \begin{data} を main data section とする。
    last_end_def = None
    for i, line in enumerate(lines):
        if line.strip().startswith(r"\end{def}"):
            last_end_def = i
    if last_end_def is None:
        search_start = 0
    else:
        data_start_idx = None
        for i in range(last_end_def + 1, len(lines)):
            if lines[i].strip().startswith(r"\begin{data}"):
                data_start_idx = i
                break
        if data_start_idx is None:
            raise ValueError(
                r"\begin{data} not found after the last \end{def} in template"
            )
        search_start = data_start_idx + 1

    # 1) field のマッチを探す (search_start 以降のみ)
    start_idx = None
    pattern = re.compile(rf"^{re.escape(field_name)}\s*:")
    for i in range(search_start, len(lines)):
        if pattern.match(lines[i].lstrip()):
            start_idx = i
            break

    if start_idx is None:
        # field が無い → \end{data} の直前に追加
        for i, line in enumerate(lines):
            if line.strip().startswith("\\end{data}"):
                lines.insert(i, new_value)
                logger.info("patch_dpm: field %s not found, inserted before \\end{data}", field_name)
                return "\n".join(lines)
        raise ValueError(
            f"\\end{{data}} not found in template; cannot insert {field_name}"
        )

    # 2) start_idx から開始ブレース ('{' or '[') の位置を探し、 brace counting で
    #    対応する閉じブレースまでの行範囲を確定
    start_line = lines[start_idx]
    open_idx = -1
    open_char = None
    for ch in ("{", "["):
        pos = start_line.find(ch)
        if pos >= 0 and (open_idx < 0 or pos < open_idx):
            open_idx = pos
            open_char = ch
    if open_idx < 0:
        raise ValueError(
            f"opening brace not found on {field_name}'s start line: {start_line!r}"
        )
    close_char = "}" if open_char == "{" else "]"

    # 3) brace counting で end line を探す
    depth = 0
    end_idx = None
    end_pos_in_line = None
    for i in range(start_idx, len(lines)):
        line = lines[i]
        scan_from = open_idx if i == start_idx else 0
        for j, ch in enumerate(line[scan_from:], start=scan_from):
            if ch == open_char:
                depth += 1
            elif ch == close_char:
                depth -= 1
                if depth == 0:
                    end_idx = i
                    end_pos_in_line = j
                    break
        if end_idx is not None:
            break

    if end_idx is None:
        raise ValueError(
            f"closing brace for {field_name} not found (template malformed?)"
        )

    # 4) start_idx の field 行頭 + end_idx の閉じブレース直後までを new_value で置換
    before = "\n".join(lines[:start_idx])
    after_in_end_line = lines[end_idx][end_pos_in_line + 1:]
    after_lines = lines[end_idx + 1:]
    after = (after_in_end_line + "\n" + "\n".join(after_lines)) if after_lines else after_in_end_line

    rebuilt = (before + ("\n" if before else "") + new_value
               + ("\n" if not new_value.endswith("\n") else "")
               + after)
    return rebuilt


# --- block builders (B 案 で patch する 5 ブロック) ----------------------

def _build_segment_models(spec: DpdSystemSpec) -> str:
    """SegmentModel[]:[ {...}, {...}, ... ] を生成 (各 segment 1 entry)。

    各 entry のフィールド構造は dpm spec に依存するため、 **最小限の値だけ書き、
    細部は user が J-OCTA で編集する前提**。 主要 fields:
      - segment 名
      - radius (4.5 が dpm-sample default)
      - mass param (25.0 が dpm-sample default、 fcews aii_val と同じ)
      - 色 / 表示 (適当な初期値、 user 編集前提)
    """
    seg_names = spec.segment_names()
    lines = ["SegmentModel[]:["]
    for name in seg_names:
        lines.append(f' {{')
        lines.append(f'  "{name}",')
        lines.append(f'  0,')
        lines.append(f'  0,')
        lines.append(f'  "false",')
        lines.append(f'  "",')
        lines.append(f'  1,')
        lines.append(f'  0.0,')
        lines.append(f'  0.0,')
        lines.append(f'  0.0,')
        lines.append(f'  1.0,')
        lines.append(f'  1.0,')
        lines.append(f'  "false",')
        lines.append(f'  "true",')
        lines.append(f'  0,')
        lines.append(f'  4.5,')
        lines.append(f'  {spec.aij.aii},')
        for _ in range(6):
            lines.append(f'  "false",')
        lines.append(f'  "",')
        lines.append(f'  "",')
        lines.append(f'  {{"Virtual"}}')
        lines.append(f'  {{80,40,{{178,0,0,255}}}}')
        lines.append(f'  []')
        lines.append(f'  {{')
        lines.append(f'   0,')
        for _ in range(5):
            lines.append(f'   []')
        lines.append(f'  }}')
        lines.append(f'  [{{"NO-H",-1.00000000000000e+99,-1.00000000000000e+99}}]')
        lines.append(f'  {{"","","",0.0,0.0,0.0,1.0}}')
        lines.append(f'  {{0.0,0.0,0.0,0.0,0.0,0.0}}')
        lines.append(f'  {{"false",0}}')
        lines.append(f' }}')
    lines.append("]")
    return "\n".join(lines)


def _build_segment_pair_models(spec: DpdSystemSpec) -> str:
    """SegmentPairModel[]:[ {...}, ... ] を生成 (各 segment 順序組合せ 1 entry)。

    fcews aij の各 pair (i != j のみ、 同種は SegmentModel に統合) を反映。
    """
    a_values = {
        (i, j): v for (i, j, v) in spec.aij.to_a_values()
    }
    seg_names = spec.segment_names()
    lines = ["SegmentPairModel[]:["]
    for i in seg_names:
        for j in seg_names:
            if i == j:
                continue
            # 対称: (i, j) or (j, i) から値を取る (なければ default aii)
            val = a_values.get((i, j), a_values.get((j, i), spec.aij.aii))
            lines.append(f' {{')
            lines.append(f'  "{i}",')
            lines.append(f'  "{j}",')
            for _ in range(8):
                lines.append(f'  "false",')
            lines.append(f'  0,')
            lines.append(f'  0.0,')
            lines.append(f'  0.0,')
            lines.append(f'  {{4.5,{val}}}')
            lines.append(f'  {{')
            lines.append(f'   0.0, 0.0, 0.0, 0.0,')
            lines.append(f'   []')
            lines.append(f'   []')
            lines.append(f'   []')
            for _ in range(7):
                lines.append(f'   "",')
            lines.append(f'  }}')
            for _ in range(5):
                lines.append(f'  0.0,')
            lines.append(f'  {{0,0,"","","",0.0,0.0,0.0}}')
            lines.append(f' }}')
    lines.append("]")
    return "\n".join(lines)


def _build_polymer_models(spec: DpdSystemSpec) -> str:
    """PolymerModel[]:[ {polymer1}, ... ] を生成 (各 monomer 1 entry)。

    cg_segmenter monomer の bond12 / particle_names を反映。
    """
    lines = ["PolymerModel[]:["]
    for mono in spec.monomers:
        lines.append(f' {{')
        lines.append(f'  "{mono.name}",')
        lines.append(f'  [')
        # particle_names を `{name, count, x, y, charge, ...}` 形式で
        for k, p_name in enumerate(mono.particle_names):
            lines.append(f'   {{"{p_name}",1,0,0,0.0,0,0,0,0,0}}')
        lines.append(f'  ]')
        # bond12 を polymer 内 bond として
        lines.append(f'  [')
        for (i, j) in mono.bond12:
            lines.append(f'   {{')
            lines.append(f'    [{{{i},{j},0,0}}]')
            lines.append(f'   }}')
        lines.append(f'  ]')
        lines.append(f'  []')
        lines.append(f'  []')
        lines.append(f'  {mono.n_particles}')
        lines.append(f' }}')
    lines.append("]")
    return "\n".join(lines)


def _build_dpd_input(spec: DpdSystemSpec) -> str:
    """DpdInput:{ ... } を生成 (DPDBond / DPDAngle / solvent / 物理パラメータ)。

    cg_segmenter の bond12 (path 1 のみ DPD bond として登録、 bond13/bond14 は
    angle 経由なので除外)、 angle13 (cognac 余角 convention) を流用。
    """
    calc = spec.calc_sett
    dt = calc.dt if calc else 0.05
    total_steps = calc.step_list[0] if calc and calc.step_list else 10000
    out_interval = calc.output_interval if calc else 100
    density = calc.density if calc else 3.0

    lines = ["DpdInput:{"]
    lines.append(" 0.0,")              # Atom_Mass
    lines.append(" 0.0,")              # Interaction_Site_Range
    lines.append(f" {dt},")            # delta_T
    lines.append(f" {total_steps},")   # Total_Steps
    lines.append(f" {out_interval},")  # Output_Interval_Steps
    lines.append(f" {density},")       # Density
    lines.append(" 4.0,")              # NumWater (default)
    lines.append(" 0.71132543032264,") # UnitLength (default from dpm-sample)
    lines.append(" 142.119094428142,") # UnitTime (default)
    lines.append(" 0.65,")             # DPD_lambda (default)
    lines.append(" 1,")                # angle_sw (on)
    lines.append(' "",')               # angle_polymer_filename
    # DpdBond[] : bond12 (path 1) のみを登録
    lines.append(" [")
    seen_bond_pairs = set()
    for mono in spec.monomers:
        for (i_idx, j_idx) in mono.bond12:
            i_name = mono.particle_names[i_idx]
            j_name = mono.particle_names[j_idx]
            key = tuple(sorted([i_name, j_name]))
            if key in seen_bond_pairs:
                continue
            seen_bond_pairs.add(key)
            # bond12h の dist を採用 (default 0.86)、 K = 100.0
            lines.append(f'  {{"{i_name}","{j_name}","Harmonic",100.0,0.86}}')
    lines.append(" ]")
    # DpdAngle[] : angle13 (cognac 余角 convention) を登録
    lines.append(" [")
    seen_angle_triples = set()
    for mono in spec.monomers:
        for k_idx, (a, b, c) in enumerate(mono.angle13):
            a_name = mono.particle_names[a]
            b_name = mono.particle_names[b]
            c_name = mono.particle_names[c]
            key = (a_name, b_name, c_name)
            if key in seen_angle_triples:
                continue
            seen_angle_triples.add(key)
            data_row = mono.angle13data[k_idx] if k_idx < len(mono.angle13data) else None
            theta0 = float(data_row[3]) if data_row else 0.0
            K = float(data_row[4]) if data_row else 5.0
            lines.append(f'  {{"{a_name}","{b_name}","{c_name}",{theta0},{K}}}')
    lines.append(" ]")
    # solvent[]: monomer のうち n_particles == 1 のものを solvent と推定
    lines.append(" [")
    for mono in spec.monomers:
        if mono.n_particles == 1:
            lines.append(f'  {{"{mono.name}",1000}}')
    lines.append(" ]")
    lines.append(" 0,")               # calc_chi_param
    lines.append(" 1,")               # calc_a_from_chi (a 値を chi から計算する)
    lines.append(" 3,")               # calc_a_type
    lines.append(' "",')              # calc_from_volume_flag
    lines.append(' "",')              # calc_from_volume
    lines.append(' "true",')          # non_bond_1_3
    lines.append(' "true",')          # position_generation_from_cognac
    lines.append(" {0,0,0,{0,0,0}0.0,0}")   # density_output (empty)
    lines.append("}")
    return "\n".join(lines)


def _build_fcews_param(spec: DpdSystemSpec) -> str:
    """FcewsParam:{ ... } を生成。 abmptools 経路では空 (FCEWS 未連携) で OK。"""
    lines = ["FcewsParam:{"]
    lines.append(' "false",')
    lines.append(' 0,')
    lines.append(' "",')
    lines.append(' "",')
    lines.append(' {"","","","","","",""}')
    lines.append(' {0,"","","false",0,0,0,"","",0,0,0,"","","","",""}')
    lines.append(' {"",0.0,"false",0,0,0,0,"","",0,0,"","false",0.0}')
    lines.append(' []')
    lines.append('}')
    return "\n".join(lines)
