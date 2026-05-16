# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.models
-----------------------
Data classes for DPD system specification (R1 = UDF / R2 = DPM 経路共通)。

cg_segmenter (R0) が `{name}_monomer` + `{name}_calc_sett` を吐き、 fcews が
`aij.dat` (Python 辞書スクリプト) を生成する。 これらを **構造化された Python
オブジェクト** として保持し、 `udf_writer` / `dpm_writer` の入力にする。

設計指針:
- Frozen dataclass、 Pydantic 風 (cg_peptide / cg_membrane と同じスタイル)
- JSON roundtrip 可能 (asdict / from_dict)
- `cg_segmenter.dpdgen_exporter` の出力と byte-equivalent な構造を再構築できる
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class AijMatrix:
    """fcews aij.dat (Python 辞書スクリプト) の構造化表現。

    aij.dat は `aij = [[seg_i, seg_j, value], ...]` または
    `chi = [[seg_i, seg_j, value], ...]` の Python list of [str, str, float]。

    Attributes
    ----------
    segments : List[str]
        参加 segment 名のリスト (登場順)。 N segments → N*(N+1)/2 ペア。
    pairs : List[Tuple[str, str, float]]
        各 segment ペア (i, j, value)。 a_ij または chi_ij。 i=j (自分相互作用) を含む。
    mode : str
        ``"a"`` (Groot-Warren a パラメータ直接) または ``"chi"`` (Flory-Huggins χ)。
    aii : float
        同種粒子間 a パラメータ (DPD default 25.0)。 chi → a 変換時に必要。
    """

    segments: List[str]
    pairs: List[Tuple[str, str, float]]
    mode: str = "a"
    aii: float = 25.0

    def __post_init__(self) -> None:
        if self.mode not in ("a", "chi"):
            raise ValueError(f"mode must be 'a' or 'chi', got {self.mode!r}")

    def to_a_values(self) -> List[Tuple[str, str, float]]:
        """chi → a 変換 (Groot-Warren: a_ij = a_ii + chi_ij / 0.306)。

        Returns
        -------
        List of (seg_i, seg_j, a_ij)。 mode='a' の場合は self.pairs をそのまま返す。
        """
        if self.mode == "a":
            return list(self.pairs)
        return [(i, j, self.aii + v / 0.306) for (i, j, v) in self.pairs]

    def get(self, seg_i: str, seg_j: str) -> Optional[float]:
        """対称 lookup: (i, j) または (j, i) の値を返す。 なければ None。"""
        for (a, b, v) in self.pairs:
            if (a == seg_i and b == seg_j) or (a == seg_j and b == seg_i):
                return v
        return None


@dataclass
class MonomerSpec:
    """cg_segmenter `{name}_monomer` ファイルの構造化表現。

    `{name}_monomer` は Python script で以下を定義する:
    - `bond12` / `bond12h` (path 1 = 隣接 bond)
    - `bond13_150` / `bond13_150h` (1-skip ring bond、 cg_segmenter で生成)
    - `bond14_150` / `bond14_150h` (2-skip ring bond、 cg_segmenter で生成)
    - `angle13` / `angle13data` (angle ポテンシャル、 cognac 余角)

    Attributes
    ----------
    name : str
        Monomer 名 (例: "chol", "wat")。 calc_sett の ratio_list / monomer-lib の
        ディレクトリ名にも使う。
    particle_names : List[str]
        各 CG 粒子のラベル (segment 名)。 長さ = n_particles。
    bond12, bond13_150, bond14_150 : List[Tuple[int, int]]
        粒子 index ペア (0-origin)。
    bond12h, bond13_150h, bond14_150h : List[List]
        各 6 要素 `[type, i, j, dist, stiff, 0]`。
    angle13 : List[Tuple[int, int, int]]
        粒子 index トリプル (中央 = b)。
    angle13data : List[List]
        各 5 要素 `[a, b, c, eq_余角, stiff]` (cognac 余角 convention)。
    """

    name: str
    particle_names: List[str]
    bond12: List[Tuple[int, int]] = field(default_factory=list)
    bond12h: List[List] = field(default_factory=list)
    bond13_150: List[Tuple[int, int]] = field(default_factory=list)
    bond13_150h: List[List] = field(default_factory=list)
    bond14_150: List[Tuple[int, int]] = field(default_factory=list)
    bond14_150h: List[List] = field(default_factory=list)
    angle13: List[Tuple[int, int, int]] = field(default_factory=list)
    angle13data: List[List] = field(default_factory=list)

    @property
    def n_particles(self) -> int:
        return len(self.particle_names)

    def all_bonds(self) -> List[List]:
        """全 bond ポテンシャル entry を結合 (bondNNh の連結)。"""
        return list(self.bond12h) + list(self.bond13_150h) + list(self.bond14_150h)


@dataclass
class CalcSett:
    """calc_sett ファイルの構造化表現。

    cg_segmenter dpdgen_exporter が生成する Python script を読み込んで構造化。

    Attributes
    ----------
    total_num_list : List[int]
        各 ratio 設定での total 粒子数。
    step_list : List[int]
        MD ステップ数。
    monomer_file : str
        monomer ファイル名 (例: "chol_monomer")。
    aij_file : str
        aij ファイル名 (例: "aij.dat")。
    box_size : Tuple[float, float, float]
        DPD box (x, y, z)。
    density : float
        粒子密度 (default 3.0)。
    dt : float
        time step (default 0.05)。
    ratio_list : List[Dict[str, Any]]
        各 monomer のモル比設定。
    name_tail_list : List[str]
        出力ファイル名 suffix。
    phys_param : Dict[str, Any]
        物理パラメータ全般 (gamma / lambda / cutoff / temperature 等)。
    default_file : Dict[str, str]
        補助ファイル名辞書 (filename / init_pos_file / aij_file / ...)。
    """

    total_num_list: List[int]
    step_list: List[int]
    monomer_file: str
    aij_file: str
    box_size: Tuple[float, float, float]
    density: float = 3.0
    dt: float = 0.05
    ratio_list: List[Dict[str, Any]] = field(default_factory=list)
    name_tail_list: List[str] = field(default_factory=list)
    phys_param: Dict[str, Any] = field(default_factory=dict)
    default_file: Dict[str, str] = field(default_factory=dict)
    output_interval: int = 100
    aii_val: float = 25.0


@dataclass
class DpdSystemSpec:
    """R1 (UDF) / R2 (DPM) 出力に渡す全体システム仕様。

    Attributes
    ----------
    monomers : List[MonomerSpec]
        参加 monomer (cholesterol、 water 等の混合系で複数指定可)。
    aij : AijMatrix
        相互作用パラメータ。
    calc_sett : Optional[CalcSett]
        R1 (UDF 出力) 時に必須、 R2 (DPM 出力) のみなら None で可。
    project_name : str
        Cognac の `ProjectName` フィールド。
    comment : str
        ヘッダー comment。
    """

    monomers: List[MonomerSpec]
    aij: AijMatrix
    calc_sett: Optional[CalcSett] = None
    project_name: str = "abmptools-cg-dpd"
    comment: str = "Generated by abmptools.cg.dpd"

    def segment_names(self) -> List[str]:
        """全 monomer の particle_names を union (登場順)。"""
        seen: Dict[str, None] = {}
        for m in self.monomers:
            for p in m.particle_names:
                seen.setdefault(p, None)
        return list(seen.keys())

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
