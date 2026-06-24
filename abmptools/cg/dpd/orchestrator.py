# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.orchestrator
-----------------------------
``CGDpdBuilder`` クラス: cg_segmenter の monomer + fcews の aij.dat + (R1 のみ
追加で calc_sett) から DPD 用ファイルを生成するパイプライン。

Routes
------
- **R1 (UDF 直接)**: ``build_udf(output_path)`` — Cognac DPD 入力 UDF
- **R2 (DPM 経由)**: ``build_dpm(template_path, output_dir, virtual_mom_template)``
  — J-OCTA GUI 用 ``*.dpm`` + ``monomer-lib/<seg>/Virtual.mom`` + ``#Message.txt``
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .aij_io import read_aij
from .calc_sett_io import read_calc_sett
from .dpm_writer import (
    patch_dpm, propagate_virtual_mom, write_message_txt,
)
from .models import AijMatrix, CalcSett, DpdSystemSpec, MonomerSpec
from .monomer_io import assign_particle_names, read_monomer
from .udf_writer import write_dpd_udf
from .udf_writer_udfm import write_dpd_udf_udfm

logger = logging.getLogger(__name__)


@dataclass
class CGDpdBuilder:
    """DPD 系構築 pipeline (R1 + R2 共通)。

    通常のコンストラクタは :meth:`from_files` を使う:

        >>> sg = CGDpdBuilder.from_files(
        ...     monomer="chol_monomer", aij="aij.dat",
        ...     calc_sett="chol_calc_sett",  # R1 (UDF) では必須、 R2 (DPM) では省略可
        ... )
        >>> sg.build_dpm(  # R2: B 案 (user template + patch)
        ...     template="empty.dpm",
        ...     output_dir="./chol_project",
        ...     virtual_mom_template="Virtual.mom",  # user が J-OCTA で 1 回作成
        ... )
    """

    spec: DpdSystemSpec

    @classmethod
    def from_files(
        cls,
        monomer: Union[str, Path],
        aij: Union[str, Path],
        calc_sett: Optional[Union[str, Path]] = None,
        *,
        particle_names: Optional[List[str]] = None,
        project_name: str = "abmptools-cg-dpd",
    ) -> "CGDpdBuilder":
        """既存の ``{name}_monomer`` + ``aij.dat`` + (任意で) ``{name}_calc_sett`` から構築。

        Parameters
        ----------
        monomer : str | Path
            cg_segmenter dpdgen_exporter が生成した monomer file。
        aij : str | Path
            fcews 統一形式の aij.dat (``aij = [['W','A', 55.0], ...]`` Python script)。
        calc_sett : str | Path, optional
            cg_segmenter が生成した calc_sett file。 R1 (UDF) で必要、 R2 (DPM) のみ
            なら省略可。
        particle_names : List[str], optional
            cg_segmenter は ``P0..Pn-1`` の汎用 label しか出さないので、 実 segment 名
            (例 ``["segA","segB","segC","segA","WAT"]``) を渡したい場合に指定。
        project_name : str
            dpm / UDF の ``ProjectName`` ヘッダ。
        """
        mono = read_monomer(monomer)
        if particle_names is not None:
            mono = assign_particle_names(mono, particle_names)

        aii_val = 25.0
        sett: Optional[CalcSett] = None
        if calc_sett is not None:
            sett = read_calc_sett(calc_sett)
            aii_val = sett.aii_val
        aij_matrix = read_aij(aij, aii=aii_val)

        spec = DpdSystemSpec(
            monomers=[mono],
            aij=aij_matrix,
            calc_sett=sett,
            project_name=project_name,
        )
        logger.info(
            "CGDpdBuilder: %d monomer(s), %d aij pair(s), %s calc_sett",
            len(spec.monomers), len(aij_matrix.pairs),
            "with" if sett else "without",
        )
        return cls(spec=spec)

    @classmethod
    def from_multi_files(
        cls,
        monomer_specs: List[Dict[str, Any]],
        aij: Union[str, Path],
        calc_sett: Optional[Union[str, Path]] = None,
        *,
        project_name: str = "abmptools-cg-dpd",
    ) -> "CGDpdBuilder":
        """複数 monomer の混合系 (cholesterol + water 等) を構築する。

        Parameters
        ----------
        monomer_specs : List[Dict[str, Any]]
            各要素 ``{"monomer_file": "path/to/xxx_monomer", "particle_names": [...]}``。
            ``particle_names`` を ``None`` or 省略すれば cg_segmenter 出力の汎用ラベル
            (P0..Pn-1) を維持。
        aij : str | Path
            混合系全体での fcews aij.dat (全 monomer の segment 名を含む必要あり)。
        calc_sett : str | Path, optional
            R1 (UDF) で必要、 R2 (DPM) のみなら省略可。
        project_name : str
            dpm / UDF の ProjectName ヘッダ。
        """
        monomers: List[MonomerSpec] = []
        for spec_dict in monomer_specs:
            if "monomer_file" not in spec_dict:
                raise ValueError(
                    f"monomer_specs entry missing 'monomer_file' key: {spec_dict}"
                )
            mono = read_monomer(spec_dict["monomer_file"])
            names = spec_dict.get("particle_names")
            if names is not None and len(names) > 0:
                mono = assign_particle_names(mono, list(names))
            monomers.append(mono)

        aii_val = 25.0
        sett: Optional[CalcSett] = None
        if calc_sett is not None:
            sett = read_calc_sett(calc_sett)
            aii_val = sett.aii_val
        aij_matrix = read_aij(aij, aii=aii_val)

        spec = DpdSystemSpec(
            monomers=monomers, aij=aij_matrix, calc_sett=sett,
            project_name=project_name,
        )

        # 整合性 warning: monomer の segment 名が aij.dat にあるか
        seg_set = set(spec.segment_names())
        aij_set = set(aij_matrix.segments)
        missing = seg_set - aij_set
        if missing:
            logger.warning(
                "CGDpdBuilder.from_multi_files: %d segment(s) not in aij.dat: %s "
                "(R1/R2 出力で default aii=%g が fallback されます)",
                len(missing), sorted(missing), aii_val,
            )

        logger.info(
            "CGDpdBuilder.from_multi_files: %d monomers, %d aij pair(s), "
            "%d unique segments",
            len(spec.monomers), len(aij_matrix.pairs), len(seg_set),
        )
        return cls(spec=spec)

    def validate(self) -> List[str]:
        """spec の整合性を検証して warning メッセージ list を返す。

        チェック項目:
        - monomer の particle 名が aij.dat の segments にあるか (missing)
        - aij.dat の segments で monomer から参照されないものがあるか (extra)

        Returns
        -------
        List[str]
            warning メッセージ (空 list なら問題なし)。
        """
        warnings: List[str] = []
        seg_set = set(self.spec.segment_names())
        aij_set = set(self.spec.aij.segments)
        missing = seg_set - aij_set
        if missing:
            warnings.append(
                f"{len(missing)} segment(s) used by monomer but NOT in aij.dat: "
                f"{sorted(missing)} — R1/R2 出力で default aii={self.spec.aij.aii} に fallback"
            )
        extra = aij_set - seg_set
        if extra:
            warnings.append(
                f"{len(extra)} segment(s) defined in aij.dat but NOT used by any monomer: "
                f"{sorted(extra)} — 余分なペアが Pair_Interaction[] / SegmentPairModel[] に残る"
            )
        return warnings

    def build_dpm(
        self,
        template: Union[str, Path],
        output_dir: Union[str, Path],
        virtual_mom_template: Optional[Union[str, Path]] = None,
        *,
        dpm_filename: str = "system.dpm",
    ) -> Path:
        """**R2**: B 案で dpm ファイル + monomer-lib/<seg>/Virtual.mom + #Message.txt 生成。

        Parameters
        ----------
        template : str | Path
            ユーザーが J-OCTA で作成した空 dpm template。
        output_dir : str | Path
            出力ディレクトリ。 以下を生成:

            - ``{output_dir}/{dpm_filename}``
            - ``{output_dir}/{name}/#Message.txt``
            - ``{output_dir}/{name}/monomer-lib/<seg>/Virtual.mom`` (virtual_mom_template
              指定時)
        virtual_mom_template : str | Path, optional
            J-OCTA Monomer Modeler で作成した Virtual.mom。 省略時は Virtual.mom 配置
            をスキップ (dpm のみ生成)。
        dpm_filename : str
            dpm ファイル名 (default ``"system.dpm"``)。
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 1) dpm patch
        dpm_path = output_dir / dpm_filename
        patch_dpm(template_path=template, output_path=dpm_path, spec=self.spec)

        # 2) monomer-lib/<seg>/Virtual.mom 配置 (virtual_mom_template 指定時)
        if virtual_mom_template is not None:
            stem = Path(dpm_filename).stem
            seg_dir = output_dir / stem / "monomer-lib"
            segs = self.spec.segment_names()
            propagate_virtual_mom(
                template_path=virtual_mom_template,
                segments=segs,
                output_dir=seg_dir,
            )
            # 3) #Message.txt
            write_message_txt(output_dir / stem)

        logger.info("CGDpdBuilder.build_dpm: -> %s", dpm_path)
        return dpm_path

    def build_udf(
        self,
        output_path: Union[str, Path],
        *,
        include_file: str = "cognac112.udf",
    ) -> Path:
        """**R1**: Cognac DPD 入力 UDF を **UDFManager** で出力 (cognac112 スキーマ準拠)。

        named-path put で組み立てるため cognac でロード可能な UDF になる。
        ``UDFManager`` (OCTA 同梱) が必須 (`include_file` の cognac クラス定義も
        UDFManager が解決できる場所に必要)。

        Parameters
        ----------
        output_path : str | Path
            出力 UDF path (例: ``chol_uin.udf``)。
        include_file : str
            冒頭 ``\\include`` で参照する Cognac class 定義 file (default
            ``"cognac112.udf"``、 J-OCTA 11.x 同梱)。
        """
        return write_dpd_udf_udfm(self.spec, output_path, include_file=include_file)
