# -*- coding: utf-8 -*-
"""``python -m abmptools.cg.dpd`` CLI."""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from .orchestrator import CGDpdBuilder


def _load_multi_monomer(json_path: str) -> list:
    """multi-monomer JSON ファイルを読んで list[dict] を返す。

    Expected format::

        [
          {"monomer_file": "chol_monomer",
           "particle_names": ["A_Chol","B_Chol","C_Chol","D_Chol","Tail_Chol"]},
          {"monomer_file": "wat_monomer", "particle_names": ["W"]}
        ]
    """
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list):
        raise ValueError(
            f"multi-monomer JSON must be a list, got {type(data).__name__}"
        )
    return data


def _build_builder(args):
    """共通 builder 構築 (single / multi monomer 切替)。"""
    if args.multi_monomer and args.monomer:
        print("ERROR: --monomer と --multi-monomer は排他です", file=sys.stderr)
        sys.exit(2)
    if not args.multi_monomer and not args.monomer:
        print("ERROR: --monomer または --multi-monomer のどちらか必須です", file=sys.stderr)
        sys.exit(2)
    if args.multi_monomer:
        mspecs = _load_multi_monomer(args.multi_monomer)
        return CGDpdBuilder.from_multi_files(
            monomer_specs=mspecs, aij=args.aij,
            calc_sett=getattr(args, "calc_sett", None),
            project_name=args.project_name,
        )
    names = args.particle_names.split(",") if args.particle_names else None
    return CGDpdBuilder.from_files(
        monomer=args.monomer, aij=args.aij,
        calc_sett=getattr(args, "calc_sett", None),
        particle_names=names, project_name=args.project_name,
    )


def main(argv=None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s [%(name)s] %(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="python -m abmptools.cg.dpd",
        description=(
            "abmptools.cg.dpd: DPD 系構築 (R1=UDF / R2=DPM)。 "
            "cg_segmenter monomer + fcews aij.dat + calc_sett から COGNAC DPD "
            "入力 UDF / J-OCTA dpm を生成。"
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # build-dpm (R2)
    p2 = sub.add_parser("build-dpm", help="R2: J-OCTA dpm + monomer-lib 生成 (B 案: user template + patch)")
    p2.add_argument("--monomer", default=None, help="cg_segmenter monomer file (single)")
    p2.add_argument("--multi-monomer", default=None,
                    help="multi-monomer 混合系: JSON file with [{monomer_file, particle_names}, ...]")
    p2.add_argument("--aij", required=True, help="fcews aij.dat (Python 辞書)")
    p2.add_argument("--calc-sett", default=None, help="calc_sett file (optional for R2)")
    p2.add_argument("--template", required=True, help="user 作成の空 dpm template")
    p2.add_argument("--virtual-mom", default=None, help="user 作成の Virtual.mom (各 segment dir に copy)")
    p2.add_argument("--output-dir", required=True, help="出力ディレクトリ")
    p2.add_argument("--dpm-filename", default="system.dpm", help="dpm ファイル名 (default system.dpm)")
    p2.add_argument("--project-name", default="abmptools-cg-dpd", help="dpm ProjectName")
    p2.add_argument("--particle-names", default=None,
                    help="comma-separated segment 名 (single monomer 時、 例 'segA,segB,WAT')")

    # verify: spec 整合性チェック (実機 dry-run 用)
    pv = sub.add_parser("verify", help="cg.dpd spec の整合性をチェック (実機 build 前 dry-run)")
    pv.add_argument("--monomer", default=None, help="single monomer file")
    pv.add_argument("--multi-monomer", default=None, help="multi-monomer JSON file")
    pv.add_argument("--aij", required=True)
    pv.add_argument("--calc-sett", default=None)
    pv.add_argument("--particle-names", default=None)
    pv.add_argument("--project-name", default="abmptools-cg-dpd")

    # build-udf (R1)
    p1 = sub.add_parser("build-udf", help="R1: Cognac DPD 入力 UDF を生成")
    p1.add_argument("--monomer", default=None, help="cg_segmenter monomer file (single)")
    p1.add_argument("--multi-monomer", default=None,
                    help="multi-monomer 混合系: JSON file with [{monomer_file, particle_names}, ...]")
    p1.add_argument("--aij", required=True)
    p1.add_argument("--calc-sett", required=True)
    p1.add_argument("--output", required=True, help="出力 UDF path (例 chol_uin.udf)")
    p1.add_argument("--include-file", default="cognac112.udf",
                    help="冒頭 \\include の Cognac class 定義 file (default cognac112.udf)")
    p1.add_argument("--particle-names", default=None,
                    help="comma-separated segment 名 (single monomer 時、 例 'segA,segB,WAT')")
    p1.add_argument("--project-name", default="abmptools-cg-dpd")

    args = parser.parse_args(argv)

    if args.cmd == "verify":
        builder = _build_builder(args)
        warnings = builder.validate()
        if not warnings:
            print(
                f"OK: {len(builder.spec.monomers)} monomer(s), "
                f"{len(builder.spec.aij.pairs)} aij pair(s), "
                f"{len(builder.spec.segment_names())} unique segment(s) — 整合性 OK"
            )
            return 0
        print(f"WARNING: {len(warnings)} issue(s) found:")
        for w in warnings:
            print(f"  - {w}")
        return 1

    if args.cmd == "build-dpm":
        builder = _build_builder(args)
        dpm_path = builder.build_dpm(
            template=args.template,
            output_dir=args.output_dir,
            virtual_mom_template=args.virtual_mom,
            dpm_filename=args.dpm_filename,
        )
        print(f"dpm: {dpm_path}")
        return 0

    if args.cmd == "build-udf":
        builder = _build_builder(args)
        out = builder.build_udf(args.output, include_file=args.include_file)
        print(f"udf: {out}")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
