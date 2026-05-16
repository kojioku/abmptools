# -*- coding: utf-8 -*-
"""``python -m abmptools.cg.dpd`` CLI."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .orchestrator import CGDpdBuilder


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
    p2.add_argument("--monomer", required=True, help="cg_segmenter monomer file")
    p2.add_argument("--aij", required=True, help="fcews aij.dat (Python 辞書)")
    p2.add_argument("--calc-sett", default=None, help="calc_sett file (optional for R2)")
    p2.add_argument("--template", required=True, help="user 作成の空 dpm template")
    p2.add_argument("--virtual-mom", default=None, help="user 作成の Virtual.mom (各 segment dir に copy)")
    p2.add_argument("--output-dir", required=True, help="出力ディレクトリ")
    p2.add_argument("--dpm-filename", default="system.dpm", help="dpm ファイル名 (default system.dpm)")
    p2.add_argument("--project-name", default="abmptools-cg-dpd", help="dpm ProjectName")
    p2.add_argument("--particle-names", default=None,
                    help="comma-separated segment 名 (例 'segA,segB,segC,segA,WAT')")

    # build-udf (R1) - placeholder for Phase 3
    p1 = sub.add_parser("build-udf", help="R1: Cognac DPD 入力 UDF 生成 (Phase 3 で実装予定)")
    p1.add_argument("--monomer", required=True)
    p1.add_argument("--aij", required=True)
    p1.add_argument("--calc-sett", required=True)
    p1.add_argument("--output", required=True)

    args = parser.parse_args(argv)

    if args.cmd == "build-dpm":
        names = args.particle_names.split(",") if args.particle_names else None
        builder = CGDpdBuilder.from_files(
            monomer=args.monomer, aij=args.aij, calc_sett=args.calc_sett,
            particle_names=names, project_name=args.project_name,
        )
        dpm_path = builder.build_dpm(
            template=args.template,
            output_dir=args.output_dir,
            virtual_mom_template=args.virtual_mom,
            dpm_filename=args.dpm_filename,
        )
        print(f"dpm: {dpm_path}")
        return 0

    if args.cmd == "build-udf":
        print("ERROR: build-udf は Phase 3 で実装予定 (まだ未実装)", file=sys.stderr)
        return 1

    return 1


if __name__ == "__main__":
    sys.exit(main())
