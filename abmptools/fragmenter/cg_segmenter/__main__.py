# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.__main__
------------------------------------------
CLI for CG segment construction.

Usage:
    python -m abmptools.fragmenter.cg_segmenter build \\
        --pdb input.pdb --output-dir ./segs --target-mw 200

    python -m abmptools.fragmenter.cg_segmenter example > config.json
"""
from __future__ import annotations

import argparse
import logging
import sys

from .models import CGSegmenterConfig
from .orchestrator import CGSegmenter


def cmd_build(args: argparse.Namespace) -> int:
    config = CGSegmenterConfig(
        pdb_path=args.pdb,
        output_dir=args.output_dir,
        target_mw=args.target_mw,
        separate_rings=not args.no_separate_rings,
        allow_atom_sharing=not args.no_atom_sharing,
        absorb_single_substituent=not args.no_absorb_substituent,
    )
    if args.h_only:
        config.hetero_cap_methyl_elements = []

    sg = CGSegmenter.from_pdb(config)
    result = sg.export()

    print(f"Wrote {result.total_atoms_with_cap} atoms across "
          f"{len(result.segments)} segments to: {config.output_dir}")
    for seg in result.segments:
        n_methyl = sum(1 for c in seg.cap_atoms if c.is_methyl_cap)
        n_h = sum(1 for c in seg.cap_atoms if not c.is_methyl_cap)
        print(
            f"  seg {seg.segment_id:>3} ({seg.kind:<26}): "
            f"{len(seg.atom_indices):>3} heavy, {n_methyl} CH3 cap, {n_h} H cap"
        )

    if result.shared_atom_pairs:
        print(f"\nShared atoms (fused-ring boundaries): {len(result.shared_atom_pairs)} pair(s)")
        for a, s1, s2 in result.shared_atom_pairs[:10]:
            print(f"  atom {a:>3}: seg {s1} & seg {s2}")
        if len(result.shared_atom_pairs) > 10:
            print(f"  ... ({len(result.shared_atom_pairs) - 10} more)")
    return 0


def cmd_example(args: argparse.Namespace) -> int:
    config = CGSegmenterConfig(pdb_path="input.pdb")
    print(config.to_json())
    return 0


def cmd_dpdgen(args: argparse.Namespace) -> int:
    """build に加えて DPDgen 入力 (monomer + calc_sett) も生成する。"""
    config = CGSegmenterConfig(
        pdb_path=args.pdb,
        output_dir=args.output_dir,
        target_mw=args.target_mw,
    )
    sg = CGSegmenter.from_pdb(config)
    result = sg.export()
    monomer_path, calc_sett_path = sg.export_dpdgen(
        output_dir=args.output_dir,
        monomer_name=args.monomer_name,
        box_size=(args.box, args.box, args.box),
        total_num=args.total_num,
        step=args.step,
        aij_file=args.aij_file,
    )
    print(f"CG segments: {len(result.segments)} -> {config.output_dir}")
    print(f"DPDgen monomer:   {monomer_path}")
    print(f"DPDgen calc_sett: {calc_sett_path}")
    print(f"  Edit ratio_list / aij_file / phys_param in calc_sett, then run:")
    print(f"    python makeudf_dpd.py -p {calc_sett_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m abmptools.fragmenter.cg_segmenter",
        description=(
            "Coarse-grained (CG) segment construction from PDB. "
            "Splits the molecule into ring / chain segments with H or CH3 caps "
            "at boundaries; allows atom sharing across fused-ring segments."
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_build = sub.add_parser("build", help="build CG segments from a PDB file")
    p_build.add_argument("--pdb", required=True, help="input PDB file")
    p_build.add_argument("--output-dir", default="./cg_segments",
                         help="output dir (per-segment PDB + XYZ + segments.json)")
    p_build.add_argument("--target-mw", type=float, default=200.0,
                         help="chain split target MW (g/mol) -- same as fragmenter")
    p_build.add_argument("--no-separate-rings", action="store_true",
                         help="do not split rings into separate segments")
    p_build.add_argument("--no-atom-sharing", action="store_true",
                         help="for fused rings, do not share atoms across segments")
    p_build.add_argument("--no-absorb-substituent", action="store_true",
                         help="do not absorb 1-heavy-atom substituents (-OH, -NH2, -F) into rings")
    p_build.add_argument("--h-only", action="store_true",
                         help="cap ALL boundaries with H (skip the hetero -> CH3 rule)")
    p_build.set_defaults(func=cmd_build)

    p_ex = sub.add_parser("example", help="emit a sample CGSegmenterConfig as JSON to stdout")
    p_ex.set_defaults(func=cmd_example)

    p_dpd = sub.add_parser(
        "dpdgen",
        help="build CG segments + DPDgen input (monomer + calc_sett) for OCTA COGNAC",
    )
    p_dpd.add_argument("--pdb", required=True, help="input PDB file")
    p_dpd.add_argument("--output-dir", default="./cg_dpdgen",
                       help="output dir (per-segment PDB/XYZ + DPDgen monomer/calc_sett)")
    p_dpd.add_argument("--target-mw", type=float, default=200.0)
    p_dpd.add_argument("--monomer-name", default="cg",
                       help="output filename prefix: {name}_monomer / {name}_calc_sett")
    p_dpd.add_argument("--box", type=float, default=10.0,
                       help="cubic box size (DPD unit) for phys_param x=y=z")
    p_dpd.add_argument("--total-num", type=int, default=100000,
                       help="total_num_list[0] in calc_sett")
    p_dpd.add_argument("--step", type=int, default=100,
                       help="step_list[0] in calc_sett")
    p_dpd.add_argument("--aij-file", default="aij.dat",
                       help="aij_file path written into calc_sett "
                            "(file itself is NOT generated; user provides via fcews-manybody)")
    p_dpd.set_defaults(func=cmd_dpdgen)

    return parser


def main(argv=None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s [%(name)s] %(message)s",
    )
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
