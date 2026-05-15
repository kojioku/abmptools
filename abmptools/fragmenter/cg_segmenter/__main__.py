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
