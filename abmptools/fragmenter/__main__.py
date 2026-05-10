# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.__main__
-----------------------------
CLI エントリポイント。

サブコマンド:
    suggest   PDB から自動提案 + レビューバンドル出力
    apply     編集 JSON を読み込んで segment_data.dat を出力
    example   example config を stdout 出力

使用例:
    python -m abmptools.fragmenter suggest \
        --pdb input.pdb --target-mw 200 --output-dir ./review
    python -m abmptools.fragmenter apply \
        --pdb input.pdb --review-dir ./review --output segment_data.dat
    python -m abmptools.fragmenter example > config.json
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .auto_split import suggest_cuts_for_groups
from .expand_to_system import export_to_system
from .grouping import group_by_smiles
from .headless_io import export_review_bundle, import_edited_review, resync_member_indices
from .models import FragmenterConfig
from .pdb_loader import load_pdb_molecules


def cmd_suggest(args: argparse.Namespace) -> int:
    config = FragmenterConfig(
        pdb_path=args.pdb,
        output_dir=args.output_dir,
        target_mw=args.target_mw,
        exclude_ring_cc=not args.no_exclude_ring,
        exclude_multibond=not args.no_exclude_multibond,
        exclude_heteroneighbor=not args.no_exclude_hetero,
        skip_protein_dna=True,
        include_residue_name=args.split_by_resname,
        include_c_heteroatom=args.include_c_heteroatom,
    )

    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded, include_residue_name=config.include_residue_name)
    suggest_cuts_for_groups(groups, loaded, config)
    export_review_bundle(config.output_dir, groups, loaded, config)

    print(f"Wrote review bundle: {config.output_dir}")
    print(f"  Groups: {len(groups)}")
    for i, g in enumerate(groups, start=1):
        print(f"  - Group {i}: SMILES={g.smiles}, n_copies={g.n_copies}, "
              f"cut_sites={len(g.cut_sites)}")
    return 0


def cmd_apply(args: argparse.Namespace) -> int:
    review_dir = args.review_dir
    output = args.output
    pdb_path = args.pdb

    config_path = Path(review_dir) / "config.json"
    if config_path.exists():
        config = FragmenterConfig.from_json(str(config_path))
        # PDB が config と異なる場合は CLI の --pdb を優先
        if pdb_path:
            config.pdb_path = pdb_path
    else:
        if not pdb_path:
            print("ERROR: --pdb is required when config.json is missing in review_dir.",
                  file=sys.stderr)
            return 1
        config = FragmenterConfig(pdb_path=pdb_path)

    loaded = load_pdb_molecules(config.pdb_path)
    edited_groups = import_edited_review(review_dir)
    resync_member_indices(edited_groups, loaded)

    result = export_to_system(config.pdb_path, edited_groups, loaded, output)

    print(f"Wrote segment_data: {output}")
    print(f"  Groups: {len(edited_groups)}")
    print(f"  Total fragments: {result.total_fragments}")
    return 0


def cmd_example(args: argparse.Namespace) -> int:
    config = FragmenterConfig(pdb_path="input.pdb")
    print(config.to_json())
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m abmptools.fragmenter",
        description=(
            "FMO automatic fragment splitter. Outputs a log2config-compatible "
            "segment_data.dat that pdb2fmo can consume."
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # suggest
    p_sg = sub.add_parser("suggest", help="auto-suggest cut sites + write review bundle")
    p_sg.add_argument("--pdb", required=True, help="input PDB file")
    p_sg.add_argument("--output-dir", default="./fragmenter_review",
                      help="output directory for review bundle (default: ./fragmenter_review)")
    p_sg.add_argument("--target-mw", type=float, default=200.0,
                      help="cut-site target molecular weight in g/mol (default: 200)")
    p_sg.add_argument("--no-exclude-ring", action="store_true",
                      help="allow cuts on ring-internal C-C bonds")
    p_sg.add_argument("--no-exclude-multibond", action="store_true",
                      help="allow cuts on double / triple bonds")
    p_sg.add_argument("--no-exclude-hetero", action="store_true",
                      help="allow cuts on C-C bonds adjacent to heteroatoms")
    p_sg.add_argument("--include-c-heteroatom", action="store_true",
                      help="also cut C-X single bonds (X = N/O/S/...). BDA = C side.")
    p_sg.add_argument("--split-by-resname", action="store_true",
                      help="separate groups by residue name even when SMILES match")
    p_sg.set_defaults(func=cmd_suggest)

    # apply
    p_ap = sub.add_parser("apply", help="apply edited review JSON -> segment_data.dat")
    p_ap.add_argument("--pdb", help="input PDB file (overrides config.json)")
    p_ap.add_argument("--review-dir", required=True, help="directory containing edited group_*.json")
    p_ap.add_argument("--output", default="segment_data.dat", help="output segment_data.dat path")
    p_ap.set_defaults(func=cmd_apply)

    # example
    p_ex = sub.add_parser("example", help="emit a sample FragmenterConfig as JSON to stdout")
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
