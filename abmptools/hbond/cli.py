"""
cli.py
------
Command-line entry point for abmptools.hbond.

Usage::

    python -m abmptools.hbond <bdf_path> [options]
"""
from __future__ import annotations

import argparse
import sys

from .analyzer import Analyzer, AnalyzerConfig
from .hbond_detector import HBondCriteria


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m abmptools.hbond",
        description="Detect carboxyl/amide H-bonds in COGNAC UDF/BDF trajectory",
    )
    parser.add_argument("bdf_path", help="Input COGNAC UDF or BDF file")
    parser.add_argument(
        "--out-prefix", "-o", default="hbond_result",
        help="Prefix for output files (default: hbond_result)"
    )
    parser.add_argument(
        "--criteria", choices=["luzar-chandler", "strict", "custom"],
        default="luzar-chandler",
        help="H-bond geometric criteria (default: luzar-chandler)"
    )
    parser.add_argument(
        "--d-da-max", type=float, default=None,
        help="(custom) max donor-acceptor distance Å"
    )
    parser.add_argument(
        "--d-ha-max", type=float, default=None,
        help="(custom) max H-acceptor distance Å"
    )
    parser.add_argument(
        "--angle-min", type=float, default=None,
        help="(custom) min D-H-A angle in degrees"
    )
    parser.add_argument(
        "--record-start", type=int, default=0,
        help="Starting record (default: 0)"
    )
    parser.add_argument(
        "--record-end", type=int, default=-1,
        help="End record (exclusive); -1 = all"
    )
    parser.add_argument(
        "--mol-name", default="IMC",
        help="Base molecule name for renamed groups (default: IMC)"
    )
    parser.add_argument(
        "--no-colorize", action="store_true",
        help="Skip writing colored BDF"
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Skip writing count plot"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Suppress per-frame output"
    )
    args = parser.parse_args(argv)

    custom = None
    if args.criteria == "custom":
        custom = HBondCriteria(
            d_da_max=args.d_da_max if args.d_da_max is not None else 3.5,
            d_ha_max=args.d_ha_max,
            angle_min=args.angle_min if args.angle_min is not None else 120.0,
        )

    config = AnalyzerConfig(
        bdf_path=args.bdf_path,
        out_prefix=args.out_prefix,
        criteria_mode=args.criteria,
        custom_criteria=custom,
        record_start=args.record_start,
        record_end=args.record_end,
        base_mol_name=args.mol_name,
        do_colorize=not args.no_colorize,
        do_plot=not args.no_plot,
        verbose=not args.quiet,
    )

    analyzer = Analyzer(config)
    analyzer.load()
    analyzer.run()
    out = analyzer.write_outputs()

    print("\nGenerated outputs:")
    for kind, path in out.items():
        print(f"  {kind}: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
