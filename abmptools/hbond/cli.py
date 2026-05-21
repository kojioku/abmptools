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
        "--classify-mode", choices=["imc", "generic"], default="imc",
        help="imc: COOH-centric 4-species (dual/chain/single/free, IMC default). "
             "generic: donor-type x acceptor-type pair stats (PVA/peptide/etc.)"
    )
    parser.add_argument(
        "--no-colorize", action="store_true",
        help="Skip writing any colored output"
    )
    parser.add_argument(
        "--colorize-mode", choices=["molname", "action", "both"],
        default="molname",
        help="molname: rename Mol_Name + Draw_Attributes (default, v1.25 legacy); "
             "action: Python action (.act) per functional group (Mol_Name kept, "
             "J-OCTA pre-render compatible); "
             "both: emit both"
    )
    parser.add_argument(
        "--no-copy-uncolored", action="store_true",
        help="Skip writing uncolored BDF copy (<prefix>.bdf), "
             "which is intended for J-OCTA pre-render compatibility"
    )
    parser.add_argument(
        "--no-write-attributes", action="store_true",
        help="Skip appending per-atom Attributes[] hbond=Dual/Single/Free/Accept "
             "tags to <prefix>.bdf"
    )
    parser.add_argument(
        "--attributes-name", default="hbond",
        help="Attribute Name to write (default: hbond)"
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Skip writing count plot"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Suppress per-frame output"
    )
    # Force field + functional group selection
    parser.add_argument(
        "--force-field", default=None,
        help="Force field name (auto-detect if not specified). "
             "Built-in: GAFF2, OPLS-AA, CHARMM36, OpenFF"
    )
    parser.add_argument(
        "--donor-groups", default=None,
        help="Comma-separated donor groups. Options: "
             "carboxyl, amide_donor, amine_donor, hydroxyl. "
             "Default: carboxyl (= COOH OH only)"
    )
    parser.add_argument(
        "--acceptor-groups", default=None,
        help="Comma-separated acceptor groups. Options: "
             "carboxyl_O, amide_O, hydroxyl_O, ether_O. "
             "Default: carboxyl_O,amide_O"
    )
    # Lifetime analysis
    parser.add_argument(
        "--no-lifetime", action="store_true",
        help="Skip lifetime / autocorrelation computation"
    )
    parser.add_argument(
        "--gap-tolerance", type=int, default=0,
        help="Frames of gap to bridge for intermittent lifetime (default 0=strict)"
    )
    parser.add_argument(
        "--dt", type=float, default=1.0,
        help="Time between consecutive records (user units, e.g. ps). "
             "Default 1.0 (= per-frame)"
    )
    parser.add_argument(
        "--autocorr-max-lag", type=int, default=None,
        help="Max lag (in records) for autocorrelation C(t). Default: N/2"
    )
    args = parser.parse_args(argv)

    custom = None
    if args.criteria == "custom":
        custom = HBondCriteria(
            d_da_max=args.d_da_max if args.d_da_max is not None else 3.5,
            d_ha_max=args.d_ha_max,
            angle_min=args.angle_min if args.angle_min is not None else 120.0,
        )

    donor_groups = (
        [s.strip() for s in args.donor_groups.split(",") if s.strip()]
        if args.donor_groups else None
    )
    acceptor_groups = (
        [s.strip() for s in args.acceptor_groups.split(",") if s.strip()]
        if args.acceptor_groups else None
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
        colorize_mode=args.colorize_mode,
        classify_mode=args.classify_mode,
        do_copy_uncolored=not args.no_copy_uncolored,
        do_write_attributes=not args.no_write_attributes,
        attributes_name=args.attributes_name,
        do_plot=not args.no_plot,
        verbose=not args.quiet,
        force_field=args.force_field,
        donor_groups=donor_groups,
        acceptor_groups=acceptor_groups,
        compute_lifetime=not args.no_lifetime,
        gap_tolerance=args.gap_tolerance,
        dt=args.dt,
        autocorr_max_lag=args.autocorr_max_lag,
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
