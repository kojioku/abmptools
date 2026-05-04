# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.cli
--------------------------
Command-line interface for the Martini 3 peptide-membrane US builder.

Subcommands (filled in incrementally across implementation steps):
    example         print an example JSON config to stdout (Step 1)
    validate        check config + dependencies + Martini 3 .itp files (Step 2)
    build           run the full 7-stage build pipeline (Step 6)
    make-windows    extract per-window start frames from pull.xtc (Step 6)
    wham            run gmx wham over completed windows (Step 5)
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from typing import List, Optional

from pathlib import Path

from . import forcefield_check, pmf, pulling
from ._subprocess import setup_logging
from .builder import MembraneCGBuilder
from .models import (
    EquilibrationCGProtocol,
    LipidMix,
    MembraneCGBuildConfig,
    PeptideMembraneSpec,
    PullingCGProtocol,
    UmbrellaCGProtocol,
)

logger = logging.getLogger(__name__)


def _default_example() -> MembraneCGBuildConfig:
    """A minimal, consistent example for ``example`` subcommand."""
    return MembraneCGBuildConfig(
        lipids=[LipidMix(resname="POPC", n_per_leaflet=64)],
        peptide=PeptideMembraneSpec(name="kgg", sequence="KGG"),
        martini_itp_dir="./ff",
        insane_d_nm=8.0,
        box_z_nm=14.0,
    )


def _load_config(path: str) -> MembraneCGBuildConfig:
    p = Path(path)
    if p.suffix.lower() in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore
        except ImportError as exc:
            raise RuntimeError(
                "PyYAML is required for YAML configs. "
                "Install via: pip install abmptools[cg]"
            ) from exc
        raw = yaml.safe_load(p.read_text())
        lipids = [LipidMix(**lp) for lp in raw.pop("lipids", [])]
        pep_raw = raw.pop("peptide", None)
        peptide = PeptideMembraneSpec(**pep_raw) if pep_raw else None
        eq = EquilibrationCGProtocol(**raw.pop("equilibration"))
        pull = PullingCGProtocol(**raw.pop("pulling"))
        umb = UmbrellaCGProtocol(**raw.pop("umbrella"))
        return MembraneCGBuildConfig(
            lipids=lipids,
            peptide=peptide,
            equilibration=eq,
            pulling=pull,
            umbrella=umb,
            **raw,
        )
    return MembraneCGBuildConfig.from_json(str(p))


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def _cmd_example(args: argparse.Namespace) -> int:
    cfg = _default_example()
    print(json.dumps(asdict(cfg), indent=2, ensure_ascii=False))
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    n_lipids = sum(lp.n_per_leaflet for lp in cfg.lipids) * 2
    print(
        f"Configuration: OK ({cfg.lipids[0].resname} bilayer, "
        f"{n_lipids} lipids total, peptide={cfg.peptide.name}, "
        f"T={cfg.equilibration.temperature_K} K, "
        f"NaCl {cfg.nacl_molar} M, "
        f"{cfg.umbrella.n_windows} umbrella windows)\n"
    )

    ff_dir = args.ff_dir or cfg.martini_itp_dir or ""
    tools = forcefield_check.check_external_tools(
        insane=cfg.insane_path,
        martinize2=cfg.martinize2_path,
        gmx=cfg.gmx_path,
        tleap=cfg.tleap_path,
    )
    files = forcefield_check.check_martini_files(ff_dir)
    ok = forcefield_check.report(tools, files, ff_dir)
    return 0 if ok else 1


def _cmd_build(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    if args.output_dir:
        cfg.output_dir = args.output_dir
    if args.ff_dir:
        cfg.martini_itp_dir = args.ff_dir
    builder = MembraneCGBuilder(cfg)
    result = builder.build()
    print(f"Build complete. Output: {result['output_dir']}")
    print(f"  final coordinates: {result['gro']}")
    print(f"  topology:          {result['top']}")
    print(f"  index:             {result['ndx']}")
    print(f"  run script:        {result['run_script']}")
    print(f"  windows: {result['n_windows']}, lipid: {result['lipid_resname']}")
    return 0


def _cmd_make_windows(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    if args.gmx_path:
        cfg.gmx_path = args.gmx_path
    out = pulling.extract_window_frames(
        pull_xtc=args.pull_xtc,
        pull_tpr=args.pull_tpr,
        pullx_xvg=args.pull_xvg,
        config=cfg,
        out_dir=args.windows_dir,
    )
    print(f"Extracted {len(out)} window frames into {args.windows_dir}")
    return 0


def _cmd_wham(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    if args.gmx_path:
        cfg.gmx_path = args.gmx_path
    result = pmf.run_wham(
        windows_dir=Path(args.windows_dir),
        analysis_dir=Path(args.analysis_dir),
        config=cfg,
        bootstrap_n=args.bootstrap_n,
        temperature_K=args.temperature,
    )
    print(f"WHAM complete. PMF: {result['pmf']}")
    print(f"  histograms: {result['histo']}")
    return 0


# ``build`` / ``make-windows`` are added in the subsequent implementation
# steps. The parser registers their stubs so argparse error messages are
# useful even before the handlers are wired.

def _cmd_not_yet_implemented(args: argparse.Namespace) -> int:
    print(
        f"Subcommand '{args.cmd}' is not yet implemented at this revision.",
        file=sys.stderr,
    )
    return 2


# ---------------------------------------------------------------------------
# Parser / main
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="abmptools.cg.membrane",
        description=(
            "Build a Martini 3 peptide-membrane umbrella-sampling system "
            "for GROMACS MD."
        ),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable DEBUG-level logging.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_ex = sub.add_parser(
        "example",
        help="Print an example JSON config to stdout.",
    )
    p_ex.set_defaults(func=_cmd_example)

    p_val = sub.add_parser(
        "validate",
        help=(
            "Validate config + check dependencies (insane / gmx / "
            "martinize2 / tleap) + Martini 3 .itp files."
        ),
    )
    p_val.add_argument(
        "--config", required=True, type=str,
        help="Path to config JSON or YAML.",
    )
    p_val.add_argument(
        "--ff-dir", type=str, default=None,
        help="Override martini_itp_dir.",
    )
    p_val.set_defaults(func=_cmd_validate)

    p_wham = sub.add_parser(
        "wham",
        help="Run gmx wham over completed windows (PMF analysis).",
    )
    p_wham.add_argument("--config", required=True, type=str)
    p_wham.add_argument(
        "--windows-dir", required=True, type=str,
        help="Path to the windows/ directory.",
    )
    p_wham.add_argument(
        "--analysis-dir", required=True, type=str,
        help="Output directory for pmf.xvg / histo.xvg.",
    )
    p_wham.add_argument(
        "--bootstrap-n", type=int, default=0,
        help="Bootstrap samples (0 disables).",
    )
    p_wham.add_argument(
        "--temperature", type=float, default=None,
        help="Override temperature in K (default: config.equilibration).",
    )
    p_wham.add_argument(
        "--gmx-path", type=str, default=None,
        help="Override gmx CLI path.",
    )
    p_wham.set_defaults(func=_cmd_wham)

    p_bld = sub.add_parser(
        "build",
        help="Run the full 7-stage build pipeline (insane + topology compose + MDPs + run.sh).",
    )
    p_bld.add_argument("--config", required=True, type=str)
    p_bld.add_argument(
        "-o", "--output-dir", type=str, default=None,
        help="Override config.output_dir.",
    )
    p_bld.add_argument(
        "--ff-dir", type=str, default=None,
        help="Override martini_itp_dir.",
    )
    p_bld.set_defaults(func=_cmd_build)

    p_mw = sub.add_parser(
        "make-windows",
        help="Extract per-window starting frames from pull.xtc.",
    )
    p_mw.add_argument("--config", required=True, type=str)
    p_mw.add_argument("--pull-tpr", required=True, type=str)
    p_mw.add_argument("--pull-xtc", required=True, type=str)
    p_mw.add_argument("--pull-xvg", required=True, type=str)
    p_mw.add_argument(
        "--windows-dir", required=True, type=str,
        help="Output directory for windows/win_NNN/start.gro files.",
    )
    p_mw.add_argument(
        "--gmx-path", type=str, default=None,
        help="Override gmx CLI path.",
    )
    p_mw.set_defaults(func=_cmd_make_windows)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    setup_logging(verbose=getattr(args, "verbose", False))
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
