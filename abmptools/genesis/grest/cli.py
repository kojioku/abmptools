# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest.cli
---------------------------
Command-line interface for the GENESIS gREST_SSCR builder.

Subcommands (filled in incrementally across implementation steps):

    example      print an example JSON config to stdout (Step 6)
    validate     check config + external tools (Step 6)
    build        run the 5-stage build pipeline (Step 7)
    analyze      run remd_convert + transition / acceptance / PMF (Step 8)
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from typing import List, Optional

from pathlib import Path

from . import analysis, forcefield_check
from ._subprocess import setup_logging
from .builder import GrestBuilder
from .models import (
    GrestBuildConfig,
    ReplicaTemperatureSpec,
    RESTSelectionSpec,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Defaults / helpers
# ---------------------------------------------------------------------------

def _default_example() -> GrestBuildConfig:
    """A POC-style minimal config that ``validate`` accepts cleanly."""
    return GrestBuildConfig(
        input_pdb="protein.pdb",
        project_name="grest_run",
        rest_selection=RESTSelectionSpec(
            mode="explicit",
            residues=["1-138"],
        ),
        replica_temperatures=ReplicaTemperatureSpec(
            mode="manual",
            temperatures=[300.0, 318.11, 337.11, 357.10],
        ),
    )


def _load_config(path: str) -> GrestBuildConfig:
    """Load a config from JSON (YAML support deferred to v1.20.x)."""
    return GrestBuildConfig.from_json(path)


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def _cmd_example(args: argparse.Namespace) -> int:
    cfg = _default_example()
    print(json.dumps(asdict(cfg), indent=2, ensure_ascii=False))
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    ok = forcefield_check.report(cfg)
    return 0 if ok else 1


def _cmd_build(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    if args.output_dir:
        cfg.output_dir = args.output_dir
    builder = GrestBuilder(cfg)
    result = builder.build()
    print(f"Build complete. Output: {result['output_dir']}")
    print(f"  prmtop:        {result['prmtop']}")
    print(f"  coor:          {result['coor']}")
    print(
        f"  REST:          {result['rest_selection_string']} "
        f"({result['n_rest_residues']} residues)"
    )
    print(
        f"  ladder:        {result['n_replicas']} replicas "
        f"({result['temperature_ladder'][0]:.2f}-"
        f"{result['temperature_ladder'][-1]:.2f} K)"
    )
    print(f"  run script:    {result['run_script']}")
    return 0


def _cmd_analyze(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    run_dir = Path(args.run_dir).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else run_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    tasks = set(args.what)
    if "all" in tasks:
        tasks = {"sort", "transition", "acceptance", "pmf"}

    if "sort" in tasks:
        inp = run_dir / "inp" / "step5_remd_convert.inp"
        if not inp.exists():
            print(
                f"sort: missing {inp} -- did you run 'build' first?",
                file=sys.stderr,
            )
            return 1
        param1 = analysis.run_remd_convert(
            cfg, inp_path=inp, workdir=run_dir / "build"
        )
        print(f"sort: lowest-T trajectory -> {param1}")

    if "transition" in tasks:
        rem_files = sorted((run_dir / "build").glob("step3_rep*.rem"))
        if not rem_files:
            print("transition: no .rem files found in build/", file=sys.stderr)
            return 1
        analysis.plot_replica_transition(
            rem_files=rem_files,
            out_png=out_dir / "replica_transition.png",
            out_csv=out_dir / "replica_transition.csv",
        )

    if "acceptance" in tasks:
        log = run_dir / "logs" / "step3.log"
        if not log.exists():
            print(f"acceptance: missing {log}", file=sys.stderr)
            return 1
        analysis.plot_acceptance_ratio(
            log_path=log,
            out_png=out_dir / "acceptance_ratio.png",
            out_csv=out_dir / "acceptance_ratio.csv",
            burn_in=100,
        )

    if "pmf" in tasks:
        if not (args.mask1 and args.mask2):
            print(
                "pmf: --mask1 and --mask2 required (e.g. "
                "'rno:96 and an:NZ').",
                file=sys.stderr,
            )
            return 1
        prmtop = run_dir / "build" / "system.prmtop"
        dcd = run_dir / "build" / "param1.dcd"
        T = args.temperature
        if T is None:
            T = cfg.replica_temperatures.temperatures[0] if (
                cfg.replica_temperatures.mode == "manual"
            ) else cfg.replica_temperatures.T_min
        distances = analysis.compute_distance_timeseries_cpptraj(
            prmtop=prmtop, dcd=dcd,
            mask1=args.mask1, mask2=args.mask2,
            workdir=out_dir,
            cpptraj_path=cfg.cpptraj_path,
        )
        analysis.compute_distance_pmf(
            distances=distances,
            T_K=T,
            out_png=out_dir / "dist_pmf.png",
            out_xvg=out_dir / "dist_pmf.xvg",
        )
    return 0


def _cmd_not_yet_implemented(args: argparse.Namespace) -> int:
    print(
        f"Subcommand '{args.cmd}' is not yet implemented at this revision. "
        "It is wired up in subsequent steps of the v1.20.0 release.",
        file=sys.stderr,
    )
    return 2


# ---------------------------------------------------------------------------
# Parser / main
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="abmptools.genesis.grest",
        description=(
            "Build a GENESIS gREST_SSCR replica-exchange MD system "
            "(AMBER ff19SB + TIP3P)."
        ),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable DEBUG-level logging.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # example
    p_ex = sub.add_parser(
        "example",
        help="Print an example JSON config to stdout.",
    )
    p_ex.set_defaults(func=_cmd_example)

    # validate
    p_val = sub.add_parser(
        "validate",
        help=(
            "Validate config + check external tools (tleap / spdyn / "
            "atdyn / remd_convert / mpirun / cpptraj)."
        ),
    )
    p_val.add_argument(
        "--config", required=True, type=str,
        help="Path to config JSON.",
    )
    p_val.set_defaults(func=_cmd_validate)

    # build (scaffold; wired up in Step 7)
    p_bld = sub.add_parser(
        "build",
        help="Run the 5-stage build pipeline (Step 7).",
    )
    p_bld.add_argument("--config", required=True, type=str)
    p_bld.add_argument(
        "-o", "--output-dir", type=str, default=None,
        help="Override config.output_dir.",
    )
    p_bld.set_defaults(func=_cmd_build)

    # analyze (scaffold; wired up in Step 8)
    p_an = sub.add_parser(
        "analyze",
        help=(
            "Run remd_convert + replica transition + acceptance ratio "
            "+ 1D distance PMF (Step 8)."
        ),
    )
    p_an.add_argument("--config", required=True, type=str)
    p_an.add_argument("--run-dir", required=True, type=str)
    p_an.add_argument(
        "--what", nargs="+",
        default=["all"],
        choices=["sort", "transition", "acceptance", "pmf", "all"],
        help="Which analyses to run.",
    )
    p_an.add_argument(
        "--out-dir", type=str, default=None,
        help="Output directory for plots / xvg / csv.",
    )
    p_an.add_argument(
        "--mask1", type=str, default=None,
        help="GENESIS-style mask for PMF group1 (e.g. 'rno:96 and an:NZ').",
    )
    p_an.add_argument(
        "--mask2", type=str, default=None,
        help="GENESIS-style mask for PMF group2.",
    )
    p_an.add_argument(
        "--temperature", type=float, default=None,
        help="Override T in K for PMF (-kT log P); default: lowest replica T.",
    )
    p_an.set_defaults(func=_cmd_analyze)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    setup_logging(verbose=getattr(args, "verbose", False))
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
