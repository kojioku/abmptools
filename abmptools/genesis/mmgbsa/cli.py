# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa.cli
----------------------------
Command-line interface for the GENESIS MM/GBSA pipeline.

7 sub-commands (mode "A"):

    example       print a default JSON config to stdout
    validate      check config + external tools + Python deps
    divide        Stage 1 only (PDB splitter)
    parameterize  Stage 2 only (acpype + tleap)
    run           Stage 3 only (atdyn GBSA)
    analyze       Stage 4 only (CSV + plot from existing logs)
    pipeline      all 4 stages + folder-mode shortcut (-i / -r / -c)
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path
from typing import List, Optional

from . import analysis, forcefield_check
from ._subprocess import setup_logging
from .builder import MMGBSAOrchestrator, synthesize_targets_from_folder
from .models import (
    MMGBSABuildConfig,
    TargetSpec,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

def _default_example() -> MMGBSABuildConfig:
    """Default JSON example: 4 POC targets at ligand_resno=201."""
    return MMGBSABuildConfig(
        project_name="poc_reproduction",
        targets=[
            TargetSpec(pdb="3772L-rename.pdb", ligand_resno=201),
            TargetSpec(pdb="9MM52-rename.pdb", ligand_resno=201),
            TargetSpec(pdb="L7759-rename.pdb", ligand_resno=201),
            TargetSpec(pdb="M99KZ-rename.pdb", ligand_resno=201),
        ],
        input_dir="./input",
    )


def _load_config(path: str) -> MMGBSABuildConfig:
    return MMGBSABuildConfig.from_json(path)


# ---------------------------------------------------------------------------
# Sub-command handlers
# ---------------------------------------------------------------------------

def _cmd_example(args: argparse.Namespace) -> int:
    cfg = _default_example()
    print(json.dumps(asdict(cfg), indent=2, ensure_ascii=False))
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    ok = forcefield_check.report(cfg)
    return 0 if ok else 1


def _apply_overrides(cfg: MMGBSABuildConfig, args: argparse.Namespace) -> None:
    if getattr(args, "output_dir", None):
        cfg.output_dir = args.output_dir


def _cmd_divide(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    _apply_overrides(cfg, args)
    orch = MMGBSAOrchestrator(cfg)
    orch.divide()
    print(f"Stage 1 complete. Output: {orch.output_dir}")
    return 0


def _cmd_parameterize(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    _apply_overrides(cfg, args)
    orch = MMGBSAOrchestrator(cfg)
    orch.divide()
    orch.parameterize()
    print(f"Stage 2 complete. Output: {orch.output_dir}")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    _apply_overrides(cfg, args)
    orch = MMGBSAOrchestrator(cfg)
    orch.divide()
    orch.parameterize()
    orch.run_gbsa()
    print(f"Stage 3 complete. Output: {orch.output_dir}")
    return 0


def _cmd_analyze(args: argparse.Namespace) -> int:
    """Stage 4 only: re-aggregate logs from an existing run dir.

    Two modes:
    - ``--config`` + ``--run-dir``: discover per-target dirs from config's
      target list under run_dir.
    - ``--run-dir`` + ``--target-dirs``: explicit list of per-target dirs
      (POC compatibility -- replays against legacy outputs).
    """
    if args.target_dirs:
        target_dirs = [Path(d) for d in args.target_dirs]
    else:
        cfg = _load_config(args.config)
        run_dir = Path(args.run_dir).resolve()
        target_dirs = []
        for t in cfg.targets:
            tname = t.name or Path(t.pdb).stem
            tdir = run_dir / tname
            if tdir.is_dir():
                target_dirs.append(tdir)
            else:
                print(f"WARN: target dir missing: {tdir}", file=sys.stderr)

    if not target_dirs:
        print("No target directories found to analyze.", file=sys.stderr)
        return 1

    out_dir = Path(args.out_dir).resolve() if args.out_dir else (
        target_dirs[0].parent / "analysis"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    result = analysis.analyze(
        target_dirs=target_dirs,
        out_csv=out_dir / "analysis_results.csv",
        out_png=out_dir / "dg_bind_plot.png",
    )
    print(f"Stage 4 complete. CSV: {result.csv_path}, PNG: {result.png_path}")
    for t in result.targets:
        print(f"  {t.name}: ΔG_bind = {t.dg_bind:+.4f} kcal/mol")
    return 0


def _cmd_pipeline(args: argparse.Namespace) -> int:
    """All 4 stages. Supports folder-mode shortcut ``-i / -r / -c``."""
    if args.config:
        cfg = _load_config(args.config)
    else:
        # Folder-mode shortcut: synthesise config from CLI flags.
        if not (args.input_dir and args.ligand_resno):
            print(
                "pipeline requires either --config or both -i/--input-dir "
                "and -r/--ligand-resno.",
                file=sys.stderr,
            )
            return 1
        targets = synthesize_targets_from_folder(
            input_dir=Path(args.input_dir),
            ligand_resno=args.ligand_resno,
            chain=args.chain,
        )
        cfg = MMGBSABuildConfig(
            targets=targets,
            input_dir=str(args.input_dir),
            project_name=args.project_name or "mmgbsa_run",
        )

    _apply_overrides(cfg, args)
    orch = MMGBSAOrchestrator(cfg)
    result = orch.run()
    print(f"Pipeline complete. Output: {result['output_dir']}")
    print(f"  succeeded: {result['n_succeeded']} / {result['n_targets']}")
    if result["n_failed"]:
        print(f"  failed: {result['n_failed']} (see logger output)")
    if result["csv_path"] is not None:
        print(f"  CSV: {result['csv_path']}")
        print(f"  PNG: {result['png_path']}")
        for t in result["targets"]:
            print(f"    {t['name']}: ΔG_bind = {t['dg_bind']:+.4f} kcal/mol")
    return 0


# ---------------------------------------------------------------------------
# Parser / main
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="abmptools.genesis.mmgbsa",
        description=(
            "Build and analyse a GENESIS MM/GBSA single-point ΔG_bind "
            "calculation (AMBER ff14SB + GAFF/GAFF2)."
        ),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable DEBUG-level logging.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # example
    p_ex = sub.add_parser("example", help="Print an example JSON config to stdout.")
    p_ex.set_defaults(func=_cmd_example)

    # validate
    p_val = sub.add_parser(
        "validate",
        help="Validate config + check external tools + Python modules.",
    )
    p_val.add_argument("--config", required=True, type=str)
    p_val.set_defaults(func=_cmd_validate)

    # divide / parameterize / run -- shared shape.
    for cmd, helpmsg, fn in [
        ("divide",       "Stage 1 only: split PDB into receptor + ligand.",
         _cmd_divide),
        ("parameterize", "Stage 2 only: acpype + tleap (3 systems per target).",
         _cmd_parameterize),
        ("run",          "Stage 3 only: atdyn GBSA single-point on 3 systems.",
         _cmd_run),
    ]:
        p = sub.add_parser(cmd, help=helpmsg)
        p.add_argument("--config", required=True, type=str)
        p.add_argument(
            "-o", "--output-dir", type=str, default=None,
            help="Override config.output_dir.",
        )
        p.set_defaults(func=fn)

    # analyze
    p_an = sub.add_parser(
        "analyze",
        help="Stage 4: aggregate logs + ΔG_bind CSV + bar plot.",
    )
    p_an.add_argument("--config", type=str, default=None)
    p_an.add_argument("--run-dir", type=str, default=None)
    p_an.add_argument(
        "--target-dirs", nargs="+", default=None,
        help="Explicit per-target directories (overrides --run-dir).",
    )
    p_an.add_argument("--out-dir", type=str, default=None)
    p_an.set_defaults(func=_cmd_analyze)

    # pipeline (config or folder-mode shortcut)
    p_pl = sub.add_parser(
        "pipeline",
        help=(
            "All 4 stages. Use --config OR folder-mode shortcut "
            "(-i / -r [-c]) for POC compatibility."
        ),
    )
    p_pl.add_argument("--config", type=str, default=None)
    p_pl.add_argument("-i", "--input-dir", type=str, default=None)
    p_pl.add_argument("-r", "--ligand-resno", type=int, default=None)
    p_pl.add_argument("-c", "--chain", type=str, default=None)
    p_pl.add_argument("--project-name", type=str, default=None)
    p_pl.add_argument(
        "-o", "--output-dir", type=str, default=None,
        help="Override config.output_dir (or set output dir for folder mode).",
    )
    p_pl.set_defaults(func=_cmd_pipeline)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    setup_logging(verbose=getattr(args, "verbose", False))
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
