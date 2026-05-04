# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide.cli
-------------------------
Command-line interface for the Martini 3 peptide CG builder.

Subcommands:
    build      run the full pipeline from a JSON (or YAML) config
    validate   check config + dependencies + Martini 3 .itp presence
    example    print a minimal JSON config to stdout
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path
from typing import List, Optional

from . import forcefield_check
from ._subprocess import setup_logging
from .builder import PeptideCGBuilder
from .models import PeptideBuildConfig, PeptideSpec

logger = logging.getLogger(__name__)


def _default_example() -> PeptideBuildConfig:
    return PeptideBuildConfig(
        peptides=[
            PeptideSpec(name="kgg", sequence="KGG", count=5),
        ],
        box_size_nm=8.0,
        martini_itp_dir="./ff",
    )


def _load_config(path: str) -> PeptideBuildConfig:
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
        peps = [PeptideSpec(**spec) for spec in raw.pop("peptides", [])]
        return PeptideBuildConfig(peptides=peps, **raw)
    return PeptideBuildConfig.from_json(str(p))


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def _cmd_example(args: argparse.Namespace) -> int:
    cfg = _default_example()
    print(json.dumps(asdict(cfg), indent=2, ensure_ascii=False))
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    n_species = len(cfg.peptides)
    n_total = sum(p.count for p in cfg.peptides)
    print(
        f"Configuration: OK ({n_species} peptide species, "
        f"{n_total} copies total, "
        f"box_size_nm={cfg.box_size_nm}, "
        f"T={cfg.temperature} K, NaCl {cfg.nacl_molar} M)\n"
    )

    ff_dir = args.ff_dir or cfg.martini_itp_dir or ""
    tools = forcefield_check.check_external_tools(
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

    builder = PeptideCGBuilder(cfg)
    result = builder.build()
    print(f"Build complete. Output: {result['output_dir']}")
    print(f"  final coordinates: {result['gro']}")
    print(f"  topology:          {result['top']}")
    print(f"  run script:        {result['run_script']}")
    return 0


# ---------------------------------------------------------------------------
# Parser / main
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="abmptools.cg.peptide",
        description="Build a Martini 3 peptide CG system for GROMACS MD.",
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
        help="Validate config + check dependencies and Martini 3 .itp files.",
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

    p_bld = sub.add_parser(
        "build",
        help="Run the full 6-stage build pipeline.",
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

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    setup_logging(verbose=getattr(args, "verbose", False))
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
