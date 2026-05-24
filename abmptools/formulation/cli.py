# -*- coding: utf-8 -*-
"""argparse CLI for abmptools.formulation.

Subcommands: ``example`` / ``validate`` / ``build`` / ``analyze`` /
``release_us``.
"""
from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
from pathlib import Path
from typing import Optional, Sequence

from .builder import FormulationBuilder
from .models import (
    BileSaltSpec,
    EnhancerSpec,
    EquilibrationProtocol,
    FormulationBuildConfig,
    PeptideSpec,
    ProductionProtocol,
    SystemSpec,
    USProtocol,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# example
# ---------------------------------------------------------------------------


def _example_config() -> FormulationBuildConfig:
    """Return a minimal example config (KGGGG + caprate + taurocholate)."""
    return FormulationBuildConfig(
        system=SystemSpec(
            peptides=[PeptideSpec(
                name="kggggg", sequence="KGGGG", n_copies=2,
                cap_n="ACE", cap_c="NME",
            )],
            enhancers=[EnhancerSpec(
                name="Na-caprate", resname="CPR",
                smiles_neutral="CCCCCCCCCC(=O)O",
                smiles_charged="CCCCCCCCCC(=O)[O-]",
                n_neutral=16, n_charged=16,
            )],
            bile_salts=[BileSaltSpec(
                name="taurocholate", resname="TCH",
                smiles="C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@H]2[C@@H]3[C@@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC(O)[C@]12C",
                n_copies=2, net_charge=-1,
            )],
            box_size_nm=10.0,
            salt_concentration_M=0.15,
        ),
        equilibration=EquilibrationProtocol(),
        production=ProductionProtocol(),
        release_us=None,
    )


def _cmd_example(args: argparse.Namespace) -> int:
    cfg = _example_config()
    if args.output:
        cfg.to_json(args.output)
        print(f"wrote example config to {args.output}")
    else:
        from dataclasses import asdict
        from .models import _to_jsonable
        print(json.dumps(_to_jsonable(cfg), indent=2, ensure_ascii=False))
    return 0


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------


def _cmd_validate(args: argparse.Namespace) -> int:
    cfg = FormulationBuildConfig.from_json(args.config)
    print(
        f"Configuration: OK ("
        f"{len(cfg.system.peptides)} peptide species × "
        f"{cfg.system.total_peptide_copies} total copies, "
        f"{len(cfg.system.enhancers)} enhancer species, "
        f"{len(cfg.system.bile_salts)} bile salt species, "
        f"box={cfg.system.box_size_nm} nm, NaCl={cfg.system.salt_concentration_M} M)"
    )
    print("External tools:")
    tools = [
        ("tleap", cfg.tleap_path),
        ("acpype", cfg.acpype_path),
        ("packmol", cfg.packmol_path),
        ("gmx", cfg.gmx_path),
    ]
    for name, path in tools:
        located = shutil.which(path)
        mark = "[OK]" if located else "[NOT FOUND]"
        print(f"  {mark}  {name:10s}  ({located or path})")
    return 0


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------


def _cmd_build(args: argparse.Namespace) -> int:
    cfg = FormulationBuildConfig.from_json(args.config)
    if args.output_dir:
        cfg.output_dir = args.output_dir
    art = FormulationBuilder(cfg).build()
    print("Build complete. Output:", art.output_dir)
    print(f"  final coordinates: {art.gro}")
    print(f"  topology:          {art.top}")
    print(f"  index:             {art.ndx}")
    print(f"  run script:        {art.run_script}")
    print(f"  atoms (system):    {art.atoms_estimate}")
    print(f"  n_peptides total:  {art.n_peptides_total}")
    return 0


# ---------------------------------------------------------------------------
# analyze / release_us (stubs that import lazily)
# ---------------------------------------------------------------------------


def _cmd_analyze(args: argparse.Namespace) -> int:
    from .analysis import run_analysis  # lazy
    run_analysis(
        traj=args.traj, top=args.top, out_dir=args.out,
        enhancer_resnames=args.enhancer_resnames.split(",") if args.enhancer_resnames else (),
        bile_salt_resnames=args.bile_salt_resnames.split(",") if args.bile_salt_resnames else (),
    )
    return 0


def _cmd_release_us(args: argparse.Namespace) -> int:
    from .umbrella_release import run_release_us  # lazy
    cfg = FormulationBuildConfig.from_json(args.config)
    if args.output_dir:
        cfg.output_dir = args.output_dir
    run_release_us(
        config=cfg,
        target_peptide_idx=args.target_peptide_idx,
    )
    return 0


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="abmptools.formulation",
        description=(
            "AA peptide-formulation mixed-solution builder "
            "(Hossain 2023-style, AMBER ff14SB + GAFF2 + TIP3P)."
        ),
    )
    p.add_argument("-v", "--verbose", action="store_true")

    sub = p.add_subparsers(dest="cmd", required=True)

    p_ex = sub.add_parser("example",
                          help="Print an example JSON config to stdout (or file).")
    p_ex.add_argument("--output", "-o", default=None,
                      help="Write JSON to this path instead of stdout.")
    p_ex.set_defaults(func=_cmd_example)

    p_v = sub.add_parser("validate",
                         help="Validate a config JSON + check tool availability.")
    p_v.add_argument("--config", "-c", required=True)
    p_v.set_defaults(func=_cmd_validate)

    p_b = sub.add_parser("build",
                         help="Run the full 7-stage build pipeline.")
    p_b.add_argument("--config", "-c", required=True)
    p_b.add_argument("--output-dir", "-o", default=None,
                     help="Override config.output_dir.")
    p_b.set_defaults(func=_cmd_build)

    p_a = sub.add_parser("analyze", help="Run aggregate / contact / SS / SASA / hbond analysis.")
    p_a.add_argument("--traj", required=True)
    p_a.add_argument("--top", required=True)
    p_a.add_argument("--out", required=True)
    p_a.add_argument("--enhancer-resnames", default="")
    p_a.add_argument("--bile-salt-resnames", default="")
    p_a.set_defaults(func=_cmd_analyze)

    p_us = sub.add_parser("release_us", help="Set up + (later) run peptide-from-aggregate US.")
    p_us.add_argument("--config", "-c", required=True)
    p_us.add_argument("--output-dir", "-o", default=None)
    p_us.add_argument("--target-peptide-idx", type=int, default=0)
    p_us.set_defaults(func=_cmd_release_us)

    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
