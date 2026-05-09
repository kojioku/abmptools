# -*- coding: utf-8 -*-
"""
abmptools.crystal.cli
----------------------
Command-line interface for the crystal-FMO subpackage.

Eight subcommands (Phase C-6 full set, mirrors mmgbsa):

    example       print a default YAML config snippet to stdout
    validate      check optional dependencies (ase / pyyaml / abinitmp)
    expand        Stage 1 only (CIF → cifout/layer<L>/{pdb,xyz}/)
    fragment      Stages 1+2 (PDB → for_abmp/*.{ajf,pdb})
    jobs          Stages 1+2+3 (+ HPC jobscripts + runbatch.sh)
    pipeline      all stages (with optional --run-local in Phase D)
    postproc      Stage 5 only (FMO logs → IFIE/PIEDA + nearest atoms)
    nearest       Standalone nearest-atom lookup on a single PDB

Invoke via:
    python -m abmptools.crystal <subcommand> --config crystal.yaml [options]
or the bundled console script:
    abmp-crystal <subcommand> ...
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
from .atom_distance import find_nearest_atoms
from .models import CrystalBuildConfig

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Default example config (also returned by `abmp-crystal example`)
# ---------------------------------------------------------------------------

_EXAMPLE_YAML = """\
# abmptools.crystal — example config
#
# Drop into a copy of `sample/crystal/csp7_smoke/` (provides the
# required UNK.ajf template) and run:
#
#   python -m abmptools.crystal pipeline --config crystal_config.yaml
#
project_name: csp7_layer5_around6
output_dir: ./out
inputs:
  - cif: XXXI-MMFF-R00001.cif
    layer: 5
    atoms_in_mol: [32]
cif_engine:
  engine: legacy   # 'legacy' | 'ase'
fragment:
  cutmode: around
  solutes: [0]
  criteria: 6.0
  molname: [UNK]
  pieda: true
  cmm: true
  template_ajf: ./UNK.ajf
fmo:
  method: MP2
  basis_set: 6-31Gdag
  memory: 6000
  is_xyz: true
hpc:
  scheduler: PJM
  queue: small
  group: hp190133
  nodes: 12
  proc_per_node: 2
  omp_threads: 24
  elapse: '24:00:00'
  abinit_dir: /data/hp190133/programs/ABINIT-MP
  binary_name: abinitmp_smp
postproc:
  enable: false
  frag_target: '1-10'
  annotate_nearest_atoms: true
"""


def _load_config(path: str) -> CrystalBuildConfig:
    """Dispatch on file extension; YAML for ``.yml/.yaml``, JSON otherwise."""
    suffix = Path(path).suffix.lower()
    if suffix in (".yml", ".yaml"):
        return CrystalBuildConfig.from_yaml(path)
    return CrystalBuildConfig.from_json(path)


def _build_orchestrator(args: argparse.Namespace):
    """Materialise the orchestrator from a CLI ``--config``."""
    from .builder import CrystalOrchestrator

    cfg = _load_config(args.config)
    if getattr(args, "output_dir", None):
        cfg.output_dir = args.output_dir
    if getattr(args, "engine", None):
        cfg.cif_engine.engine = args.engine
    if getattr(args, "run_local", False):
        cfg.run_local = True
    config_dir = str(Path(args.config).resolve().parent)
    return CrystalOrchestrator(cfg, config_dir=config_dir)


# ---------------------------------------------------------------------------
# Sub-command handlers
# ---------------------------------------------------------------------------

def _cmd_example(args: argparse.Namespace) -> int:
    sys.stdout.write(_EXAMPLE_YAML)
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    ok = forcefield_check.report(config_path=args.config)
    if args.config:
        try:
            cfg = _load_config(args.config)
            print(
                f"\nConfig valid: {len(cfg.inputs)} input(s), engine="
                f"{cfg.cif_engine.engine}, scheduler={cfg.hpc.scheduler}"
            )
        except Exception as exc:
            print(f"\nConfig parse error: {exc}", file=sys.stderr)
            ok = False
    return 0 if ok else 1


def _cmd_expand(args: argparse.Namespace) -> int:
    orch = _build_orchestrator(args)
    orch.expand_cif()
    print(f"Stage 1 complete. Output: {orch.output_dir}")
    return 0


def _cmd_fragment(args: argparse.Namespace) -> int:
    orch = _build_orchestrator(args)
    orch.expand_cif()
    orch.generate_fmo()
    print(f"Stages 1+2 complete. Output: {orch.output_dir}")
    return 0


def _cmd_jobs(args: argparse.Namespace) -> int:
    orch = _build_orchestrator(args)
    orch.expand_cif()
    orch.generate_fmo()
    orch.write_jobs()
    print(f"Stages 1+2+3 complete. Output: {orch.output_dir}")
    print(f"  jobscripts: {len(orch.job_scripts)}")
    return 0


def _cmd_pipeline(args: argparse.Namespace) -> int:
    orch = _build_orchestrator(args)
    summary = orch.run()
    print(f"Pipeline complete. Output: {summary['output_dir']}")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def _cmd_postproc(args: argparse.Namespace) -> int:
    """Stage 5 standalone — runs IFIE/PIEDA + (optional) nearest atoms."""
    from .postproc import run_postprocess

    orch = _build_orchestrator(args)
    spec = orch.config.postproc
    output_dir = orch.output_dir / spec.output_dir
    log_paths = [Path(p) for p in args.logs]
    pdb = Path(args.pdb) if args.pdb else None
    result = run_postprocess(
        spec=spec,
        output_dir=output_dir,
        log_files=log_paths,
        pdb_for_annotation=pdb,
        z_prime=args.z_prime,
    )
    print(f"Stage 5 complete. csvs={len(result['csvs'])}, "
          f"annotated={len(result['annotated'])}")
    for c in result["csvs"]:
        print(f"  {c}")
    return 0


def _cmd_nearest(args: argparse.Namespace) -> int:
    """Standalone: print the n nearest atoms to a residue centre."""
    rows = find_nearest_atoms(
        pdb_path=args.pdb,
        center_res_seq=args.center,
        n_neighbors=args.n,
    )
    for r in rows:
        print(
            f"  {r.serial:5d} {r.element:2s} ({r.name:<4s}) "
            f"res={r.res_seq} d={r.distance:.4f} Å"
        )
    return 0


# ---------------------------------------------------------------------------
# Parser / main
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="abmptools.crystal",
        description=(
            "Crystal-FMO workflow: CIF → supercell PDB → FMO inputs → "
            "(HPC|local) ABINIT-MP → IFIE/PIEDA. "
            "Phase C: stages 1-3 + postproc subcommand. "
            "Phase D adds --run-local & tutorials."
        ),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable DEBUG-level logging.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # example / validate
    p_ex = sub.add_parser("example",
                          help="Print a placeholder YAML config to stdout.")
    p_ex.set_defaults(func=_cmd_example)

    p_val = sub.add_parser("validate",
                           help="Report optional dependency status.")
    p_val.add_argument("--config", type=str, default=None)
    p_val.set_defaults(func=_cmd_validate)

    # Shared shape for stage runners that take --config.
    for cmd, helpmsg, fn in [
        ("expand",   "Stage 1 only: CIF -> cifout/layer<L>/{pdb,xyz}/",
         _cmd_expand),
        ("fragment", "Stages 1+2: PDB -> for_abmp/*.{ajf,pdb}",
         _cmd_fragment),
        ("jobs",     "Stages 1+2+3 (+ HPC jobscripts + runbatch.sh)",
         _cmd_jobs),
    ]:
        p = sub.add_parser(cmd, help=helpmsg)
        p.add_argument("--config", required=True, type=str)
        p.add_argument("-o", "--output-dir", type=str, default=None,
                       help="Override config.output_dir.")
        p.add_argument("--engine", choices=["legacy", "ase"], default=None,
                       help="Override config.cif_engine.engine.")
        p.set_defaults(func=fn)

    # pipeline (with --run-local Phase D placeholder).
    p_pl = sub.add_parser(
        "pipeline",
        help="All stages. --run-local invokes ABINIT-MP locally (Phase D).",
    )
    p_pl.add_argument("--config", required=True, type=str)
    p_pl.add_argument("-o", "--output-dir", type=str, default=None)
    p_pl.add_argument("--engine", choices=["legacy", "ase"], default=None)
    p_pl.add_argument("--run-local", action="store_true",
                      help="(Phase D) Invoke ABINIT-MP locally per AJF.")
    p_pl.set_defaults(func=_cmd_pipeline)

    # postproc (logs -> CSV + nearest atoms)
    p_pp = sub.add_parser(
        "postproc",
        help="Stage 5 only: FMO logs -> IFIE/PIEDA CSV + nearest atoms.",
    )
    p_pp.add_argument("--config", required=True, type=str)
    p_pp.add_argument("--logs", nargs="+", required=True,
                      help="ABINIT-MP log files to parse.")
    p_pp.add_argument("--pdb", type=str, default=None,
                      help="PDB for nearest-atom annotation.")
    p_pp.add_argument("--z-prime", type=int, default=1,
                      help="Z' value for the -zp flag of getifiepieda.")
    p_pp.set_defaults(func=_cmd_postproc)

    # nearest (standalone)
    p_nr = sub.add_parser("nearest",
                          help="Standalone nearest-atom lookup on a PDB.")
    p_nr.add_argument("--pdb", required=True, type=str)
    p_nr.add_argument("--center", required=True, type=int,
                      help="Centre residue sequence number.")
    p_nr.add_argument("-n", type=int, default=3,
                      help="Number of nearest atoms (default 3).")
    p_nr.set_defaults(func=_cmd_nearest)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    setup_logging(verbose=getattr(args, "verbose", False))
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
