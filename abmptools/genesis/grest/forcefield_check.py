# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest.forcefield_check
----------------------------------------
External tool resolution for the gREST_SSCR pipeline.

Required tools (subprocess only, never bundled):

- ``tleap``        -- AmberTools, AMBER ``prmtop`` + ``coor`` generation.
- ``spdyn``        -- GENESIS, equilibration + production replica MD.
- ``atdyn``        -- GENESIS, minimisation.
- ``remd_convert`` -- GENESIS, post-MD parameter sort.
- ``mpirun``       -- OpenMPI / MPICH, parallel replica execution.

Optional:

- ``cpptraj``      -- AmberTools, ``around``-mode REST residue resolution.

The :func:`report` function prints a human-readable status block;
returns ``True`` only if all required tools resolve and the config
itself is buildable.

GENESIS build instructions (LGPL-3.0-or-later):
https://github.com/genesis-release-r-ccs/genesis

POC-known build caveat (icx removes ftello64/fseeko64): patch
``fileio_data_.c`` to use ``ftello`` / ``fseeko`` (drop the ``64``
suffix); then ``./configure CC=icx && make``.
"""
from __future__ import annotations

import logging
import shutil
from typing import List

# Re-use ToolStatus from cg.peptide to avoid duplicating the dataclass.
from abmptools.cg.peptide.forcefield_check import ToolStatus

from .models import GrestBuildConfig
from .replica_temperatures import generate_ladder, ladder_ratios
from .rest_selection import (
    RESTSelectionResult,
    parse_explicit_residues,
    format_genesis_selection,
)

logger = logging.getLogger(__name__)


GENESIS_REPO_URL = "https://github.com/genesis-release-r-ccs/genesis"
GENESIS_LICENSE = "LGPL-3.0-or-later"

ICX_PATCH_HINT = (
    "POC build caveat: under icx, fileio_data_.c uses ftello64 / fseeko64\n"
    "  which clang doesn't auto-import. Patch with sed:\n"
    "    sed -i 's/ftello64/ftello/g; s/fseeko64/fseeko/g' \\\n"
    "        analysis/src/lib/fileio_data_.c\n"
    "  then: ./configure CC=icx && make"
)


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

def check_external_tools(
    tleap: str = "tleap",
    spdyn: str = "spdyn",
    atdyn: str = "atdyn",
    remd_convert: str = "remd_convert",
    mpirun: str = "mpirun",
    cpptraj: str = "cpptraj",
) -> List[ToolStatus]:
    """Resolve each required + optional external tool to an absolute path.

    The order in the returned list matches the printout order in
    :func:`report`.
    """
    spec = [
        ("tleap", tleap,
         "AmberTools -- AMBER prmtop + coor generation",
         True),
        ("spdyn", spdyn,
         f"GENESIS -- replica MD ({GENESIS_LICENSE})",
         True),
        ("atdyn", atdyn,
         f"GENESIS -- minimisation ({GENESIS_LICENSE})",
         True),
        ("remd_convert", remd_convert,
         "GENESIS -- post-MD parameter sort",
         True),
        ("mpirun", mpirun,
         "OpenMPI / MPICH -- parallel replica execution",
         True),
        ("cpptraj", cpptraj,
         "AmberTools -- around-mode REST residue resolution (optional)",
         False),
    ]
    out: List[ToolStatus] = []
    for label, exe, purpose, required in spec:
        resolved = shutil.which(exe)
        out.append(ToolStatus(
            name=label, purpose=purpose, found=resolved is not None,
            path=resolved, required=required,
        ))
    return out


# ---------------------------------------------------------------------------
# Pre-build sanity checks
# ---------------------------------------------------------------------------

def preview_rest_selection(cfg: GrestBuildConfig) -> RESTSelectionResult:
    """Resolve the REST selection in pure Python where possible.

    For ``mode="explicit"`` this returns the parsed residue list and
    string. For ``mode="around"`` we cannot probe the prmtop yet (it
    is built later by :mod:`system_builder`), so we return an empty
    placeholder result -- ``report`` displays the centre + radius
    instead.
    """
    spec = cfg.rest_selection
    if spec.mode == "explicit":
        residues = parse_explicit_residues(spec.residues)
        return RESTSelectionResult(
            residues=residues,
            n_residues=len(residues),
            selection_string=format_genesis_selection(residues),
        )
    return RESTSelectionResult(
        residues=[], n_residues=0,
        selection_string=f"(around {spec.center} radius {spec.radius_A:.2f} Å)",
    )


# ---------------------------------------------------------------------------
# Human-readable report
# ---------------------------------------------------------------------------

def report(cfg: GrestBuildConfig) -> bool:
    """Print a status block for :class:`GrestBuildConfig`.

    Returns ``True`` if the configuration is buildable: all REQUIRED
    tools are found AND the cross-field invariants in ``__post_init__``
    have already passed (the caller is expected to have constructed
    the config successfully). On missing tools, falsy return.
    """
    # Header.
    rt = cfg.replica_temperatures
    if rt.mode == "manual":
        n_rep = len(rt.temperatures)
        ladder = generate_ladder(rt)
    else:
        n_rep = rt.n_replicas
        ladder = generate_ladder(rt)
    sel = preview_rest_selection(cfg)
    print(
        f"Configuration: OK (project={cfg.project_name}, "
        f"AMBER {cfg.ff_protein.split('.')[-1]} + "
        f"{cfg.ff_water.split('.')[-1]}, "
        f"{n_rep} replicas, REST={cfg.rest_selection.mode} "
        f"({sel.n_residues} residues), "
        f"T=[{ladder[0]:.2f}, ..., {ladder[-1]:.2f}] K)"
    )

    # Tools.
    tools = check_external_tools(
        tleap=cfg.tleap_path,
        spdyn=cfg.spdyn_path,
        atdyn=cfg.atdyn_path,
        remd_convert=cfg.remd_convert_path,
        mpirun=cfg.mpirun_path,
        cpptraj=cfg.cpptraj_path,
    )
    print("\nExternal tools:")
    required_ok = True
    for t in tools:
        if t.found:
            print(f"  [OK]    {t.name:<13} ({t.path})")
        elif t.required:
            required_ok = False
            print(f"  [MISS]  {t.name:<13} (required) -- {t.purpose}")
        else:
            print(f"  [WARN]  {t.name:<13} -- {t.purpose}")

    # REST preview.
    print("\nREST selection preview:")
    print(f"  mode       = {cfg.rest_selection.mode}")
    if cfg.rest_selection.mode == "explicit":
        print(f"  selection  = {sel.selection_string}")
        print(f"  n_residues = {sel.n_residues}")
    else:
        print(
            f"  centre     = {cfg.rest_selection.center}, "
            f"radius     = {cfg.rest_selection.radius_A:.2f} Å"
        )
        print("  (resolved via cpptraj at build time)")
    print(f"  param_type = {' '.join(cfg.rest_selection.param_types)}  "
          f"(GENESIS Setup_Remd_Solute_Tempering tokens)")

    # Replica ladder preview + ratios.
    print(f"\nReplica temperature ladder ({rt.mode}/{rt.method}):")
    ratios = ladder_ratios(ladder)
    for i, t in enumerate(ladder, start=1):
        if i == 1:
            print(f"  rep{i:02d} -> {t:>8.3f} K     ratio T_i/T_0 = 1.000")
        else:
            print(
                f"  rep{i:02d} -> {t:>8.3f} K     ratio T_i/T_0 = "
                f"{t/ladder[0]:.4f}  "
                f"adj T_i/T_{i-2} = {ratios[i-2]:.4f}"
            )

    # MPI mapping preview.
    n_total = n_rep * cfg.mpi_processes_per_replica
    n_cores = n_total * cfg.omp_num_threads
    print("\nMPI mapping preview:")
    print(
        f"  total processes = {n_rep} replicas x "
        f"{cfg.mpi_processes_per_replica} MPI/replica = {n_total}"
    )
    print(f"  threads/proc    = {cfg.omp_num_threads} (OMP_NUM_THREADS)")
    print(f"  total cores     = {n_cores}")

    # Build hints if required tools missing.
    if not required_ok:
        print(f"\nGENESIS build instructions: {GENESIS_REPO_URL}")
        print(ICX_PATCH_HINT)

    return required_ok


__all__ = [
    "GENESIS_REPO_URL",
    "GENESIS_LICENSE",
    "ICX_PATCH_HINT",
    "check_external_tools",
    "preview_rest_selection",
    "report",
]
