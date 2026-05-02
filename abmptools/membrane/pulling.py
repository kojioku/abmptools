# -*- coding: utf-8 -*-
"""
abmptools.membrane.pulling
--------------------------
Reaction-coordinate generation: peptide pulled along z through the bilayer.

The output of stage 4 is a single MD trajectory (``pull.xtc`` + ``pullx.xvg``)
plus snapshots at uniform z-spacing that seed each US window.

Why pulling
-----------
Each US window needs an initial structure where the peptide centre-of-mass
sits near the window's reference z. The cheapest way to get those is a
single steered-MD run with constant-velocity pulling:

    pull-coord1-init  = <z at start of pulling>
    pull-coord1-rate  = <nm/ps> (signed; negative pulls peptide toward bilayer)
    pull-coord1-k     = <large k>

Then we extract frames at z = z_min, z_min + dz, … z_max from the
trajectory. The pulling MDP shares the NPT-semiisotropic chassis from
:mod:`mdp_us_protocol` and adds a pull-code block.
"""
from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from .models import MembraneConfig
from .mdp_us_protocol import render_npt_mdp, render_pull_block
from .bilayer import _resolve_amberhome

logger = logging.getLogger(__name__)


def write_pulling_mdp(
    *, config: MembraneConfig, pull_dir: str,
    pull_init_nm: float,
) -> str:
    """Write pull.mdp combining NPT-semiisotropic + pull-code block.

    Parameters
    ----------
    config : MembraneConfig
        Build configuration.
    pull_dir : str
        Output directory for pull.mdp.
    pull_init_nm : float
        Initial value of the pull coordinate at the start of the pulling
        run, in nm. Should match the peptide-bilayer COM-z separation at
        the end of the npt-equilibration step (use
        :func:`estimate_initial_pull_coord` to compute this from the
        equilibrated .gro).
    """
    pd = Path(pull_dir).resolve()
    pd.mkdir(parents=True, exist_ok=True)

    # Override the NPT body's nsteps with the pulling-stage nsteps.
    p = config.pulling
    npt_body = render_npt_mdp(config=config, pcoupl="c-rescale")
    npt_body = _override_mdp_field(npt_body, "nsteps", str(p.nsteps))
    npt_body = _override_mdp_field(npt_body, "dt", f"{config.equilibration.dt_ps:.3f}")
    npt_body = _override_mdp_field(npt_body, "nstxout-compressed",
                                   str(p.nstxout_compressed))

    pull_body = render_pull_block(
        config=config,
        pull_init_nm=pull_init_nm,
        pull_rate_nm_per_ps=p.pull_rate_nm_per_ps,
        pull_k=p.pull_force_constant,
        nst_pull_xf=max(50, p.nstxout_compressed // 10),
    )

    out = (
        "; pull.mdp — constant-velocity steered MD\n"
        "; (NPT-semiisotropic chassis + pull-code block)\n"
        + npt_body.replace("; npt.mdp — semiisotropic NPT\n", "", 1)
        + "\n"
        + pull_body
    )
    out_path = pd / "pull.mdp"
    out_path.write_text(out)
    return str(out_path)


def estimate_initial_pull_coord(
    *, gro_path: str, ndx_path: str,
    pull_group1: str, pull_group2: str,
) -> float:
    """Compute peptide_COM_z - bilayer_COM_z from a .gro + .ndx.

    Uses *unweighted* mean (atom count) — an approximation to the
    mass-weighted COM that GROMACS uses internally; usually within
    ~0.1 nm for realistic systems and good enough as the *starting*
    value for the pulling run (the actual coord is then evolved by
    GROMACS using the proper mass-weighted COM).
    """
    pos = _read_gro_positions(gro_path)
    groups = _read_ndx_groups(ndx_path)
    if pull_group1 not in groups:
        raise KeyError(
            f"index group {pull_group1!r} not found in {ndx_path}; "
            f"available: {sorted(groups)}"
        )
    if pull_group2 not in groups:
        raise KeyError(
            f"index group {pull_group2!r} not found in {ndx_path}; "
            f"available: {sorted(groups)}"
        )
    z1 = sum(pos[i - 1][2] for i in groups[pull_group1]) / len(groups[pull_group1])
    z2 = sum(pos[i - 1][2] for i in groups[pull_group2]) / len(groups[pull_group2])
    return z2 - z1


def extract_window_frames(
    *, pull_xtc: str, pull_tpr: str, pullx_xvg: str,
    config: MembraneConfig, out_dir: str,
) -> Dict[int, str]:
    """Extract per-window starting GRO files from the pulling trajectory.

    For each window *i* (i ∈ [0, n_windows)):
      1. Compute target z = z_min + i * dz.
      2. Find the time in *pullx_xvg* where the pull-coord is closest to
         the target.
      3. Run ``gmx trjconv -dump <time>`` to write
         ``<out_dir>/win_{i:03d}/start.gro``.

    Returns a dict mapping window index → starting .gro path.
    """
    od = Path(out_dir).resolve()
    od.mkdir(parents=True, exist_ok=True)
    times, zs = parse_pullx_xvg(pullx_xvg)
    if not times:
        raise RuntimeError(f"no data rows in {pullx_xvg}")

    out: Dict[int, str] = {}
    for i in range(config.umbrella.n_windows):
        z_target = (
            config.umbrella.z_min_nm + i * config.umbrella.window_spacing_nm
        )
        # find time with z closest to target
        best_idx = min(range(len(zs)), key=lambda k: abs(zs[k] - z_target))
        t_dump = times[best_idx]
        win_dir = od / f"win_{i:03d}"
        win_dir.mkdir(parents=True, exist_ok=True)
        out_gro = str(win_dir / "start.gro")
        _gmx_trjconv_dump(
            tpr=pull_tpr, xtc=pull_xtc, t_ps=t_dump, out_gro=out_gro,
            gmx_path=config.gmx_path,
        )
        out[i] = out_gro
        logger.info(
            "window %3d  z_target=%+5.2f  matched z=%+5.2f at t=%.1f ps",
            i, z_target, zs[best_idx], t_dump,
        )
    return out


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def parse_pullx_xvg(path: str) -> Tuple[List[float], List[float]]:
    """Parse a GROMACS pullx.xvg → (times[ps], pull_coord_values[nm]).

    Skips lines starting with ``@`` or ``#``. For multi-coord pull the
    *first* data column after time is returned (i.e. coord 1).
    """
    times: List[float] = []
    zs: List[float] = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("@", "#")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                times.append(float(parts[0]))
                zs.append(float(parts[1]))
            except ValueError:
                continue
    return times, zs


def _read_gro_positions(path: str) -> List[Tuple[float, float, float]]:
    """Return list of (x, y, z) in nm for each atom (1-based aligned)."""
    out: List[Tuple[float, float, float]] = []
    with open(path) as f:
        f.readline()
        n_atoms = int(f.readline())
        for _ in range(n_atoms):
            line = f.readline()
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            out.append((x, y, z))
    return out


def _read_ndx_groups(path: str) -> Dict[str, List[int]]:
    """Parse a .ndx file → {group_name: [atom_indices_1based]}."""
    groups: Dict[str, List[int]] = {}
    current: str = ""
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("[") and s.endswith("]"):
                current = s[1:-1].strip()
                groups.setdefault(current, [])
                continue
            if current:
                for tok in s.split():
                    try:
                        groups[current].append(int(tok))
                    except ValueError:
                        pass
    return groups


def _gmx_trjconv_dump(
    *, tpr: str, xtc: str, t_ps: float, out_gro: str, gmx_path: str = "gmx",
) -> None:
    """Run ``gmx trjconv -dump <t> -o <out_gro>`` (selecting System group)."""
    env = os.environ.copy()
    amberhome = _resolve_amberhome(gmx_path)
    if amberhome:
        env.setdefault("AMBERHOME", amberhome)
    cmd = [
        gmx_path, "trjconv",
        "-s", tpr, "-f", xtc,
        "-dump", f"{t_ps:.3f}",
        "-o", out_gro,
    ]
    # gmx trjconv asks for a group selection from stdin; "0" = System
    result = subprocess.run(cmd, input="0\n", capture_output=True,
                            text=True, env=env)
    if result.returncode != 0 or not Path(out_gro).is_file():
        raise RuntimeError(
            f"gmx trjconv failed (rc={result.returncode}) for t={t_ps}.\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}"
        )


def _override_mdp_field(mdp_text: str, field: str, value: str) -> str:
    """Replace a ``field = value`` line in MDP text (matches at line start)."""
    out_lines = []
    replaced = False
    for line in mdp_text.splitlines():
        # match e.g. "nsteps          = 5000"
        stripped = line.lstrip()
        if stripped.startswith(f"{field} ") or stripped.startswith(f"{field}="):
            indent = line[:len(line) - len(stripped)]
            new_line = f"{indent}{field:<15s} = {value}"
            out_lines.append(new_line)
            replaced = True
        else:
            out_lines.append(line)
    if not replaced:
        out_lines.append(f"{field} = {value}")
    return "\n".join(out_lines) + ("\n" if mdp_text.endswith("\n") else "")


# ---------------------------------------------------------------------------
# CLI entry point (used by run.sh)
# ---------------------------------------------------------------------------

def _main() -> None:
    """``python -m abmptools.membrane.pulling`` — extract window frames."""
    import argparse
    from .models import MembraneConfig

    parser = argparse.ArgumentParser(
        description="Extract per-window starting frames from a pulling run.",
    )
    parser.add_argument("--pull-tpr", required=True)
    parser.add_argument("--pull-xtc", required=True)
    parser.add_argument("--pull-xvg", required=True)
    parser.add_argument("--config", required=True,
                        help="Path to MembraneConfig JSON.")
    parser.add_argument("--windows-dir", required=True)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(name)s | %(message)s")
    cfg = MembraneConfig.from_json(args.config)
    extract_window_frames(
        pull_tpr=args.pull_tpr, pull_xtc=args.pull_xtc,
        pullx_xvg=args.pull_xvg,
        config=cfg, out_dir=args.windows_dir,
    )


if __name__ == "__main__":
    _main()
