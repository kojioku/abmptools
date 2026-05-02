# -*- coding: utf-8 -*-
"""
abmptools.membrane.bilayer
--------------------------
Bilayer + peptide construction via *packmol-memgen* (AmberTools).

Why packmol-memgen
------------------
packmol-memgen (Schott-Verdugo & Gohlke, 2019) is the most convenient
license-OK route for assembling a bilayer + water + ions + peptide PDB
without CHARMM-GUI:

- Bundled with AmberTools — fully commercial OK
- Knows AMBER Lipid21 / Lipid17 residue conventions out-of-the-box
  (head-tail split for POPC/POPE/POPG/etc.)
- Optional ``--charmm`` flag converts atom names via the bundled
  ``charmmlipid2amber.py`` (Madej, also AmberTools — license OK)
- Handles peptide insertion via ``--solute <pep.pdb>``
- Salt + neutralisation via ``--saltcon`` / ``--solute_charge``

This module produces the **bilayer + peptide PDB only**. Force-field
assignment (AMBER tleap or CHARMM36 pdb2gmx) happens downstream in
:mod:`abmptools.membrane.parameterize_amber` /
:mod:`abmptools.membrane.parameterize_charmm`.

Environment requirement
-----------------------
``packmol-memgen`` requires ``AMBERHOME`` to be set so it can find
its bundled tleap / leap data files. We accept this either via the
already-set environment variable or via a *config*-derived fallback
(``Path(packmol_memgen_path).parents[1]`` for conda-installed copies).
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

from .models import MembraneConfig, PeptideSpec

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sequence → 3-letter conversion
# ---------------------------------------------------------------------------

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIE", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def one_to_three(seq: str) -> List[str]:
    """Convert a one-letter amino-acid code to a list of 3-letter codes.

    Histidine is mapped to ``HIE`` (epsilon-protonated, the AMBER default
    for neutral pH). Other tautomers must be specified by editing the
    resulting list.
    """
    out = []
    for ch in seq.upper():
        if ch not in _ONE_TO_THREE:
            raise ValueError(
                f"Unknown amino-acid one-letter code: {ch!r}. "
                f"Allowed: {sorted(_ONE_TO_THREE)}"
            )
        out.append(_ONE_TO_THREE[ch])
    return out


def estimate_peptide_charge(seq: str) -> int:
    """Estimate net charge of a peptide at neutral pH.

    Uses standard ionisation states:
      - ASP (D), GLU (E): -1
      - LYS (K), ARG (R): +1
      - HIS (H): 0 (assumed neutral HIE; user can override sequence-level)
      - All others (including caps): 0

    For non-standard protonation (HIP, CYM, TYM, …) supply ``pdb_path``
    instead of ``sequence`` and pass the explicit charge through tleap.
    """
    pos = sum(1 for c in seq.upper() if c in ("K", "R"))
    neg = sum(1 for c in seq.upper() if c in ("D", "E"))
    return pos - neg


# ---------------------------------------------------------------------------
# AMBERHOME resolution (packmol-memgen requirement)
# ---------------------------------------------------------------------------

def _resolve_amberhome(packmol_memgen_path: str) -> Optional[str]:
    """Resolve AMBERHOME for invoking packmol-memgen.

    Order of precedence:
      1. ``$AMBERHOME`` environment variable, if set.
      2. The conda env root containing the binary (``../..`` from the
         binary path, e.g. ``$ENV/bin/packmol-memgen`` → ``$ENV``).
      3. ``None`` (caller must surface the missing env).
    """
    env = os.environ.get("AMBERHOME")
    if env:
        return env
    resolved = shutil.which(packmol_memgen_path) or packmol_memgen_path
    p = Path(resolved).resolve()
    if p.is_file() and p.parent.name == "bin":
        candidate = p.parent.parent
        if (candidate / "dat" / "leap").is_dir():
            return str(candidate)
    return None


# ---------------------------------------------------------------------------
# Stage entry points
# ---------------------------------------------------------------------------

def write_peptide_pdb(
    *, config: MembraneConfig, out_path: str,
) -> str:
    """Build an initial peptide PDB from sequence (or copy ``pdb_path``).

    For *sequence* input, uses tleap with the configured protein force
    field to produce an extended-chain capped peptide. The output is
    written to *out_path* and returned.

    Parameters
    ----------
    config : MembraneConfig
        Build configuration (uses ``config.peptide`` and
        ``config.amber_protein_ff`` / ``config.tleap_path``).
    out_path : str
        Destination PDB path.

    Returns
    -------
    str
        Absolute path to the written PDB.
    """
    pep = config.peptide
    if pep is None:
        raise ValueError("config.peptide is None")

    out_abs = str(Path(out_path).resolve())

    if pep.pdb_path:
        shutil.copyfile(pep.pdb_path, out_abs)
        logger.info("peptide PDB copied from %s → %s", pep.pdb_path, out_abs)
        return out_abs

    work_dir = Path(out_abs).parent
    work_dir.mkdir(parents=True, exist_ok=True)
    tleap_in = work_dir / "build_peptide.tleap"
    tleap_in.write_text(_render_peptide_tleap(pep, config.amber_protein_ff,
                                              out_abs))

    env = os.environ.copy()
    amberhome = _resolve_amberhome(config.tleap_path)
    if amberhome:
        env["AMBERHOME"] = amberhome

    cmd = [config.tleap_path, "-f", str(tleap_in)]
    logger.info("running tleap to build peptide: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=work_dir, env=env,
                            capture_output=True, text=True)
    if result.returncode != 0 or not Path(out_abs).is_file():
        raise RuntimeError(
            f"tleap failed (rc={result.returncode}) when building peptide.\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}"
        )
    return out_abs


def _gcd_list(values: List[int]) -> int:
    from math import gcd
    out = values[0]
    for v in values[1:]:
        out = gcd(out, v)
    return out or 1


def estimate_distxy_angstrom(config: MembraneConfig,
                             apl_angstrom2: float = 65.0) -> float:
    """Estimate lipid-patch xy edge length (Å) from per-leaflet counts.

    Uses ``area = sum(n_per_leaflet * APL)``; default APL = 65 Å² is a
    middle-of-the-road value for liquid-disordered phospholipids near
    physiological temperature. Override via *apl_angstrom2* for
    saturated lipids (~50 Å²) or cholesterol-rich mixtures.
    """
    total_area = sum(l.n_per_leaflet * apl_angstrom2 for l in config.lipids)
    return total_area ** 0.5


def assemble_packmol_memgen_cmd(
    *, config: MembraneConfig, peptide_pdb: str, output_pdb: str,
    apl_angstrom2: float = 65.0,
) -> List[str]:
    """Compose the packmol-memgen argv from a config object.

    Pure function — no side effects. Useful for unit-testing the command
    construction without invoking the binary.
    """
    if not config.lipids:
        raise ValueError("config.lipids must not be empty")

    lipid_str = ":".join(l.resname for l in config.lipids)
    # packmol-memgen --ratio is a molar ratio. Reduce n_per_leaflet by
    # gcd so packmol-memgen can pack the actual count freely (the absolute
    # count is set by --distxy_fix * APL).
    counts = [l.n_per_leaflet for l in config.lipids]
    g = _gcd_list(counts)
    ratio_str = ":".join(str(n // g) for n in counts)

    distxy_a = estimate_distxy_angstrom(config, apl_angstrom2=apl_angstrom2)

    pep = config.peptide
    if pep is None:
        raise ValueError("config.peptide is None")
    solute_charge = (
        estimate_peptide_charge(pep.sequence) if pep.sequence else 0
    )

    cmd: List[str] = [
        config.packmol_memgen_path,
        "--lipids", lipid_str,
        "--ratio", ratio_str,
        "--solute", peptide_pdb,
        "--solute_con", str(pep.n_copies),
        "--solute_charge", str(solute_charge),
        "--output", output_pdb,
        "--ffprot", _strip_leaprc(config.amber_protein_ff),
        "--fflip", _fflip_from_leaprc(config.amber_lipid_ff),
        "--ffwat", _ffwat_from_leaprc(config.amber_water_ff),
        # Lipid-patch lateral size derived from target n_per_leaflet * APL.
        "--distxy_fix", f"{distxy_a:.2f}",
        # Buffer between lipid edge and box (xy) — peptide placed in water.
        "--dist", f"{config.distance_to_lipid_nm * 10.0:.2f}",
        # Water thickness above and below the membrane (z).
        "--dist_wat", f"{config.water_thickness_nm * 10.0:.2f}",
        "--salt",
        "--saltcon", str(config.ions.salt_concentration_M),
        "--salt_c", config.ions.cation,
        "--salt_a", config.ions.anion,
        "--overwrite",
        "--noprogress",
    ]

    if config.backend == "charmm36":
        # Bundle's charmmlipid2amber.py runs in reverse to convert lipid
        # residue / atom names to CHARMM convention (POPC → split into
        # CHARMM POPC residue with atom names like P, O11, etc.).
        # Note: the solute peptide is *passed through* unchanged by
        # packmol-memgen, so we still need to translate it ourselves —
        # see parameterize_charmm._translate_pdb_amber_to_charmm.
        cmd.append("--charmm")

    if config.seed is not None:
        cmd += ["--leapline", f"# seed={config.seed}"]

    return cmd


def build_bilayer(*, config: MembraneConfig, build_dir: str) -> str:
    """Build a bilayer (+ peptide) PDB and return its absolute path.

    Parameters
    ----------
    config : MembraneConfig
        Build configuration.
    build_dir : str
        Directory in which to invoke packmol-memgen.

    Returns
    -------
    str
        Absolute path to the assembled PDB
        (``{build_dir}/bilayer_peptide.pdb``).
    """
    bd = Path(build_dir).resolve()
    bd.mkdir(parents=True, exist_ok=True)

    peptide_pdb = str(bd / "peptide.pdb")
    write_peptide_pdb(config=config, out_path=peptide_pdb)

    output_pdb = str(bd / "bilayer_peptide.pdb")
    cmd = assemble_packmol_memgen_cmd(
        config=config, peptide_pdb=peptide_pdb, output_pdb=output_pdb,
    )

    env = os.environ.copy()
    amberhome = _resolve_amberhome(config.packmol_memgen_path)
    if amberhome:
        env["AMBERHOME"] = amberhome
    else:
        logger.warning(
            "AMBERHOME could not be resolved; packmol-memgen may fail. "
            "Set $AMBERHOME or use a conda-installed packmol-memgen."
        )

    logger.info("running packmol-memgen: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=bd, env=env,
                            capture_output=True, text=True)
    log_path = bd / "packmol_memgen.log"
    log_path.write_text(
        f"=== argv ===\n{' '.join(cmd)}\n\n"
        f"=== stdout ===\n{result.stdout}\n"
        f"=== stderr ===\n{result.stderr}\n"
        f"=== returncode ===\n{result.returncode}\n"
    )
    if result.returncode != 0 or not Path(output_pdb).is_file():
        raise RuntimeError(
            f"packmol-memgen failed (rc={result.returncode}). "
            f"See {log_path} for details."
        )
    logger.info("bilayer PDB written: %s", output_pdb)
    return output_pdb


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _render_peptide_tleap(
    pep: PeptideSpec, leaprc: str, output_pdb: str,
) -> str:
    """Render a tleap input that builds a capped extended peptide chain."""
    res3 = one_to_three(pep.sequence)
    if pep.cap_n:
        res3 = [pep.cap_n] + res3
    if pep.cap_c:
        res3 = res3 + [pep.cap_c]
    seq_block = " ".join(res3)
    return (
        f"# Auto-generated by abmptools.membrane.bilayer\n"
        f"source {leaprc}\n"
        f"peptide = sequence {{ {seq_block} }}\n"
        f"savepdb peptide {output_pdb}\n"
        f"quit\n"
    )


def _strip_leaprc(leaprc_str: str) -> str:
    """Map ``leaprc.protein.ff19SB`` → ``ff19SB`` for packmol-memgen --ffprot."""
    return leaprc_str.split(".")[-1]


def _fflip_from_leaprc(leaprc_str: str) -> str:
    """Map ``leaprc.lipid21`` → ``lipid21`` for packmol-memgen --fflip."""
    return leaprc_str.split(".")[-1]


def _ffwat_from_leaprc(leaprc_str: str) -> str:
    """Map ``leaprc.water.tip3p`` → ``tip3p`` for packmol-memgen --ffwat."""
    return leaprc_str.split(".")[-1]
