# -*- coding: utf-8 -*-
"""
abmptools.membrane.parameterize_charmm
--------------------------------------
CHARMM36 backend: CHARMM36m + CHARMM36 lipids → GROMACS top/gro.

Pipeline
--------
1. Stage the CHARMM36 force-field directory into the build dir
   (``charmm36-jul2022.ff/`` etc.). Source: Klauda lab's GROMACS port,
   specified by :attr:`MembraneConfig.charmm_ff_dir`. **Not downloaded
   by this package** — must already be present (license/network policy
   left to the user).

2. Translate the AMBER-named PDB produced by tleap +
   packmol-memgen-with-``--charmm`` into a fully-CHARMM-named PDB.
   Lipids and water/ions are already CHARMM-named by packmol-memgen's
   ``--charmm`` flag (uses bundled charmmlipid2amber.py in reverse);
   only the *peptide solute* needs a name pass.

3. Run ``gmx pdb2gmx -f bilayer_peptide.pdb -ff charmm36-jul2022
   -water tip3p -ter`` to obtain system.top + system.gro.

4. Generate a GROMACS index (.ndx) using
   :func:`parameterize_amber.write_index_from_gro` (the residue-class
   tables already include CHARMM aliases such as HSD/HSE/HSP, SOD/CLA,
   TIP3).

License
-------
Uses only **CHARMM36 force-field parameter values** (free for academic
and industrial use per MacKerell lab's stated license). Does **NOT**:

- call CGenFF Web server (commercial use forbidden by Silcsbio license)
- consume CHARMM-GUI generated files (commercial use requires subscription)

Therefore peptide + standard lipid systems are fully commercial-OK on
this route, provided no novel small molecules requiring CGenFF appear.
For novel small molecules, switch to backend="amber" and use GAFF2.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Tuple

from .models import IonSpec, MembraneConfig
from .bilayer import _resolve_amberhome
from .parameterize_amber import write_index_from_gro

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# AMBER → CHARMM translation tables
# ---------------------------------------------------------------------------

#: Residue-name map (AMBER → CHARMM36).
#:
#: Standard 20 amino acids share names between AMBER and CHARMM36
#: (ALA, GLY, LYS, ...) so they don't appear here. Only protomers,
#: caps, and solvent / ions differ.
AMBER_TO_CHARMM_RESNAME: Dict[str, str] = {
    # histidine tautomers
    "HID": "HSD",
    "HIE": "HSE",
    "HIP": "HSP",
    # protonation states
    "ASH": "ASPP",
    "GLH": "GLUP",
    "LYN": "LSN",
    # disulfide / thiolate cysteines
    "CYX": "CYS",   # disulfide handled via pdb2gmx LINK / -ss
    "CYM": "CYS",   # rare; pdb2gmx may need patch
    # caps
    "NME": "CT3",
    # water
    "WAT": "TIP3",
    "HOH": "TIP3",
    # ions (also handled via translate_ion_name())
    "Na+": "SOD",
    "Cl-": "CLA",
    "K+":  "POT",
    "Mg2+": "MG",
    "Ca2+": "CAL",
}

#: Atom-name overrides per residue (AMBER → CHARMM).
#:
#: Outer key is the **CHARMM** residue name (post-rename), inner key is
#: the AMBER atom name, value is the CHARMM atom name. Atoms not listed
#: pass through unchanged. Most standard residues only need ``H → HN``
#: which is handled by :data:`UNIVERSAL_BACKBONE_ATOM_MAP`.
PER_RESIDUE_ATOM_MAP: Dict[str, Dict[str, str]] = {
    "ACE": {       # acetyl cap (residue name unchanged)
        "H1":  "HY1",
        "H2":  "HY2",
        "H3":  "HY3",
        "CH3": "CAY",
        "C":   "CY",
        "O":   "OY",
    },
    "CT3": {       # N-methylamide cap (was AMBER NME)
        "CH3":  "CAT",
        "HH31": "HT1",
        "HH32": "HT2",
        "HH33": "HT3",
    },
    # Monatomic ions: in CHARMM the atom name matches the residue name.
    "SOD": {"Na+": "SOD", "NA+": "SOD", "NA": "SOD"},
    "CLA": {"Cl-": "CLA", "CL-": "CLA", "CL": "CLA"},
    "POT": {"K+":  "POT", "K":   "POT"},
    "MG":  {"Mg2+": "MG", "MG2+": "MG"},
    "CAL": {"Ca2+": "CAL", "CA2+": "CAL"},
    # Water: AMBER WAT (O, H1, H2) → CHARMM TIP3 (OH2, H1, H2).
    "TIP3": {"O": "OH2"},
}

#: Universal backbone atom replacement applied to *all* protein residues
#: (caps excluded by virtue of having different residue names).
#: AMBER's amide hydrogen ``H`` becomes ``HN`` in CHARMM.
UNIVERSAL_BACKBONE_ATOM_MAP: Dict[str, str] = {
    "H": "HN",
}

#: Residues for which :data:`UNIVERSAL_BACKBONE_ATOM_MAP` does NOT apply
#: (caps with their own atom-name conventions).
SKIP_BACKBONE_RENAME: frozenset = frozenset({"ACE", "CT3", "TIP3", "SOD",
                                             "CLA", "POT", "MG", "CAL"})


# ---------------------------------------------------------------------------
# Public translation helpers
# ---------------------------------------------------------------------------

def translate_residue_name(amber_resname: str) -> str:
    """Map an AMBER residue name to its CHARMM36 equivalent.

    Returns the input unchanged when no mapping is defined (most standard
    amino acids share names; e.g., ``"ALA" -> "ALA"``).
    """
    return AMBER_TO_CHARMM_RESNAME.get(amber_resname.strip(),
                                       amber_resname.strip())


def translate_atom_name(charmm_resname: str, amber_atom: str) -> str:
    """Map an AMBER atom name to its CHARMM36 equivalent.

    *charmm_resname* is the **post-rename** residue name (e.g. ``CT3``,
    not ``NME``) so caps are resolved correctly. Returns the input
    unchanged when no mapping applies.
    """
    aa = amber_atom.strip()
    rr = charmm_resname.strip()
    # 1) per-residue overrides take precedence
    per = PER_RESIDUE_ATOM_MAP.get(rr, {})
    if aa in per:
        return per[aa]
    # 2) universal backbone H → HN (skip for caps / non-protein)
    if rr not in SKIP_BACKBONE_RENAME and aa in UNIVERSAL_BACKBONE_ATOM_MAP:
        return UNIVERSAL_BACKBONE_ATOM_MAP[aa]
    return aa


def translate_ion_name(amber_ion: str) -> str:
    """Map an AMBER ion name (``Na+``, ``Cl-``, ...) to CHARMM36 (``SOD``, ``CLA``, ...).

    If the input doesn't match any AMBER convention it is returned
    unchanged, allowing users to pass CHARMM names verbatim
    (e.g. ``IonSpec(cation='SOD')`` is respected as-is).
    """
    return AMBER_TO_CHARMM_RESNAME.get(amber_ion.strip(), amber_ion.strip())


def translate_ions_for_charmm(ions: IonSpec) -> IonSpec:
    """Return a copy of *ions* with cation/anion translated to CHARMM names."""
    return IonSpec(
        cation=translate_ion_name(ions.cation),
        anion=translate_ion_name(ions.anion),
        salt_concentration_M=ions.salt_concentration_M,
        neutralize=ions.neutralize,
    )


# ---------------------------------------------------------------------------
# PDB translator (whole-file pass)
# ---------------------------------------------------------------------------

def translate_pdb_amber_to_charmm(*, input_pdb: str, output_pdb: str) -> str:
    """Rewrite *input_pdb* with CHARMM36 residue / atom names → *output_pdb*.

    Handles ATOM and HETATM records. Other lines (HEADER, CONECT, END,
    REMARK, ...) are passed through unchanged. PDB columns are preserved
    using fixed-column substitution to avoid breaking downstream parsers.
    """
    out_lines = []
    with open(input_pdb) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                # Special case: residue names like "Na+" / "Cl-" only fit
                # the 4-character residue field if present; AMBER uses
                # 4-char "Na+ " / "Cl- " with the space.
                resname_4 = line[17:21].strip()
                # Translate residue name (try 4-char first, fall back to 3-char)
                new_res = translate_residue_name(resname_4) \
                          if resname_4 in AMBER_TO_CHARMM_RESNAME \
                          else translate_residue_name(resname)
                # Translate atom name (using the new residue name)
                new_atom = translate_atom_name(new_res, atom_name)

                # Reformat atom name to 4 chars (PDB cols 13-16) and
                # residue name to 3 chars (PDB cols 18-20).
                # PDB atom-name column has special right/left alignment
                # rules; here we use a conservative left-aligned form.
                atom_field = _format_atom_name(new_atom)
                # CHARMM residue names are <= 4 chars; pdb residue field
                # is 3 chars (col 17-19) plus 1-char chain code (col 20)
                # which packmol-memgen leaves blank or 'A'. We write 4
                # chars into 17-20 if needed (CHARMM TIP3 etc.) and
                # leave chain blank.
                resname_field = (new_res + "    ")[:4]
                line = (
                    line[:12] + atom_field + " "
                    + resname_field + line[21:]
                )
            out_lines.append(line)

    Path(output_pdb).write_text("".join(out_lines))
    return str(Path(output_pdb).resolve())


def _format_atom_name(atom_name: str) -> str:
    """Format an atom name for PDB cols 13–16 (4 chars).

    Convention: 1- or 2-char element name names usually right-justify
    starting at col 14 (so "CA" → " CA "), but we keep things simple and
    left-justify at col 13 to avoid ambiguity when the atom name is
    > 3 chars (e.g. ``HT1`` from CT3 cap). pdb2gmx is forgiving about
    the leading-space convention.
    """
    s = atom_name.strip()
    return (s + "    ")[:4]


# ---------------------------------------------------------------------------
# Force-field staging
# ---------------------------------------------------------------------------

def stage_forcefield(*, charmm_ff_dir: str, build_dir: str) -> str:
    """Copy / symlink the CHARMM36 .ff directory into *build_dir*.

    pdb2gmx looks for the force field in ``./<ff_name>.ff`` relative to
    its current working directory. We mirror the source dir into
    ``build_dir/<basename>``.

    Returns the in-build path (``build_dir/<basename>``) that pdb2gmx's
    ``-ff <basename without .ff>`` will resolve.
    """
    src = Path(charmm_ff_dir).resolve()
    if not src.is_dir():
        raise FileNotFoundError(
            f"charmm_ff_dir does not exist: {src}. "
            f"See docs/membrane.md (CHARMM36 GROMACS port acquisition) "
            f"for the license-OK download route."
        )

    dst = Path(build_dir).resolve() / src.name
    if dst.exists():
        # Idempotent: assume the existing copy is fresh enough.
        return str(dst)
    # Prefer symlink (cheap, mirrors any updates) if same filesystem.
    try:
        os.symlink(src, dst)
    except OSError:
        # Fallback: full copy for cross-fs / Windows-without-symlink.
        shutil.copytree(src, dst)
    return str(dst)


def ff_name_from_dir(charmm_ff_dir: str) -> str:
    """Derive the pdb2gmx ``-ff`` argument from a force-field directory.

    GROMACS expects ``-ff <name>`` where the directory is named
    ``<name>.ff`` next to the working dir. So
    ``charmm36-jul2022.ff/`` → ``charmm36-jul2022``.
    """
    base = Path(charmm_ff_dir).name
    return base[:-3] if base.endswith(".ff") else base


# ---------------------------------------------------------------------------
# pdb2gmx wrapper
# ---------------------------------------------------------------------------

def run_pdb2gmx(
    *, input_pdb: str, ff_name: str, build_dir: str,
    water_model: str = "tip3p",
    n_term: str = "0",   # default: ACE-cap retained (none)
    c_term: str = "0",   # default: CT3-cap retained (none)
    gmx_path: str = "gmx",
) -> Tuple[str, str]:
    """Run ``gmx pdb2gmx`` and return (top_path, gro_path).

    *n_term* / *c_term* are the answers fed to pdb2gmx's interactive
    terminal-selection prompts (one numeric token per line):

    - ``"0"`` = first menu entry, typically ``NH3+`` for N-term
      (or ``"None"`` if the residue is already a cap like ACE).
    - With ACE / CT3 already present in the PDB, pdb2gmx detects them
      and may auto-skip the prompt; if not, ``"0"`` (None) is the
      safe default.
    """
    bd = Path(build_dir).resolve()
    bd.mkdir(parents=True, exist_ok=True)

    top_path = str(bd / "system.top")
    gro_path = str(bd / "system.gro")

    cmd = [
        gmx_path, "pdb2gmx",
        "-f", input_pdb,
        "-o", gro_path,
        "-p", top_path,
        "-ff", ff_name,
        "-water", water_model,
        "-ter",
        "-merge", "no",
    ]

    env = os.environ.copy()
    amberhome = _resolve_amberhome(gmx_path)
    if amberhome:
        env.setdefault("AMBERHOME", amberhome)

    logger.info("running pdb2gmx: %s", " ".join(cmd))
    stdin = f"{n_term}\n{c_term}\n"
    result = subprocess.run(
        cmd, cwd=bd, env=env, input=stdin,
        capture_output=True, text=True,
    )
    log_path = bd / "system_pdb2gmx.log"
    log_path.write_text(
        f"=== argv ===\n{' '.join(cmd)}\n\n"
        f"=== stdin ===\n{stdin}\n"
        f"=== stdout ===\n{result.stdout}\n"
        f"=== stderr ===\n{result.stderr}\n"
        f"=== returncode ===\n{result.returncode}\n"
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"gmx pdb2gmx failed (rc={result.returncode}). See {log_path}."
        )
    if not Path(top_path).is_file() or not Path(gro_path).is_file():
        raise RuntimeError(
            f"pdb2gmx returned 0 but did not produce {top_path} / {gro_path}. "
            f"Check {log_path}."
        )
    return top_path, gro_path


# ---------------------------------------------------------------------------
# Top-level entry point (matches parameterize_amber.parameterize signature)
# ---------------------------------------------------------------------------

def parameterize(
    *, config: MembraneConfig, input_pdb: str, build_dir: str,
) -> Tuple[str, str, str]:
    """Run CHARMM36 parameterisation and emit GROMACS files.

    Parameters
    ----------
    config : MembraneConfig
        Build configuration. Must have ``backend="charmm36"`` and
        a non-empty ``charmm_ff_dir``.
    input_pdb : str
        Bilayer + peptide PDB from stage 1 (packmol-memgen output with
        ``--charmm`` flag — lipids & water & ions already in CHARMM
        names; peptide may still carry AMBER names from tleap).
    build_dir : str
        Directory to stage the force field, run pdb2gmx, and write
        the GROMACS files.

    Returns
    -------
    (top, gro, ndx) : tuple[str, str, str]
        Absolute paths to system.top, system.gro, system.ndx.
    """
    if not config.charmm_ff_dir:
        raise ValueError(
            "MembraneConfig.charmm_ff_dir is empty. "
            "See docs/membrane.md (CHARMM36 GROMACS port acquisition)."
        )

    bd = Path(build_dir).resolve()
    bd.mkdir(parents=True, exist_ok=True)

    # 1. Stage the .ff dir so pdb2gmx -ff <name> resolves it.
    staged = stage_forcefield(charmm_ff_dir=config.charmm_ff_dir,
                              build_dir=str(bd))
    ff_name = ff_name_from_dir(staged)

    # 2. Translate the bilayer+peptide PDB to fully-CHARMM names.
    #    (packmol-memgen --charmm has already translated lipids & water
    #     & ions; this catches the peptide which is passed through.)
    translated_pdb = str(bd / "bilayer_peptide_charmm.pdb")
    translate_pdb_amber_to_charmm(
        input_pdb=input_pdb, output_pdb=translated_pdb,
    )

    # 3. pdb2gmx → top + gro.
    top_path, gro_path = run_pdb2gmx(
        input_pdb=translated_pdb,
        ff_name=ff_name,
        build_dir=str(bd),
        water_model="tip3p",
        gmx_path=config.gmx_path,
    )

    # 4. Index file (residue tables already cover CHARMM aliases).
    ndx_path = str(bd / "system.ndx")
    write_index_from_gro(gro_path=gro_path, ndx_path=ndx_path)

    return top_path, gro_path, ndx_path
