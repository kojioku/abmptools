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

#: Residue-name map (AMBER → CHARMM36 GROMACS port, Klauda lab).
#:
#: The Klauda lab port keeps the *AMBER* names for many residues that
#: classic CHARMM (top_all36_prot.rtf) renames, so this table is
#: smaller than one might expect:
#:
#: - **NME** stays as ``NME`` (not ``CT3``). The port defines
#:   ``[ NME ]`` directly in ``aminoacids.rtp``.
#: - **CYM** stays as ``CYM`` (the port has ``[ CYM ]`` directly).
#: - **HIP** stays as ``HIP`` (the port has both ``[ HIP ]`` and
#:   ``[ HSP ]``; HIP is acceptable verbatim).
#: - **ACE** stays as ``ACE`` but the *atom* names need translation
#:   (see PER_RESIDUE_ATOM_MAP).
#:
#: Standard 20 amino acids share names between AMBER and CHARMM36
#: (ALA, GLY, LYS, ...) so they don't appear here.
AMBER_TO_CHARMM_RESNAME: Dict[str, str] = {
    # histidine tautomers (HID/HIE → HSD/HSE; HIP shared, kept as-is)
    "HID": "HSD",
    "HIE": "HSE",
    # protonation states (ASPP/GLUP/LSN exist in port, ASH/GLH/LYN don't)
    "ASH": "ASPP",
    "GLH": "GLUP",
    "LYN": "LSN",
    # disulfide cysteine: AMBER CYX → CHARMM CYS (with -ss / LINK patch)
    "CYX": "CYS",
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
    # Acetyl N-terminal cap. The Klauda port's ACE expects
    # CH3 / HH31 / HH32 / HH33 / C / O — **NOT** the classic
    # CAY / HY1..HY3 / CY / OY. Only the methyl hydrogens need
    # renaming from AMBER's H1/H2/H3 to HH31/HH32/HH33.
    "ACE": {
        "H1": "HH31",
        "H2": "HH32",
        "H3": "HH33",
        # CH3, C, O — pass through (port keeps these AMBER-style names)
    },
    # AMBER tleap-built NME has atoms: N, H, **C**, **H1, H2, H3**.
    # The CHARMM36 port NME expects: N, HN, **CH3**, **HH31, HH32, HH33**.
    # Both are 6-atom residues representing the same chemistry; only
    # the methyl C / H names differ. (H → HN is handled by the universal
    # backbone map.)
    "NME": {
        "C":  "CH3",
        "H1": "HH31",
        "H2": "HH32",
        "H3": "HH33",
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

#: Residues for which :data:`UNIVERSAL_BACKBONE_ATOM_MAP` does NOT apply.
#: ACE has no backbone amide H (it's the N-terminal cap providing the C=O).
#: Ions / water also have no protein backbone. NME stays in scope
#: because its amide H must become HN (handled by the universal map).
SKIP_BACKBONE_RENAME: frozenset = frozenset({"ACE", "TIP3", "SOD",
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
                # Read the full 4-char residue field (cols 18-21).
                # CHARMM residues like POPC / TIP3 use all 4 chars;
                # 3-char names (ALA / CHL) leave col 21 blank, and the
                # ``.strip()`` below normalises both forms.
                # AMBER's "Na+ " (3-char + sign) also fits here.
                resname_4 = line[17:21].strip()
                # ``translate_residue_name`` returns its input verbatim
                # when no mapping is defined (POPC → POPC, ALA → ALA),
                # so 4-char names pass through unchanged. The earlier
                # if/else fallback to ``line[17:20]`` was a bug — it
                # truncated 4-char unmapped names like POPC to "POP",
                # which pdb2gmx then rejected with
                # "Residue 'POP' not found in residue topology database".
                new_res = translate_residue_name(resname_4)
                # Translate atom name (using the new residue name)
                new_atom = translate_atom_name(new_res, atom_name)

                # Reformat atom name to 4 chars (PDB cols 13-16) and
                # residue name to 4 chars (PDB cols 18-21). For 3-char
                # residues, col 21 is naturally a trailing space; for
                # 4-char residues (POPC / TIP3 / ...) the chain ID at
                # col 22 is preserved by ``line[21:]`` after the slice.
                atom_field = _format_atom_name(new_atom)
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

#: Hardcoded "None" terminus indices per known CHARMM36 GROMACS port.
#:
#: pdb2gmx prints an interactive menu listing terminal-patch options
#: (``NH3+``, ``COO-``, ``None``, ...). The menu order depends on which
#: ``.n.tdb`` / ``.c.tdb`` files the FF loads. We hardcode the known
#: ports because programmatic probing via ``subprocess.run(... input=...)``
#: turned out to be unreliable: pdb2gmx 2021.3 spins forever on the
#: bad-index input instead of erroring and exiting (observed 99% CPU
#: for 1+ hour during development).
KNOWN_CHARMM_TERMINUS_NONE: Dict[str, Tuple[int, int]] = {
    # Klauda lab GROMACS port, Feb 2026 update, all 4 variants:
    #   N-term: 0:NH3+ 1:NH2 2:HYD1 3:MET1 4:5TER 5:5MET 6:5PHO 7:5POM 8:None
    #   C-term: 0:COO- 1:COOH 2:CT2 3:CT1 4:HYD2 5:MET2 6:3TER 7:None
    "charmm36-feb2026_cgenff-5.0":         (8, 7),
    "charmm36-feb2026_cgenff-4.6":         (8, 7),
    "charmm36-feb2026_ljpme_cgenff-5.0":   (8, 7),
    "charmm36-feb2026_ljpme_cgenff-4.6":   (8, 7),
    # Older July 2022 port (same .tdb structure):
    "charmm36-jul2022":                    (8, 7),
    "charmm36_ljpme-jul2022":              (8, 7),
}


def _resolve_terminus_none_indices(ff_name: str) -> Tuple[int, int]:
    """Return ``(n_term_idx, c_term_idx)`` for the ``None`` patch.

    Looks up :data:`KNOWN_CHARMM_TERMINUS_NONE` first; falls back to
    ``(8, 7)`` (the layout used by all current Klauda ports). Users
    with custom CHARMM ports should pass ``n_term`` / ``c_term``
    explicitly to :func:`run_pdb2gmx`.
    """
    return KNOWN_CHARMM_TERMINUS_NONE.get(ff_name, (8, 7))


def _strip_water_spurious_angles(itp_path: str) -> int:
    """Remove pdb2gmx's auto-generated spurious O-H-H angles in water itps.

    The CHARMM port's ``solvent.rtp`` lists three "bonds" for TIP3:
    ``OH2-H1``, ``OH2-H2``, **and** ``H1-H2`` (the third is the rigid-water
    SETTLE constraint, listed as a bond for topology purposes). pdb2gmx
    interprets this as a fully-bonded 3-cycle and writes three angles
    per water:

      - ``H1-OH2-H2`` (real H-O-H angle, defined in ffbonded.itp)
      - ``OH2-H1-H2`` (spurious — H is the middle atom)
      - ``OH2-H2-H1`` (spurious — symmetric)

    The CHARMM ffbonded.itp only defines the ``HT-OT-HT`` angle type, so
    the spurious O-H-H entries trigger "No default U-B types" at grompp,
    one error per spurious angle per water.

    This function rewrites the itp in place, dropping any ``[ angles ]``
    line whose middle atom (the second of three indices) is **not** an
    OT-typed (water oxygen) atom. Lipid / protein chains never have OT
    atoms, so their angles are passed through unchanged — the function
    is safe to call on every per-chain itp.

    Returns the number of angle lines removed.
    """
    path = Path(itp_path)
    text = path.read_text().splitlines(keepends=True)

    # Pass 1: locate [ atoms ] block, build idx → type map for OT atoms only
    in_atoms = False
    is_ot: Dict[int, bool] = {}
    for line in text:
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            in_atoms = (s == "[ atoms ]")
            continue
        if not in_atoms or not s or s.startswith(";"):
            continue
        toks = s.split()
        # [ atoms ] format: idx type resnr resname atname cgnr q m
        try:
            idx = int(toks[0])
        except (ValueError, IndexError):
            continue
        atom_type = toks[1] if len(toks) >= 2 else ""
        is_ot[idx] = (atom_type == "OT")

    if not any(is_ot.values()):
        return 0  # no water in this chain — nothing to strip

    # Pass 2: rewrite, filtering [ angles ] section
    out: List[str] = []
    in_angles = False
    dropped = 0
    for line in text:
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            in_angles = (s == "[ angles ]")
            out.append(line)
            continue
        if not in_angles or not s or s.startswith(";"):
            out.append(line)
            continue
        toks = s.split()
        if len(toks) < 4:
            out.append(line)
            continue
        try:
            ai, aj, ak = int(toks[0]), int(toks[1]), int(toks[2])
        except ValueError:
            out.append(line)
            continue
        # Keep the angle if middle atom is an OT (oxygen of water).
        # In water, only H-O-H (middle=O) has a valid CHARMM angle type;
        # O-H-H (middle=H) does not.
        if is_ot.get(aj, False):
            out.append(line)
        else:
            dropped += 1

    path.write_text("".join(out))
    return dropped


def run_pdb2gmx(
    *, input_pdb: str, ff_name: str, build_dir: str,
    water_model: str = "tip3p",
    n_term: Optional[str] = None,
    c_term: Optional[str] = None,
    gmx_path: str = "gmx",
) -> Tuple[str, str]:
    """Run ``gmx pdb2gmx`` and return (top_path, gro_path).

    *n_term* / *c_term* are the answers fed to pdb2gmx's interactive
    terminal-selection prompts (one numeric token per line). When left
    as ``None`` (default), the indices of the ``"None"`` (no-op)
    terminus entries are auto-discovered by probing pdb2gmx's menu —
    this is the typical use case for ACE-NME-capped peptides where the
    caps are already present in the PDB and no further patching is
    desired.

    The probe is needed because different CHARMM36 ports load different
    .tdb files (aminoacids / nucleic acids / lipids), shifting the menu
    order. ``charmm36-feb2026_cgenff-*`` for example places ``None`` at
    index 8 (N-term) and 7 (C-term), not 0.
    """
    bd = Path(build_dir).resolve()
    bd.mkdir(parents=True, exist_ok=True)

    top_path = str(bd / "system.top")
    gro_path = str(bd / "system.gro")

    env = os.environ.copy()
    amberhome = _resolve_amberhome(gmx_path)
    if amberhome:
        env.setdefault("AMBERHOME", amberhome)

    # Resolve None indices for any unspecified terminus from the
    # KNOWN_CHARMM_TERMINUS_NONE lookup table.
    if n_term is None or c_term is None:
        n_idx, c_idx = _resolve_terminus_none_indices(ff_name)
        if n_term is None:
            n_term = str(n_idx)
        if c_term is None:
            c_term = str(c_idx)
        logger.info(
            "terminus 'None' indices for %s: N=%s C=%s",
            ff_name, n_term, c_term,
        )

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

    # Post-process: strip pdb2gmx's spurious water angles from per-chain
    # itps. CHARMM TIP3's rigid SETTLE constraint listed as an H1-H2
    # bond causes pdb2gmx to write extra O-H-H angles that ffbonded.itp
    # cannot resolve (one "No default U-B types" error each).
    total_dropped = 0
    for chain_itp in sorted(bd.glob("system_Other_chain_*.itp")):
        d = _strip_water_spurious_angles(str(chain_itp))
        total_dropped += d
        if d:
            logger.info(
                "stripped %d spurious water O-H-H angles from %s",
                d, chain_itp.name,
            )
    if total_dropped:
        logger.info(
            "removed %d spurious water angles in total (CHARMM TIP3 fix)",
            total_dropped,
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
