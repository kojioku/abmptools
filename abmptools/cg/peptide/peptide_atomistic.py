# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide.peptide_atomistic
---------------------------------------
Atomistic peptide PDB generation for downstream Martini 3 CG mapping.

Two backends:

1. **tleap (推奨)**
       AmberTools tleap で ``leaprc.protein.ff14SB`` の sequence builder
       を使い、フル sidechain の atomistic PDB を生成。

2. **Extended backbone (fallback)**
       tleap が不在のとき python だけで beta-strand 風の伸び切り構造を
       書き出す。N / CA / C / O + CB のみ (GLY は CB なし)。

       警告: sidechain heavy atom が欠落するため、芳香族残基 (W/F/Y)
       で martinize2 が CG bead 座標を NaN にする可能性あり。研究品質
       では tleap を推奨する。

Public API: :func:`build_atomistic_pdb`。
"""
from __future__ import annotations

import logging
import shutil
import tempfile
from pathlib import Path
from typing import Dict, List

from ._subprocess import ensure_dir, run_command, write_text

logger = logging.getLogger(__name__)


# 1-letter -> 3-letter (tleap-compatible)
AA_3LETTER: Dict[str, str] = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}

#: Aromatic residues whose Martini 3 sidechain mapping (multi-bead ring)
#: is most likely to produce NaN bead coordinates when fed an extended
#: backbone PDB lacking sidechain heavy atoms.
ARTIFACT_PRONE_RESIDUES = set("WFY")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_atomistic_pdb(
    sequence: str,
    name: str,
    output_dir: Path,
    *,
    prefer_tleap: bool = True,
    tleap_path: str = "tleap",
) -> Path:
    """Generate an atomistic peptide PDB suitable for martinize2 input.

    Tries tleap first when ``prefer_tleap=True`` (default); falls back to
    a simple extended backbone if tleap is not on PATH. The fallback path
    logs a WARNING; sequences containing aromatic residues (W/F/Y) trigger
    an additional WARNING about possible NaN bead coordinates after CG
    mapping.

    Parameters
    ----------
    sequence
        One-letter amino acid sequence (uppercase; validated upstream by
        :class:`PeptideSpec`).
    name
        Molecule name; used for output filename.
    output_dir
        Destination directory (created if missing).
    prefer_tleap
        Try tleap before falling back. Default True.
    tleap_path
        tleap executable name or absolute path.

    Returns
    -------
    Path
        Path to the written ``<name>_atomistic.pdb``.
    """
    output_dir = ensure_dir(Path(output_dir)).resolve()

    if prefer_tleap and shutil.which(tleap_path):
        return _build_with_tleap(sequence, name, output_dir, tleap_path)

    if prefer_tleap:
        logger.warning(
            "tleap (%r) not found on PATH. Falling back to extended-backbone "
            "atomistic structure. 研究品質では tleap (AmberTools) を推奨。 "
            "conda 経由インストール: mamba install -c conda-forge ambertools",
            tleap_path,
        )
    else:
        logger.info(
            "tleap disabled by caller (prefer_tleap=False); "
            "using extended-backbone fallback."
        )

    artifact_residues = set(sequence) & ARTIFACT_PRONE_RESIDUES
    if artifact_residues:
        logger.warning(
            "Sequence contains aromatic residues sensitive to atomistic "
            "completeness (%s). Extended-backbone fallback omits sidechain "
            "heavy atoms; martinize2 may produce NaN CG bead coordinates "
            "for these residues. Install tleap to avoid this artifact.",
            "".join(sorted(artifact_residues)),
        )

    return _build_extended_backbone_pdb(sequence, name, output_dir)


# ---------------------------------------------------------------------------
# tleap backend
# ---------------------------------------------------------------------------

def _build_with_tleap(
    sequence: str,
    name: str,
    output_dir: Path,
    tleap_path: str,
) -> Path:
    """Run tleap (ff14SB) to produce a full-atom peptide PDB."""
    tleap_resnames = " ".join(AA_3LETTER[aa] for aa in sequence)
    pdb_path = output_dir / f"{name}_atomistic.pdb"

    leap_input = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".leap", delete=False, dir=output_dir,
        ) as f:
            f.write("source leaprc.protein.ff14SB\n")
            f.write(f"seq = sequence {{ {tleap_resnames} }}\n")
            f.write(f"savepdb seq {pdb_path}\n")
            f.write("quit\n")
            leap_input = Path(f.name)

        run_command([tleap_path, "-f", str(leap_input)], cwd=output_dir)
    finally:
        if leap_input is not None:
            leap_input.unlink(missing_ok=True)

    if not pdb_path.exists():
        raise RuntimeError(f"tleap failed to generate {pdb_path}")

    logger.info(
        "Built atomistic peptide PDB (tleap): %s (%d residues)",
        pdb_path, len(sequence),
    )
    return pdb_path


# ---------------------------------------------------------------------------
# Extended backbone fallback
# ---------------------------------------------------------------------------

def _build_extended_backbone_pdb(
    sequence: str, name: str, output_dir: Path,
) -> Path:
    """Write a beta-strand-like extended atomistic PDB (no full sidechain)."""
    atoms = _build_extended_backbone(sequence)
    content = _atoms_to_pdb(atoms, title=f"{name} ({sequence})")
    pdb_path = output_dir / f"{name}_atomistic.pdb"
    write_text(pdb_path, content)
    logger.info(
        "Built atomistic peptide PDB (extended-backbone fallback): "
        "%s (%d residues)",
        pdb_path, len(sequence),
    )
    return pdb_path


def _build_extended_backbone(sequence: str) -> List[Dict]:
    """Generate backbone atoms in extended (beta-strand-like) conformation.

    Coordinates in Angstrom. Atom set: N, CA, C, O + CB (CB skipped for GLY).
    """
    atoms: List[Dict] = []
    serial = 1
    ca_spacing = 3.80  # CA-CA distance, Angstrom

    for i, aa in enumerate(sequence):
        resname = AA_3LETTER[aa]
        resid = i + 1
        x_base = i * ca_spacing

        atoms.append(_atom(serial, "N",  resname, resid,
                           x_base - 1.20, 0.0, 0.0));  serial += 1
        atoms.append(_atom(serial, "CA", resname, resid,
                           x_base,        0.0, 0.0));  serial += 1
        atoms.append(_atom(serial, "C",  resname, resid,
                           x_base + 0.90, 0.0, 0.0));  serial += 1
        atoms.append(_atom(serial, "O",  resname, resid,
                           x_base + 0.90, 1.24, 0.0)); serial += 1

        if aa != "G":
            atoms.append(_atom(serial, "CB", resname, resid,
                               x_base, -1.53, 0.0));   serial += 1

    return atoms


def _atom(serial, name, resname, resid, x, y, z) -> Dict:
    return {
        "serial": serial, "name": name, "resname": resname, "resid": resid,
        "x": x, "y": y, "z": z,
    }


def _atoms_to_pdb(atoms: List[Dict], title: str = "Peptide") -> str:
    lines = [f"TITLE     {title}", "MODEL     1"]
    for atom in atoms:
        lines.append(
            f"ATOM  {atom['serial']:5d} {atom['name']:<4s} "
            f"{atom['resname']:>3s}  {atom['resid']:4d}    "
            f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
            f"  1.00  0.00           {atom['name'][0]:>1s}"
        )
    lines.append("ENDMDL")
    lines.append("END")
    return "\n".join(lines) + "\n"
