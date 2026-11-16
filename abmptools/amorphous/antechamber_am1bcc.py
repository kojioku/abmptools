# -*- coding: utf-8 -*-
"""AM1-BCC charges through antechamber on Windows, where OpenFF cannot start it.

OpenFF's ``AmberToolsToolkitWrapper.assign_partial_charges`` runs
``subprocess.check_output(["antechamber", ...])`` with ``-pf yes``. Two things stop
that on Windows with the native Windows build of AmberTools:

- that build keeps the program as ``bin\\wrapped_progs\\antechamber.exe`` behind
  ``bin\\antechamber.bat``; Windows starts a bare ``antechamber`` only as an ``.exe``
  on PATH, so the ``.bat`` is never reached and the call fails;
- ``-pf yes`` makes antechamber remove its scratch files by running ``rm -f`` through
  the shell, and Windows has no ``rm``: antechamber stops ("Cannot properly run rm -f").

:func:`assign_am1bcc` does what that method does (one RDKit conformer with the same
RMS cutoff, the same two antechamber commands, the charges read back in atom order,
normalised to the formal charge) with the ``.exe`` started directly, ``AMBERHOME`` set
and ``bin`` on PATH (where the build keeps its DLLs and ``sqm``), and ``-pf n``. The
molecule then carries its charges and Interchange is told to keep them
(``charge_from_molecules``), as for ``charge_method='nagl'`` / ``'gasteiger'``.
Linux and macOS keep OpenFF's own path.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def needed() -> bool:
    """Whether AM1-BCC must go through :func:`assign_am1bcc` here (Windows)."""
    return sys.platform == "win32"


def antechamber_command() -> Optional[Tuple[List[str], Dict[str, str]]]:
    """``([antechamber.exe], env)`` for AmberTools' real antechamber, or ``None``.

    ``AMBERHOME``'s ``bin\\wrapped_progs\\antechamber.exe`` first, then the one beside an
    ``antechamber.bat`` found on PATH, then an ``antechamber(.exe)`` on PATH as it is.
    """
    roots: List[Path] = []
    if os.environ.get("AMBERHOME"):
        roots.append(Path(os.environ["AMBERHOME"]))
    found = shutil.which("antechamber")
    if found and Path(found).suffix.lower() in (".bat", ".cmd"):
        roots.append(Path(found).resolve().parent.parent)
    for root in roots:
        for name in ("antechamber.exe", "antechamber"):
            exe = root / "bin" / "wrapped_progs" / name
            if exe.is_file():
                env = dict(os.environ)
                env["AMBERHOME"] = str(root)
                env["PATH"] = str(root / "bin") + os.pathsep + env.get("PATH", "")
                return [str(exe)], env
    if found and Path(found).suffix.lower() not in (".bat", ".cmd"):
        return [found], dict(os.environ)
    return None


def _run(argv: List[str], cwd: str, env: Dict[str, str]) -> None:
    r = subprocess.run(argv, cwd=cwd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        tail = "\n".join(((r.stdout or "") + (r.stderr or "")).strip().splitlines()[-8:])
        raise RuntimeError("antechamber failed (exit %d): %s\n%s"
                           % (r.returncode, " ".join(argv), tail))


def assign_am1bcc(mol: Any) -> None:
    """Set AM1-BCC partial charges on the OpenFF Molecule *mol* with antechamber / sqm."""
    import numpy as np
    from openff.toolkit.topology import Molecule
    from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
    from openff.units import unit

    cmd = antechamber_command()
    if cmd is None:
        raise RuntimeError(
            "AM1-BCC charges need AmberTools' antechamber, which was not found "
            "(AMBERHOME, or antechamber on PATH). Install AmberTools, or use "
            "charge_method='nagl' / 'gasteiger'.")
    antechamber, env = cmd

    rdkit = RDKitToolkitWrapper()
    mol_copy = Molecule(mol)
    # as AmberToolsToolkitWrapper does for am1bcc (rec_confs = 1)
    mol_copy.generate_conformers(n_conformers=1, rms_cutoff=0.25 * unit.angstrom,
                                 toolkit_registry=rdkit)
    net_charge = mol_copy.total_charge.m_as(unit.elementary_charge)
    with tempfile.TemporaryDirectory(prefix="abmptools_am1bcc_") as tmpdir:
        rdkit.to_file(mol_copy, os.path.join(tmpdir, "molecule.sdf"), file_format="sdf")
        _run(antechamber + ["-i", "molecule.sdf", "-fi", "sdf", "-o", "charged.mol2",
                            "-fo", "mol2", "-pf", "n", "-dr", "n", "-c", "bcc",
                            "-nc", str(net_charge)], tmpdir, env)
        _run(antechamber + ["-dr", "n", "-i", "charged.mol2", "-fi", "mol2",
                            "-o", "charges2.mol2", "-fo", "mol2", "-c", "wc",
                            "-cf", "charges.txt", "-pf", "n"], tmpdir, env)
        path = os.path.join(tmpdir, "charges.txt")
        if not os.path.exists(path):
            raise RuntimeError("antechamber / sqm wrote no charges for %s" % mol.to_smiles())
        with open(path, encoding="ascii") as fh:
            tokens = fh.read().split()
    if len(tokens) != mol.n_atoms:
        raise RuntimeError("antechamber returned %d charges for %d atoms (%s)"
                           % (len(tokens), mol.n_atoms, mol.to_smiles()))
    charges = np.array([float(t) for t in tokens], dtype=np.float64)
    mol.partial_charges = unit.Quantity(charges, unit.elementary_charge)
    mol._normalize_partial_charges()
