# -*- coding: utf-8 -*-
"""Small-molecule GAFF2/AM1-BCC parameterization for formulation builds.

Wraps :func:`abmptools.genesis.mmgbsa.ligand_parameterize.run_acpype`
(GPL-3.0 subprocess only) and adds:

1. SMILES → 3D PDB seeding via RDKit (lazy import).
2. Per-species ``resname`` rewrite (3- or 4-letter tag).
3. Discovery of acpype's GROMACS outputs (``_GMX.itp`` + ``_GMX.gro``)
   in addition to the ``frcmod`` + ``mol2`` returned by the upstream
   helper, so the formulation builder can assemble ``topol.top`` and
   pack monomer PDBs from a single result object.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from ..genesis.mmgbsa.ligand_parameterize import (
    AcpypeResult,
    run_acpype as _genesis_run_acpype,
)
from ..genesis.mmgbsa.models import LigandParameterization

logger = logging.getLogger(__name__)


@dataclass
class SmallMoleculeResult:
    """All artifacts of one parameterised small-molecule species."""

    resname: str
    monomer_pdb: Path           # one-molecule PDB (input to packmol)
    acpype_dir: Path
    frcmod: Path
    mol2: Path
    itp: Path                   # GROMACS topology fragment
    gro: Path                   # one-molecule GROMACS coordinate
    n_atoms: int
    net_charge: int


def smiles_to_pdb(
    smiles: str,
    *,
    resname: str,
    out_pdb: Path,
    seed: int = 42,
) -> Path:
    """Embed *smiles* in 3D (MMFF94) and write a single-molecule PDB.

    RDKit is imported lazily — it is part of ``abmptools[formulation]``
    extras only.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ImportError(
            "rdkit is required to seed small molecules from SMILES. "
            "Install via `pip install abmptools[formulation]` or "
            "`pip install rdkit-pypi`."
        ) from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(
            f"RDKit ETKDGv3 embedding failed for SMILES {smiles!r}"
        )
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("MMFF optimize failed (continuing): %s", exc)

    # Rewrite residue name + atom name on every atom so the resulting
    # PDB stays within standard column widths (resname = 3 chars,
    # atom name = 4 chars including leading space).
    rname_3 = resname.strip().ljust(3)[:3]
    elem_counts: dict = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            info = Chem.AtomPDBResidueInfo()
        sym = atom.GetSymbol()
        elem_counts[sym] = elem_counts.get(sym, 0) + 1
        # PDB atom-name field is 4 chars; prefix with space so the
        # first letter lands in column 14 (heavy-atom convention).
        aname = f" {sym}{elem_counts[sym]}"[:4].ljust(4)
        info.SetName(aname)
        info.SetResidueName(rname_3)
        info.SetResidueNumber(1)
        info.SetChainId("A")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)

    out_pdb = Path(out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.PDBWriter(str(out_pdb))
    writer.write(mol)
    writer.close()
    return out_pdb


def _locate_gmx_outputs(acpype_dir: Path) -> tuple[Path, Path]:
    """Locate ``*_GMX.itp`` + ``*_GMX.gro`` produced by acpype."""
    itps = sorted(acpype_dir.glob("*_GMX.itp"))
    gros = sorted(acpype_dir.glob("*_GMX.gro"))
    if not itps:
        raise FileNotFoundError(f"No *_GMX.itp in {acpype_dir}")
    if not gros:
        raise FileNotFoundError(f"No *_GMX.gro in {acpype_dir}")
    return itps[0], gros[0]


def _locate_gaff_mol2(acpype_dir: Path, atom_type: str = "gaff2") -> Path:
    """Locate the GAFF/GAFF2-atom-type mol2 (handles both BCC and Gasteiger).

    Upstream ``genesis.mmgbsa.ligand_parameterize._locate_outputs`` hard-codes
    the ``*_bcc_*.mol2`` pattern, so Gasteiger-charge runs (which emit
    ``*_gas_*.mol2``) trip it. We re-scan with a broader pattern and pick
    whichever charge variant exists.
    """
    suffix = "_gaff2.mol2" if atom_type == "gaff2" else "_gaff.mol2"
    candidates = sorted(acpype_dir.glob(f"*_*{suffix}"))
    # Prefer BCC over GAS over OPLS / others.
    preferred_prefixes = ("bcc_", "gas_", "user_")
    for prefix in preferred_prefixes:
        for c in candidates:
            if f"_{prefix}" in c.name:
                return c
    if candidates:
        return candidates[0]
    raise FileNotFoundError(
        f"No *{suffix} found in {acpype_dir}; acpype may have failed."
    )


def _count_atoms_in_pdb(pdb_path: Path) -> int:
    n = 0
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("ATOM ", "HETATM")):
            n += 1
    return n


def parameterize_small_molecule(
    *,
    name: str,
    resname: str,
    smiles: str = "",
    pdb_path: str = "",
    net_charge: int = 0,
    workdir: Path,
    config: Optional[LigandParameterization] = None,
    acpype_path: str = "acpype",
    seed: int = 42,
) -> SmallMoleculeResult:
    """Generate GAFF2/AM1-BCC parameters for one species.

    Provide either *smiles* (preferred; reproducible from spec) or a
    *pdb_path* (when 3D coordinates are needed exactly).
    """
    if not smiles and not pdb_path:
        raise ValueError(
            f"parameterize_small_molecule[{name}] requires smiles or pdb_path."
        )
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    monomer_pdb = workdir / f"{resname}.pdb"
    if smiles:
        smiles_to_pdb(
            smiles, resname=resname, out_pdb=monomer_pdb, seed=seed,
        )
    else:
        src = Path(pdb_path).resolve()
        monomer_pdb.write_text(src.read_text())

    if config is None:
        config = LigandParameterization(
            charge_method="bcc",
            net_charge=net_charge,
            atom_type="gaff2",
            skip_if_cached=True,
        )

    acp: AcpypeResult = _genesis_run_acpype(
        ligand_pdb=monomer_pdb,
        workdir=workdir,
        config=config,
        acpype_path=acpype_path,
    )
    itp, gro = _locate_gmx_outputs(acp.acpype_dir)
    # Re-locate mol2 with the broader pattern so Gasteiger / user-charge
    # variants (``*_gas_*.mol2`` / ``*_user_*.mol2``) work too.
    mol2 = _locate_gaff_mol2(acp.acpype_dir, atom_type=config.atom_type)
    n_atoms = _count_atoms_in_pdb(monomer_pdb)
    return SmallMoleculeResult(
        resname=resname,
        monomer_pdb=monomer_pdb,
        acpype_dir=acp.acpype_dir,
        frcmod=acp.frcmod,
        mol2=mol2,
        itp=itp,
        gro=gro,
        n_atoms=n_atoms,
        net_charge=net_charge,
    )


__all__ = [
    "SmallMoleculeResult",
    "parameterize_small_molecule",
    "smiles_to_pdb",
]
