# -*- coding: utf-8 -*-
"""tests/cg_segmenter/conftest.py — PDB fixtures generated via RDKit."""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("rdkit")


def _build_pdb(smiles: str, out_path: Path, seed: int = 42) -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ret = AllChem.EmbedMolecule(mol, randomSeed=seed)
    if ret == -1:
        AllChem.EmbedMolecule(mol, randomSeed=seed + 100, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    Chem.MolToPDBFile(mol, str(out_path))
    return out_path


@pytest.fixture(scope="session")
def propane_pdb(tmp_path_factory):
    return _build_pdb("CCC", tmp_path_factory.mktemp("cg_prop") / "propane.pdb")


@pytest.fixture(scope="session")
def octane_pdb(tmp_path_factory):
    return _build_pdb("CCCCCCCC", tmp_path_factory.mktemp("cg_oct") / "octane.pdb")


@pytest.fixture(scope="session")
def methyl_acetate_pdb(tmp_path_factory):
    return _build_pdb("CC(=O)OC", tmp_path_factory.mktemp("cg_me_ace") / "me_ace.pdb")


@pytest.fixture(scope="session")
def benzene_pdb(tmp_path_factory):
    return _build_pdb("c1ccccc1", tmp_path_factory.mktemp("cg_benz") / "benzene.pdb")


@pytest.fixture(scope="session")
def naphthalene_pdb(tmp_path_factory):
    return _build_pdb("c1ccc2ccccc2c1", tmp_path_factory.mktemp("cg_naph") / "naph.pdb")


@pytest.fixture(scope="session")
def cyclohexanol_pdb(tmp_path_factory):
    return _build_pdb("OC1CCCCC1", tmp_path_factory.mktemp("cg_choh") / "cyclohexanol.pdb")


@pytest.fixture(scope="session")
def cyclohexyl_octyl_pdb(tmp_path_factory):
    return _build_pdb("CCCCCCCCC1CCCCC1", tmp_path_factory.mktemp("cg_co") / "ch_oct.pdb")
