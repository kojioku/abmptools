# -*- coding: utf-8 -*-
"""tests/fragmenter/conftest.py

PDB fixtures (rdkit で生成、tmp_path_factory にキャッシュ)。
"""
from __future__ import annotations

from pathlib import Path

import pytest

# rdkit が無い場合はこのテストモジュール全体を skip
pytest.importorskip("rdkit")


def _build_pdb(smiles: str, out_path: Path, n_copies: int = 1, seed: int = 42) -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    combined = mol
    for _ in range(n_copies - 1):
        combined = Chem.CombineMols(combined, mol)
    Chem.MolToPDBFile(combined, str(out_path))
    return out_path


@pytest.fixture(scope="session")
def propane_pdb(tmp_path_factory):
    p = tmp_path_factory.mktemp("frag_prop") / "propane.pdb"
    return _build_pdb("CCC", p, n_copies=1)


@pytest.fixture(scope="session")
def octane_pdb(tmp_path_factory):
    p = tmp_path_factory.mktemp("frag_oct") / "octane.pdb"
    return _build_pdb("CCCCCCCC", p, n_copies=1)


@pytest.fixture(scope="session")
def propane_acetone_mix_pdb(tmp_path_factory):
    """propane x3 + acetone x2 を 1 PDB に combine。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    propane = Chem.AddHs(Chem.MolFromSmiles("CCC"))
    acetone = Chem.AddHs(Chem.MolFromSmiles("CC(=O)C"))
    AllChem.EmbedMolecule(propane, randomSeed=42)
    AllChem.EmbedMolecule(acetone, randomSeed=43)
    combined = propane
    for _ in range(2):
        combined = Chem.CombineMols(combined, propane)
    for _ in range(2):
        combined = Chem.CombineMols(combined, acetone)

    p = tmp_path_factory.mktemp("frag_mix") / "mix.pdb"
    Chem.MolToPDBFile(combined, str(p))
    return p


@pytest.fixture(scope="session")
def pe_n10_n11_mix_pdb(tmp_path_factory):
    """PE N=10 + PE N=11 を 1 PDB に combine (γ 経路テスト用)。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    pe10 = Chem.AddHs(Chem.MolFromSmiles("CCCCCCCCCC"))      # N=10
    pe11 = Chem.AddHs(Chem.MolFromSmiles("CCCCCCCCCCC"))     # N=11
    AllChem.EmbedMolecule(pe10, randomSeed=42)
    AllChem.EmbedMolecule(pe11, randomSeed=43)
    combined = Chem.CombineMols(pe10, pe11)

    p = tmp_path_factory.mktemp("frag_pe") / "pe10_11.pdb"
    Chem.MolToPDBFile(combined, str(p))
    return p
