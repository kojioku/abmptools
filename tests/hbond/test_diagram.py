"""Tests for the 2D H-bond-site diagram (abmptools.hbond.diagram)."""
from types import SimpleNamespace

import pytest

from abmptools.hbond.diagram import (
    element_from_name, build_mol_from_topology, draw_hbond_diagram,
)


def test_element_from_name():
    assert element_from_name("Cl") == "Cl"
    assert element_from_name("CL") == "Cl"      # case-insensitive trailing
    assert element_from_name("cl") == "Cl"
    assert element_from_name("Br1") == "Br"
    assert element_from_name("C10") == "C"      # element + index
    assert element_from_name("HO") == "H"       # H named after its partner
    assert element_from_name("O") == "O"
    assert element_from_name("N4") == "N"
    assert element_from_name("") == ""


def _mock_topology_from_rdkit(mol):
    """Build a MoleculeTopology-like object + coords from an embedded RDKit mol."""
    conf = mol.GetConformer()
    atoms = [SimpleNamespace(atom_id=i, atom_name=a.GetSymbol(),
                             atom_type=a.GetSymbol())
             for i, a in enumerate(mol.GetAtoms())]
    bonds = [SimpleNamespace(type_name="", atom1=b.GetBeginAtomIdx(),
                             atom2=b.GetEndAtomIdx(), order=1.0)
             for b in mol.GetBonds()]
    topo = SimpleNamespace(mol_name="LIG", atoms=atoms, bonds=bonds)
    coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    return topo, coords


def test_build_mol_from_topology_perceives_bond_orders():
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors

    # acetamide CC(=O)N: single-bonded topology + 3D coords must round-trip to
    # the correct formula with the C=O double bond perceived from geometry.
    ref = Chem.AddHs(Chem.MolFromSmiles("CC(=O)N"))
    assert AllChem.EmbedMolecule(ref, randomSeed=1) == 0
    topo, coords = _mock_topology_from_rdkit(ref)

    mol, reason = build_mol_from_topology(topo, coords, charge=0)
    assert mol is not None, reason
    assert rdMolDescriptors.CalcMolFormula(mol) == "C2H5NO"
    # one C=O double bond perceived
    assert sum(1 for b in mol.GetBonds()
               if b.GetBondTypeAsDouble() == 2.0) == 1


def test_build_mol_from_topology_atom_count_mismatch():
    topo = SimpleNamespace(mol_name="X",
                           atoms=[SimpleNamespace(atom_name="C", atom_type="C")],
                           bonds=[])
    mol, reason = build_mol_from_topology(topo, coords=[], charge=0)
    assert mol is None and "mismatch" in reason


def test_draw_hbond_diagram_writes_files(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    ref = Chem.AddHs(Chem.MolFromSmiles("CC(=O)N"))
    assert AllChem.EmbedMolecule(ref, randomSeed=1) == 0
    topo, coords = _mock_topology_from_rdkit(ref)
    mol, reason = build_mol_from_topology(topo, coords, charge=0)
    assert mol is not None, reason

    # acetamide: amide N is the donor heavy atom, carbonyl O the acceptor.
    n_idx, c_idx, o_idx = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[NX3][CX3]=[OX1]"))[0]
    png = str(tmp_path / "lig_diagram.png")
    svg = str(tmp_path / "lig_diagram.svg")
    written = draw_hbond_diagram(mol, [n_idx], [o_idx], png, svg,
                                 legend="LIG test")
    assert set(written) == {png, svg}
    import os
    assert os.path.getsize(png) > 0 and os.path.getsize(svg) > 0


def test_draw_hbond_diagram_both_role(tmp_path):
    # an atom present in BOTH donor and acceptor lists (e.g. an OH oxygen) must
    # render without error via the donor+acceptor (magenta) branch.
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    ref = Chem.AddHs(Chem.MolFromSmiles("CC(=O)N"))
    assert AllChem.EmbedMolecule(ref, randomSeed=1) == 0
    topo, coords = _mock_topology_from_rdkit(ref)
    mol, reason = build_mol_from_topology(topo, coords, charge=0)
    assert mol is not None, reason
    n_idx, c_idx, o_idx = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[NX3][CX3]=[OX1]"))[0]
    png = str(tmp_path / "both.png")
    written = draw_hbond_diagram(mol, [n_idx, o_idx], [o_idx], png)  # O in both
    import os
    assert written == [png] and os.path.getsize(png) > 0
