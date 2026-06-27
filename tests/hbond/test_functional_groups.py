"""
Test functional group detection on the IMC BDF.
"""
import pytest

from abmptools.hbond.bdf_reader import BDFTrajectory
from abmptools.hbond.functional_groups import (
    detect_amides,
    detect_carboxyls,
    detect_hydroxyls,
    summarize_groups,
)

IMC_BDF = "/home/okuwaki/llm-project/SI/IMC_result450.0_out_rec900.bdf"


@pytest.fixture(scope="module")
def imc_traj():
    traj = BDFTrajectory(IMC_BDF)
    traj.load_topology()
    return traj


def test_load_topology(imc_traj):
    assert len(imc_traj.molecules) == 125
    assert len(imc_traj.molecules[0].atoms) == 41
    assert len(imc_traj.molecules[0].bonds) == 43


def test_detect_carboxyls(imc_traj):
    carboxyls = detect_carboxyls(imc_traj.molecules)
    assert len(carboxyls) == 125
    cg0 = carboxyls[0]
    assert cg0.mol_index == 0
    assert cg0.c_atom == 18
    assert cg0.o_atom == 4
    assert cg0.oh_atom == 3
    assert cg0.ho_atom == 40


def test_detect_amides(imc_traj):
    amides = detect_amides(imc_traj.molecules)
    assert len(amides) == 125
    ag0 = amides[0]
    assert ag0.mol_index == 0
    assert ag0.c_atom == 14
    assert ag0.o_atom == 2
    assert ag0.n_atom == 5
    assert ag0.tert is True


def test_detect_hydroxyls(imc_traj):
    """Only the carboxyl OH should be found; with exclude_carboxyl=True → 0."""
    hyds_excl = detect_hydroxyls(imc_traj.molecules, exclude_carboxyl=True)
    hyds_all = detect_hydroxyls(imc_traj.molecules, exclude_carboxyl=False)
    assert len(hyds_excl) == 0
    assert len(hyds_all) == 125


def test_summarize_groups(imc_traj):
    s = summarize_groups(imc_traj.molecules)
    assert s["n_carboxyls"] == 125
    assert s["n_amides"] == 125
    assert s["n_mols_with_carboxyl"] == 125
    assert s["n_mols_with_amide"] == 125


# --- ether detection (element + connectivity, FF-agnostic) -------------------

from abmptools.hbond.bdf_reader import AtomInfo, BondInfo, MoleculeTopology
from abmptools.hbond.functional_groups import detect_ethers


def _mol(atom_names, bonds):
    atoms = [AtomInfo(atom_id=i, atom_name=nm, atom_type=nm)
             for i, nm in enumerate(atom_names)]
    bond_objs = [BondInfo(type_name="", atom1=a, atom2=b, order=1.0)
                 for a, b in bonds]
    return MoleculeTopology(mol_name="LIG", atoms=atoms, bonds=bond_objs)


def test_detect_ethers_dimethyl_ether():
    # CH3-O-CH3: O bonded to two carbons, no H -> one ether
    names = ["C", "O", "C", "H", "H", "H", "H", "H", "H"]
    bonds = [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5), (2, 6), (2, 7), (2, 8)]
    ethers = detect_ethers([_mol(names, bonds)])
    assert len(ethers) == 1
    assert ethers[0].mol_index == 0
    assert ethers[0].o_atom == 1


def test_detect_ethers_excludes_hydroxyl():
    # methanol C-O-H: the O has an H neighbour -> not an ether
    names = ["C", "O", "H", "H", "H", "H"]
    bonds = [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5)]
    assert detect_ethers([_mol(names, bonds)]) == []


def test_detect_ethers_excludes_carbonyl():
    # acetone C(=O): the carbonyl O has a single heavy neighbour -> not an ether
    names = ["C", "O", "C", "C", "H", "H", "H", "H", "H", "H"]
    bonds = [(0, 1), (0, 2), (0, 3),
             (2, 4), (2, 5), (2, 6), (3, 7), (3, 8), (3, 9)]
    assert detect_ethers([_mol(names, bonds)]) == []
