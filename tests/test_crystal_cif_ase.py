# -*- coding: utf-8 -*-
"""Unit tests for ``abmptools.crystal.cif_engine_ase``.

The ASE backend is opt-in (``CIFEngineConfig.engine='ase'``); these
tests exercise the supercell + molecule-detection path against the
csp7 R00001 fixture. ASE's coordinate origin differs from the legacy
parser (``Atoms.repeat`` does not centre), so byte-equivalence with
the Phase B fixture is **not** expected — the tests check structural
invariants (atom counts, molecule-size partitioning) instead.
"""
from __future__ import annotations

import os

import pytest

# ASE is an optional dependency; skip the module if it's not installed.
ase = pytest.importorskip("ase", reason="ase not installed")

TESTS_DIR = os.path.dirname(__file__)
REF_MAIN = os.path.join(TESTS_DIR, "regression", "reference", "main")
CRYSTAL_FIXTURE = os.path.join(REF_MAIN, "crystal_csp7")
R00001_CIF = os.path.join(CRYSTAL_FIXTURE, "R00001", "XXXI-MMFF-R00001.cif")
R00002_CIF = os.path.join(CRYSTAL_FIXTURE, "R00002", "XXXI-MMFF-R00002.cif")
R00004_CIF = os.path.join(CRYSTAL_FIXTURE, "R00004", "XXXI-MMFF-R00004.cif")


# ---------------------------------------------------------------------------
# read_cif_to_atoms
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_read_cif_returns_expanded_asymmetric_unit():
    """csp7 R00001 should expand to 128 atoms (32 in asym unit × 4 symm)."""
    from abmptools.crystal.cif_engine_ase import read_cif_to_atoms

    atoms = read_cif_to_atoms(R00001_CIF)
    # P21/n has 4 symmetry operations; asymmetric unit holds 32 atoms.
    assert len(atoms) == 128, f"got {len(atoms)} atoms, expected 128"


@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_read_cif_cell_parameters_match_csp7():
    """The CIF's cell parameters should round-trip through ASE."""
    from abmptools.crystal.cif_engine_ase import read_cif_to_atoms

    atoms = read_cif_to_atoms(R00001_CIF)
    a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
    # Approximate cell from CIF (5.7264, 12.8467, 18.5441, 90, 91.733, 90).
    assert a == pytest.approx(5.7264, abs=1e-3)
    assert b == pytest.approx(12.8467, abs=1e-3)
    assert c == pytest.approx(18.5441, abs=1e-3)
    assert alpha == pytest.approx(90.0, abs=1e-3)
    assert beta == pytest.approx(91.733, abs=1e-2)
    assert gamma == pytest.approx(90.0, abs=1e-3)


# ---------------------------------------------------------------------------
# expand_supercell
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_expand_supercell_layer5_atom_count():
    """Layer 5 -> 5^3 cells * 128 atoms/cell = 16000 atoms (matches legacy)."""
    from abmptools.crystal.cif_engine_ase import (
        expand_supercell,
        read_cif_to_atoms,
    )

    atoms = read_cif_to_atoms(R00001_CIF)
    super_atoms = expand_supercell(atoms, layer=5)
    assert len(super_atoms) == 128 * 125  # = 16000


def test_expand_supercell_invalid_layer():
    from abmptools.crystal.cif_engine_ase import expand_supercell

    # Build a minimal Atoms manually so we don't need a CIF.
    from ase import Atoms
    atoms = Atoms("H2", positions=[[0, 0, 0], [0.74, 0, 0]])
    atoms.set_cell([10, 10, 10])
    with pytest.raises(ValueError, match="layer must be"):
        expand_supercell(atoms, layer=0)


# ---------------------------------------------------------------------------
# detect_molecules
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_detect_molecules_layer5_partition():
    """csp7 R00001 layer 5: 16000 atoms -> 500 molecules of 32 atoms each."""
    from abmptools.crystal.cif_engine_ase import (
        detect_molecules,
        expand_supercell,
        read_cif_to_atoms,
    )

    atoms = read_cif_to_atoms(R00001_CIF)
    super_atoms = expand_supercell(atoms, layer=5)
    molecules = detect_molecules(super_atoms, atoms_in_mol=[32])
    assert len(molecules) == 500
    # All molecules must have 32 atoms.
    assert all(len(m) == 32 for m in molecules)
    # Total atom count matches.
    assert sum(len(m) for m in molecules) == 16000


@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_detect_molecules_wrong_atoms_in_mol_raises():
    """An incorrect atoms_in_mol should fail with a clear message."""
    from abmptools.crystal.cif_engine_ase import (
        detect_molecules,
        expand_supercell,
        read_cif_to_atoms,
    )

    atoms = read_cif_to_atoms(R00001_CIF)
    super_atoms = expand_supercell(atoms, layer=2)  # smaller for speed
    with pytest.raises(ValueError, match="do not match atoms_in_mol"):
        detect_molecules(super_atoms, atoms_in_mol=[31])  # wrong count


# ---------------------------------------------------------------------------
# run_ase end-to-end
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.parametrize(
    "cif_path,layer,n_mol_expected",
    [
        # csp7 R00001: P21/n (Z=4) -> 128 atoms/cell * 125 cells = 16000 atoms = 500 mol
        (R00001_CIF, 5, 500),
        # csp7 R00002: Z=8 -> 256 atoms/cell * 125 = 32000 atoms = 1000 mol
        (R00002_CIF, 5, 1000),
        # csp7 R00004: P-1 (Z=2) -> 64 atoms/cell * 125 = 8000 atoms = 250 mol
        (R00004_CIF, 5, 250),
    ],
)
def test_run_ase_layer5_atom_partition(cif_path, layer, n_mol_expected):
    """Top-level orchestrator: csp7 small/mid/large structures via ASE.

    Atom-count expectation depends on the space group and Z value:

    - R00001 (P21/n, Z=4): 32 atoms/asym × 4 symm = 128 atoms/cell -> 500 mol
    - R00002 (Z=8): 32 × 8 = 256 atoms/cell -> 1000 mol
    - R00004 (P-1, Z=2): 32 × 2 = 64 atoms/cell -> 250 mol
    """
    if not os.path.isfile(cif_path):
        pytest.skip(f"fixture missing: {cif_path}")

    from abmptools.crystal.cif_engine_ase import run_ase

    super_atoms, molecules = run_ase(
        cif=cif_path,
        layer=layer,
        atoms_in_mol=[32],
    )
    assert len(super_atoms) == 32 * n_mol_expected
    assert len(molecules) == n_mol_expected


# ---------------------------------------------------------------------------
# write_pdb_for_abmp / write_xyz_for_abmp / run_ase_pipeline (Phase D-2)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_write_pdb_for_abmp_layout(tmp_path):
    """The HETATM records must match the legacy column layout
    (HETATM serial atom_name UNK res_seq x y z occ temp element)."""
    from abmptools.crystal.cif_engine_ase import (
        detect_molecules,
        expand_supercell,
        read_cif_to_atoms,
        write_pdb_for_abmp,
    )
    atoms = read_cif_to_atoms(R00001_CIF)
    super_atoms = expand_supercell(atoms, layer=2)  # 8 cells = 32 mol
    mols = detect_molecules(super_atoms, atoms_in_mol=[32])
    pdb = tmp_path / "ase_out.pdb"
    write_pdb_for_abmp(super_atoms, mols, molname="UNK", output_pdb=str(pdb))

    text = pdb.read_text()
    hetatm = [l for l in text.splitlines() if l.startswith("HETATM")]
    assert len(hetatm) == 32 * 32  # 8 cells * 4 mol/cell * 32 atoms

    line0 = hetatm[0]
    # Column-by-column sanity (PDB v3.3):
    assert line0[0:6] == "HETATM"
    # Serial number (right-justified 5-char field).
    assert int(line0[6:11].strip()) == 1
    # Residue name UNK at columns 17-19.
    assert line0[17:20] == "UNK"
    # Residue sequence (1-based, right-justified 4 chars at 22-25).
    assert int(line0[22:26].strip()) == 1
    # x is %8.3f, columns 30-37.
    float(line0[30:38])  # raises if malformed


@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_write_xyz_for_abmp_full_precision(tmp_path):
    """The XYZ writer must preserve the ASE coordinates verbatim
    (full precision, no %8.3f truncation)."""
    from abmptools.crystal.cif_engine_ase import (
        detect_molecules,
        expand_supercell,
        read_cif_to_atoms,
        write_xyz_for_abmp,
    )
    atoms = read_cif_to_atoms(R00001_CIF)
    super_atoms = expand_supercell(atoms, layer=2)
    mols = detect_molecules(super_atoms, atoms_in_mol=[32])
    xyz_path = tmp_path / "ase_out.xyz"
    write_xyz_for_abmp(super_atoms, mols, output_xyz=str(xyz_path))
    text = xyz_path.read_text().splitlines()
    assert int(text[0]) == 32 * 32
    # Third line is the first atom; must contain a coordinate with > 4
    # decimals (full precision, not %.3f).
    first_atom = text[2].split()
    assert len(first_atom) == 4  # element + 3 coords
    # Heuristic: at least one coordinate should have > 4 decimal digits.
    decimal_lengths = [
        len(c.split(".")[1]) if "." in c else 0 for c in first_atom[1:]
    ]
    assert max(decimal_lengths) > 4, (
        f"XYZ should preserve full precision; got decimals "
        f"{decimal_lengths} in {first_atom!r}"
    )


@pytest.mark.slow
@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_run_ase_pipeline_emits_pdb_and_xyz(tmp_path):
    """Full Phase D-2 pipeline: CIF → cifout/layer<L>/{pdb,xyz}/."""
    import shutil as _shutil
    from abmptools.crystal.cif_engine_ase import run_ase_pipeline

    _shutil.copy(R00001_CIF, tmp_path / "XXXI-MMFF-R00001.cif")
    layer_dir = run_ase_pipeline(
        cif="XXXI-MMFF-R00001.cif",
        layer=5,
        atoms_in_mol=[32],
        cwd=str(tmp_path),
    )
    assert layer_dir.is_dir()
    pdb = layer_dir / "pdb" / "XXXI-MMFF-R00001layer5Zp1.pdb"
    xyz = layer_dir / "xyz" / "XXXI-MMFF-R00001layer5Zp1.xyz"
    assert pdb.is_file() and xyz.is_file()
    # 16000 atoms expected (P21/n Z=4, layer 5 = 125 cells).
    assert int(xyz.read_text().splitlines()[0]) == 16000
    hetatm = [l for l in pdb.read_text().splitlines() if l.startswith("HETATM")]
    assert len(hetatm) == 16000


@pytest.mark.slow
@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_orchestrator_engine_ase_csp7_r00001(tmp_path):
    """``CrystalOrchestrator`` with ``engine='ase'`` runs the full
    pipeline (CIF → for_abmp/*.{ajf,pdb}) without raising. The output
    is *not* byte-equivalent to the legacy fixture (different supercell
    origin and atom ordering), but it must be a valid AJF (Natom > 0
    matching the &XYZ block size)."""
    import re
    import shutil as _shutil
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFEngineConfig,
        CIFInputSpec,
        CrystalBuildConfig,
        FMOMethod,
        FragmentTemplate,
        HPCJobSpec,
    )

    UNK_AJF = os.path.abspath(os.path.join(
        TESTS_DIR, os.pardir, "sample", "crystal", "csp7_smoke", "UNK.ajf"
    ))
    if not os.path.isfile(UNK_AJF):
        pytest.skip("UNK.ajf template missing")

    _shutil.copy(R00001_CIF, tmp_path / "XXXI-MMFF-R00001.cif")
    cfg = CrystalBuildConfig(
        project_name="ase_engine_smoke",
        output_dir=str(tmp_path / "out"),
        inputs=[CIFInputSpec(
            cif=str(tmp_path / "XXXI-MMFF-R00001.cif"),
            layer=5, atoms_in_mol=[32],
        )],
        cif_engine=CIFEngineConfig(engine="ase"),
        fragment=FragmentTemplate(
            cutmode="around", solutes=[0], criteria=6.0,
            molname=["UNK"], pieda=True, cmm=True, getmode="rfile",
            template_ajf=UNK_AJF,
        ),
        fmo=FMOMethod(
            method="MP2", basis_set="6-31Gdag", memory=6000,
            cpfflag=True, abinit_ver="rev23", npro=1, is_xyz=True,
        ),
        hpc=HPCJobSpec(
            scheduler="local", abinit_dir="/tmp", binary_name="dummy",
            nodes=1, proc_per_node=1, omp_threads=1,
        ),
    )
    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))
    summary = orch.run()
    assert summary["engine"] == "ase"
    assert summary["n_ajf_total"] == 1

    # The for_abmp/*.ajf is the FMO calculation input we want to inspect;
    # the template UNK.ajf staged in the pdb/ directory must be skipped.
    ajf = next((tmp_path / "out").rglob("for_abmp/*.ajf"))
    body = ajf.read_text()
    natom = int(re.search(r"Natom=(\d+)", body).group(1))
    assert natom > 0, "Natom should be > 0 in ASE-emitted AJF"
    xyz_block = body.split("&XYZ")[1].split("/")[0].strip().splitlines()
    assert len(xyz_block) == natom, (
        f"&XYZ block has {len(xyz_block)} atoms but Natom={natom}"
    )


# ---------------------------------------------------------------------------
# unwrap_molecules — PBC boundary fix for Atoms.repeat
# ---------------------------------------------------------------------------


def test_unwrap_molecules_keeps_molecules_intact():
    """A diatomic ``A-B`` straddling the cell boundary must, after
    unwrap, sit contiguously in real space (max intramolecular
    distance = bond length, not 1 cell - bond length).
    """
    pytest.importorskip("ase")
    from ase import Atoms
    from abmptools.crystal.cif_engine_ase import (
        detect_molecules,
        expand_supercell,
        unwrap_molecules,
    )

    # Two N atoms across the x-boundary: 0.95 and 0.05 of a 5 Å cell.
    # Without unwrap: nearest-image distance is 1.0 Å (bond), but in
    # real-space coordinates atoms are 4.0 Å apart.
    atoms = Atoms(
        symbols="NN",
        positions=[(4.75, 2.5, 2.5), (0.25, 2.5, 2.5)],  # spans the x face
        cell=[5.0, 5.0, 5.0],
        pbc=True,
    )
    super_atoms = expand_supercell(atoms, layer=2)
    molecules = detect_molecules(super_atoms, atoms_in_mol=[2], bond_tolerance=0.4)
    # 2x2x2 supercell -> 8 N2 molecules.
    assert len(molecules) == 8

    pos_before = super_atoms.get_positions().copy()
    # Pre-unwrap: at least one molecule has its two atoms at opposite
    # faces of the supercell.
    spans_before = []
    for mol in molecules:
        d = abs(pos_before[mol[0]] - pos_before[mol[1]])
        spans_before.append(max(d))
    assert max(spans_before) > 4.0, "test setup expected pre-unwrap split"

    unwrap_molecules(super_atoms, molecules, bond_tolerance=0.4)
    pos_after = super_atoms.get_positions()

    # Post-unwrap: every molecule's two atoms are within ~1 Å (bond),
    # because find_mic returned the minimum-image vector.
    for mol in molecules:
        import numpy as np
        d = np.linalg.norm(pos_after[mol[0]] - pos_after[mol[1]])
        assert d < 1.5, f"unwrap failed for molecule {mol}: dist={d:.3f}"


@pytest.mark.skipif(not os.path.isfile(R00001_CIF), reason="R00001 fixture missing")
def test_run_ase_unwraps_molecules_for_csp7():
    """run_ase calls unwrap_molecules under the hood; csp7 layer 2
    should leave each 32-atom molecule contiguous (max C-C / C-O
    intramolecular distance < ~12 Å, not the 18+ Å that boundary
    splits would produce in a 9 Å unit cell)."""
    pytest.importorskip("ase")
    import numpy as np
    from abmptools.crystal.cif_engine_ase import run_ase

    super_atoms, molecules = run_ase(
        cif=R00001_CIF, layer=2, atoms_in_mol=[32],
    )
    pos = super_atoms.get_positions()
    max_intramol = 0.0
    for mol in molecules:
        coords = pos[mol]
        # naive O(N^2) span calculation; cheap for 32 atoms.
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                d = np.linalg.norm(coords[i] - coords[j])
                if d > max_intramol:
                    max_intramol = d
    # csp7 R00001 molecules are about 6 Å end-to-end; allowing 15 Å
    # rules out a wraparound contribution (would otherwise add the
    # ~10 Å cell vector).
    assert max_intramol < 15.0, (
        f"intramolecular distance {max_intramol:.2f} Å suggests boundary split"
    )
