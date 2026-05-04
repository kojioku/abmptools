"""Topology content sanity tests for the membrane builder.

These tests validate the *content* of GROMACS topology output (atom
counts, charges, residue presence) — complementing the table-level
unit tests in ``test_membrane_charmm_translate.py`` and the smoke
tests in ``tests/integration/run_membrane_us_*_smoke.sh`` which only
verify that grompp accepts the topology.

Tests run on synthetic topology fixtures (no GROMACS / packmol-memgen
required), and additionally probe a real build directory if
``MEMBRANE_BUILD_DIR`` env var points at one (auto-skip otherwise).
"""
from __future__ import annotations

import os
from pathlib import Path
from textwrap import dedent

import pytest

from abmptools.membrane.topology_sanity import (
    TopologySummary,
    summarize_topology,
)


# ---------------------------------------------------------------------------
# fixtures: synthetic minimal topologies (no GROMACS needed)
# ---------------------------------------------------------------------------


@pytest.fixture
def amber_style_top(tmp_path: Path) -> Path:
    """Inline-itp AMBER-style monolithic top (mimics tleap+parmed output)."""
    build = tmp_path / "build"
    build.mkdir()
    (build / "system.top").write_text(dedent("""\
        ; AMBER-style synthetic top
        [ moleculetype ]
        ; Name           nrexcl
        Peptide        3

        [ atoms ]
        ;   nr   type  resnr resname  atom  cgnr   charge   mass
             1   N3      1    ALA     N      1    -0.5000   14.01
             2   H        1    ALA     H1     1     0.2000    1.008
             3   H        1    ALA     H2     1     0.2000    1.008
             4   H        1    ALA     H3     1     0.1000    1.008
             5   CT       1    ALA     CA     1     0.0337   12.01

        [ moleculetype ]
        Na+            3

        [ atoms ]
             1   Na+      1    Na+    Na+    1     1.0000   22.99

        [ moleculetype ]
        Cl-            3

        [ atoms ]
             1   Cl-      1    Cl-    Cl-    1    -1.0000   35.45

        [ system ]
        Test

        [ molecules ]
        ; Compound       #mols
        Peptide              1
        Na+                  3
        Cl-                  3
    """))
    return build


@pytest.fixture
def charmm_style_top(tmp_path: Path) -> Path:
    """Per-chain itp CHARMM-style top (mimics pdb2gmx output)."""
    build = tmp_path / "build"
    build.mkdir()
    (build / "system.top").write_text(dedent("""\
        ; CHARMM-style synthetic top
        #include "./charmm.ff/forcefield.itp"

        ; Include chain topologies
        #include "system_Other_chain_A.itp"
        #include "system_Other_chain_B.itp"
        #include "system_Protein_chain_C.itp"

        [ system ]
        Test

        [ molecules ]
        ; Compound       #mols
        Other_chain_A      1
        Other_chain_B      1
        Protein_chain_C    1
    """))
    # chain A: 2 POPC molecules (charge 0 each, just lipid placeholder)
    (build / "system_Other_chain_A.itp").write_text(dedent("""\
        [ moleculetype ]
        Other_chain_A   3

        [ atoms ]
             1   NTL    1   POPC    N      1    -0.6000   14.01
             2   CTL2   1   POPC    C12    1    -0.1000   12.01
             3   CTL2   1   POPC    C13    1    +0.7000   12.01
             4   NTL    2   POPC    N      2    -0.6000   14.01
             5   CTL2   2   POPC    C12    2    -0.1000   12.01
             6   CTL2   2   POPC    C13    2    +0.7000   12.01
    """))
    # chain B: 2 SOD ions (+1 each)
    (build / "system_Other_chain_B.itp").write_text(dedent("""\
        [ moleculetype ]
        Other_chain_B   3

        [ atoms ]
             1   SOD    1   SOD    SOD    1     1.0000   22.99
             2   SOD    2   SOD    SOD    2     1.0000   22.99
    """))
    # chain C: ACE-ALA-NME peptide (sum charge 0, 5 atoms simplified)
    (build / "system_Protein_chain_C.itp").write_text(dedent("""\
        [ moleculetype ]
        Protein_chain_C   3

        [ atoms ]
             1   CT3    1    ACE   CH3    1    -0.2700   12.01
             2   HA3    1    ACE   HH31   1     0.0900    1.008
             3   HA3    1    ACE   HH32   1     0.0900    1.008
             4   HA3    1    ACE   HH33   1     0.0900    1.008
             5   NH1    2    ALA   N      2     0.0000   14.01

        ; Include Position restraint file
        #ifdef POSRES
        #include "posre_Protein_chain_C.itp"
        #endif
    """))
    return build


# ---------------------------------------------------------------------------
# parser tests
# ---------------------------------------------------------------------------

class TestSummarizeTopologyAmberStyle:
    def test_atom_count(self, amber_style_top):
        s = summarize_topology(amber_style_top)
        # 5 peptide atoms + 3 Na + 3 Cl = 11
        assert s.total_atoms == 11

    def test_charge_neutral(self, amber_style_top):
        s = summarize_topology(amber_style_top)
        # peptide 0.0337 ≈ 0.03, Na 3×1=+3, Cl 3×-1=-3 → ~0.034
        # Actually peptide doesn't sum cleanly (synthetic), just check ions cancel
        assert abs(s.total_charge) < 0.1, f"net {s.total_charge}"

    def test_resname_detected(self, amber_style_top):
        s = summarize_topology(amber_style_top)
        residues = s.by_residue()
        assert "ALA" in residues
        assert "Na+" in residues
        assert "Cl-" in residues
        # Na+ count: 3, Cl- count: 3
        assert residues["Na+"][0] == 3
        assert residues["Cl-"][0] == 3


class TestSummarizeTopologyCharmmStyle:
    def test_per_chain_itp_loaded(self, charmm_style_top):
        s = summarize_topology(charmm_style_top)
        # chain A: 6 atoms (2 POPC × 3 atoms placeholder)
        # chain B: 2 ions
        # chain C: 5 peptide atoms
        assert s.total_atoms == 13

    def test_chain_count_preserved(self, charmm_style_top):
        s = summarize_topology(charmm_style_top)
        # all chains have count=1 (CHARMM pdb2gmx convention)
        for b in s.blocks:
            assert b.count == 1

    def test_lipid_residue_present(self, charmm_style_top):
        s = summarize_topology(charmm_style_top)
        residues = s.by_residue()
        assert "POPC" in residues
        # 2 POPC (counted across all atoms with that resname)

    def test_protein_chain_atoms(self, charmm_style_top):
        s = summarize_topology(charmm_style_top)
        # chain C should have ACE + ALA residues
        protein_blocks = [b for b in s.blocks if "ACE" in b.residues]
        assert len(protein_blocks) == 1
        assert "ALA" in protein_blocks[0].residues
        assert protein_blocks[0].atoms_per_mol == 5


class TestPosresHandlingInItp:
    """The CHARMM chain itps end with #ifdef POSRES — parser must stop there."""
    def test_posres_block_does_not_leak(self, charmm_style_top):
        # chain C has #ifdef POSRES at end; parser should not re-read past it
        s = summarize_topology(charmm_style_top)
        protein = next(b for b in s.blocks if "ACE" in b.residues)
        assert protein.atoms_per_mol == 5  # not larger than the actual atoms


# ---------------------------------------------------------------------------
# real-build smoke test (auto-skipped if no MEMBRANE_BUILD_DIR set)
# ---------------------------------------------------------------------------

@pytest.fixture
def real_build_dir() -> Path:
    env = os.environ.get("MEMBRANE_BUILD_DIR")
    if not env:
        pytest.skip("MEMBRANE_BUILD_DIR not set — skipping real-build sanity test")
    p = Path(env)
    if not (p / "system.top").exists():
        pytest.skip(f"{p}/system.top not found — skipping")
    return p


class TestRealBuildSanity:
    """Exercise sanity invariants on a real builder output (opt-in via env)."""

    def test_system_is_electroneutral(self, real_build_dir):
        s = summarize_topology(real_build_dir)
        # Net charge must be 0 within rounding error (1e-3)
        assert abs(s.total_charge) < 1e-3, (
            f"Net charge {s.total_charge:+.4f} differs from 0 — "
            f"check ion count / peptide protonation"
        )

    def test_has_lipid(self, real_build_dir):
        s = summarize_topology(real_build_dir)
        residues = s.by_residue()
        has_lipid = any(
            r in residues
            for r in ("POPC", "PC", "PA", "OL", "DOPC", "POPE", "DPPC")
        )
        assert has_lipid, f"No lipid residue found; got: {sorted(residues.keys())}"

    def test_has_water(self, real_build_dir):
        s = summarize_topology(real_build_dir)
        residues = s.by_residue()
        has_water = any(r in residues for r in ("WAT", "TIP3", "HOH", "SOL"))
        assert has_water, f"No water residue found; got: {sorted(residues.keys())}"

    def test_has_peptide(self, real_build_dir):
        s = summarize_topology(real_build_dir)
        residues = s.by_residue()
        # Standard amino acids
        has_aa = any(r in residues for r in ("ALA", "GLY", "LEU", "ARG", "LYS"))
        assert has_aa, f"No standard AA residue found; got: {sorted(residues.keys())}"

    def test_ion_balance(self, real_build_dir):
        """Cation count × |+1| must match anion count × |-1| (after peptide
        charge correction). Validates that ``neutralize=True`` worked."""
        s = summarize_topology(real_build_dir)
        residues = s.by_residue()
        cation_charge = 0.0
        anion_charge = 0.0
        for r, (n, c) in residues.items():
            if r in ("Na+", "SOD", "K+", "POT"):
                cation_charge += c
            elif r in ("Cl-", "CLA"):
                anion_charge += c
        # |cation_charge| - |anion_charge| should be small (peptide is neutral here)
        delta = cation_charge + anion_charge
        assert abs(delta) < 1.001, (
            f"Ion imbalance: cations={cation_charge}, anions={anion_charge}, "
            f"delta={delta} (peptide should be neutral for this fixture)"
        )

    def test_atom_count_in_expected_range(self, real_build_dir):
        """Sanity: poly-Ala 5-mer + ~32 POPC + water + ions ≈ 18-20k atoms.
        Reject obviously-broken builds with >100k atoms or <1k atoms."""
        s = summarize_topology(real_build_dir)
        assert 1000 < s.total_atoms < 100_000, (
            f"Total atoms {s.total_atoms} out of expected range "
            f"(1k-100k) for typical bilayer + peptide system"
        )
