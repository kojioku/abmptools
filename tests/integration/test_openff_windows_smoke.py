# -*- coding: utf-8 -*-
"""Windows-native OpenFF route smoke test (Phase 2-C).

formulation の OpenFF route (`force_field_route="openff"`) が **全 OS** で
動くことの最小検証。 tleap/acpype (AmberTools、 Windows native 不可) を一切
使わず、 PDBFixer + OpenFF ``Topology.from_pdb`` + ff14SB SMIRNOFF だけで
peptide を parametrize できることを確認する。

CI では ``windows-latest`` runner 上でこのファイルだけを走らせ、
「Windows native build 可能」を実機で確証する (.github/workflows/windows.yml)。

heavy deps (openff-toolkit / openff-interchange / openff-amber-ff-ports /
pdbfixer / PeptideBuilder) が無い環境では skip。 packmol / gmx は使わない
(parametrization 層だけを検証 — それらは外部ツールで別途用意される)。
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

# --- gating: OpenFF route の Windows native 依存が全部揃っているか ---
pytest.importorskip("openff.toolkit", reason="openff-toolkit not installed")
pytest.importorskip("openff.interchange", reason="openff-interchange not installed")
pytest.importorskip("openff.amber_ff_ports", reason="openff-amber-ff-ports not installed")
pytest.importorskip("pdbfixer", reason="pdbfixer not installed")
pytest.importorskip("openmm", reason="openmm not installed")
pytest.importorskip("PeptideBuilder", reason="PeptideBuilder not installed")

pytestmark = pytest.mark.slow

FF14SB = "ff14sb_off_impropers_0.0.4.offxml"
SAGE = "openff_unconstrained-2.1.0.offxml"
TIP3P = "tip3p.offxml"


def test_sequence_to_ff14sb_interchange_windows_native(tmp_path):
    """sequence → PeptideBuilder 3D → PDBFixer → Topology.from_pdb → ff14SB
    SMIRNOFF Interchange。 AmberTools 非依存の全経路を 1 本で検証。
    """
    from abmptools.formulation.peptide_atomistic_openff import (
        _build_peptide_from_sequence_openff,
        _pdbfix_protein,
    )
    from openff.toolkit import ForceField, Topology
    from openff.interchange import Interchange
    from openff.units import unit
    import numpy as np

    # 1. sequence → extended chain (PeptideBuilder、 natural L-AA)
    raw_pdb = _build_peptide_from_sequence_openff(
        sequence="GAG", name="tri", output_dir=tmp_path,
    )
    assert Path(raw_pdb).is_file()

    # 2. PDBFixer: explicit H 付加 (OpenFF 必須)
    fixed_pdb = _pdbfix_protein(str(raw_pdb), str(tmp_path / "tri_fixed.pdb"))
    assert Path(fixed_pdb).is_file()

    # 3. Topology.from_pdb (Molecule.from_polymer_pdb でなく)
    top = Topology.from_pdb(fixed_pdb)
    assert top.n_molecules == 1
    mol = top.molecule(0)
    assert mol.n_atoms > 0

    # 4. ff14SB SMIRNOFF + Sage で Interchange (library charges、 sqm 不要)
    single = Topology()
    single.add_molecule(mol)
    single.box_vectors = np.eye(3) * 4.0 * unit.nanometer
    ff = ForceField(FF14SB, SAGE, TIP3P)
    ic = Interchange.from_smirnoff(force_field=ff, topology=single)
    assert ic.topology.n_atoms == mol.n_atoms
    # library charges が乗っていること (全 0 でない)
    charges = ic["Electrostatics"].charges
    assert len(charges) == mol.n_atoms
    assert any(abs(c.m) > 1e-6 for c in charges.values())


def test_pdbfixer_strips_water_adds_hydrogens(tmp_path):
    """PDBFixer が water 除去 + H 付加することを最小確認 (Windows native 前処理)。"""
    from abmptools.formulation.peptide_atomistic_openff import (
        _build_peptide_from_sequence_openff,
        _pdbfix_protein,
    )

    raw_pdb = _build_peptide_from_sequence_openff(
        sequence="AA", name="di", output_dir=tmp_path,
    )
    fixed = _pdbfix_protein(str(raw_pdb), str(tmp_path / "di_fixed.pdb"))
    text = Path(fixed).read_text()
    # H atom が付加されている (PeptideBuilder の bare 出力に H が増える)
    n_h = sum(1 for ln in text.splitlines()
              if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip().startswith("H"))
    assert n_h > 0, "PDBFixer should have added explicit hydrogens"
    # HOH が残っていない
    assert "HOH" not in text
