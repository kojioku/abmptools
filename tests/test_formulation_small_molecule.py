# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.small_molecule."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from abmptools.formulation import small_molecule as sm


def _make_acpype_artifacts(acp_dir: Path, basename: str) -> None:
    """Lay down minimal acpype outputs for `_locate_gmx_outputs` + upstream."""
    acp_dir.mkdir(parents=True, exist_ok=True)
    (acp_dir / f"{basename}_AC.frcmod").write_text("frcmod")
    (acp_dir / f"{basename}_bcc_gaff2.mol2").write_text("mol2")
    (acp_dir / f"{basename}_GMX.itp").write_text("[ moleculetype ]\n")
    (acp_dir / f"{basename}_GMX.gro").write_text("LIG\n  3\n")


def test_locate_gmx_outputs_finds_itp_and_gro(tmp_path):
    _make_acpype_artifacts(tmp_path / "x.acpype", "x")
    itp, gro = sm._locate_gmx_outputs(tmp_path / "x.acpype")
    assert itp.name == "x_GMX.itp"
    assert gro.name == "x_GMX.gro"


def test_locate_gmx_outputs_missing_itp_raises(tmp_path):
    d = tmp_path / "y.acpype"
    d.mkdir()
    with pytest.raises(FileNotFoundError, match="No \\*_GMX.itp"):
        sm._locate_gmx_outputs(d)


def test_smiles_to_pdb_with_rdkit(tmp_path):
    """Real RDKit test (not mocked) — RDKit is in abmptoolsenv."""
    out = tmp_path / "etoh.pdb"
    sm.smiles_to_pdb("CCO", resname="ETO", out_pdb=out, seed=7)
    assert out.is_file()
    text = out.read_text()
    # 1-letter resname padded to 3 chars
    assert "ETO" in text


def test_smiles_to_pdb_invalid(tmp_path):
    with pytest.raises(ValueError, match="failed to parse SMILES"):
        sm.smiles_to_pdb("not-a-smiles!!!", resname="XXX",
                          out_pdb=tmp_path / "x.pdb")


def test_parameterize_small_molecule_requires_input(tmp_path):
    with pytest.raises(ValueError, match="smiles or pdb_path"):
        sm.parameterize_small_molecule(
            name="x", resname="XXX", workdir=tmp_path,
        )


def test_parameterize_small_molecule_calls_acpype_with_smiles(tmp_path):
    fake_acpype_dir = tmp_path / "CPR.acpype"
    fake_acpype_dir.mkdir()
    _make_acpype_artifacts(fake_acpype_dir, "CPR")

    with patch.object(sm, "_genesis_run_acpype") as mock_acp:
        from abmptools.genesis.mmgbsa.ligand_parameterize import AcpypeResult
        mock_acp.return_value = AcpypeResult(
            acpype_dir=fake_acpype_dir,
            frcmod=fake_acpype_dir / "CPR_AC.frcmod",
            mol2=fake_acpype_dir / "CPR_bcc_gaff2.mol2",
        )
        result = sm.parameterize_small_molecule(
            name="caprate", resname="CPR",
            smiles="CCCCCCCCCC(=O)O", net_charge=0,
            workdir=tmp_path,
        )
        # acpype was called
        assert mock_acp.called
        # monomer pdb was actually written by smiles_to_pdb
        assert result.monomer_pdb.is_file()
        # discovered GMX outputs
        assert result.itp.name == "CPR_GMX.itp"
        assert result.gro.name == "CPR_GMX.gro"
        assert result.resname == "CPR"
        assert result.net_charge == 0


def test_parameterize_small_molecule_from_pdb_skips_rdkit(tmp_path):
    """When pdb_path is given, skip SMILES→3D embedding."""
    src = tmp_path / "src.pdb"
    src.write_text("HETATM    1  C   XYZ A   1       0.000   0.000   0.000\nEND\n")

    fake_acp = tmp_path / "XYZ.acpype"
    fake_acp.mkdir()
    _make_acpype_artifacts(fake_acp, "XYZ")

    with patch.object(sm, "_genesis_run_acpype") as mock_acp, \
         patch.object(sm, "smiles_to_pdb") as mock_smiles:
        from abmptools.genesis.mmgbsa.ligand_parameterize import AcpypeResult
        mock_acp.return_value = AcpypeResult(
            acpype_dir=fake_acp,
            frcmod=fake_acp / "XYZ_AC.frcmod",
            mol2=fake_acp / "XYZ_bcc_gaff2.mol2",
        )
        result = sm.parameterize_small_molecule(
            name="xyz", resname="XYZ", pdb_path=str(src),
            workdir=tmp_path,
        )
        mock_smiles.assert_not_called()
        assert result.monomer_pdb.is_file()
