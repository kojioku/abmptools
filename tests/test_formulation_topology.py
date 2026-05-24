# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.topology."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from abmptools.formulation import topology as ftop
from abmptools.formulation.models import (
    BileSaltSpec,
    EnhancerSpec,
    FormulationBuildConfig,
    PeptideSpec,
    SystemSpec,
)
from abmptools.formulation.small_molecule import SmallMoleculeResult


def _make_config() -> FormulationBuildConfig:
    return FormulationBuildConfig(
        system=SystemSpec(
            peptides=[PeptideSpec(name="p", sequence="AAAAA", n_copies=2)],
            enhancers=[EnhancerSpec(
                name="caprate", resname="CPR",
                smiles_neutral="CCCCCCCCCC(=O)O", n_neutral=4,
            )],
            bile_salts=[BileSaltSpec(
                name="tch", resname="TCH",
                smiles="O", n_copies=2, net_charge=-1,
            )],
        ),
    )


def _make_sm_result(tmp_path: Path, resname: str) -> SmallMoleculeResult:
    d = tmp_path / f"{resname}.acpype"
    d.mkdir(parents=True, exist_ok=True)
    frcmod = d / f"{resname}_AC.frcmod"
    mol2 = d / f"{resname}_bcc_gaff2.mol2"
    itp = d / f"{resname}_GMX.itp"
    gro = d / f"{resname}_GMX.gro"
    for p in (frcmod, mol2, itp, gro):
        p.write_text("")
    pdb = tmp_path / f"{resname}.pdb"
    pdb.write_text("HETATM    1  C   "
                   f"{resname:>3s} A   1       0.000   0.000   0.000\nEND\n")
    return SmallMoleculeResult(
        resname=resname,
        monomer_pdb=pdb,
        acpype_dir=d, frcmod=frcmod, mol2=mol2, itp=itp, gro=gro,
        n_atoms=1, net_charge=0,
    )


# ---------------------------------------------------------------------------
# render_unified_tleap_input
# ---------------------------------------------------------------------------


def test_tleap_input_sources_ff14sb(tmp_path):
    cfg = _make_config()
    sm = _make_sm_result(tmp_path, "CPR")
    text = ftop.render_unified_tleap_input(
        config=cfg,
        mixture_pdb="/tmp/mix.pdb",
        small_molecules=[sm],
        solvatebox_margin_nm=0.5,
        target_nacl_molarity=0.15,
        prmtop_out="/tmp/system.prmtop",
        inpcrd_out="/tmp/system.inpcrd",
    )
    assert "source leaprc.protein.ff14SB" in text
    assert "source leaprc.water.tip3p" in text


def test_tleap_input_loads_each_ligand(tmp_path):
    cfg = _make_config()
    sm_cpr = _make_sm_result(tmp_path, "CPR")
    sm_tch = _make_sm_result(tmp_path, "TCH")
    text = ftop.render_unified_tleap_input(
        config=cfg, mixture_pdb="/tmp/mix.pdb",
        small_molecules=[sm_cpr, sm_tch],
        solvatebox_margin_nm=0.5, target_nacl_molarity=0.15,
        prmtop_out="/tmp/x", inpcrd_out="/tmp/y",
    )
    assert f"loadAmberParams {sm_cpr.frcmod}" in text
    assert f"CPR = loadmol2 {sm_cpr.mol2}" in text
    assert f"loadAmberParams {sm_tch.frcmod}" in text
    assert f"TCH = loadmol2 {sm_tch.mol2}" in text


def test_tleap_input_loadpdb_mixture(tmp_path):
    cfg = _make_config()
    text = ftop.render_unified_tleap_input(
        config=cfg, mixture_pdb="/tmp/mix.pdb",
        small_molecules=[],
        solvatebox_margin_nm=0.5, target_nacl_molarity=0.15,
        prmtop_out="/tmp/x", inpcrd_out="/tmp/y",
    )
    assert "sys = loadpdb /tmp/mix.pdb" in text


def test_tleap_input_solvatebox_in_angstrom(tmp_path):
    cfg = _make_config()
    text = ftop.render_unified_tleap_input(
        config=cfg, mixture_pdb="/tmp/x.pdb", small_molecules=[],
        solvatebox_margin_nm=1.0, target_nacl_molarity=0.0,
        prmtop_out="/tmp/x", inpcrd_out="/tmp/y",
    )
    # 1.0 nm = 10.00 Å
    assert "solvatebox sys TIP3PBOX 10.00" in text


def test_tleap_input_neutralize_default(tmp_path):
    cfg = _make_config()
    text = ftop.render_unified_tleap_input(
        config=cfg, mixture_pdb="/tmp/x.pdb", small_molecules=[],
        solvatebox_margin_nm=0.5, target_nacl_molarity=0.15,
        prmtop_out="/tmp/x", inpcrd_out="/tmp/y",
    )
    assert "addions sys Na+ 0" in text
    assert "addions sys Cl- 0" in text


def test_tleap_input_save_paths(tmp_path):
    cfg = _make_config()
    text = ftop.render_unified_tleap_input(
        config=cfg, mixture_pdb="/tmp/x.pdb", small_molecules=[],
        solvatebox_margin_nm=0.5, target_nacl_molarity=0.0,
        prmtop_out="/build/sys.prmtop", inpcrd_out="/build/sys.inpcrd",
    )
    assert "saveAmberParm sys /build/sys.prmtop /build/sys.inpcrd" in text


def test_tleap_input_extra_bond_directives(tmp_path):
    cfg = _make_config()
    text = ftop.render_unified_tleap_input(
        config=cfg, mixture_pdb="/tmp/x.pdb", small_molecules=[],
        solvatebox_margin_nm=0.5, target_nacl_molarity=0.0,
        prmtop_out="/tmp/p", inpcrd_out="/tmp/c",
        extra_bonds=["bond sys.6.SG sys.11.SG"],
    )
    assert "bond sys.6.SG sys.11.SG" in text


# ---------------------------------------------------------------------------
# build_topology orchestration (mocked subprocess + parmed)
# ---------------------------------------------------------------------------


def test_build_topology_runs_tleap_and_parmed(tmp_path):
    cfg = _make_config()
    cfg.output_dir = str(tmp_path)
    sm = _make_sm_result(tmp_path, "CPR")
    mix = tmp_path / "mix.pdb"
    mix.write_text("END\n")

    def fake_run(argv, **kwargs):
        # Touch prmtop + inpcrd to simulate tleap success
        wd = Path(kwargs["cwd"])
        (wd / "system.prmtop").write_text("PRMTOP")
        (wd / "system.inpcrd").write_text("INPCRD")
        from subprocess import CompletedProcess
        return CompletedProcess(argv, returncode=0, stdout="", stderr="")

    with patch("abmptools.formulation.topology.run_command", side_effect=fake_run), \
         patch("abmptools.formulation.topology.amber_to_gromacs") as mock_a2g:
        # Simulate parmed writing files
        def fake_a2g(prmtop, inpcrd, top_out, gro_out):
            Path(top_out).write_text("[ defaults ]\n")
            Path(gro_out).write_text("0\n")
        mock_a2g.side_effect = fake_a2g

        art = ftop.build_topology(
            config=cfg, mixture_pdb=str(mix),
            small_molecules=[sm], workdir=str(tmp_path / "build"),
        )
        assert Path(art.prmtop).is_file()
        assert Path(art.top).is_file()
        assert Path(art.gro).is_file()
        assert Path(art.tleap_input).is_file()
        assert Path(art.tleap_log).is_file()


def test_build_topology_raises_if_tleap_missing_outputs(tmp_path):
    cfg = _make_config()
    sm = _make_sm_result(tmp_path, "CPR")
    mix = tmp_path / "mix.pdb"
    mix.write_text("END\n")

    def fake_run_nofile(argv, **kwargs):
        from subprocess import CompletedProcess
        return CompletedProcess(argv, returncode=0, stdout="", stderr="")

    with patch("abmptools.formulation.topology.run_command", side_effect=fake_run_nofile):
        with pytest.raises(RuntimeError, match="missing"):
            ftop.build_topology(
                config=cfg, mixture_pdb=str(mix), small_molecules=[sm],
                workdir=str(tmp_path / "b"),
            )
