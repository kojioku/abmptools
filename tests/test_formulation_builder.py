# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.builder (subprocess mocked)."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from abmptools.formulation.builder import (
    FormulationArtifacts,
    FormulationBuilder,
    render_run_script,
)
from abmptools.formulation.models import (
    BileSaltSpec,
    EnhancerSpec,
    EquilibrationProtocol,
    FormulationBuildConfig,
    PeptideSpec,
    ProductionProtocol,
    SystemSpec,
)
from abmptools.formulation.small_molecule import SmallMoleculeResult
from abmptools.formulation.topology import TopologyArtifacts


def _make_config(output_dir: Path) -> FormulationBuildConfig:
    return FormulationBuildConfig(
        system=SystemSpec(
            peptides=[PeptideSpec(
                name="ala5", sequence="AAAAA", n_copies=2,
                cap_n="ACE", cap_c="NME",
            )],
            enhancers=[EnhancerSpec(
                name="caprate", resname="CPR",
                smiles_neutral="CCCCCCCCCC(=O)O",
                smiles_charged="CCCCCCCCCC(=O)[O-]",
                n_neutral=4, n_charged=4,
            )],
            bile_salts=[BileSaltSpec(
                name="tch", resname="TCH",
                smiles="O", n_copies=2, net_charge=-1,
            )],
            box_size_nm=5.0,
        ),
        equilibration=EquilibrationProtocol(em_steps=100, nvt_nsteps=100, npt_nsteps=100),
        production=ProductionProtocol(nsteps=200),
        output_dir=str(output_dir),
    )


def _make_sm_result(workdir: Path, resname: str) -> SmallMoleculeResult:
    d = workdir / f"{resname}.acpype"
    d.mkdir(parents=True, exist_ok=True)
    for ext in (
        f"{resname}_AC.frcmod", f"{resname}_bcc_gaff2.mol2",
        f"{resname}_GMX.itp", f"{resname}_GMX.gro",
    ):
        (d / ext).write_text("")
    pdb = workdir / f"{resname}.pdb"
    pdb.write_text(
        f"HETATM    1  C   {resname:>3s} A   1       0.000   0.000   0.000\nEND\n"
    )
    return SmallMoleculeResult(
        resname=resname, monomer_pdb=pdb,
        acpype_dir=d,
        frcmod=d / f"{resname}_AC.frcmod",
        mol2=d / f"{resname}_bcc_gaff2.mol2",
        itp=d / f"{resname}_GMX.itp",
        gro=d / f"{resname}_GMX.gro",
        n_atoms=1, net_charge=0,
    )


# ---------------------------------------------------------------------------
# render_run_script
# ---------------------------------------------------------------------------


def test_render_run_script_includes_em_nvt_npt_prod():
    text = render_run_script(gmx_path="gmx")
    assert "md/em.mdp" in text
    assert "md/nvt.mdp" in text
    assert "md/npt.mdp" in text
    assert "md/prod.mdp" in text


def test_render_run_script_uses_custom_gmx():
    text = render_run_script(gmx_path="gmx_mpi")
    assert "GMX=\"${GMX:-gmx_mpi}\"" in text


# ---------------------------------------------------------------------------
# 7-stage orchestrator (full mock)
# ---------------------------------------------------------------------------


def test_builder_build_runs_all_stages(tmp_path):
    """End-to-end with every external tool mocked.

    Verifies that the artifacts dict contains every expected key
    and that file paths exist after build().
    """
    cfg = _make_config(tmp_path / "out")
    b = FormulationBuilder(cfg)

    sm_cprN = _make_sm_result(tmp_path, "CPRN")
    sm_cprC = _make_sm_result(tmp_path, "CPRC")
    sm_tch = _make_sm_result(tmp_path, "TCH")

    fake_top = TopologyArtifacts(
        prmtop=str(tmp_path / "fake.prmtop"),
        inpcrd=str(tmp_path / "fake.inpcrd"),
        top=str(tmp_path / "fake.top"),
        gro=str(tmp_path / "fake.gro"),
        tleap_input=str(tmp_path / "fake.tleap"),
        tleap_log=str(tmp_path / "fake.log"),
    )
    # Plausible .gro: title + n_atoms + box
    Path(fake_top.top).write_text("[ defaults ]\n")
    Path(fake_top.gro).write_text(
        "FAKE\n   3\n"
        "    1ALA      N    1   0.100   0.100   0.100\n"
        "    1ALA     CA    2   0.150   0.100   0.100\n"
        "    2HOH     OW    3   0.500   0.500   0.500\n"
        "  5.000  5.000  5.000\n"
    )

    with patch(
        "abmptools.formulation.builder.build_peptide_from_sequence"
    ) as mock_tleap, patch(
        "abmptools.formulation.builder.parameterize_small_molecule"
    ) as mock_acp, patch(
        "abmptools.formulation.builder.pack_mixed_solution"
    ) as mock_packmol, patch(
        "abmptools.formulation.builder.build_topology"
    ) as mock_top:
        # peptide tleap → write a stub PDB
        def fake_tleap(**kw):
            Path(kw["out_pdb"]).parent.mkdir(parents=True, exist_ok=True)
            Path(kw["out_pdb"]).write_text(
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000\nEND\n"
            )
            return kw["out_pdb"]
        mock_tleap.side_effect = fake_tleap

        # acpype results: 3 (CPR neutral + CPR charged + TCH)
        mock_acp.side_effect = [sm_cprN, sm_cprC, sm_tch]
        mock_packmol.side_effect = lambda **kw: Path(kw["out_pdb"]).write_text("MIX")
        mock_top.return_value = fake_top

        art = b.build()

        # 7-stage order: peptide → small_mol × 3 → packmol → topology
        assert mock_tleap.called
        assert mock_acp.call_count == 3
        assert mock_packmol.called
        assert mock_top.called

        # All artifacts present
        assert Path(art.output_dir).is_dir()
        assert Path(art.top).is_file()
        assert Path(art.gro).is_file()
        assert Path(art.ndx).is_file()
        assert Path(art.run_script).is_file()
        assert Path(art.config_json).is_file()
        assert set(art.mdp_files.keys()) == {"em", "nvt", "npt", "prod"}
        for p in art.mdp_files.values():
            assert Path(p).is_file()

        # n_peptides_total = 2
        assert art.n_peptides_total == 2
        # atoms_estimate from fake .gro = 3
        assert art.atoms_estimate == 3


def test_builder_stage1_uses_existing_pdb_path(tmp_path):
    """If PeptideSpec has pdb_path, normalize the source instead of tleap."""
    src = tmp_path / "src.pdb"
    src.write_text(
        "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n"
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000\n"
        "END\n"
    )
    cfg = _make_config(tmp_path / "out")
    cfg.system.peptides = [PeptideSpec(
        name="insulin", pdb_path=str(src), n_copies=1,
    )]
    b = FormulationBuilder(cfg)
    b._stage0_ff_staging()
    pdbs = b._stage1_peptide_atomistic()
    text = Path(pdbs[0]).read_text()
    assert "CRYST1" not in text
    assert "ALA" in text


def test_builder_run_script_chmod_executable(tmp_path):
    """run.sh has +x permission bits set after stage 7."""
    cfg = _make_config(tmp_path / "out")
    b = FormulationBuilder(cfg)
    b._stage0_ff_staging()
    run_path = b._stage7_run_script()
    mode = run_path.stat().st_mode
    assert mode & 0o111  # any execute bit
