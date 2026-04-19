# -*- coding: utf-8 -*-
"""
Integration smoke test for AmorphousBuilder.build().

Runs the real pipeline on a tiny system (methane x 10, fixed box) to
verify that the produced artifacts are structurally valid. No MD is
executed — only the pre-MD outputs from build() are checked.

Skipped automatically when any of the heavy scientific dependencies
are missing (OpenFF, RDKit, Packmol binary, AmberTools `sqm`).
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pytest

# -- gating ---------------------------------------------------------------

_openff = pytest.importorskip(
    "openff.toolkit", reason="openff-toolkit not installed")
pytest.importorskip(
    "openff.interchange", reason="openff-interchange not installed")
pytest.importorskip("rdkit", reason="rdkit not installed")

_packmol = shutil.which("packmol")
_sqm = shutil.which("sqm")  # AmberTools provides this for AM1-BCC charges

pytestmark = [
    pytest.mark.skipif(_packmol is None,
                       reason="packmol binary not available on PATH"),
    pytest.mark.skipif(_sqm is None,
                       reason="AmberTools `sqm` not available on PATH "
                              "(required for AM1-BCC charges)"),
    pytest.mark.slow,
]


# -- fixtures -------------------------------------------------------------

@pytest.fixture(scope="module")
def build_result(tmp_path_factory):
    """Run a real but tiny build once; every test reads from the result."""
    from abmptools.amorphous import (
        AmorphousBuilder, BuildConfig, ComponentSpec,
    )

    outdir = tmp_path_factory.mktemp("methane_smoke")
    cfg = BuildConfig(
        components=[
            ComponentSpec(smiles="C", name="methane", n_mol=10),
        ],
        box_size_nm=2.0,           # fixed to make the test deterministic-ish
        temperature=300.0,
        T_high=400.0,
        seed=42,
        output_dir=str(outdir),
    )
    builder = AmorphousBuilder(cfg)
    return builder.build()


# -- result-dict shape ----------------------------------------------------

def test_build_returns_all_documented_keys(build_result):
    for key in ("output_dir", "gro", "top", "ndx", "mdp_files",
                "run_script", "wrap_script", "config_json",
                "box_nm", "counts"):
        assert key in build_result, f"missing key: {key!r}"


def test_counts_match_requested(build_result):
    assert build_result["counts"] == [10]


def test_box_size_matches_request(build_result):
    # box was fixed at 2.0 nm on input.
    assert build_result["box_nm"] == pytest.approx(2.0, rel=1e-6)


# -- files on disk --------------------------------------------------------

def test_expected_artifacts_exist(build_result):
    must_exist = [
        build_result["gro"],
        build_result["top"],
        build_result["ndx"],
        build_result["config_json"],
        build_result["run_script"],
        build_result["wrap_script"],
    ] + list(build_result["mdp_files"])
    for p in must_exist:
        assert Path(p).is_file(), f"not found: {p}"


def test_five_mdp_files_with_canonical_names(build_result):
    names = [Path(p).name for p in build_result["mdp_files"]]
    assert names == [
        "01_em.mdp",
        "02_nvt_highT.mdp",
        "03_npt_highT.mdp",
        "04_anneal.mdp",
        "05_npt_final.mdp",
    ]


def test_scripts_are_executable(build_result):
    assert os.access(build_result["run_script"], os.X_OK)
    assert os.access(build_result["wrap_script"], os.X_OK)


# -- file-content sanity checks ------------------------------------------

def test_gro_has_expected_atom_count_and_box(build_result):
    """methane CH4 = 5 atoms * 10 molecules = 50 atoms, box 2.0 nm."""
    lines = Path(build_result["gro"]).read_text().splitlines()
    assert int(lines[1].strip()) == 50
    # last line is the box vectors (nm)
    bx = [float(x) for x in lines[-1].split()[:3]]
    assert all(abs(b - 2.0) < 1e-6 for b in bx)


def test_top_has_required_sections(build_result):
    text = Path(build_result["top"]).read_text()
    for section in ("[ defaults ]", "[ atomtypes ]",
                    "[ moleculetype ]", "[ system ]", "[ molecules ]"):
        assert section in text, f"topology missing section {section}"


def test_ndx_has_system_and_component_groups(build_result):
    text = Path(build_result["ndx"]).read_text()
    assert "[ System ]" in text
    assert "[ methane ]" in text


def test_config_json_roundtrips(build_result):
    data = json.loads(Path(build_result["config_json"]).read_text())
    assert data["components"][0]["smiles"] == "C"
    assert data["components"][0]["n_mol"] == 10
    assert data["box_size_nm"] == pytest.approx(2.0)
    assert data["temperature"] == pytest.approx(300.0)
    assert data["seed"] == 42


def test_run_all_sh_invokes_all_stages(build_result):
    text = Path(build_result["run_script"]).read_text()
    for stage in ("01_em", "02_nvt_highT", "03_npt_highT",
                  "04_anneal", "05_npt_final"):
        assert stage in text


def test_wrap_pbc_sh_uses_gmx_trjconv_with_pbc_mol(build_result):
    text = Path(build_result["wrap_script"]).read_text()
    assert "gmx trjconv" in text
    assert "-pbc mol" in text
    assert "-ur compact" in text
