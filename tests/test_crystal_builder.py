# -*- coding: utf-8 -*-
"""End-to-end tests for ``abmptools.crystal.builder.CrystalOrchestrator``.

Stages 1+2+3 (expand_cif → generate_fmo → write_jobs) on csp7 R00001
must reproduce the Phase B fixture byte-for-byte. ``run_abinit`` and
``postprocess`` raise ``NotImplementedError`` at this point — they
land in Phase D and Phase C-6 respectively.
"""
from __future__ import annotations

import os

import pytest

from .test_regression import _compare_files, _compare_output_dir  # noqa: F401

TESTS_DIR = os.path.dirname(__file__)
REF_MAIN = os.path.join(TESTS_DIR, "regression", "reference", "main")
CRYSTAL_FIXTURE = os.path.join(REF_MAIN, "crystal_csp7")
SAMPLE_CSP7_SMOKE = os.path.join(
    TESTS_DIR, os.pardir, "sample", "crystal", "csp7_smoke",
)
UNK_AJF_TEMPLATE = os.path.abspath(
    os.path.join(SAMPLE_CSP7_SMOKE, "UNK.ajf")
)

_HAS_FIXTURE = os.path.isdir(os.path.join(CRYSTAL_FIXTURE, "R00001"))
_HAS_TEMPLATE = os.path.isfile(UNK_AJF_TEMPLATE)


# ---------------------------------------------------------------------------
# Builder unit pieces (no external pipeline run)
# ---------------------------------------------------------------------------

def test_emit_input_param_keys_match_setfmo():
    """The synthesised input_param must carry the keys
    :func:`abmptools.setfmo.setrfmoparam` reads."""
    from abmptools.crystal.builder import _emit_input_param
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FMOMethod,
        FragmentTemplate,
    )

    cfg = CrystalBuildConfig(
        inputs=[CIFInputSpec(cif="x.cif")],
        fragment=FragmentTemplate(
            cutmode="around", solutes=[0], criteria=6.0,
            molname=["UNK"], pieda=True, cmm=True,
        ),
        fmo=FMOMethod(method="MP2", basis_set="6-31Gdag", memory=6000),
    )
    text = _emit_input_param(cfg)
    # Required keys verified against setfmo.setrfmoparam.
    for key in (
        "'cutmode'", "'solutes'", "'criteria'", "'tgtpos'", "'molname'",
        "'getmode'", "'piedaflag'", "'cmmflag'", "'solv_flag'",
        "'ajf_method'", "'ajf_basis_set'", "'cpfflag'", "'abinit_ver'",
        "'memory'", "'npro'",
    ):
        assert key in text, f"missing key in input_param: {key}"
    # Sanity: the literal must parse as Python.
    namespace: dict = {}
    exec(text, namespace)
    assert namespace["param"]["cutmode"] == "around"
    assert namespace["param"]["piedaflag"] is True


def test_emit_segment_data_basic():
    from abmptools.crystal.builder import _emit_segment_data
    from abmptools.crystal.models import CIFInputSpec

    spec = CIFInputSpec(cif="x.cif", atoms_in_mol=[32])
    text = _emit_segment_data(spec, molname="UNK")
    namespace: dict = {}
    exec(text, namespace)
    seg = namespace["seg_data"][0]
    assert seg["name"] == "UNK"
    assert seg["atom"] == [32]
    assert seg["seg_info"] == [list(range(1, 33))]


# ---------------------------------------------------------------------------
# Full pipeline (slow): R00001 only — validates byte-equivalence with
# the Phase B fixture for stages 1+2.
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.skipif(
    not (_HAS_FIXTURE and _HAS_TEMPLATE),
    reason="crystal_csp7 fixture or csp7_smoke UNK.ajf missing",
)
def test_orchestrator_csp7_r00001_byte_equivalent(tmp_path):
    """``CrystalOrchestrator.run`` must reproduce the Phase B
    R00001 fixture exactly through stages 1+2; stage 3 must produce
    one jobscript per input AJF."""
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FMOMethod,
        FragmentTemplate,
        HPCJobSpec,
    )

    fixture_dir = os.path.join(CRYSTAL_FIXTURE, "R00001")
    cif_src = os.path.join(fixture_dir, "XXXI-MMFF-R00001.cif")
    cif_dst = tmp_path / "XXXI-MMFF-R00001.cif"
    import shutil
    shutil.copy(cif_src, cif_dst)

    cfg = CrystalBuildConfig(
        project_name="csp7_r00001_smoke",
        output_dir=str(tmp_path / "out"),
        inputs=[
            CIFInputSpec(cif=str(cif_dst), layer=5, atoms_in_mol=[32]),
        ],
        fragment=FragmentTemplate(
            cutmode="around", solutes=[0], criteria=6.0,
            molname=["UNK"], pieda=True, cmm=True, getmode="rfile",
            template_ajf=UNK_AJF_TEMPLATE,
        ),
        fmo=FMOMethod(
            method="MP2", basis_set="6-31Gdag", memory=6000,
            cpfflag=True, abinit_ver="rev23", npro=1, is_xyz=True,
        ),
        hpc=HPCJobSpec(
            scheduler="PJM", queue="small", group="hp190133",
            nodes=12, proc_per_node=2, omp_threads=24,
            elapse="24:00:00",
            abinit_dir="/data/hp190133/programs/ABINIT-MP",
        ),
    )

    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))
    summary = orch.run()
    assert summary["n_inputs"] == 1
    assert summary["n_ajf_total"] == 1
    assert summary["n_jobscripts"] == 1

    # Stages 1+2: byte-equivalence with the Phase B fixture.
    out_for_abmp = (
        tmp_path / "out" / "XXXI-MMFF-R00001"
        / "cifout" / "layer5" / "pdb" / "for_abmp"
    )
    ref_for_abmp = os.path.join(fixture_dir, "for_abmp")
    _compare_output_dir(str(out_for_abmp), ref_for_abmp)

    # Stage 3: one PJM jobscript was rendered next to the AJF.
    scripts = list(out_for_abmp.glob("*.sh"))
    assert any("12n-2p-24t" in s.name for s in scripts), (
        f"expected a 12n-2p-24t PJM jobscript among {[s.name for s in scripts]}"
    )
    runbatch = out_for_abmp / "runbatch.sh"
    assert runbatch.is_file(), "runbatch.sh not emitted"


@pytest.mark.skipif(
    not _HAS_TEMPLATE, reason="csp7_smoke UNK.ajf missing",
)
def test_orchestrator_postprocess_stub_raises(tmp_path):
    """``postprocess`` is still a stub at the orchestrator level; the
    standalone Stage-5 driver in :mod:`abmptools.crystal.postproc` is
    the supported entry point until full integration in Phase D-3."""
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FragmentTemplate,
        PostProcessSpec,
    )

    cfg = CrystalBuildConfig(
        project_name="stub",
        output_dir=str(tmp_path / "out"),
        inputs=[CIFInputSpec(cif="never_called.cif")],
        fragment=FragmentTemplate(template_ajf=UNK_AJF_TEMPLATE),
        postproc=PostProcessSpec(enable=True),
    )
    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))
    with pytest.raises(NotImplementedError, match="Phase C-6"):
        orch.postprocess()


# ---------------------------------------------------------------------------
# run_abinit (Phase D-1) — binary resolution + missing-binary errors
# ---------------------------------------------------------------------------

def test_resolve_abinit_binary_missing(tmp_path, monkeypatch):
    """Unknown binary + empty abinit_dir raises FileNotFoundError."""
    import abmptools.crystal.builder as builder_mod
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FragmentTemplate,
        HPCJobSpec,
    )

    cfg = CrystalBuildConfig(
        project_name="missing_abinit",
        output_dir=str(tmp_path / "out"),
        inputs=[CIFInputSpec(cif="never_called.cif")],
        fragment=FragmentTemplate(template_ajf=UNK_AJF_TEMPLATE),
        hpc=HPCJobSpec(
            scheduler="local",
            abinit_dir="",
            binary_name="abinitmp_does_not_exist_xyz",
        ),
    )
    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))
    monkeypatch.setattr(builder_mod.shutil, "which", lambda name: None)
    with pytest.raises(FileNotFoundError, match="abinitmp binary not found"):
        orch._resolve_abinit_binary()


def test_resolve_abinit_binary_via_path(tmp_path, monkeypatch):
    """``shutil.which(binary_name)`` is the second resolution step."""
    import abmptools.crystal.builder as builder_mod
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FragmentTemplate,
        HPCJobSpec,
    )

    fake = tmp_path / "fake_abinitmp"
    fake.write_text("#!/bin/bash\necho fake\n")
    fake.chmod(0o755)

    cfg = CrystalBuildConfig(
        project_name="path_abinit",
        output_dir=str(tmp_path / "out"),
        inputs=[CIFInputSpec(cif="never_called.cif")],
        fragment=FragmentTemplate(template_ajf=UNK_AJF_TEMPLATE),
        hpc=HPCJobSpec(
            scheduler="local",
            abinit_dir="",
            binary_name="fake_abinitmp",
        ),
    )
    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))
    monkeypatch.setattr(
        builder_mod.shutil, "which",
        lambda name: str(fake) if name == "fake_abinitmp" else None,
    )
    resolved = orch._resolve_abinit_binary()
    assert resolved == fake


def test_resolve_abinit_binary_via_abinit_dir(tmp_path, monkeypatch):
    """``abinit_dir/binary_name`` takes precedence over PATH lookup."""
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FragmentTemplate,
        HPCJobSpec,
    )

    bin_dir = tmp_path / "abinitmp_bins"
    bin_dir.mkdir()
    bin_path = bin_dir / "abinitmp_smp"
    bin_path.write_text("#!/bin/bash\necho fake\n")
    bin_path.chmod(0o755)

    cfg = CrystalBuildConfig(
        project_name="dir_abinit",
        output_dir=str(tmp_path / "out"),
        inputs=[CIFInputSpec(cif="never_called.cif")],
        fragment=FragmentTemplate(template_ajf=UNK_AJF_TEMPLATE),
        hpc=HPCJobSpec(
            scheduler="local",
            abinit_dir=str(bin_dir),
            binary_name="abinitmp_smp",
        ),
    )
    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))
    resolved = orch._resolve_abinit_binary()
    assert resolved == bin_path


def test_run_abinit_with_fake_binary(tmp_path, monkeypatch):
    """``run_abinit`` should invoke the binary, capture log/err, and
    return the log paths. Uses a stub binary that echoes its stdin to
    stdout so we don't hit a real ABINIT-MP run."""
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FragmentTemplate,
        HPCJobSpec,
    )

    # Fake binary: cat-like, exit 0.
    fake = tmp_path / "fake_abinitmp"
    fake.write_text("#!/bin/bash\ncat\n")
    fake.chmod(0o755)

    # Skeleton run that bypasses generate_fmo: hand-construct fmo_results.
    cfg = CrystalBuildConfig(
        project_name="fake_run",
        output_dir=str(tmp_path / "out"),
        inputs=[CIFInputSpec(cif="dummy.cif", name="dummy")],
        fragment=FragmentTemplate(template_ajf=UNK_AJF_TEMPLATE),
        hpc=HPCJobSpec(
            scheduler="local",
            abinit_dir=str(tmp_path),
            binary_name="fake_abinitmp",
        ),
    )
    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))

    # Inject a fake fmo_results pointing at a real AJF on disk.
    for_abmp_dir = tmp_path / "for_abmp"
    for_abmp_dir.mkdir()
    fake_ajf = for_abmp_dir / "stub.ajf"
    fake_ajf.write_text("STUB AJF CONTENT\n")
    orch.fmo_results = {
        "dummy": {
            "ajf_files": [fake_ajf],
            "pdb_files": [],
            "for_abmp_dir": for_abmp_dir,
        },
    }

    logs = orch.run_abinit(timeout_s=30)
    assert len(logs) == 1
    assert logs[0].is_file()
    # The cat-like fake binary copies STDIN to STDOUT, so the log is
    # the AJF content itself.
    assert "STUB AJF CONTENT" in logs[0].read_text()


def test_run_abinit_propagates_failure(tmp_path):
    """Non-zero exit raises by default (fail_fast=True)."""
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FragmentTemplate,
        HPCJobSpec,
    )

    fake = tmp_path / "boom"
    fake.write_text("#!/bin/bash\nexit 7\n")
    fake.chmod(0o755)

    cfg = CrystalBuildConfig(
        project_name="boom_run",
        output_dir=str(tmp_path / "out"),
        inputs=[CIFInputSpec(cif="dummy.cif", name="d")],
        fragment=FragmentTemplate(template_ajf=UNK_AJF_TEMPLATE),
        hpc=HPCJobSpec(
            scheduler="local",
            abinit_dir=str(tmp_path),
            binary_name="boom",
        ),
        fail_fast=True,
    )
    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))

    for_abmp_dir = tmp_path / "for_abmp_boom"
    for_abmp_dir.mkdir()
    ajf = for_abmp_dir / "x.ajf"
    ajf.write_text("body\n")
    orch.fmo_results = {
        "d": {
            "ajf_files": [ajf],
            "pdb_files": [],
            "for_abmp_dir": for_abmp_dir,
        },
    }

    with pytest.raises(RuntimeError, match="exited non-zero"):
        orch.run_abinit(timeout_s=10)


@pytest.mark.skipif(
    not _HAS_TEMPLATE, reason="csp7_smoke UNK.ajf missing",
)
def test_orchestrator_template_ajf_required(tmp_path):
    """``FragmentTemplate.template_ajf=None`` should raise a clear error."""
    from abmptools.crystal.builder import CrystalOrchestrator
    from abmptools.crystal.models import (
        CIFInputSpec,
        CrystalBuildConfig,
        FragmentTemplate,
    )

    cif_src = os.path.join(CRYSTAL_FIXTURE, "R00001", "XXXI-MMFF-R00001.cif")
    if not os.path.isfile(cif_src):
        pytest.skip("R00001 fixture missing")
    cif_dst = tmp_path / "XXXI-MMFF-R00001.cif"
    import shutil
    shutil.copy(cif_src, cif_dst)

    cfg = CrystalBuildConfig(
        project_name="missing_template",
        output_dir=str(tmp_path / "out"),
        inputs=[CIFInputSpec(cif=str(cif_dst), atoms_in_mol=[32])],
        fragment=FragmentTemplate(template_ajf=None),
    )
    orch = CrystalOrchestrator(cfg, config_dir=str(tmp_path))
    with pytest.raises(ValueError, match="template_ajf is required"):
        orch.run()
