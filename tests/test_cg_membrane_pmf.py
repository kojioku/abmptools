# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.pmf.run_wham (mocked gmx wham)."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from abmptools.cg.membrane import pmf
from abmptools.cg.membrane.models import (
    LipidMix,
    MembraneCGBuildConfig,
    PeptideMembraneSpec,
)


def _cfg():
    return MembraneCGBuildConfig(
        lipids=[LipidMix()],
        peptide=PeptideMembraneSpec(),
        martini_itp_dir="./ff",
    )


def _make_dummy_windows(tmp_path, n=3):
    """Stub windows/win_NNN/{window.tpr,pullx.xvg} for run_wham."""
    wd = tmp_path / "windows"
    wd.mkdir()
    for i in range(n):
        d = wd / f"win_{i:03d}"
        d.mkdir()
        (d / "window.tpr").write_text("tpr stub\n")
        (d / "pullx.xvg").write_text("0 -1.0\n")
    return wd


def test_run_wham_invokes_aa_implementation(tmp_path, monkeypatch):
    """The CG run_wham should delegate to abmptools.membrane.pmf.run_wham."""
    wd = _make_dummy_windows(tmp_path, n=3)
    ad = tmp_path / "analysis"

    def fake_aa_run_wham(*, windows_dir, analysis_dir, config, **kw):
        # Make sure we got the CG config (duck-typed)
        assert hasattr(config, "gmx_path")
        assert hasattr(config, "equilibration")
        assert config.equilibration.temperature_K == 310.0
        Path(analysis_dir).mkdir(parents=True, exist_ok=True)
        return {"pmf": Path(analysis_dir) / "pmf.xvg",
                "histo": Path(analysis_dir) / "histo.xvg",
                "tpr_list": [],
                "pullx_list": []}

    monkeypatch.setattr(pmf, "_aa_run_wham", fake_aa_run_wham)
    result = pmf.run_wham(
        windows_dir=wd, analysis_dir=ad, config=_cfg(),
    )
    assert "pmf" in result
    assert "histo" in result


def test_run_wham_passes_bootstrap_n(tmp_path, monkeypatch):
    wd = _make_dummy_windows(tmp_path, n=3)
    captured = {}

    def fake_aa_run_wham(*, windows_dir, analysis_dir, config, bootstrap_n, temperature_K):
        captured["bootstrap_n"] = bootstrap_n
        captured["temperature_K"] = temperature_K
        Path(analysis_dir).mkdir(parents=True, exist_ok=True)
        return {"pmf": "x", "histo": "y", "tpr_list": [], "pullx_list": []}

    monkeypatch.setattr(pmf, "_aa_run_wham", fake_aa_run_wham)
    pmf.run_wham(
        windows_dir=wd,
        analysis_dir=tmp_path / "ana",
        config=_cfg(),
        bootstrap_n=50,
        temperature_K=298.15,
    )
    assert captured["bootstrap_n"] == 50
    assert captured["temperature_K"] == 298.15


def test_run_wham_default_temperature_passes_none(tmp_path, monkeypatch):
    wd = _make_dummy_windows(tmp_path, n=3)
    captured = {}

    def fake_aa_run_wham(*, windows_dir, analysis_dir, config, bootstrap_n, temperature_K):
        captured["temperature_K"] = temperature_K
        Path(analysis_dir).mkdir(parents=True, exist_ok=True)
        return {"pmf": "x", "histo": "y", "tpr_list": [], "pullx_list": []}

    monkeypatch.setattr(pmf, "_aa_run_wham", fake_aa_run_wham)
    pmf.run_wham(
        windows_dir=wd,
        analysis_dir=tmp_path / "ana",
        config=_cfg(),
    )
    # When None, AA's run_wham fills it from config.equilibration.temperature_K.
    assert captured["temperature_K"] is None


# ---------------------------------------------------------------------------
# Real AA call with subprocess.run mocked (end-to-end through the wrapper)
# ---------------------------------------------------------------------------

def _fake_gmx_wham_run(captured_argv):
    """Build a fake subprocess.run that materializes gmx wham outputs."""
    def fake_run(cmd, **kw):
        captured_argv.append(list(cmd))
        # AA's run_wham checks that ``-o pmf.xvg`` and ``-hist histo.xvg``
        # files exist after the call. Materialize them.
        for i, tok in enumerate(cmd):
            if tok in ("-o", "-hist") and i + 1 < len(cmd):
                Path(cmd[i + 1]).parent.mkdir(parents=True, exist_ok=True)
                Path(cmd[i + 1]).write_text("# stub xvg\n0.0  0.0\n")
        return MagicMock(returncode=0, stdout="", stderr="")
    return fake_run


def test_full_chain_with_mocked_subprocess(tmp_path, monkeypatch):
    """Verify the full call chain -> AA -> subprocess.run with -temp set."""
    wd = _make_dummy_windows(tmp_path, n=2)
    captured_argv = []

    monkeypatch.setattr(
        "abmptools.membrane.pmf.subprocess.run",
        _fake_gmx_wham_run(captured_argv),
    )
    pmf.run_wham(
        windows_dir=wd,
        analysis_dir=tmp_path / "ana",
        config=_cfg(),
    )
    assert captured_argv, "gmx wham was never called"
    argv = captured_argv[0]
    assert "wham" in argv
    assert "-temp" in argv
    # Default 310.0 K from EquilibrationCGProtocol
    temp_idx = argv.index("-temp")
    assert float(argv[temp_idx + 1]) == 310.0


def test_full_chain_uses_custom_gmx_path(tmp_path, monkeypatch):
    wd = _make_dummy_windows(tmp_path, n=2)
    cfg = _cfg()
    cfg.gmx_path = "/opt/bin/my_gmx"
    captured_argv = []

    monkeypatch.setattr(
        "abmptools.membrane.pmf.subprocess.run",
        _fake_gmx_wham_run(captured_argv),
    )
    pmf.run_wham(
        windows_dir=wd,
        analysis_dir=tmp_path / "ana",
        config=cfg,
    )
    argv = captured_argv[0]
    assert argv[0] == "/opt/bin/my_gmx"
