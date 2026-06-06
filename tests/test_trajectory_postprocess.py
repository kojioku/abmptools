"""Unit tests for :mod:`abmptools.trajectory.postprocess`.

gmx subprocess を monkeypatch して、 引数組立と output path 命名を検証。
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import List, Tuple

import pytest

from abmptools.trajectory import (
    GmxError,
    gmx_energy,
    nojump,
    run_trjconv,
    thin,
    thin_and_nojump,
    wrap_pbc,
)
from abmptools.trajectory import postprocess as pp


@pytest.fixture()
def gmx_recorder(monkeypatch):
    """``subprocess.run`` を mock し、 (cmd, stdin) を記録する。"""
    records: List[Tuple[List[str], bytes]] = []

    class _Result:
        returncode = 0
        stdout = b""
        stderr = b""

    def _fake_run(cmd, input=None, capture_output=False, check=False, **kw):
        records.append((list(cmd), input or b""))
        # output file を空で作って下流の Path.resolve() を成功させる
        try:
            out_idx = cmd.index("-o")
            Path(cmd[out_idx + 1]).touch()
        except (ValueError, IndexError):
            pass
        return _Result()

    monkeypatch.setattr(subprocess, "run", _fake_run)
    # _resolve_gmx を素通しに
    monkeypatch.setattr(pp, "_resolve_gmx", lambda gmx: gmx)
    return records


@pytest.fixture()
def fake_inputs(tmp_path):
    traj = tmp_path / "prod.xtc"
    tpr = tmp_path / "prod.tpr"
    traj.write_bytes(b"fake xtc")
    tpr.write_bytes(b"fake tpr")
    return traj, tpr


def test_thin_and_nojump_default_output_name(gmx_recorder, fake_inputs):
    traj, tpr = fake_inputs
    out = thin_and_nojump(trajectory=traj, tpr=tpr, skip=10)
    assert out.name == "prod_nojump_skip10.xtc"

    cmd, stdin = gmx_recorder[0]
    assert "-pbc" in cmd and "nojump" in cmd
    assert "-skip" in cmd and "10" in cmd
    assert stdin == b"System\n"


def test_thin_and_nojump_custom_skip(gmx_recorder, fake_inputs):
    traj, tpr = fake_inputs
    out = thin_and_nojump(trajectory=traj, tpr=tpr, skip=5)
    assert out.name == "prod_nojump_skip5.xtc"
    cmd, _ = gmx_recorder[0]
    assert "5" in cmd


def test_nojump_no_skip_flag(gmx_recorder, fake_inputs):
    traj, tpr = fake_inputs
    out = nojump(trajectory=traj, tpr=tpr)
    assert out.name == "prod_nojump.xtc"
    cmd, _ = gmx_recorder[0]
    assert "-skip" not in cmd
    assert "nojump" in cmd


def test_thin_no_pbc_flag(gmx_recorder, fake_inputs):
    traj, tpr = fake_inputs
    out = thin(trajectory=traj, tpr=tpr, skip=20)
    assert out.name == "prod_skip20.xtc"
    cmd, _ = gmx_recorder[0]
    assert "-pbc" not in cmd
    assert "-skip" in cmd and "20" in cmd


def test_wrap_pbc_default(gmx_recorder, fake_inputs):
    traj, tpr = fake_inputs
    out = wrap_pbc(trajectory=traj, tpr=tpr)
    assert out.name == "prod_pbc.xtc"
    cmd, stdin = gmx_recorder[0]
    assert "-pbc" in cmd and "mol" in cmd
    assert "-ur" in cmd and "compact" in cmd
    assert "-center" not in cmd
    assert stdin == b"System\n"


def test_wrap_pbc_with_center(gmx_recorder, fake_inputs):
    traj, tpr = fake_inputs
    wrap_pbc(trajectory=traj, tpr=tpr, center="Peptide", group="System")
    cmd, stdin = gmx_recorder[0]
    assert "-center" in cmd
    # center group が 1 行目、 output group が 2 行目
    assert stdin == b"Peptide\nSystem\n"


def test_ndx_flag_added(gmx_recorder, fake_inputs, tmp_path):
    traj, tpr = fake_inputs
    ndx = tmp_path / "system.ndx"
    ndx.write_text("[ System ]\n1\n")
    thin_and_nojump(trajectory=traj, tpr=tpr, ndx=ndx)
    cmd, _ = gmx_recorder[0]
    assert "-n" in cmd
    n_idx = cmd.index("-n")
    assert cmd[n_idx + 1] == str(ndx)


def test_custom_output_path(gmx_recorder, fake_inputs, tmp_path):
    traj, tpr = fake_inputs
    custom = tmp_path / "sub" / "myout.xtc"
    out = thin_and_nojump(trajectory=traj, tpr=tpr, output=custom)
    assert out == custom.resolve()
    # parent dir 自動作成
    assert custom.parent.is_dir()


def test_missing_trajectory_raises(tmp_path):
    tpr = tmp_path / "prod.tpr"
    tpr.write_bytes(b"x")
    with pytest.raises(FileNotFoundError, match="trajectory not found"):
        nojump(trajectory=tmp_path / "missing.xtc", tpr=tpr)


def test_missing_reference_raises(tmp_path):
    traj = tmp_path / "prod.xtc"
    traj.write_bytes(b"x")
    with pytest.raises(FileNotFoundError, match="reference structure not found"):
        nojump(trajectory=traj, tpr=tmp_path / "missing.tpr")


def test_gmx_error_attaches_stderr(monkeypatch, fake_inputs):
    traj, tpr = fake_inputs

    class _Result:
        returncode = 99
        stdout = b"out blah"
        stderr = b"err blah"

    monkeypatch.setattr(subprocess, "run",
                        lambda *a, **kw: _Result())
    monkeypatch.setattr(pp, "_resolve_gmx", lambda gmx: gmx)

    with pytest.raises(GmxError) as exc:
        nojump(trajectory=traj, tpr=tpr)
    assert exc.value.returncode == 99
    assert "err blah" in str(exc.value)
    assert "out blah" in str(exc.value)


def test_resolve_gmx_missing(monkeypatch, fake_inputs):
    """gmx が PATH に無いと FileNotFoundError."""
    traj, tpr = fake_inputs
    monkeypatch.setattr("shutil.which", lambda _: None)
    with pytest.raises(FileNotFoundError, match="not found in PATH"):
        nojump(trajectory=traj, tpr=tpr, gmx="nonexistent-gmx-binary")


def test_run_trjconv_sequence_group(gmx_recorder, fake_inputs):
    """group が sequence なら複数行で stdin に渡る."""
    traj, tpr = fake_inputs
    out = tmp_out = fake_inputs[0].parent / "out.xtc"
    run_trjconv(
        trajectory=traj, reference=tpr, output=out,
        group=["GroupA", "GroupB"],
        extra_args=("-pbc", "mol"),
    )
    cmd, stdin = gmx_recorder[0]
    assert stdin == b"GroupA\nGroupB\n"


def test_gmx_energy_default_terms(gmx_recorder, tmp_path):
    edr = tmp_path / "prod.edr"
    edr.write_bytes(b"fake edr")
    out = tmp_path / "prod_energy.xvg"
    gmx_energy(edr=edr, output=out)
    cmd, stdin = gmx_recorder[0]
    assert "energy" in cmd
    assert "-f" in cmd and str(edr) in cmd
    # default terms = 1..50 → stdin に 50 行
    lines = stdin.decode().splitlines()
    assert lines[0] == "1"
    assert lines[-1] == "50"
    assert len(lines) == 50


def test_gmx_energy_custom_terms(gmx_recorder, tmp_path):
    edr = tmp_path / "prod.edr"
    edr.write_bytes(b"x")
    out = tmp_path / "out.xvg"
    gmx_energy(edr=edr, output=out, terms=[1, 5, 10, "Potential"])
    _, stdin = gmx_recorder[0]
    assert stdin == b"1\n5\n10\nPotential\n"


def test_gmx_energy_missing_edr_raises(tmp_path):
    out = tmp_path / "out.xvg"
    with pytest.raises(FileNotFoundError, match="edr not found"):
        gmx_energy(edr=tmp_path / "missing.edr", output=out)


def test_cli_energy_entry(monkeypatch, tmp_path, capsys):
    from abmptools.trajectory.cli import main

    edr = tmp_path / "prod.edr"
    edr.write_bytes(b"x")
    out = tmp_path / "prod_energy.xvg"

    def _fake_run(cmd, input=None, capture_output=False, check=False, **kw):
        out_idx = cmd.index("-o")
        Path(cmd[out_idx + 1]).touch()

        class R:
            returncode = 0
            stdout = b""
            stderr = b""

        return R()

    monkeypatch.setattr(subprocess, "run", _fake_run)
    monkeypatch.setattr(pp, "_resolve_gmx", lambda gmx: gmx)

    rc = main([
        "energy",
        "--edr", str(edr),
        "--out", str(out),
        "--terms-max", "30",
    ])
    assert rc == 0
    out_msg = capsys.readouterr().out
    assert "saved:" in out_msg


def test_cli_thin_nojump_entry(monkeypatch, fake_inputs, capsys):
    """``python -m abmptools.trajectory thin_nojump ...`` の動作確認."""
    from abmptools.trajectory.cli import main

    traj, tpr = fake_inputs

    captured: list = []

    def _fake_run(cmd, input=None, capture_output=False, check=False, **kw):
        captured.append(list(cmd))
        out_idx = cmd.index("-o")
        Path(cmd[out_idx + 1]).touch()

        class R:
            returncode = 0
            stdout = b""
            stderr = b""

        return R()

    monkeypatch.setattr(subprocess, "run", _fake_run)
    monkeypatch.setattr(pp, "_resolve_gmx", lambda gmx: gmx)

    rc = main([
        "thin_nojump",
        "--traj", str(traj),
        "--tpr", str(tpr),
        "--skip", "5",
    ])
    assert rc == 0
    out_msg = capsys.readouterr().out
    assert "saved:" in out_msg
    assert "_nojump_skip5.xtc" in out_msg
