# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.umbrella."""
from __future__ import annotations

from abmptools.cg.membrane import umbrella
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


# ---------------------------------------------------------------------------
# write_window_mdps
# ---------------------------------------------------------------------------

class TestWriteWindowMdps:
    def test_creates_n_windows_dirs(self, tmp_path):
        result = umbrella.write_window_mdps(
            config=_cfg(), windows_dir=tmp_path,
        )
        assert len(result) == 13  # default UmbrellaCGProtocol
        for i in range(13):
            wd = tmp_path / f"win_{i:03d}"
            assert wd.exists()
            assert (wd / "window.mdp").exists()

    def test_window_mdp_z_value_in_text(self, tmp_path):
        umbrella.write_window_mdps(
            config=_cfg(), windows_dir=tmp_path,
        )
        # Window 0: z = -1.5; Window 12: z = +1.5; Window 6: z = 0.0
        win0 = (tmp_path / "win_000" / "window.mdp").read_text()
        assert "z = -1.500" in win0
        win12 = (tmp_path / "win_012" / "window.mdp").read_text()
        assert "z = +1.500" in win12
        win6 = (tmp_path / "win_006" / "window.mdp").read_text()
        assert "z = +0.000" in win6

    def test_window_mdp_static_pull(self, tmp_path):
        umbrella.write_window_mdps(
            config=_cfg(), windows_dir=tmp_path,
        )
        # All windows have rate=0
        for i in range(13):
            mdp = (tmp_path / f"win_{i:03d}" / "window.mdp").read_text()
            assert "pull-coord1-rate           = 0.000000" in mdp

    def test_window_mdp_init_progression(self, tmp_path):
        umbrella.write_window_mdps(
            config=_cfg(), windows_dir=tmp_path,
        )
        # init goes from -1.5 to +1.5 in 0.25 steps
        for i in range(13):
            z = -1.5 + i * 0.25
            mdp = (tmp_path / f"win_{i:03d}" / "window.mdp").read_text()
            assert f"pull-coord1-init           = {z:.4f}" in mdp

    def test_window_mdp_pbc_atom_propagates(self, tmp_path):
        umbrella.write_window_mdps(
            config=_cfg(), windows_dir=tmp_path,
            pbc_atom_g1=999,
        )
        win0 = (tmp_path / "win_000" / "window.mdp").read_text()
        assert "pull-group1-pbcatom        = 999" in win0


# ---------------------------------------------------------------------------
# write_run_script
# ---------------------------------------------------------------------------

class TestWriteRunScript:
    def _setup(self, tmp_path):
        # Build minimal stub paths inside tmp_path for relative-path math.
        (tmp_path / "topol.top").write_text("; stub\n")
        (tmp_path / "system_ions.gro").write_text("; stub\n")
        (tmp_path / "index.ndx").write_text("[ System ]\n")
        equil = tmp_path / "mdp"
        equil.mkdir()
        em = equil / "em.mdp"
        em.write_text("integrator = steep\n")
        nvt = equil / "nvt.mdp"
        nvt.write_text("integrator = md\n")
        npt = equil / "npt.mdp"
        npt.write_text("integrator = md\n")
        pull_d = tmp_path / "pull"
        pull_d.mkdir()
        pull_mdp = pull_d / "pull.mdp"
        pull_mdp.write_text("integrator = md\n")
        win_dir = tmp_path / "windows"
        win_dir.mkdir()
        win_mdps = {}
        for i in range(13):
            wd = win_dir / f"win_{i:03d}"
            wd.mkdir()
            mdp = wd / "window.mdp"
            mdp.write_text("integrator = md\n")
            win_mdps[i] = mdp
        return {
            "top": tmp_path / "topol.top",
            "gro": tmp_path / "system_ions.gro",
            "ndx": tmp_path / "index.ndx",
            "equil_mdps": {"em": em, "nvt": nvt, "npt": npt},
            "pull_mdp": pull_mdp,
            "window_mdps": win_mdps,
        }

    def test_writes_executable_script(self, tmp_path):
        s = self._setup(tmp_path)
        script = umbrella.write_run_script(
            config=_cfg(), output_dir=tmp_path, **s,
        )
        assert script.exists()
        assert script.name == "run.sh"
        # Check x permission
        assert script.stat().st_mode & 0o100

    def test_run_script_includes_all_stages(self, tmp_path):
        s = self._setup(tmp_path)
        script = umbrella.write_run_script(
            config=_cfg(), output_dir=tmp_path, **s,
        )
        text = script.read_text()
        assert "Stage 1: energy minimisation" in text
        assert "Stage 2: NVT" in text
        assert "Stage 3: NPT" in text
        assert "Stage 4: pulling" in text
        assert "Stage 5: extract per-window starting frames" in text
        assert "Stage 6: per-window MD" in text
        assert "Stage 7: PMF analysis" in text

    def test_run_script_uses_cg_membrane_python_invocations(self, tmp_path):
        s = self._setup(tmp_path)
        script = umbrella.write_run_script(
            config=_cfg(), output_dir=tmp_path, **s,
        )
        text = script.read_text()
        assert "python -m abmptools.cg.membrane make-windows" in text
        assert "python -m abmptools.cg.membrane wham" in text

    def test_run_script_window_loop_count(self, tmp_path):
        s = self._setup(tmp_path)
        script = umbrella.write_run_script(
            config=_cfg(), output_dir=tmp_path, **s,
        )
        text = script.read_text()
        # 13 windows -> seq 0 12
        assert "seq -f '%03g' 0 12" in text

    def test_run_script_uses_configured_gmx_path(self, tmp_path):
        s = self._setup(tmp_path)
        cfg = _cfg()
        cfg.gmx_path = "/opt/bin/gmx-2024"
        script = umbrella.write_run_script(
            config=cfg, output_dir=tmp_path, **s,
        )
        text = script.read_text()
        assert "/opt/bin/gmx-2024" in text
