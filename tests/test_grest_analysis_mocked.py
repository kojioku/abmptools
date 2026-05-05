# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.grest.analysis."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

# Ensure non-GUI backend is used inside CI.
import matplotlib
matplotlib.use("Agg")

from abmptools.genesis.grest.analysis import (
    KB,
    _genesis_to_cpptraj_mask,
    compute_distance_pmf,
    parse_acceptance_log,
    parse_rem_file,
    parse_rem_files,
    plot_acceptance_ratio,
    plot_replica_transition,
    run_remd_convert,
)
from abmptools.genesis.grest.models import (
    GrestBuildConfig,
    ReplicaTemperatureSpec,
    RESTSelectionSpec,
)


# ---------------------------------------------------------------------------
# .rem file parsing
# ---------------------------------------------------------------------------

class TestParseRemFile:
    def test_basic(self, tmp_path):
        f = tmp_path / "rep1.rem"
        f.write_text(
            "# step replica parameter\n"
            "0 1 1\n"
            "500 1 1\n"
            "1000 1 2\n"
            "1500 1 1\n"
        )
        arr = parse_rem_file(f)
        assert arr.shape == (4, 3)
        assert arr[0].tolist() == [0, 1, 1]
        assert arr[2].tolist() == [1000, 1, 2]

    def test_skips_comments_and_blanks(self, tmp_path):
        f = tmp_path / "x.rem"
        f.write_text(
            "# header\n\n# more\n0 1 1\n   \n500 1 2\n"
        )
        arr = parse_rem_file(f)
        assert arr.shape == (2, 3)

    def test_empty_returns_empty(self, tmp_path):
        f = tmp_path / "empty.rem"
        f.write_text("# only comments\n\n")
        arr = parse_rem_file(f)
        assert arr.shape == (0, 3)


class TestParseRemFiles:
    def test_per_replica_files(self, tmp_path):
        # 4 replicas, 3 steps each. Replica id is constant per file.
        for rep in range(1, 5):
            (tmp_path / f"rep{rep}.rem").write_text(
                f"0 {rep} {rep}\n500 {rep} {rep}\n1000 {rep} {(rep % 4) + 1}\n"
            )
        files = sorted(tmp_path.glob("rep*.rem"))
        d = parse_rem_files(files)
        assert sorted(d.keys()) == [1, 2, 3, 4]
        for rep, arr in d.items():
            assert arr.shape == (3, 2)
            assert arr[0, 0] == 0


# ---------------------------------------------------------------------------
# replica transition plot
# ---------------------------------------------------------------------------

class TestPlotReplicaTransition:
    def test_writes_png_and_csv(self, tmp_path):
        for rep in range(1, 5):
            (tmp_path / f"rep{rep}.rem").write_text(
                f"0 {rep} {rep}\n500 {rep} {rep}\n1000 {rep} {(rep % 4) + 1}\n"
            )
        rem_files = sorted(tmp_path.glob("rep*.rem"))
        out_png = tmp_path / "transition.png"
        out_csv = tmp_path / "transition.csv"
        series = plot_replica_transition(rem_files, out_png, out_csv)
        assert out_png.exists()
        assert out_csv.exists()
        # Each replica had 3 data points.
        for rep, arr in series.items():
            assert arr.shape == (3, 2)
        # CSV header sane.
        first_line = out_csv.read_text().splitlines()[0]
        assert first_line == "step,replica,parameter"

    def test_no_data_raises(self, tmp_path):
        (tmp_path / "empty.rem").write_text("# nothing\n")
        with pytest.raises(ValueError, match="No usable data"):
            plot_replica_transition([tmp_path / "empty.rem"], tmp_path / "x.png")


# ---------------------------------------------------------------------------
# Acceptance log parsing + plot
# ---------------------------------------------------------------------------

POC_LOG_EXCERPT = """\
[STEP5] Perform Replica-Exchange MD Simulation

REMD> Step:         80   Dimension:    1   ExchangePattern:    2
  Replica      ExchangeTrial             AcceptanceRatio      Before       After
        1          1 >     0   N           0 /         0     300.000     300.000
        2          2 >     3   A           1 /         1     318.110     337.110
        3          3 >     2   A           1 /         1     337.110     318.110
        4          4 >     0   N           0 /         0     357.100     357.100

  Parameter    :    300.000   337.110   318.110   357.100
  RepIDtoParmID:          1         3         2         4
  ParmIDtoRepID:          1         3         2         4

REMD> Step:        160   Dimension:    1   ExchangePattern:    1
        1          1 >     2   A           1 /         1     300.000     318.110
        2          2 >     1   A           1 /         1     318.110     300.000
        3          3 >     4   A           1 /         1     337.110     357.100
        4          4 >     3   A           1 /         1     357.100     337.110
"""


class TestParseAcceptanceLog:
    def test_parses_poc_excerpt(self, tmp_path):
        log = tmp_path / "step3.log"
        log.write_text(POC_LOG_EXCERPT)
        events = parse_acceptance_log(log)
        # 8 events (4 replica entries x 2 patterns).
        assert len(events) == 8
        # First event: rep 1 trial with rep 0 (boundary), N (not accepted).
        assert events[0].step == 80
        assert events[0].replica_a == 1
        assert events[0].replica_b == 0
        assert events[0].accepted is False
        # Third event: rep 3 swap with rep 2, accepted.
        assert events[2].replica_a == 3
        assert events[2].replica_b == 2
        assert events[2].accepted is True
        # Pattern 1 events from step 160.
        assert events[4].step == 160
        assert events[4].pattern == 1
        assert events[4].accepted is True


class TestPlotAcceptanceRatio:
    def test_writes_png_and_csv(self, tmp_path):
        # Repeat the POC excerpt enough times to clear burn_in.
        log = tmp_path / "step3.log"
        log.write_text(POC_LOG_EXCERPT * 30)
        out_png = tmp_path / "acc.png"
        out_csv = tmp_path / "acc.csv"
        series = plot_acceptance_ratio(
            log, out_png, burn_in=10, out_csv=out_csv
        )
        assert out_png.exists()
        assert out_csv.exists()
        # Series has at least the 1<->2 and 3<->4 pairs.
        keys = set(series.keys())
        assert (1, 2) in keys or (2, 3) in keys

    def test_burn_in_drops_all_raises(self, tmp_path):
        log = tmp_path / "tiny.log"
        log.write_text(POC_LOG_EXCERPT)  # ~8 events
        with pytest.raises(ValueError, match="No exchange events"):
            plot_acceptance_ratio(log, tmp_path / "x.png", burn_in=1000)


# ---------------------------------------------------------------------------
# distance PMF
# ---------------------------------------------------------------------------

class TestComputeDistancePmf:
    def test_writes_xvg_and_png(self, tmp_path):
        # Synthetic Gaussian-like distribution.
        rng = np.random.default_rng(seed=42)
        distances = rng.normal(loc=3.0, scale=0.3, size=10_000)
        out_png = tmp_path / "pmf.png"
        out_xvg = tmp_path / "pmf.xvg"
        result = compute_distance_pmf(
            distances=distances, T_K=300.0,
            out_png=out_png, out_xvg=out_xvg, n_bins=40,
        )
        assert out_png.exists()
        assert out_xvg.exists()
        assert result.shape == (40, 2)
        # Min PMF == 0 by construction.
        assert result[:, 1].min() == pytest.approx(0.0, abs=1e-9)
        # Centre of distribution should be around 3.0 A and have lowest PMF.
        i_min = int(np.argmin(result[:, 1]))
        assert 2.7 < result[i_min, 0] < 3.3

    def test_empty_raises(self, tmp_path):
        with pytest.raises(ValueError, match="empty distance"):
            compute_distance_pmf(
                np.array([]), 300.0,
                tmp_path / "x.png", tmp_path / "x.xvg",
            )


# ---------------------------------------------------------------------------
# Genesis -> cpptraj mask conversion
# ---------------------------------------------------------------------------

class TestGenesisToCpptrajMask:
    def test_residue_with_atom(self):
        assert _genesis_to_cpptraj_mask("rno:96 and an:NZ") == ":96@NZ"

    def test_residue_only(self):
        assert _genesis_to_cpptraj_mask("rno:96") == ":96"

    def test_no_residue_raises(self):
        with pytest.raises(ValueError, match="rno:"):
            _genesis_to_cpptraj_mask("an:CA")


# ---------------------------------------------------------------------------
# remd_convert orchestration
# ---------------------------------------------------------------------------

class TestRunRemdConvert:
    def test_invokes_and_finds_output(self, tmp_path):
        cfg = GrestBuildConfig(
            input_pdb="/tmp/p.pdb",
            rest_selection=RESTSelectionSpec(mode="explicit", residues=["1"]),
            replica_temperatures=ReplicaTemperatureSpec(
                mode="manual", temperatures=[300.0, 357.10]
            ),
        )
        inp = tmp_path / "step5.inp"
        inp.write_text("# fake")
        workdir = tmp_path / "wd"

        def fake_run(cmd, cwd=None, capture=True, **kwargs):
            (Path(cwd) / "param1.dcd").write_text("# fake dcd")
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "", "")

        with patch(
            "abmptools.genesis.grest.analysis.run_command",
            side_effect=fake_run,
        ):
            out = run_remd_convert(cfg, inp, workdir)
        assert out == workdir / "param1.dcd"
        assert out.exists()

    def test_missing_output_raises(self, tmp_path):
        cfg = GrestBuildConfig(
            input_pdb="/tmp/p.pdb",
            rest_selection=RESTSelectionSpec(mode="explicit", residues=["1"]),
            replica_temperatures=ReplicaTemperatureSpec(
                mode="manual", temperatures=[300.0, 357.10]
            ),
        )
        inp = tmp_path / "step5.inp"
        inp.write_text("# fake")
        workdir = tmp_path / "wd"

        def fake_run_no_output(cmd, cwd=None, capture=True, **kwargs):
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "", "")

        from abmptools.genesis.grest._subprocess import CommandError
        with patch(
            "abmptools.genesis.grest.analysis.run_command",
            side_effect=fake_run_no_output,
        ):
            with pytest.raises(CommandError, match="missing"):
                run_remd_convert(cfg, inp, workdir)
