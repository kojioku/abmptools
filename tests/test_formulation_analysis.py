# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.analysis."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def test_per_frame_clusters_isolated():
    pytest.importorskip("networkx")
    from abmptools.formulation.analysis.aggregate_transition import per_frame_clusters
    coms = np.array([[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]])
    labels = per_frame_clusters(coms, cutoff_nm=0.5)
    assert labels[0] != labels[1]


def test_per_frame_clusters_aggregated():
    pytest.importorskip("networkx")
    from abmptools.formulation.analysis.aggregate_transition import per_frame_clusters
    coms = np.array([[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]])
    labels = per_frame_clusters(coms, cutoff_nm=0.5)
    # All three are within cutoff (chain 0-0.3-0.6, each link <= 0.5)
    assert len(set(labels)) == 1


def test_per_frame_clusters_partial():
    pytest.importorskip("networkx")
    from abmptools.formulation.analysis.aggregate_transition import per_frame_clusters
    coms = np.array([[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [10.0, 10.0, 10.0]])
    labels = per_frame_clusters(coms, cutoff_nm=0.5)
    # Two clusters: {0, 1} and {2}
    assert len(set(labels)) == 2
    assert labels[0] == labels[1]
    assert labels[2] != labels[0]


def test_cluster_size_per_peptide():
    from abmptools.formulation.analysis.aggregate_transition import cluster_size_per_peptide
    labels = [0, 0, 1, 0, 2]
    sizes = cluster_size_per_peptide(labels)
    assert sizes == [3, 3, 1, 3, 1]


def test_parse_dssp_xpm_synthetic(tmp_path):
    from abmptools.formulation.analysis.secondary_structure import parse_dssp_xpm
    xpm = tmp_path / "ss.xpm"
    xpm.write_text(
        "static char *gromacs_xpm[] = {\n"
        '"HHHH",\n'  # residue 1: 100% H
        '"EEEE",\n'  # residue 2: 100% E
        '"CCCC",\n'  # residue 3: 100% coil
        "};\n"
    )
    fr = parse_dssp_xpm(str(xpm))
    assert fr["H_fraction"].shape == (4,)
    # Each frame: 1H + 1E + 1C → fractions 1/3, 1/3, 1/3
    assert abs(fr["H_fraction"][0] - 1.0 / 3.0) < 1e-9
    assert abs(fr["E_fraction"][0] - 1.0 / 3.0) < 1e-9
    assert abs(fr["C_fraction"][0] - 1.0 / 3.0) < 1e-9


def test_parse_dssp_xpm_empty(tmp_path):
    from abmptools.formulation.analysis.secondary_structure import parse_dssp_xpm
    xpm = tmp_path / "ss.xpm"
    xpm.write_text("header but no data\n")
    fr = parse_dssp_xpm(str(xpm))
    assert fr["H_fraction"].size == 0


def test_plot_aggregate_size_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    from abmptools.formulation.analysis.plots import plot_aggregate_size
    csv = tmp_path / "agg.csv"
    csv.write_text(
        "frame,max_cluster_size,n_clusters\n"
        "0,3,2\n"
        "1,4,1\n"
        "2,4,1\n"
    )
    png = tmp_path / "agg.png"
    plot_aggregate_size(str(csv), str(png))
    assert png.is_file()


def test_plot_contact_heatmap_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    from abmptools.formulation.analysis.plots import plot_contact_heatmap
    npy = tmp_path / "cm.npy"
    np.save(npy, np.array([[0.1, 0.2], [0.3, 0.5]]))
    png = tmp_path / "cm.png"
    plot_contact_heatmap(
        str(npy), str(png),
        row_labels=["R1", "R2"], col_labels=["CPR", "TCH"],
    )
    assert png.is_file()


def test_run_gmx_dssp_invokes_subprocess(tmp_path):
    from unittest.mock import patch
    from abmptools.formulation.analysis.secondary_structure import run_gmx_dssp
    from subprocess import CompletedProcess

    with patch(
        "abmptools.formulation.analysis.secondary_structure.run_command"
    ) as mock:
        mock.return_value = CompletedProcess(
            ["gmx", "dssp"], 0, "", ""
        )
        out = tmp_path / "ss.xpm"
        run_gmx_dssp(traj="t.xtc", tpr="t.tpr", out_xpm=str(out))
        assert mock.called
        argv = mock.call_args.args[0]
        assert "dssp" in argv


def test_run_gmx_sasa_invokes_subprocess(tmp_path):
    from unittest.mock import patch
    from abmptools.formulation.analysis.sasa import run_gmx_sasa
    from subprocess import CompletedProcess

    with patch("abmptools.formulation.analysis.sasa.run_command") as mock:
        mock.return_value = CompletedProcess(["gmx", "sasa"], 0, "", "")
        run_gmx_sasa(traj="t.xtc", tpr="t.tpr", out_xvg=str(tmp_path / "s.xvg"))
        assert mock.called


def test_run_gmx_hbond_invokes_subprocess(tmp_path):
    from unittest.mock import patch
    from abmptools.formulation.analysis.hbond import run_gmx_hbond
    from subprocess import CompletedProcess

    with patch("abmptools.formulation.analysis.hbond.run_command") as mock:
        mock.return_value = CompletedProcess(["gmx", "hbond"], 0, "", "")
        run_gmx_hbond(
            traj="t.xtc", tpr="t.tpr",
            out_xvg=str(tmp_path / "h.xvg"),
            donor_group="Peptide", acceptor_group="Solvent",
        )
        assert mock.called
        # stdin contains the group selection
        kwargs = mock.call_args.kwargs
        assert "Peptide" in kwargs["input_text"]
