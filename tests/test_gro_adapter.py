# -*- coding: utf-8 -*-
"""Tests for abmptools.gro2udf.gro_adapter."""
import math

import pytest

from abmptools.gro2udf.gro_parser import GROFrame, GROAtomLine
from abmptools.gro2udf.gro_adapter import _norm, _dot, GROAdapter


# ---------------------------------------------------------------------------
# Helper math functions
# ---------------------------------------------------------------------------

def test_norm():
    assert _norm([3, 4, 0]) == pytest.approx(5.0)


def test_dot_orthogonal():
    assert _dot([1, 0, 0], [0, 1, 0]) == pytest.approx(0.0)


def test_dot_general():
    assert _dot([1, 2, 3], [4, 5, 6]) == pytest.approx(32.0)


# ---------------------------------------------------------------------------
# Molecule boundary detection
# ---------------------------------------------------------------------------

def test_molecule_boundary_same_residue():
    """Two atoms with the same residue_col belong to 1 molecule."""
    frame = GROFrame(
        title="test", step=0, time=0.0,
        atoms=[
            GROAtomLine(residue_col="    1", x=0.1, y=0.2, z=0.3, vx=0, vy=0, vz=0),
            GROAtomLine(residue_col="    1", x=0.4, y=0.5, z=0.6, vx=0, vy=0, vz=0),
        ],
        box_vals=[1.0, 2.0, 3.0],
    )
    adapter = GROAdapter()
    positions, _ = adapter.to_positions_and_cell(frame)
    # Both atoms should have mol_id == 0
    assert positions[0].mol_id == 0
    assert positions[1].mol_id == 0
    assert positions[1].atom_id == 1


def test_molecule_boundary_different_residue():
    """Changing residue_col starts a new molecule."""
    frame = GROFrame(
        title="test", step=0, time=0.0,
        atoms=[
            GROAtomLine(residue_col="    1", x=0.1, y=0.2, z=0.3, vx=0, vy=0, vz=0),
            GROAtomLine(residue_col="    2", x=0.4, y=0.5, z=0.6, vx=0, vy=0, vz=0),
        ],
        box_vals=[1.0, 2.0, 3.0],
    )
    adapter = GROAdapter()
    positions, _ = adapter.to_positions_and_cell(frame)
    assert positions[0].mol_id == 0
    assert positions[1].mol_id == 1
    assert positions[1].atom_id == 0


# ---------------------------------------------------------------------------
# Cell geometry
# ---------------------------------------------------------------------------

def test_cell_from_3_value_box():
    """3-value box produces a rectangular cell (alpha=beta=gamma=90)."""
    frame = GROFrame(
        title="test", step=0, time=0.0,
        atoms=[
            GROAtomLine(residue_col="    1", x=0.1, y=0.2, z=0.3, vx=0, vy=0, vz=0),
        ],
        box_vals=[1.0, 2.0, 3.0],
    )
    adapter = GROAdapter()
    _, cell = adapter.to_positions_and_cell(frame)
    assert cell.a == pytest.approx(1.0)
    assert cell.b == pytest.approx(2.0)
    assert cell.c == pytest.approx(3.0)
    assert cell.alpha == pytest.approx(90.0)
    assert cell.beta == pytest.approx(90.0)
    assert cell.gamma == pytest.approx(90.0)


def test_cell_from_9_value_box():
    """9-value box produces a CellGeometry with computed angles."""
    # Orthogonal box: a=(1,0,0), b=(0,2,0), c=(0,0,3)
    # GROMACS order: v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y
    box_vals = [1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    frame = GROFrame(
        title="test", step=0, time=0.0,
        atoms=[
            GROAtomLine(residue_col="    1", x=0.1, y=0.2, z=0.3, vx=0, vy=0, vz=0),
        ],
        box_vals=box_vals,
    )
    adapter = GROAdapter()
    _, cell = adapter.to_positions_and_cell(frame)
    assert cell.a == pytest.approx(1.0)
    assert cell.b == pytest.approx(2.0)
    assert cell.c == pytest.approx(3.0)
    assert cell.alpha == pytest.approx(90.0)
    assert cell.beta == pytest.approx(90.0)
    assert cell.gamma == pytest.approx(90.0)


# ---------------------------------------------------------------------------
# Positions preserve coordinates
# ---------------------------------------------------------------------------

def test_positions_preserve_coordinates():
    frame = GROFrame(
        title="test", step=0, time=0.0,
        atoms=[
            GROAtomLine(residue_col="    1", x=0.1, y=0.2, z=0.3, vx=0.01, vy=0.02, vz=0.03),
        ],
        box_vals=[1.0, 2.0, 3.0],
    )
    adapter = GROAdapter()
    positions, _ = adapter.to_positions_and_cell(frame)
    pos = positions[0]
    assert pos.x == pytest.approx(0.1)
    assert pos.y == pytest.approx(0.2)
    assert pos.z == pytest.approx(0.3)
    assert pos.vx == pytest.approx(0.01)
    assert pos.vy == pytest.approx(0.02)
    assert pos.vz == pytest.approx(0.03)
