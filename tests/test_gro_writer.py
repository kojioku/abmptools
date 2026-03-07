# -*- coding: utf-8 -*-
"""Tests for abmptools.udf2gro.gromacs.writers.gro_writer."""
import math

import pytest

from abmptools.core.system_model import (
    AtomPosition,
    CellGeometry,
    SystemModel,
)
from abmptools.udf2gro.gromacs.writers.gro_writer import GroWriter


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _minimal_model(
    n_atoms: int = 3,
    cell: CellGeometry | None = None,
) -> SystemModel:
    """Return a SystemModel with *n_atoms* dummy positions."""
    positions = [
        AtomPosition(
            mol_id=1,
            mol_name_short="MOL",
            atom_gro_name=f"A{i+1:04d}",
            atom_id=i + 1,
            x=0.1 * (i + 1),
            y=0.2 * (i + 1),
            z=0.3 * (i + 1),
            vx=0.0,
            vy=0.0,
            vz=0.0,
        )
        for i in range(n_atoms)
    ]
    if cell is None:
        cell = CellGeometry(a=3.0, b=3.0, c=3.0)

    return SystemModel(
        title="test system",
        udf_path="/tmp/test.udf",
        comb_rule=2,
        flags14=0,
        fudgeLJ=0.5,
        fudgeQQ=0.8333,
        calcQQ=0,
        atom_positions=positions,
        cell=cell,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestGroWriter:

    def test_write_creates_file_with_correct_line_count(self, tmp_path):
        """Header + atom-count line + N atoms + box line = N+3 lines."""
        n_atoms = 3
        model = _minimal_model(n_atoms=n_atoms)
        filepath = str(tmp_path / "out.gro")

        GroWriter().write(model, filepath)

        with open(filepath) as f:
            lines = f.read().splitlines()

        # title, atom_count, 3 atom lines, box line
        assert len(lines) == n_atoms + 3

    def test_atom_count_line_matches_positions(self, tmp_path):
        """Second line of the .gro file must equal the number of atoms."""
        n_atoms = 5
        model = _minimal_model(n_atoms=n_atoms)
        filepath = str(tmp_path / "out.gro")

        GroWriter().write(model, filepath)

        with open(filepath) as f:
            lines = f.read().splitlines()

        assert int(lines[1]) == n_atoms

    def test_rectangular_cell_line_has_3_values(self, tmp_path):
        """A rectangular cell should produce a box line with exactly 3 values."""
        cell = CellGeometry(a=2.0, b=3.0, c=4.0)
        model = _minimal_model(cell=cell)
        filepath = str(tmp_path / "out.gro")

        GroWriter().write(model, filepath)

        with open(filepath) as f:
            lines = f.read().splitlines()

        box_line = lines[-1]
        values = box_line.split()
        assert len(values) == 3
        assert float(values[0]) == pytest.approx(2.0)
        assert float(values[1]) == pytest.approx(3.0)
        assert float(values[2]) == pytest.approx(4.0)

    def test_triclinic_cell_line_has_9_values(self, tmp_path):
        """A non-rectangular cell (alpha=80) should produce 9 box vector values."""
        cell = CellGeometry(a=2.0, b=3.0, c=4.0, alpha=80.0)
        model = _minimal_model(cell=cell)
        filepath = str(tmp_path / "out.gro")

        GroWriter().write(model, filepath)

        with open(filepath) as f:
            lines = f.read().splitlines()

        box_line = lines[-1]
        values = box_line.split()
        assert len(values) == 9

    def test_cell_line_static_rectangular(self):
        """Directly call _cell_line for a rectangular cell."""
        cell = CellGeometry(a=1.5, b=2.5, c=3.5)
        result = GroWriter._cell_line(cell)
        values = result.split()
        assert len(values) == 3
        assert float(values[0]) == pytest.approx(1.5)
        assert float(values[1]) == pytest.approx(2.5)
        assert float(values[2]) == pytest.approx(3.5)

    def test_file_content_roundtrip(self, tmp_path):
        """Write and read back: title, atom count, coordinates are recoverable."""
        model = _minimal_model(n_atoms=2)
        filepath = str(tmp_path / "out.gro")

        GroWriter().write(model, filepath)

        with open(filepath) as f:
            lines = f.read().splitlines()

        # Title
        assert lines[0] == "test system"
        # Atom count
        assert int(lines[1]) == 2
        # First atom line: check x coordinate (0.1) is present
        atom_line = lines[2]
        # x field is at chars 20-28 in gro format
        x_val = float(atom_line[20:28])
        assert x_val == pytest.approx(0.1, abs=1e-3)
