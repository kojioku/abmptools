# -*- coding: utf-8 -*-
"""Tests for abmptools.amorphous.ndx_writer module."""
import os

import pytest

from abmptools.amorphous.ndx_writer import write_ndx


class TestWriteNdx:
    """Tests for the write_ndx function."""

    def test_single_component_groups(self, tmp_path):
        """Single component produces System group and one component group."""
        out = tmp_path / "index.ndx"
        write_ndx(["Water"], [3], [2], str(out))
        content = out.read_text()
        assert "[ System ]" in content
        assert "[ Water ]" in content

    def test_two_components_sequential_ids(self, tmp_path):
        """Two components have sequential, non-overlapping atom IDs."""
        out = tmp_path / "index.ndx"
        write_ndx(["A", "B"], [2, 3], [1, 1], str(out))
        content = out.read_text()

        # Parse groups
        groups = {}
        current = None
        for line in content.splitlines():
            line = line.strip()
            if line.startswith("["):
                current = line.strip("[] ")
                groups[current] = []
            elif current and line:
                groups[current].extend(int(x) for x in line.split())

        # A: 1 mol * 2 atoms -> [1, 2]; B: 1 mol * 3 atoms -> [3, 4, 5]
        assert groups["A"] == [1, 2]
        assert groups["B"] == [3, 4, 5]

    def test_output_file_exists(self, tmp_path):
        """Output file is created on disk."""
        out = tmp_path / "result.ndx"
        write_ndx(["X"], [5], [3], str(out))
        assert out.exists()

    def test_system_group_contains_all_ids(self, tmp_path):
        """System group contains every atom ID from all components."""
        out = tmp_path / "index.ndx"
        # 2 components: 4 atoms/mol * 3 mol + 2 atoms/mol * 5 mol = 12 + 10 = 22
        write_ndx(["P", "Q"], [4, 2], [3, 5], str(out))
        content = out.read_text()

        groups = {}
        current = None
        for line in content.splitlines():
            line = line.strip()
            if line.startswith("["):
                current = line.strip("[] ")
                groups[current] = []
            elif current and line:
                groups[current].extend(int(x) for x in line.split())

        total_atoms = 4 * 3 + 2 * 5
        assert len(groups["System"]) == total_atoms
        assert groups["System"] == list(range(1, total_atoms + 1))

    def test_returns_absolute_path(self, tmp_path):
        """write_ndx returns an absolute path string."""
        out = tmp_path / "abs.ndx"
        result = write_ndx(["M"], [1], [1], str(out))
        assert os.path.isabs(result)
        assert result == str(out.resolve())
