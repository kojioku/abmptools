# -*- coding: utf-8 -*-
"""Tests for abmptools.gro2udf.gro_parser."""
import pytest

from abmptools.gro2udf.gro_parser import GROParser


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

SINGLE_FRAME_GRO = """\
test.udf t= 0.00000 step= 0
3
    1RES     A1    1   0.100   0.200   0.300  0.0010  0.0020  0.0030
    1RES     A2    2   0.400   0.500   0.600  0.0040  0.0050  0.0060
    1RES     A3    3   0.700   0.800   0.900  0.0070  0.0080  0.0090
   1.000   2.000   3.000
"""

TWO_FRAME_GRO = """\
frame1 t= 0.00000 step= 0
1
    1RES     A1    1   0.100   0.200   0.300  0.0010  0.0020  0.0030
   1.000   2.000   3.000
frame2 t= 1.00000 step= 100
1
    1RES     A1    1   0.400   0.500   0.600  0.0040  0.0050  0.0060
   4.000   5.000   6.000
"""

NO_VELOCITY_GRO = """\
novelo t= 0.00000 step= 0
1
    1RES     A1    1   0.100   0.200   0.300
   1.000   2.000   3.000
"""


@pytest.fixture
def parser():
    return GROParser()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestParseFramesSingleFrame:
    """Parse a minimal single-frame .gro file."""

    def test_single_frame_basics(self, parser, tmp_path):
        gro = tmp_path / "single.gro"
        gro.write_text(SINGLE_FRAME_GRO)

        frames = list(parser.parse_frames(str(gro)))
        assert len(frames) == 1

        f = frames[0]
        assert f.title.startswith("test.udf")
        assert f.step == 0
        assert f.time == pytest.approx(0.0)
        assert len(f.atoms) == 3
        assert f.box_vals == pytest.approx([1.0, 2.0, 3.0])


class TestParseHeader:
    """Tests for GROParser._parse_header."""

    def test_standard_header(self, parser):
        step, time = parser._parse_header("test.udf t= 0.00000 step= 0\n")
        assert step == 0
        assert time == pytest.approx(0.0)

    def test_header_no_t_field(self, parser):
        """When 't=' is absent, time defaults to 0.0."""
        step, time = parser._parse_header("just a title step= 5\n")
        assert step == 5
        assert time == pytest.approx(0.0)


class TestParseAtomLine:
    """Tests for GROParser._parse_atom_line."""

    def test_atom_line_with_velocities(self, parser):
        line = "    1RES     A1    1   0.100   0.200   0.300  0.0010  0.0020  0.0030\n"
        atom = parser._parse_atom_line(line)
        assert atom.residue_col == "    1"
        assert atom.x == pytest.approx(0.100)
        assert atom.y == pytest.approx(0.200)
        assert atom.z == pytest.approx(0.300)
        assert atom.vx == pytest.approx(0.0010)
        assert atom.vy == pytest.approx(0.0020)
        assert atom.vz == pytest.approx(0.0030)

    def test_atom_line_without_velocities(self, parser):
        line = "    1RES     A1    1   0.100   0.200   0.300\n"
        atom = parser._parse_atom_line(line)
        assert atom.x == pytest.approx(0.100)
        assert atom.vx == pytest.approx(0.0)
        assert atom.vy == pytest.approx(0.0)
        assert atom.vz == pytest.approx(0.0)


class TestMultiFrame:
    """Parse a multi-frame .gro file."""

    def test_two_frames(self, parser, tmp_path):
        gro = tmp_path / "two.gro"
        gro.write_text(TWO_FRAME_GRO)

        frames = list(parser.parse_frames(str(gro)))
        assert len(frames) == 2
        assert frames[0].step == 0
        assert frames[1].step == 100
        assert frames[1].time == pytest.approx(1.0)


class TestEmptyFile:
    """An empty .gro file should yield no frames."""

    def test_empty_file(self, parser, tmp_path):
        gro = tmp_path / "empty.gro"
        gro.write_text("")

        frames = list(parser.parse_frames(str(gro)))
        assert frames == []
