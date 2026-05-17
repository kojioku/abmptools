"""Tests for the force-field-agnostic functional tag mapping."""
import pytest

from abmptools.hbond.func_tags import (
    BUILTIN_MAPPINGS, CHARMM36, GAFF2, OPLS_AA,
    detect_force_field, get_mapping, tag_atoms,
    TAG_CARBONYL_C, TAG_CARBONYL_O, TAG_HYDROXYL_O, TAG_HYDROXYL_H,
    TAG_AMIDE_N, TAG_AMIDE_H,
)


def test_gaff2_basic_tags():
    assert GAFF2.get_tag("c") == TAG_CARBONYL_C
    assert GAFF2.get_tag("oh") == TAG_HYDROXYL_O
    assert GAFF2.get_tag("ho") == TAG_HYDROXYL_H
    assert GAFF2.get_tag("o") == TAG_CARBONYL_O
    assert GAFF2.get_tag("n") == TAG_AMIDE_N
    assert GAFF2.get_tag("hn") == TAG_AMIDE_H
    assert GAFF2.get_tag("unknown_type") is None


def test_charmm36_basic_tags():
    assert CHARMM36.get_tag("CC") == TAG_CARBONYL_C
    assert CHARMM36.get_tag("OH1") == TAG_HYDROXYL_O
    assert CHARMM36.get_tag("NH1") == TAG_AMIDE_N
    assert CHARMM36.get_tag("H") == TAG_AMIDE_H


def test_opls_aa_basic_tags():
    assert OPLS_AA.get_tag("opls_267") == TAG_CARBONYL_C
    assert OPLS_AA.get_tag("opls_268") == TAG_HYDROXYL_O
    assert OPLS_AA.get_tag("opls_270") == TAG_HYDROXYL_H


def test_get_mapping_case_insensitive():
    assert get_mapping("GAFF2").force_field == "GAFF2"
    assert get_mapping("gaff2").force_field == "GAFF2"
    assert get_mapping("OPLS-AA").force_field == "OPLS-AA"
    assert get_mapping("opls").force_field == "OPLS-AA"
    assert get_mapping("CHARMM36").force_field == "CHARMM36"
    with pytest.raises(KeyError):
        get_mapping("UnknownFF99")


def test_detect_force_field_gaff2():
    """GAFF2-style atom types (lowercase) should be detected as GAFF2."""
    types = ["c", "c3", "ca", "oh", "ho", "o", "os", "n", "hn", "hc"]
    assert detect_force_field(types) == "GAFF2"


def test_detect_force_field_charmm36():
    types = ["CT1", "CT2", "OH1", "NH1", "H", "OB", "CC"]
    assert detect_force_field(types) == "CHARMM36"


def test_detect_force_field_opls():
    types = ["opls_267", "opls_268", "opls_269", "opls_270"]
    assert detect_force_field(types) == "OPLS-AA"


def test_tag_atoms_with_unknowns():
    """Unmapped atom types should yield None tags."""
    result = tag_atoms(["c", "unknown", "oh"], GAFF2)
    assert result == [TAG_CARBONYL_C, None, TAG_HYDROXYL_O]


def test_add_mapping_validation():
    """Adding a tag should reject unknown tag strings."""
    mapping = GAFF2  # not a copy; just probe
    with pytest.raises(ValueError):
        mapping.add_mapping("X9", "not_a_real_tag")
