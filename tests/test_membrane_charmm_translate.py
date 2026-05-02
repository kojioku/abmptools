# -*- coding: utf-8 -*-
"""
tests/test_membrane_charmm_translate.py
---------------------------------------
Pure unit tests for the AMBER-to-CHARMM translation helpers in
``abmptools.membrane.parameterize_charmm``. None of these touch
GROMACS, tleap, or the CHARMM36 .ff dir, so they run quickly and
independently of CHARMM36 availability.
"""
from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import pytest

from abmptools.membrane.models import IonSpec
from abmptools.membrane.parameterize_charmm import (
    AMBER_TO_CHARMM_RESNAME,
    PER_RESIDUE_ATOM_MAP,
    SKIP_BACKBONE_RENAME,
    UNIVERSAL_BACKBONE_ATOM_MAP,
    ff_name_from_dir,
    translate_atom_name,
    translate_ion_name,
    translate_ions_for_charmm,
    translate_pdb_amber_to_charmm,
    translate_residue_name,
)


# ---------------------------------------------------------------------------
# residue name translation
# ---------------------------------------------------------------------------

class TestTranslateResidueName:
    def test_histidine_tautomers(self):
        assert translate_residue_name("HID") == "HSD"
        assert translate_residue_name("HIE") == "HSE"
        assert translate_residue_name("HIP") == "HSP"

    def test_protomers(self):
        assert translate_residue_name("ASH") == "ASPP"
        assert translate_residue_name("GLH") == "GLUP"
        assert translate_residue_name("LYN") == "LSN"

    def test_caps(self):
        # ACE keeps its name (only atom names change)
        assert translate_residue_name("ACE") == "ACE"
        # NME becomes CT3
        assert translate_residue_name("NME") == "CT3"

    def test_water(self):
        assert translate_residue_name("WAT") == "TIP3"
        assert translate_residue_name("HOH") == "TIP3"

    def test_ions(self):
        assert translate_residue_name("Na+") == "SOD"
        assert translate_residue_name("Cl-") == "CLA"
        assert translate_residue_name("K+") == "POT"

    def test_pass_through_for_unknown(self):
        # Standard 20 amino acids share names → no change
        for aa in ("ALA", "GLY", "LEU", "PHE", "ARG", "ASP", "GLU"):
            assert translate_residue_name(aa) == aa
        # Unknown residue names are returned unchanged
        assert translate_residue_name("XYZ") == "XYZ"

    def test_strips_whitespace(self):
        assert translate_residue_name(" HIE ") == "HSE"


# ---------------------------------------------------------------------------
# atom name translation
# ---------------------------------------------------------------------------

class TestTranslateAtomName:
    def test_universal_backbone_h_to_hn(self):
        # All standard amino acids: backbone H → HN
        for aa in ("ALA", "GLY", "LEU", "PHE", "HSE", "ASP", "LYS"):
            assert translate_atom_name(aa, "H") == "HN"

    def test_ace_atoms(self):
        assert translate_atom_name("ACE", "H1") == "HY1"
        assert translate_atom_name("ACE", "H2") == "HY2"
        assert translate_atom_name("ACE", "H3") == "HY3"
        assert translate_atom_name("ACE", "CH3") == "CAY"
        assert translate_atom_name("ACE", "C") == "CY"
        assert translate_atom_name("ACE", "O") == "OY"

    def test_ct3_atoms(self):
        # Note: NME → CT3 residue rename happens before atom translation.
        assert translate_atom_name("CT3", "CH3") == "CAT"
        assert translate_atom_name("CT3", "HH31") == "HT1"
        assert translate_atom_name("CT3", "HH32") == "HT2"
        assert translate_atom_name("CT3", "HH33") == "HT3"
        # Universal H→HN also applies to CT3's amide H? No — CT3 IS a cap;
        # it's in SKIP_BACKBONE_RENAME so the universal map is bypassed.
        # The CT3-specific map doesn't list H, so it passes through.
        assert translate_atom_name("CT3", "H") == "H"

    def test_caps_dont_get_universal_h_to_hn(self):
        # ACE's H1/H2/H3 are caught by per-residue map (HY1/HY2/HY3),
        # but if we somehow saw a literal "H" on ACE, it would NOT
        # become HN because ACE is in SKIP_BACKBONE_RENAME.
        assert translate_atom_name("ACE", "H") == "H"

    def test_ions(self):
        assert translate_atom_name("SOD", "Na+") == "SOD"
        assert translate_atom_name("SOD", "NA") == "SOD"
        assert translate_atom_name("CLA", "Cl-") == "CLA"
        assert translate_atom_name("POT", "K+") == "POT"
        assert translate_atom_name("MG", "Mg2+") == "MG"
        assert translate_atom_name("CAL", "Ca2+") == "CAL"

    def test_water_o_to_oh2(self):
        assert translate_atom_name("TIP3", "O") == "OH2"
        # H1 / H2 stay the same in TIP3
        assert translate_atom_name("TIP3", "H1") == "H1"
        assert translate_atom_name("TIP3", "H2") == "H2"

    def test_pass_through_for_other_atoms(self):
        # CA, CB, etc. are unchanged
        assert translate_atom_name("ALA", "CA") == "CA"
        assert translate_atom_name("ALA", "CB") == "CB"
        assert translate_atom_name("ALA", "HB1") == "HB1"
        assert translate_atom_name("ARG", "NE") == "NE"


# ---------------------------------------------------------------------------
# ion translation
# ---------------------------------------------------------------------------

class TestTranslateIon:
    def test_translate_ion_name(self):
        assert translate_ion_name("Na+") == "SOD"
        assert translate_ion_name("Cl-") == "CLA"
        assert translate_ion_name("K+") == "POT"
        assert translate_ion_name("Mg2+") == "MG"
        assert translate_ion_name("Ca2+") == "CAL"

    def test_already_charmm_passes_through(self):
        # Explicit CHARMM names are respected (override path).
        assert translate_ion_name("SOD") == "SOD"
        assert translate_ion_name("CLA") == "CLA"
        assert translate_ion_name("POT") == "POT"

    def test_translate_ions_for_charmm_keeps_other_fields(self):
        original = IonSpec(cation="Na+", anion="Cl-",
                           salt_concentration_M=0.20, neutralize=True)
        translated = translate_ions_for_charmm(original)
        assert translated.cation == "SOD"
        assert translated.anion == "CLA"
        assert translated.salt_concentration_M == 0.20
        assert translated.neutralize is True

    def test_translate_ions_with_explicit_charmm(self):
        # Already-CHARMM IonSpec is unchanged (override path).
        original = IonSpec(cation="POT", anion="CLA")
        translated = translate_ions_for_charmm(original)
        assert translated.cation == "POT"
        assert translated.anion == "CLA"


# ---------------------------------------------------------------------------
# whole-PDB translator
# ---------------------------------------------------------------------------

@pytest.fixture
def amber_peptide_pdb(tmp_path: Path) -> Path:
    """A minimal AMBER-style PDB: ACE-ALA-NME with sodium ion."""
    pdb = tmp_path / "input.pdb"
    pdb.write_text(dedent("""\
        ATOM      1  H1  ACE     1       2.000   1.000   0.000  1.00  0.00
        ATOM      2  CH3 ACE     1       2.000   2.090   0.000  1.00  0.00
        ATOM      3  C   ACE     1       3.427   2.641   0.000  1.00  0.00
        ATOM      4  O   ACE     1       4.391   1.877   0.000  1.00  0.00
        ATOM      5  N   ALA     2       3.555   3.970   0.000  1.00  0.00
        ATOM      6  H   ALA     2       2.733   4.556   0.000  1.00  0.00
        ATOM      7  CA  ALA     2       4.853   4.614   0.000  1.00  0.00
        ATOM      8  CB  ALA     2       4.853   5.450   1.260  1.00  0.00
        ATOM      9  C   ALA     2       6.063   3.700   0.000  1.00  0.00
        ATOM     10  O   ALA     2       6.063   2.481   0.000  1.00  0.00
        ATOM     11  N   NME     3       7.176   4.430   0.000  1.00  0.00
        ATOM     12  H   NME     3       7.176   5.448   0.000  1.00  0.00
        ATOM     13  CH3 NME     3       8.421   3.728   0.000  1.00  0.00
        ATOM     14 HH31 NME     3       8.421   3.111   0.890  1.00  0.00
        ATOM     15  Na+ Na+     4      10.000  10.000   0.000  1.00  0.00
    """))
    return pdb


class TestTranslatePdbAmberToCharmm:
    def test_full_translation_pipeline(self, amber_peptide_pdb, tmp_path):
        out_pdb = tmp_path / "output.pdb"
        translate_pdb_amber_to_charmm(
            input_pdb=str(amber_peptide_pdb),
            output_pdb=str(out_pdb),
        )
        text = out_pdb.read_text()

        # ACE atoms translated (atom name field is 4-char left-justified,
        # so "HY1 " "CAY " "CY  " "OY  " appear with trailing spaces).
        assert "HY1 " in text
        assert "CAY " in text
        assert "CY  " in text
        assert "OY  " in text
        # All four ACE atoms keep their residue tag
        assert text.count("ACE") == 4

        # ALA: backbone amide H → HN (only one ALA, has 1 amide H)
        assert "HN  " in text

        # NME → CT3 + atom translations
        assert "CT3" in text
        assert "CAT " in text
        assert "HT1 " in text

        # Na+ → SOD (residue name and atom name both)
        assert "SOD" in text
        # The original "Na+" should be gone everywhere
        assert "Na+" not in text

    def test_non_atom_lines_pass_through(self, tmp_path):
        # HEADER / REMARK / END / blank lines are not modified
        src = tmp_path / "src.pdb"
        src.write_text(dedent("""\
            HEADER    PEPTIDE                                01-JAN-26
            REMARK    Generated by abmptools tests
            ATOM      1  H1  ACE     1       2.000   1.000   0.000  1.00  0.00
            END
        """))
        dst = tmp_path / "dst.pdb"
        translate_pdb_amber_to_charmm(input_pdb=str(src), output_pdb=str(dst))
        out = dst.read_text()
        assert "HEADER    PEPTIDE" in out
        assert "REMARK    Generated by abmptools tests" in out
        assert "END" in out
        # ATOM line is translated (only H1 in this minimal case → HY1)
        assert "HY1" in out
        # H1 (the AMBER name) should be gone
        assert " H1 " not in out


# ---------------------------------------------------------------------------
# force-field path helpers
# ---------------------------------------------------------------------------

class TestFfNameFromDir:
    def test_strip_dot_ff_suffix(self):
        assert ff_name_from_dir("/path/to/charmm36-jul2022.ff") == "charmm36-jul2022"
        assert ff_name_from_dir("./charmm36-feb2021.ff") == "charmm36-feb2021"

    def test_no_suffix_passes_through(self):
        # In case someone passes the bare name (unusual)
        assert ff_name_from_dir("/path/to/charmm36-jul2022") == "charmm36-jul2022"


# ---------------------------------------------------------------------------
# table integrity (catch obvious typos)
# ---------------------------------------------------------------------------

class TestTableIntegrity:
    def test_residue_map_no_self_loops(self):
        # Mapping a residue to itself is suspicious (we should just omit it)
        # except for ACE which we keep on the map because the *atoms* change.
        for amber, charmm in AMBER_TO_CHARMM_RESNAME.items():
            if amber == charmm and amber != "ACE":
                pytest.fail(f"{amber!r} maps to itself; remove the entry")

    def test_universal_backbone_only_h(self):
        # Currently we only need H → HN. If we add more, extend tests.
        assert "H" in UNIVERSAL_BACKBONE_ATOM_MAP
        assert UNIVERSAL_BACKBONE_ATOM_MAP["H"] == "HN"

    def test_skip_list_contains_caps_and_ions(self):
        for r in ("ACE", "CT3", "TIP3", "SOD", "CLA"):
            assert r in SKIP_BACKBONE_RENAME, (
                f"{r!r} should skip universal backbone H→HN renaming"
            )
