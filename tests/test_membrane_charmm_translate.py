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
        # HID / HIE → HSD / HSE
        assert translate_residue_name("HID") == "HSD"
        assert translate_residue_name("HIE") == "HSE"
        # HIP is shared between AMBER and the Klauda CHARMM36 port
        # (both [ HIP ] and [ HSP ] defined); we keep HIP verbatim.
        assert translate_residue_name("HIP") == "HIP"

    def test_protomers(self):
        # AMBER → CHARMM port: ASH/GLH/LYN need translation
        assert translate_residue_name("ASH") == "ASPP"
        assert translate_residue_name("GLH") == "GLUP"
        assert translate_residue_name("LYN") == "LSN"
        # CYM is preserved as-is (port has [ CYM ] directly)
        assert translate_residue_name("CYM") == "CYM"

    def test_caps(self):
        # Both ACE and NME keep their AMBER names in the Klauda
        # GROMACS port (only atom names change for ACE; NME atoms also
        # need only the universal H→HN backbone rename).
        assert translate_residue_name("ACE") == "ACE"
        assert translate_residue_name("NME") == "NME"

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
        # Klauda CHARMM36 port ACE atoms: CH3 / HH31 / HH32 / HH33 / C / O.
        # AMBER tleap-built ACE has H1/H2/H3 instead of HH3?, others same.
        assert translate_atom_name("ACE", "H1") == "HH31"
        assert translate_atom_name("ACE", "H2") == "HH32"
        assert translate_atom_name("ACE", "H3") == "HH33"
        # CH3, C, O — pass through (port keeps the AMBER-style names)
        assert translate_atom_name("ACE", "CH3") == "CH3"
        assert translate_atom_name("ACE", "C") == "C"
        assert translate_atom_name("ACE", "O") == "O"

    def test_nme_atoms(self):
        # NME (kept as 'NME', not 'CT3') needs only the universal
        # backbone H→HN rename. CH3 / HH31 / HH32 / HH33 / N pass through.
        assert translate_atom_name("NME", "H") == "HN"
        assert translate_atom_name("NME", "CH3") == "CH3"
        assert translate_atom_name("NME", "HH31") == "HH31"
        assert translate_atom_name("NME", "HH32") == "HH32"
        assert translate_atom_name("NME", "N") == "N"

    def test_ace_h_passes_through(self):
        # ACE has no backbone amide H (it's a cap providing C=O), but
        # if a literal "H" appeared on ACE it would NOT become HN
        # because ACE is in SKIP_BACKBONE_RENAME.
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
        ATOM      2  H2  ACE     1       1.500   2.500   0.866  1.00  0.00
        ATOM      3  H3  ACE     1       1.500   2.500  -0.866  1.00  0.00
        ATOM      4  CH3 ACE     1       2.000   2.090   0.000  1.00  0.00
        ATOM      5  C   ACE     1       3.427   2.641   0.000  1.00  0.00
        ATOM      6  O   ACE     1       4.391   1.877   0.000  1.00  0.00
        ATOM      7  N   ALA     2       3.555   3.970   0.000  1.00  0.00
        ATOM      8  H   ALA     2       2.733   4.556   0.000  1.00  0.00
        ATOM      9  CA  ALA     2       4.853   4.614   0.000  1.00  0.00
        ATOM     10  CB  ALA     2       4.853   5.450   1.260  1.00  0.00
        ATOM     11  C   ALA     2       6.063   3.700   0.000  1.00  0.00
        ATOM     12  O   ALA     2       6.063   2.481   0.000  1.00  0.00
        ATOM     13  N   NME     3       7.176   4.430   0.000  1.00  0.00
        ATOM     14  H   NME     3       7.176   5.448   0.000  1.00  0.00
        ATOM     15  CH3 NME     3       8.421   3.728   0.000  1.00  0.00
        ATOM     16 HH31 NME     3       8.421   3.111   0.890  1.00  0.00
        ATOM     17  Na+ Na+     4      10.000  10.000   0.000  1.00  0.00
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

        # ACE: H1/H2/H3 → HH31/HH32/HH33; CH3/C/O pass through
        # (Klauda CHARMM36 port keeps AMBER-style ACE atom names).
        assert "HH31" in text
        assert "HH32" in text
        assert "HH33" in text
        # ACE residue name kept (only atoms change)
        # 6 ACE atoms in the fixture: H1/H2/H3/CH3/C/O
        assert text.count("ACE") == 6
        # Old (incorrect) classic-CHARMM names should NOT appear
        assert "CAY" not in text
        assert "HY1" not in text

        # ALA: backbone amide H → HN (only one ALA, has 1 amide H)
        assert "HN  " in text

        # NME stays NME (residue name unchanged, only H→HN applied)
        assert text.count("NME") >= 1
        # Old (incorrect) NME→CT3 rename must NOT have happened
        assert "CT3" not in text
        assert "CAT" not in text
        assert " HT1" not in text

        # Na+ → SOD (residue name and atom name both)
        assert "SOD" in text
        # The original "Na+" should be gone everywhere
        assert "Na+" not in text

    def test_4char_residue_not_truncated(self, tmp_path):
        """Regression: POPC / TIP3 / etc. must not be truncated to POP / TIP.

        packmol-memgen --charmm output uses 4-char residue names that fill
        cols 18-21 (extending into the alt-loc / chain field). The
        translator must preserve the 4-char width even for residues that
        are NOT in the AMBER→CHARMM rename map (like POPC, which
        passes through verbatim). Earlier code fell back to a 3-char
        ``line[17:20].strip()`` read which produced "POP", causing
        pdb2gmx to fail with "Residue 'POP' not found in residue
        topology database".
        """
        src = tmp_path / "src.pdb"
        # packmol-memgen --charmm style: 4-char POPC + chain 'A'
        src.write_text(
            "ATOM      1  N   POPCA   1       3.997  -9.117 -18.250  1.00  0.00\n"
            "ATOM      2  C12 POPCA   1       2.786  -9.916 -17.797  1.00  0.00\n"
            "ATOM      3  OW  TIP3A   2       0.000   0.000   0.000  1.00  0.00\n"
        )
        dst = tmp_path / "dst.pdb"
        translate_pdb_amber_to_charmm(input_pdb=str(src), output_pdb=str(dst))
        out = dst.read_text()
        # 4-char residue names preserved; POP must NOT appear, POPC must
        assert "POPC" in out
        assert "POP " not in out  # POP followed by space = truncated form
        # TIP3 likewise (already a CHARMM name)
        assert "TIP3" in out
        assert "TIP " not in out

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
        # ATOM line is translated (H1 → HH31 per Klauda port convention)
        assert "HH31" in out
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
        # ACE has no backbone amide H, ions / water have no protein
        # backbone, so they all skip the universal H→HN rename.
        # NME is intentionally NOT in this list — its amide H must
        # become HN, and the universal map handles it.
        for r in ("ACE", "TIP3", "SOD", "CLA"):
            assert r in SKIP_BACKBONE_RENAME, (
                f"{r!r} should skip universal backbone H→HN renaming"
            )
        assert "NME" not in SKIP_BACKBONE_RENAME, (
            "NME must use the universal H→HN rename"
        )
