"""Enhanced PIEDA (ABINIT-MP Ver.2 Rev.8, &LRD DISP='ON' + ES_RESP='YES').

The enhanced table inserts two columns into PIEDA::

    IJ-PAIR  ES        EX  CT+mix  DI(MP2)        q(I=>J)     Ver.1
    IJ-PAIR  ES(RESP)  ES  EX      CT+mix  DI(LRD)  Erest  q(I=>J)   Ver.2 Rev.8

``DI(LRD)`` is the dispersion proper and ``Erest`` the remaining correlation;
their sum is the old ``DI``. Reading such a log with the old fixed layout
returned every value one slot out without raising anything, which is what
these tests exist to prevent.
"""
from __future__ import annotations

import os

import pytest

from abmptools.anlfmo import anlfmo, pieda_columns_from_header

HARTREE = 627.5095
EXCERPT = os.path.join(os.path.dirname(__file__), "data",
                       "enhanced_pieda_excerpt.log")

STANDARD_HEADER = ("           IJ-PAIR       ES             EX             "
                   "CT+mix         DI(MP2)        q(I=>J)  ")
ENHANCED_HEADER = ("                 IJ-PAIR       ES(RESP)       ES        "
                   "     EX             CT+mix        DI(LRD)        Erest  "
                   "        q(I=>J)  ")


class TestPiedaColumnsFromHeader:
    def test_standard_layout(self):
        assert pieda_columns_from_header(STANDARD_HEADER) == [
            "I", "J", "ES", "EX", "CT-mix", "DI(MP2)", "q(I=>J)"]

    def test_enhanced_layout(self):
        assert pieda_columns_from_header(ENHANCED_HEADER) == [
            "I", "J", "ES(RESP)", "ES", "EX", "CT-mix", "DI(LRD)", "Erest",
            "q(I=>J)"]

    def test_hf_layout_keeps_the_old_column_name(self):
        """An HF log writes a bare ``DI``; it used to arrive as DI(MP2)."""
        assert pieda_columns_from_header("IJ-PAIR ES EX CT+mix DI q(I=>J)") == [
            "I", "J", "ES", "EX", "CT-mix", "DI(MP2)", "q(I=>J)"]

    def test_pb_layout(self):
        assert "Solv(ES)" in pieda_columns_from_header(
            "IJ-PAIR ES EX CT+mix Solv(ES) DI(MP2) q(I=>J)")

    def test_an_unknown_field_raises_instead_of_shifting(self):
        with pytest.raises(ValueError, match="unknown PIEDA header"):
            pieda_columns_from_header("IJ-PAIR ES EX CT+mix DI(XYZ) q(I=>J)")

    def test_a_header_without_ij_pair_raises(self):
        with pytest.raises(ValueError, match="does not start with IJ-PAIR"):
            pieda_columns_from_header("ES EX CT+mix DI(MP2) q(I=>J)")


class TestReadEnhancedLog:
    @pytest.fixture()
    def parsed(self):
        obj = anlfmo()
        rows = obj.read_pieda(EXCERPT)
        return obj, obj.getpiedadf(rows)

    def test_layout_is_taken_from_the_log(self, parsed):
        obj, df = parsed
        assert obj.pcolumn == ["I", "J", "ES(RESP)", "ES", "EX", "CT-mix",
                               "DI(LRD)", "Erest", "q(I=>J)"]
        assert list(df.columns) == obj.pcolumn
        assert len(df) == 6

    def test_values_land_in_the_right_columns(self, parsed):
        _, df = parsed
        row = df[(df["I"] == 2) & (df["J"] == 1)].iloc[0]
        assert row["ES(RESP)"] == pytest.approx(-2.141227)
        assert row["ES"] == pytest.approx(-3.479361)
        assert row["DI(LRD)"] == pytest.approx(-7.363263)
        assert row["Erest"] == pytest.approx(-1.692311)
        assert row["q(I=>J)"] == pytest.approx(-0.000562)

    def test_es_ex_ct_reproduce_hf_ifie(self, parsed):
        """The Hartree-Fock part of the pair energy must add back up."""
        _, df = parsed
        row = df[(df["I"] == 2) & (df["J"] == 1)].iloc[0]
        hf_ifie_hartree = -0.003251
        assert row["ES"] + row["EX"] + row["CT-mix"] == pytest.approx(
            hf_ifie_hartree * HARTREE, abs=2e-3)

    def test_dispersion_plus_rest_reproduce_the_correlation(self, parsed):
        """DI(LRD) + Erest is the old DI: the whole MP2 correlation."""
        _, df = parsed
        row = df[(df["I"] == 2) & (df["J"] == 1)].iloc[0]
        mp2_corr_hartree = -0.014431
        assert row["DI(LRD)"] + row["Erest"] == pytest.approx(
            mp2_corr_hartree * HARTREE, abs=2e-3)

    def test_di_column_finds_the_lrd_term(self, parsed):
        obj, df = parsed
        assert obj.di_column(df) == "DI(LRD)"

    def test_sum_terms_include_the_new_components(self, parsed):
        obj, df = parsed
        terms = obj.sum_terms(df)
        assert "DI(LRD)" in terms and "Erest" in terms and "ES(RESP)" in terms
        assert "I" not in terms and "J" not in terms


class TestEnhancedPiedaAjf:
    """The ajf side: &LRD DISP='ON' and &ANALYSIS ES_RESP='YES'."""

    def _obj(self, **kw):
        from abmptools.abinit_io import abinit_io
        obj = abinit_io()
        obj.abinit_ver = kw.pop("ver", "v2rev8")
        for k, v in kw.items():
            setattr(obj, k, v)
        return obj

    def test_lrd_is_off_by_default(self):
        assert self._obj().disp is False
        assert self._obj().es_resp is False

    def test_es_resp_needs_a_version_that_has_it(self):
        """V1DD2024 has no ES_RESP; asking for it must not write it."""
        from abmptools.abinit_io import abinit_io
        assert abinit_io().es_resp is False

    @pytest.mark.parametrize("ver,expected", [("v2rev8", True),
                                              ("rev23", False)])
    def test_es_resp_version_gate(self, ver, expected):
        obj = self._obj(ver=ver, es_resp=True)
        assert (obj.es_resp and obj.abinit_ver in ("v2rev4", "v2rev8")) is expected
