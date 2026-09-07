# -*- coding: utf-8 -*-
"""
gro2udf の guard: 変換できない .top を黙って通さないこと。

きっかけは Martini 2 の .top で、非結合が全部 [ nonbond_params ] にあり
[ atomtypes ] は c6 = c12 = 0 なので、**Pair_Interaction が 0 本の UDF が
無警告で出てきた**。
"""
import logging

import pytest

from abmptools.gro2udf.guard import (
    UnsupportedTopFeatureError,
    check_top,
    raise_if_unsupported,
    scan_sections,
)
from abmptools.gro2udf.top_parser import TopParser


MARTINI_TOP = """\
[ defaults ]
1  1  no  1.0  0.0

[ atomtypes ]
  EO   45.0 0.000 A 0.0 0.0
  P4   72.0 0.000 A 0.0 0.0

[ nonbond_params ]
  EO  EO  1 6.44779031E-02 4.07588234E-04
  EO  P4  1 1.50909015E-01 1.62668076E-03
  P4  P4  1 2.15584307E-01 2.32382966E-03

[ moleculetype ]
  TEST  1

[ atoms ]
  1 EO 1 TEST B1 1 0.0 45.0
  2 EO 1 TEST B2 2 0.0 45.0
  3 EO 1 TEST B3 3 0.0 45.0

[ bonds ]
  1 2 1 0.33 7000
  2 3 1 0.33 7000

[ angles ]
  1 2 3 2 130.0 50.0

[ system ]
  test

[ molecules ]
  TEST 1
"""

PLAIN_TOP = """\
[ defaults ]
1  2  yes  0.5  0.83333333

[ atomtypes ]
  c3  12.01 0.0 A 3.39967e-01 4.57730e-01
  hc   1.008 0.0 A 2.64953e-01 6.56888e-02

[ moleculetype ]
  MOL  3

[ atoms ]
  1 c3 1 MOL C1 1 -0.09 12.01
  2 hc 1 MOL H1 2  0.03  1.008

[ bonds ]
  1 2 1 0.1092 284512.0

[ pairs ]
  1 2 1

[ system ]
  plain

[ molecules ]
  MOL 1
"""


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


class TestScanSections:
    def test_counts_only_data_lines(self):
        counts = scan_sections(MARTINI_TOP.splitlines(keepends=True))
        assert counts["nonbond_params"] == 3
        assert counts["atomtypes"] == 2
        assert counts["angles"] == 1

    def test_comments_and_blanks_are_not_counted(self):
        counts = scan_sections([
            "[ nonbond_params ]\n", "; a comment\n", "\n", "   \n",
            "#include \"x.itp\"\n",
        ])
        assert counts["nonbond_params"] == 0

    def test_section_header_spacing_is_tolerated(self):
        counts = scan_sections(["[nonbond_params]\n", "EO EO 1 0.1 0.2\n"])
        assert counts["nonbond_params"] == 1


class TestMartiniIsRejected:
    """Martini は fatal 3 件 (comb-rule 1 / nonbond_params / angle funct 2)。"""

    @pytest.fixture
    def raw(self, tmp_path):
        return TopParser().parse(_write(tmp_path, "m.top", MARTINI_TOP))

    def test_parser_records_sections(self, raw):
        assert raw.sections["nonbond_params"] == 3

    def test_all_three_are_fatal(self, raw):
        fatal, warnings = check_top(raw, raw.sections)
        assert len(fatal) == 3
        joined = "\n".join(fatal)
        assert "comb-rule 1" in joined
        assert "[ nonbond_params ]" in joined
        assert "angle funct 2" in joined

    def test_raises(self, raw):
        with pytest.raises(UnsupportedTopFeatureError) as exc:
            raise_if_unsupported(raw, raw.sections)
        # 何が起きるはずだったかを本文に含めること
        assert "Pair_Interaction" in str(exc.value)
        assert "--allow-unsupported" in str(exc.value)

    def test_allow_unsupported_converts_but_warns(self, raw, caplog):
        with caplog.at_level(logging.WARNING):
            fatal, _ = raise_if_unsupported(raw, raw.sections,
                                            allow_unsupported=True)
        assert len(fatal) == 3
        assert "WRONG" in caplog.text


class TestPlainTopStillPasses:
    """通常の AA top を止めないこと。ここが退行すると既存運用が壊れる。"""

    @pytest.fixture
    def raw(self, tmp_path):
        return TopParser().parse(_write(tmp_path, "p.top", PLAIN_TOP))

    def test_no_findings_at_all(self, raw):
        fatal, warnings = check_top(raw, raw.sections)
        assert fatal == []
        assert warnings == []

    def test_pairs_is_not_reported(self, raw):
        """[ pairs ] は COGNAC が Scale_1_4_Pair で再現するので損失ではない。"""
        assert raw.sections["pairs"] == 1
        _, warnings = check_top(raw, raw.sections)
        assert not any("pairs" in w for w in warnings)


class TestSeveritySplit:
    """「間違って書かれる」= fatal、「落ちるだけ」= warning。"""

    def _raw(self, tmp_path, text, name):
        return TopParser().parse(_write(tmp_path, name, text))

    def test_constraints_are_fatal(self, tmp_path):
        text = PLAIN_TOP.replace(
            "[ pairs ]\n  1 2 1\n", "[ constraints ]\n  1 2 1 0.47\n")
        raw = self._raw(tmp_path, text, "c.top")
        fatal, _ = check_top(raw, raw.sections)
        assert any("[ constraints ]" in f for f in fatal)

    def test_settles_only_warns(self, tmp_path):
        text = PLAIN_TOP.replace(
            "[ pairs ]\n  1 2 1\n", "[ settles ]\n  1 1 0.09572 0.15139\n")
        raw = self._raw(tmp_path, text, "s.top")
        fatal, warnings = check_top(raw, raw.sections)
        assert fatal == []
        assert any("[ settles ]" in w for w in warnings)

    def test_unsupported_dihedral_funct_only_warns(self, tmp_path):
        """funct 2 の improper は書かれずに捨てられる — 損失だが取り違えではない。"""
        text = PLAIN_TOP.replace(
            "[ pairs ]\n  1 2 1\n",
            "[ dihedrals ]\n  1 2 1 2 2 180.0 4.6\n")
        raw = self._raw(tmp_path, text, "d.top")
        fatal, warnings = check_top(raw, raw.sections, raw.dihedral_functs)
        assert fatal == []
        assert any("dihedral funct 2" in w for w in warnings)
        assert raw.dihedral_functs == [2]

    def test_empty_section_is_ignored(self, tmp_path):
        text = PLAIN_TOP.replace("[ pairs ]\n  1 2 1\n",
                                 "[ nonbond_params ]\n; nothing here\n")
        raw = self._raw(tmp_path, text, "e.top")
        fatal, warnings = check_top(raw, raw.sections)
        assert fatal == []
        assert warnings == []


class TestIncludedItpIsSeen:
    """力場は .itp 側にある。include を展開した後で数えていること。"""

    def test_nonbond_params_inside_an_itp_is_caught(self, tmp_path):
        (tmp_path / "ff.itp").write_text(
            "[ nonbond_params ]\n  EO EO 1 0.06 0.0004\n")
        text = MARTINI_TOP.replace(
            "[ nonbond_params ]\n"
            "  EO  EO  1 6.44779031E-02 4.07588234E-04\n"
            "  EO  P4  1 1.50909015E-01 1.62668076E-03\n"
            "  P4  P4  1 2.15584307E-01 2.32382966E-03\n",
            '#include "ff.itp"\n')
        assert "[ nonbond_params ]" not in text
        raw = TopParser().parse(_write(tmp_path, "i.top", text))
        assert raw.sections["nonbond_params"] == 1
