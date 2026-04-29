"""Tests for abmptools.getifiepieda module."""

import sys
import types
import pytest


# ---------------------------------------------------------------------------
# 1. Test that the module can be imported
# ---------------------------------------------------------------------------
class TestModuleImport:
    def test_import_module(self):
        """Verify abmptools.getifiepieda is importable."""
        import abmptools.getifiepieda as mod
        assert hasattr(mod, "get_args")
        assert hasattr(mod, "setupmode")


# ---------------------------------------------------------------------------
# 2. Tests for get_args()
# ---------------------------------------------------------------------------
class TestGetArgs:
    """Test get_args with mocked sys.argv."""

    def test_frag_mode_single(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys, "argv",
            ["getifiepieda", "--frag", "10", "-i", "test.log"],
        )
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.frag == ["10"]
        assert args.input == "test.log"

    def test_frag_mode_two_ranges(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys, "argv",
            ["getifiepieda", "--frag", "1-10", "101-200", "-i", "test.log"],
        )
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.frag == ["1-10", "101-200"]

    def test_mol_mode(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys, "argv",
            ["getifiepieda", "--mol", "1-10", "-i", "test.log"],
        )
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.mol == "1-10"

    def test_defaults(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys, "argv",
            ["getifiepieda", "--frag", "1", "-i", "test.log"],
        )
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.zp == 0
        assert args.pynp == 1
        assert args.exclude == []
        assert args.nof90so is True
        assert args.noresinfo is True
        assert args.dimene == [2, 1]
        assert args.momene == [1]
        assert args.is_momdim is False
        assert args.dist is None
        assert args.molname is None
        assert args.multi is None
        assert args.tfmatrix is None
        assert args.ffmatrix is None
        assert args.fraginmol is None
        assert args.dimeres is False

    def test_dist_option(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys, "argv",
            ["getifiepieda", "--frag", "10", "-d", "0.8", "-i", "test.log"],
        )
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.dist == 0.8

    def test_inputx_option(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys, "argv",
            ["getifiepieda", "--frag", "1", "-ix", "file-xxx-tail.log"],
        )
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.inputx == "file-xxx-tail.log"
        assert args.input is None

    def test_exclude_option(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys, "argv",
            ["getifiepieda", "--frag", "1", "-i", "t.log", "--exclude", "5", "10"],
        )
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.exclude == [5, 10]

    def test_fraginmol_option(self, monkeypatch, capsys):
        monkeypatch.setattr(
            sys, "argv",
            [
                "getifiepieda",
                "--fraginmol", "1-10", "2", "WAT", "1",
                "-i", "test.log",
            ],
        )
        from abmptools.getifiepieda import get_args
        args = get_args()
        assert args.fraginmol == ["1-10", "2", "WAT", "1"]


# ---------------------------------------------------------------------------
# 3. Tests for setupmode()
# ---------------------------------------------------------------------------
class TestSetupmode:
    """Test setupmode by injecting a mock args into the module namespace
    and passing a simple mock aobj."""

    @staticmethod
    def _make_args(**overrides):
        """Return a namespace that mirrors argparse defaults."""
        defaults = dict(
            input=None,
            inputx=None,
            frag=None,
            mol=None,
            molname=None,
            fraginmol=None,
            ffmatrix=None,
            multi=None,
            tfmatrix=None,
            time=None,
            dist=None,
            pynp=1,
            exclude=[],
            nof90so=True,
            noresinfo=True,
            dimene=[2, 1],
            momene=[1],
            is_momdim=False,
            dimeres=False,
            zp=0,
        )
        defaults.update(overrides)
        return types.SimpleNamespace(**defaults)

    @staticmethod
    def _make_aobj():
        """Return a minimal mock analysis object."""
        return types.SimpleNamespace()

    def _call_setupmode(self, args_overrides):
        """Build a fake args namespace and call setupmode(aobj, args).

        ``setupmode`` was refactored to take ``args`` as a positional
        parameter rather than reading a module-level global, so tests
        now pass ``fake_args`` directly. The legacy
        ``mod.args = fake_args`` patching still happens as a defence
        in case any helper inside ``setupmode`` falls back to the
        module-level reference (none currently do, but the assignment
        is cheap and keeps the test resilient if that changes).
        """
        import abmptools.getifiepieda as mod
        fake_args = self._make_args(**args_overrides)
        mod.args = fake_args
        aobj = self._make_aobj()
        result = mod.setupmode(aobj, fake_args)
        return aobj, result

    # -- frag mode, single target ----
    def test_frag_mode_single_target(self, capsys):
        aobj, (tgt1, tgt2) = self._call_setupmode(
            dict(input="test.log", frag=["10"])
        )
        assert aobj.anlmode == "frag"
        assert tgt1 == "10"
        assert tgt2 == []
        assert aobj.ilog_head == "test.log"

    # -- frag mode, two targets ----
    def test_frag_mode_two_targets(self, capsys):
        aobj, (tgt1, tgt2) = self._call_setupmode(
            dict(input="test.log", frag=["1-10", "101-200"])
        )
        assert aobj.anlmode == "frag"
        assert aobj.tgt2type == "frag"
        assert tgt1 == "1-10"
        assert tgt2 == "101-200"

    # -- mol mode ----
    def test_mol_mode(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", mol="1-10")
        )
        assert aobj.anlmode == "mol"
        assert aobj.tgtmolid == "1-10"

    # -- fraginmol mode ----
    def test_fraginmol_mode(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", fraginmol=["5", "2", "WAT", "1"])
        )
        assert aobj.anlmode == "fraginmol"
        assert aobj.tgtmolid == 5
        assert aobj.tgt1_lofrag == 2
        assert aobj.tgt2molname == "WAT"
        assert aobj.tgt2_lofrag == 1

    def test_fraginmol_mode_string_molid(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", fraginmol=["ABC", "2", "WAT", "1"])
        )
        assert aobj.tgtmolid == "ABC"

    # -- dist mode ----
    def test_dist_mode(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", frag=["10"], dist=0.8)
        )
        assert aobj.tgt2type == "dist"
        assert aobj.dist == 0.8

    # -- molname mode ----
    def test_molname_mode(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", frag=["10"], molname="WAT")
        )
        assert aobj.tgt2type == "molname"
        assert aobj.tgt2molname == "WAT"

    # -- molname + dist (molname takes precedence for tgt2type) ----
    def test_molname_with_dist(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", frag=["10"], molname="WAT", dist=0.8)
        )
        assert aobj.tgt2type == "molname"
        assert aobj.dist == 0.8

    # -- ffmatrix mode ----
    def test_ffmatrix_mode(self, capsys):
        aobj, (tgt1, tgt2) = self._call_setupmode(
            dict(input="test.log", ffmatrix=["1-100", "101-200"])
        )
        assert aobj.anlmode == "frag"
        assert aobj.tgt2type == "frag"
        assert aobj.matrixtype == "frags-frags"
        assert tgt1 == "1-100"
        assert tgt2 == "101-200"

    # -- dimeres mode ----
    def test_dimeres_flag(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", frag=["10"], dimeres=True)
        )
        assert aobj.tgt2type == "dimer-es"

    # -- is_momdim ----
    def test_is_momdim_true(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", frag=["10"], is_momdim=True)
        )
        assert aobj.is_momdimene is True

    def test_is_momdim_false(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(input="test.log", frag=["10"], is_momdim=False)
        )
        assert aobj.is_momdimene is False

    # -- common attributes always set ----
    def test_common_attributes(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(
                input="test.log",
                frag=["10"],
                exclude=[5, 10],
                nof90so=False,
                pynp=4,
                noresinfo=False,
                dimene=[3, 2],
                momene=[2],
            )
        )
        assert aobj.exceptfrag == [5, 10]
        assert aobj.f90soflag is False
        assert aobj.pynp == 4
        assert aobj.addresinfo is False
        assert aobj.dimfrag1 == 3
        assert aobj.dimfrag2 == 2
        assert aobj.momfrag == 2

    # -- inputx mode ----
    def test_inputx_splits_on_xxx(self, capsys):
        aobj, _ = self._call_setupmode(
            dict(inputx="head-part-xxx-tail.log", frag=["1"])
        )
        assert aobj.ilog_head == "head-part-"
        assert aobj.ilog_tail == "-tail.log"

    # -- multi mode ----
    def test_multi_mode_two_targets(self, capsys):
        aobj, (tgt1, tgt2) = self._call_setupmode(
            dict(
                input='["head", "tail"]',
                multi=["1-100", "101-200"],
                time=["100", "3100", "1000"],
            )
        )
        assert aobj.anlmode == "multi"
        assert aobj.start == 100
        assert aobj.end == 3100
        assert aobj.interval == 1000
        assert aobj.ilog_head == "head"
        assert aobj.ilog_tail == "tail"
        assert tgt1 == "1-100"
        assert tgt2 == "101-200"

    def test_multi_mode_single_target(self, capsys):
        aobj, (tgt1, tgt2) = self._call_setupmode(
            dict(
                input='["head", "tail"]',
                multi=["10"],
                time=["100", "3100", "1000"],
            )
        )
        assert aobj.anlmode == "multi"
        assert tgt1 == "10"
        assert tgt2 == []

    # -- tfmatrix mode ----
    def test_tfmatrix_mode(self, capsys):
        aobj, (tgt1, tgt2) = self._call_setupmode(
            dict(
                input='["head", "tail"]',
                tfmatrix=["1-100", "101-200"],
                time=["100", "3100", "1000"],
            )
        )
        assert aobj.anlmode == "multi"
        assert aobj.tgt2type == "frag"
        assert aobj.matrixtype == "times-frags"
        assert aobj.start == 100
        assert aobj.end == 3100
        assert aobj.interval == 1000
        assert tgt1 == "1-100"
        assert tgt2 == "101-200"
