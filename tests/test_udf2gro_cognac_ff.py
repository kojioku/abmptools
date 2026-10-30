# -*- coding: utf-8 -*-
"""udf2gro: force-field decisions that used to be made silently.

- the combination rule comes from the UDF's own pair table; a pair that
  follows neither rule goes to [ nonbond_params ], and a pair missing
  between two types that both have Lennard-Jones stops the conversion
  (COGNAC leaves it without interaction; a rule would add one);
- fudgeLJ comes from the pairs that are 1-4 in the molecules, not from the
  first pair of the list; other scales go to [ pairtypes ];
- an AMBER torsion whose atoms are not a bonded chain is an improper
  (function 4);
- the phase of an AMBER torsion is decided in one place;
- the writer emits [ nonbond_params ], [ pairtypes ], [ settles ],
  residue columns and one .itp per molecule type on request.

No UDFManager needed: the helpers are pure, the writer takes a SystemModel.
"""
import math

import pytest

from abmptools.core.system_model import (
    AtomRecord, AtomType, BondRecord, MoleculeTopology, PairRecord, SystemModel)
from abmptools.udf2gro import udf_adapter as U
from abmptools.udf2gro.gromacs.writers.top_writer import TopWriter


def _types(**lj):
    return [AtomType(n, 12.0, s, e) for n, (s, e) in lj.items()]


def _pairs(table):
    out = {}
    for (a, b), (s, e, scale) in table.items():
        out[(a, b)] = out[(b, a)] = (s, e, scale)
    return out


def _lb(a, b):
    return 0.5 * (a[0] + b[0]), math.sqrt(a[1] * b[1])


A, B, C = (0.34, 0.45), (0.30, 0.60), (0.25, 0.10)


def _full(rule_cross, scale=0.5):
    lj = {"a": A, "b": B, "c": C}
    t = {(n, n): (*v, scale) for n, v in lj.items()}
    for x, y in (("a", "b"), ("a", "c"), ("b", "c")):
        t[(x, y)] = (*rule_cross(lj[x], lj[y]), scale)
    return _types(**lj), _pairs(t)


class TestCombinationRule:

    def test_lorentz_berthelot_pairs_give_rule_2(self):
        types, pairs = _full(_lb)
        assert U.UdfAdapter._combination_rule(types, pairs, 2) == (2, [])

    def test_geometric_pairs_give_rule_3(self):
        types, pairs = _full(lambda a, b: (math.sqrt(a[0] * b[0]), math.sqrt(a[1] * b[1])))
        assert U.UdfAdapter._combination_rule(types, pairs, 2) == (3, [])

    def test_an_odd_pair_becomes_a_nonbond_param(self):
        types, pairs = _full(_lb)
        pairs[("a", "b")] = pairs[("b", "a")] = (0.40, 0.11, 0.5)
        rule, rows = U.UdfAdapter._combination_rule(types, pairs, 2)
        assert rule == 2 and rows == [("a", "b", 0.40, 0.11)]

    def test_a_missing_pair_between_two_lj_types_stops(self):
        types, pairs = _full(_lb)
        del pairs[("a", "c")], pairs[("c", "a")]
        with pytest.raises(RuntimeError, match="a-c"):
            U.UdfAdapter._combination_rule(types, pairs, 2)

    def test_a_type_without_lj_needs_no_pair(self):
        """Water's H and GAFF's ho have no Lennard-Jones and no pair at all."""
        types, pairs = _full(_lb)
        types.append(AtomType("HW", 1.008, 0.0, 0.0))
        assert U.UdfAdapter._combination_rule(types, pairs, 2) == (2, [])

    def test_no_cross_pair_keeps_the_name_based_rule(self):
        types = _types(a=A, b=B)
        pairs = _pairs({("a", "a"): (*A, 0.5), ("b", "b"): (*B, 0.5)})
        assert U.UdfAdapter._combination_rule(types, pairs, 3) == (3, [])


def _topo(pairs_of_types):
    """One molecule whose 1-4 pairs join the given types."""
    atoms, pairs = [], []
    for k, (t1, t2) in enumerate(pairs_of_types):
        atoms += [AtomRecord(2 * k + 1, t1, t1, 0.0), AtomRecord(2 * k + 2, t2, t2, 0.0)]
        pairs.append(PairRecord(2 * k + 1, 2 * k + 2))
    return MoleculeTopology("M", "M", 3, atoms=atoms, pairs=pairs)


class TestFudge:

    def test_scale_comes_from_the_1_4_pairs_not_the_first_pair(self):
        """Water and ions carry 1.0, the solute 0.5; the list may start with either."""
        pairs = _pairs({("OW", "OW"): (0.31, 0.6, 1.0), ("c3", "c3"): (0.34, 0.45, 0.5),
                        ("c3", "hc"): (0.30, 0.2, 0.5)})
        fudge, rows = U.UdfAdapter._fudge_lj([_topo([("c3", "c3"), ("c3", "hc")])],
                                              [], pairs, 1.0, 2)
        assert fudge == 0.5 and rows == []

    def test_other_scales_become_pairtypes(self):
        pairs = _pairs({("c3", "c3"): (0.34, 0.45, 0.5), ("c3", "hc"): (0.30, 0.2, 0.5),
                        ("hc", "hc"): (0.26, 0.06, 1.0)})
        fudge, rows = U.UdfAdapter._fudge_lj(
            [_topo([("c3", "c3"), ("c3", "hc"), ("hc", "hc")])], [], pairs, 1.0, 2)
        assert fudge == 0.5
        assert rows == [("hc", "hc", 0.26, pytest.approx(0.06))]   # epsilon x 1.0

    def test_no_1_4_pair_keeps_the_fallback(self):
        assert U.UdfAdapter._fudge_lj([_topo([])], [], {}, 0.7, 2) == (0.7, [])


def test_proper_and_improper_from_the_bonds():
    bonded = {(1, 2), (2, 1), (2, 3), (3, 2), (3, 4), (4, 3), (3, 5), (5, 3), (3, 6), (6, 3)}
    assert U._is_proper(1, 2, 3, 4, bonded)
    assert not U._is_proper(4, 5, 3, 6, bonded)      # 3 is bonded to the other three
    assert not U._is_proper(2, 4, 3, 5, bonded)


def test_amber_phase_is_passed_through_whatever_trans_is_0_says():
    """The one place that decides; see its docstring for why."""
    for flag in (0, 1):
        assert U.amber_phase_to_gromacs(0.0, flag) == 0.0
        assert U.amber_phase_to_gromacs(180.0, flag) == 180.0


def _model(**kw):
    wat = MoleculeTopology("WAT", "WAT", 2, atoms=[
        AtomRecord(1, "OW", "OW", -0.8, resnr=5, resname="SOL"),
        AtomRecord(2, "HW", "HW1", 0.4), AtomRecord(3, "HW", "HW2", 0.4)],
        settles=(0.09572, 0.15136))
    return SystemModel(title="t", udf_path="t", comb_rule=2, flags14=1, fudgeLJ=0.5,
                       fudgeQQ=0.8333, calcQQ=1,
                       atom_types=[AtomType("OW", 16.0, 0.315, 0.636), AtomType("HW", 1.008, 0, 0)],
                       mol_topologies=[wat], mol_sequence=[("WAT", 2)], **kw)


def test_writer_emits_the_new_sections(tmp_path):
    m = _model(nonbond_params=[("OW", "HW", 0.2, 0.01)], pair_types=[("OW", "OW", 0.315, 0.318)])
    TopWriter().write(m, str(tmp_path / "t.top"))
    text = (tmp_path / "t.top").read_text()
    assert text.index("[ atomtypes ]") < text.index("[ nonbond_params ]") < text.index("[ pairtypes ]")
    assert "OW       HW       1" in text
    assert "[ settles ]" in text and "[ exclusions ]\n1  2  3\n2  1  3\n3  1  2" in text
    row = [l for l in text.splitlines() if " OW " in l and "SOL" in l]
    assert row and row[0].split()[2:4] == ["5", "SOL"]          # resnr / residue


def test_writer_splits_one_itp_per_molecule_type(tmp_path):
    TopWriter().write(_model(), str(tmp_path / "t.top"), split_itp=True)
    top = (tmp_path / "t.top").read_text()
    assert '#include "WAT.itp"' in top and "[ moleculetype ]" not in top
    assert top.index('#include "WAT.itp"') < top.index("[ molecules ]")
    assert "[ settles ]" in (tmp_path / "WAT.itp").read_text()
