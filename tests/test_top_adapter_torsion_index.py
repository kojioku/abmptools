# -*- coding: utf-8 -*-
"""Regression tests for the torsion type index mapping in gro2udf.

``TopAdapter._build_torsion_type_specs`` expands an Amber multi-term dihedral
into one ``TorsionTypeSpec`` per Fourier term, so the expanded list is longer
than ``TopRawData.torsion_types_from_mol``.  Molecule torsion records reference
the *unexpanded* index, and indexing the expanded list with it drifts by one
position per multi-term torsion already seen.  The drift is silent: names still
look well formed, so every later dihedral quietly carries another dihedral's
parameters, and the trailing references dangle.

These tests pin the mapping by rebuilding each molecule torsion's parameters
from the name it references and comparing them against the source .top.
"""
import pytest

from abmptools.gro2udf.top_parser import TopParser
from abmptools.gro2udf.top_adapter import TopAdapter


# A 6-atom chain.  The dihedral block deliberately interleaves single-term and
# multi-term entries so that any index drift shows up on the entries that
# follow a multi-term one.
MULTITERM_TOP = """\
[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
1   2   yes   0.5   0.8333

[ atomtypes ]
; name mass charge ptype sigma epsilon
c3   12.011  0.0  A  0.339967  0.457730
ca   12.011  0.0  A  0.339967  0.359824
ha    1.008  0.0  A  0.259964  0.062760

[ moleculetype ]
; Name nrexcl
MOL  3

[ atoms ]
; nr type resnr residu atom cgnr charge mass
1  c3  1  MOL  C1  1  -0.10  12.011
2  c3  1  MOL  C2  2  -0.10  12.011
3  ca  1  MOL  C3  3  -0.05  12.011
4  ca  1  MOL  C4  4  -0.05  12.011
5  c3  1  MOL  C5  5  -0.10  12.011
6  ha  1  MOL  H6  6   0.10   1.008

[ bonds ]
; ai aj funct b0 kb
1  2  1  0.1535  253634.0
2  3  1  0.1510  265265.0
3  4  1  0.1387  392459.0
4  5  1  0.1510  265265.0
5  6  1  0.1092  284512.0

[ dihedrals ]
; ai aj ak al funct phase k mult
1  2  3  4  9  0.0    0.6508  3
2  3  4  5  9  180.0  15.1670  2
2  3  4  5  9  0.0     1.2552  1
2  3  4  5  9  180.0   0.4184  3
3  4  5  6  9  0.0     0.8368  3
3  4  5  6  9  180.0   2.0920  2
2  4  3  5  4  180.0   4.6024  2
1  2  3  5  9  0.0     0.4211  3

[ molecules ]
MOL  2
"""

#: quad -> ordered list of (phase, k, mult) terms, transcribed from the block
#: above.  The dict is the ground truth the adapter output is checked against.
EXPECTED = {
    (1, 2, 3, 4): [(0.0, 0.6508, 3)],
    (2, 3, 4, 5): [(180.0, 15.1670, 2), (0.0, 1.2552, 1), (180.0, 0.4184, 3)],
    (3, 4, 5, 6): [(0.0, 0.8368, 3), (180.0, 2.0920, 2)],
    (2, 4, 3, 5): [(180.0, 4.6024, 2)],
    (1, 2, 3, 5): [(0.0, 0.4211, 3)],
}


@pytest.fixture
def built(tmp_path):
    """Parse MULTITERM_TOP and run the two adapter stages under test."""
    top = tmp_path / "multiterm.top"
    top.write_text(MULTITERM_TOP)
    raw = TopParser().parse(str(top))

    adapter = TopAdapter()
    tspecs = adapter._build_torsion_type_specs(raw.torsion_types_from_mol)
    mass_dict = adapter._build_mass_dict(raw)
    mol_specs = adapter._build_mol_specs(raw, mass_dict, [], [], tspecs)
    return raw, tspecs, mol_specs[0]


class TestTorsionTypeSpecExpansion:

    def test_multi_term_is_expanded_one_spec_per_term(self, built):
        raw, tspecs, _ = built
        # 5 unexpanded types -> 1 + 3 + 2 + 1 + 1 = 8 expanded specs
        assert len(raw.torsion_types_from_mol) == 5
        assert len(tspecs) == 8

    def test_every_spec_records_its_source_index(self, built):
        raw, tspecs, _ = built
        assert [ts.src_index for ts in tspecs] == [0, 1, 1, 1, 2, 2, 3, 4]

    def test_expanded_names_are_unique(self, built):
        _, tspecs, _ = built
        names = [ts.name for ts in tspecs]
        assert len(set(names)) == len(names)


class TestMoleculeTorsionNaming:

    def test_every_reference_resolves(self, built):
        """No molecule torsion may name a potential that was never defined."""
        _, tspecs, mol = built
        defined = {ts.name for ts in tspecs}
        dangling = [t.potential_name for t in mol.torsions
                    if t.potential_name not in defined]
        assert dangling == []

    def test_one_record_per_fourier_term(self, built):
        _, _, mol = built
        assert len(mol.torsions) == sum(len(v) for v in EXPECTED.values())

    def test_parameters_survive_the_round_trip(self, built):
        """Rebuild each quad's terms from the referenced specs.

        This is the assertion that fails when the unexpanded index is used to
        address the expanded list: names still resolve for the early entries,
        but the parameters attached to everything after the first multi-term
        dihedral belong to a different dihedral.
        """
        _, tspecs, mol = built
        by_name = {ts.name: ts for ts in tspecs}

        got = {}
        for t in mol.torsions:
            quad = (t.atom1, t.atom2, t.atom3, t.atom4)
            got.setdefault(quad, []).append(
                tuple(by_name[t.potential_name].params))

        assert set(got) == set(EXPECTED)
        for quad, terms in EXPECTED.items():
            assert len(got[quad]) == len(terms), quad
            for want, have in zip(terms, got[quad]):
                assert have[0] == pytest.approx(want[0]), quad
                assert have[1] == pytest.approx(want[1]), quad
                assert int(have[2]) == want[2], quad

    def test_improper_keeps_its_own_parameters(self, built):
        """funct 4 goes through the same branch and must not be shifted."""
        _, tspecs, mol = built
        by_name = {ts.name: ts for ts in tspecs}
        imp = [t for t in mol.torsions
               if (t.atom1, t.atom2, t.atom3, t.atom4) == (2, 4, 3, 5)]
        assert len(imp) == 1
        params = by_name[imp[0].potential_name].params
        assert params[0] == pytest.approx(180.0)
        assert params[1] == pytest.approx(4.6024)
