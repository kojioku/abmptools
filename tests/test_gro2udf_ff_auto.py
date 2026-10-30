# -*- coding: utf-8 -*-
"""``--ff auto``: name the force field of a .top for ``FF=n``, or say why not.

Representative topologies, small enough to write here: an AMBER (ff14SB)
protein, a GAFF2 and a GAFF small molecule, OPLS, and the cases where no
single force field fits -- a protein and a GAFF2 ligand together, types no
table knows, GAFF-family types that do not tell GAFF from GAFF2, [ defaults ]
that contradict the types, CHARMM (no FF=n number). Water and ions never
count. No UDFManager needed: the decision is made on the parsed topology.
"""
import textwrap

import pytest

from abmptools.gro2udf import ff_detect
from abmptools.gro2udf.top_adapter import TopAdapter
from abmptools.gro2udf.top_parser import TopParser

AMBER_DEF = "[ defaults ]\n1 2 yes 0.5 0.8333\n\n"

# (name, mass, sigma nm, epsilon kJ/mol)
AMBER = [("N", 14.01, 0.325, 0.71128), ("H", 1.008, 0.106908, 0.065689),
         ("CX", 12.01, 0.339967, 0.45773), ("HP", 1.008, 0.195998, 0.065689),
         ("C", 12.01, 0.339967, 0.359824), ("O", 16.0, 0.295992, 0.87864)]
GAFF2 = [("c3", 12.01, 0.339771, 0.451035), ("ca", 12.01, 0.331521, 0.413379),
         ("ha", 1.008, 0.262548, 0.067362), ("oh", 16.0, 0.324287, 0.389112),
         ("ho", 1.008, 0.0, 0.0), ("hc", 1.008, 0.260018, 0.087027)]
GAFF = [("c3", 12.01, 0.339967, 0.457730), ("ca", 12.01, 0.339967, 0.359824),
        ("ha", 1.008, 0.259964, 0.06276), ("oh", 16.0, 0.306647, 0.880314),
        ("ho", 1.008, 0.0, 0.0), ("hc", 1.008, 0.264953, 0.065689)]
OPLS = [("opls_135", 12.011, 0.35, 0.276144), ("opls_140", 1.008, 0.25, 0.12552),
        ("opls_154", 15.999, 0.312, 0.71128), ("opls_155", 1.008, 0.0, 0.0),
        ("opls_157", 12.011, 0.35, 0.276144), ("opls_156", 1.008, 0.25, 0.12552)]
WATER = [("OW", 16.0, 0.315061, 0.636386), ("HW", 1.008, 0.0, 0.0)]


def _top(types, molecules, defaults=AMBER_DEF, includes=()):
    s = "".join('#include "%s"\n' % i for i in includes) + "\n" + defaults
    s += "[ atomtypes ]\n" + "".join(
        "%s %s 0.0 A %s %s\n" % (n, m, sg, ep) for n, m, sg, ep in types) + "\n"
    for name, atoms in molecules:
        s += "[ moleculetype ]\n%s 3\n\n[ atoms ]\n" % name
        s += "".join("%d %s 1 %s A%d %d 0.0 %s\n" % (k + 1, t, name[:3], k + 1, k + 1, m)
                     for k, (t, m) in enumerate(atoms)) + "\n"
        if len(atoms) > 1:
            s += "[ bonds ]\n" + "".join("%d %d 1 0.15 300000.0\n" % (k + 1, k + 2) for k in range(len(atoms) - 1)) + "\n"
    s += "[ system ]\nt\n\n[ molecules ]\n" + "".join("%s 1\n" % n for n, _ in molecules)
    return s


def _gro(n):
    lines = ["t", "%5d" % n] + ["%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (1, "X", "A", k + 1, 0.1 * k, 0.0, 0.0)
                                for k in range(n)]
    return "\n".join(lines + ["   3.0 3.0 3.0"]) + "\n"


def _decide(tmp_path, types, molecules, **kw):
    top = tmp_path / "t.top"
    top.write_text(_top(types, molecules, **kw))
    n = sum(len(a) for _, a in molecules)
    (tmp_path / "t.gro").write_text(_gro(n))
    model = TopAdapter().build(TopParser().parse(str(top)), str(tmp_path / "t.gro"))
    return ff_detect.detect(model, ff_detect.include_names(str(top)))


def _mol(name, types):
    return (name, [(t[0], t[1]) for t in types])


WATER_MOL = ("SOL", [("OW", 16.0), ("HW", 1.008), ("HW", 1.008)])


@pytest.mark.parametrize("types, mols, kw, want", [
    (AMBER + WATER, [_mol("Protein", AMBER), WATER_MOL],
     {"includes": ["amber14sb.ff/forcefield.itp"]}, "amber"),
    (AMBER + WATER, [_mol("Protein", AMBER), WATER_MOL], {}, "amber"),          # 型だけでも
    (GAFF2 + WATER, [_mol("LIG", GAFF2), WATER_MOL], {}, "gaff2"),
    (GAFF + WATER, [_mol("LIG", GAFF), WATER_MOL], {}, "gaff"),
    (OPLS, [_mol("MOL", OPLS)], {"defaults": "[ defaults ]\n1 3 yes 0.5 0.5\n\n"}, "oplsaa"),
])
def test_a_single_force_field_is_named(tmp_path, types, mols, kw, want):
    d = _decide(tmp_path, types, mols, **kw)
    assert d.name == want, d.reasons


@pytest.mark.parametrize("types, mols, kw, reason", [
    (AMBER + GAFF2, [_mol("Protein", AMBER), _mol("LIG", GAFF2)], {}, "more than one"),
    ([("MOL0_1", 12.01, 0.34, 0.45), ("MOL0_2", 1.008, 0.26, 0.07)],
     [("MOL", [("MOL0_1", 12.01), ("MOL0_2", 1.008)])], {}, "nothing in the topology"),
    ([("c3", 12.01, 0.3399, 0.4577), ("ho", 1.008, 0.0, 0.0)],
     [("MOL", [("c3", 12.01), ("ho", 1.008)])], {}, "tells GAFF from GAFF2"),
    (GAFF2, [_mol("LIG", GAFF2)], {"defaults": "[ defaults ]\n1 2 yes 1.0 1.0\n\n"}, "[ defaults ]"),
    (AMBER, [_mol("Protein", AMBER)], {"includes": ["charmm36-jul2022.ff/forcefield.itp"]}, "more than one"),
    ([("CT1", 12.011, 0.356359, 0.08368), ("HA", 1.008, 0.235197, 0.092048)],
     [("PRO", [("CT1", 12.011), ("HA", 1.008)])],
     {"includes": ["charmm36-jul2022.ff/forcefield.itp"], "defaults": "[ defaults ]\n1 2 yes 1.0 1.0\n\n"},
     "charmm"),
])
def test_no_single_force_field_writes_nothing_and_says_why(tmp_path, types, mols, kw, reason):
    d = _decide(tmp_path, types, mols, **kw)
    assert d.name is None
    assert any(reason in r for r in d.reasons), d.reasons


def test_water_and_ions_do_not_count(tmp_path):
    """Water's types (and an ion's) come with every force field."""
    ion = ("NA", [("Na", 22.99)])
    d = _decide(tmp_path, GAFF2 + WATER + [("Na", 22.99, 0.243928, 0.365846)],
                [_mol("LIG", GAFF2), WATER_MOL, ion])
    assert d.name == "gaff2" and "OW" not in str(d.families) and "Na" not in str(d.families)


def test_renamed_upper_case_types_are_not_amber(tmp_path):
    """ParmEd renames types (H1, C1 ...); a name alone does not make them AMBER."""
    d = _decide(tmp_path, [("C1", 12.01, 0.336, 0.36), ("H1", 1.008, 0.21, 0.07)],
                [("MOL", [("C1", 12.01), ("H1", 1.008)])])
    assert d.name is None and "amber" not in d.families


def test_the_cli_line_says_what_was_written():
    from abmptools.gro2udf.cli import _force_field_line
    yes = ff_detect.FFDecision("gaff2", {"gaff2": ["ca"]}, [])
    assert _force_field_line(yes, "auto") == "Force field: gaff2 (FF=3), decided from the topology"
    assert _force_field_line(yes, "gaff2").endswith("as given")
    no = ff_detect.FFDecision(None, {}, ["more than one force field: amber, gaff2"])
    line = _force_field_line(no, "auto")
    assert "not written" in line and "more than one" in line and "--ff" in line


def test_a_named_force_field_is_taken_as_given():
    from abmptools.gro2udf.top_exporter import TopExporter
    e = TopExporter()
    assert e._decide_force_field(None, "gaff2", ()).name == "gaff2"
    assert e._decide_force_field(None, "3", ()).name == "3"
    assert e._decide_force_field(None, "", ()).name is None
    with pytest.raises(ValueError):
        e._decide_force_field(None, "charmm36", ())


def test_cli_default_is_auto():
    from abmptools.gro2udf.cli import _from_top_parser
    args = _from_top_parser().parse_args(["a.top", "b.gro"])
    assert args.force_field == "auto"
