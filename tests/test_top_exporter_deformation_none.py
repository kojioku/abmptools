# -*- coding: utf-8 -*-
"""gro2udf must not leave ``Deformation.Method`` as the string ``"None"``.

The COGNAC schema offers ``"None"`` for this select and COGNAC itself writes it
(its bundled sample UDFs carry it).  The GROMACS exporters around it, however,
test the field with ``len(method) > 0`` and then match it against the list of
deformations they support::

    deform = _udf_.get("Simulation_Conditions.Dynamics_Conditions.Deformation.Method")
    list_supported_deform = ["Cell_Deformation"]
    if len(deform) > 0:
        if deform in list_supported_deform:
            ...
        else:
            raise GROMACS_ConvertError(
                "Error!! deformation type '{}' is not supported.".format(deform))

-- J-OCTA 11.1 ``python/Export_GROMACS.py``, lines 1722-1737

So the four-character string ``"None"`` reads as "a deformation I don't know"
and aborts the conversion.  J-OCTA writes ``""`` in its own output, and the
bundled cognac11.2 template already used ``""``; only the cognac10.1 template
carried ``"None"``, which is exactly the template the ``--cognac-version 101``
route (OCTA8.4 / J-OCTA) selects.

``abmptools.udf2gro`` accepts both spellings since 2.13.4, but J-OCTA's
converter is not ours to patch, so the UDF we write has to use ``""``.
"""
import pytest

from abmptools.gro2udf.top_exporter import TopExporter


DEFORM = "Simulation_Conditions.Dynamics_Conditions.Deformation.Method"

#: The guard as it appears in J-OCTA's Export_GROMACS.py.
SUPPORTED = ["Cell_Deformation"]


def joctas_guard(method):
    """Return True if J-OCTA's exporter would accept *method*."""
    if len(method) > 0 and method not in SUPPORTED:
        return False
    return True


class FakeUDF:
    """Minimal stand-in recording put()s, so no OCTA install is needed."""

    def __init__(self, initial):
        self.values = dict(initial)

    def jump(self, _rec):
        pass

    def get(self, loc, *args, **kwargs):
        return self.values.get(loc, "")

    def put(self, value, loc, *args, **kwargs):
        self.values[loc] = value


@pytest.mark.parametrize("template_value", ["None", ""])
def test_no_deformation_is_written_as_empty(template_value):
    """Both spellings of "no deformation" must come out as ""."""
    u = FakeUDF({DEFORM: template_value})
    TopExporter._normalize_deformation_method(u)
    assert u.get(DEFORM) == ""


@pytest.mark.parametrize("template_value", ["Cell_Deformation", "Lees_Edwards"])
def test_real_deformation_is_left_alone(template_value):
    """A template that deliberately sets a deformation keeps it."""
    u = FakeUDF({DEFORM: template_value})
    TopExporter._normalize_deformation_method(u)
    assert u.get(DEFORM) == template_value


def test_the_guard_this_protects_against():
    """Pin the behaviour that motivated the normalisation."""
    assert joctas_guard("None") is False           # what used to be written
    assert joctas_guard("") is True                # what is written now
    assert joctas_guard("Cell_Deformation") is True


def test_bundled_templates_do_not_ship_none():
    """Neither bundled template may reintroduce the value at its source."""
    import os
    import re
    from abmptools.gro2udf import cli

    for path in (cli._BUILTIN_TEMPLATE, cli._BUILTIN_TEMPLATE_COGNAC101):
        assert os.path.exists(path), path
        with open(path) as f:
            body = f.read()
        # The Deformation block is written as a nested record; "None" must not
        # appear as a bare quoted token anywhere in Dynamics_Conditions.
        assert not re.search(r'"None"', body), (
            "%s still carries a quoted None" % os.path.basename(path))
