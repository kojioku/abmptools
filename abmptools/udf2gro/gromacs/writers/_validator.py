# -*- coding: utf-8 -*-
"""
_validator.py
--------------
Shared pre-flight checks for the GROMACS writers.

The one responsibility today is refusing to serialise a SystemModel whose
``ensemble_family`` is ``"cognac_only"`` (e.g. ``NPT_Andersen_Kremer_Grest``)
into GROMACS text formats.  Such a system can only be produced as a COGNAC
UDF; forcing a .gro/.top/.mdp would silently lose or mangle the ensemble
definition.

Used by :class:`GroWriter`, :class:`TopWriter`, :class:`MdpWriter`, and
:class:`ItpWriter` in the same package.
"""
from __future__ import annotations

from ....core.system_model import SystemModel


def raise_if_cognac_only(model: SystemModel, kind: str) -> None:
    """Raise ``ValueError`` if ``model`` is tagged as COGNAC-only.

    Parameters
    ----------
    model : SystemModel
        The intermediate representation about to be serialised.
    kind : str
        Short label of the target format, embedded into the error message
        (``"gro"``, ``"top"``, ``"mdp"``, ``"itp"``).
    """
    if getattr(model, "ensemble_family", "gromacs_ok") != "cognac_only":
        return

    algo = ""
    if model.sim_params is not None:
        algo = getattr(model.sim_params, "algorithm", "")

    raise ValueError(
        f"Cannot write a GROMACS .{kind} file: the SystemModel is marked "
        f"ensemble_family='cognac_only'"
        + (f" (algorithm={algo!r})" if algo else "")
        + ". Use `abmptools.gro2udf` with a COGNAC UDF template instead."
    )
