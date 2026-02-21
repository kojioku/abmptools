# -*- coding: utf-8 -*-
"""
abmptools.udf2gro
-----------------
Converts COGNAC-UDF files to GROMACS format (gro/top/mdp/ndx).

Public API::

    from abmptools.udf2gro import Exporter
    Exporter().export("system.udf", "output")

CLI usage (via __main__.py)::

    python -m abmptools.udf2gro system.udf output
"""
from .exporter import Exporter
from ..core.system_model import SystemModel

__all__ = ["Exporter", "SystemModel"]
