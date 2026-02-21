# -*- coding: utf-8 -*-
"""
abmptools.gro2udf
-----------------
Converts GROMACS .gro files back to COGNAC-UDF format.

Public API::

    from abmptools.gro2udf import Exporter
    Exporter().export("test.udf", "output.gro")
    # writes: test_groout.udf  in current directory

CLI usage (via __main__.py or standalone script)::

    python -m abmptools.gro2udf test.udf output.gro
    python gro2udf.py test.udf output.gro   # legacy wrapper
"""
from .exporter import Exporter

__all__ = ["Exporter"]
