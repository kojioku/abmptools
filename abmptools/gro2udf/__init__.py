# -*- coding: utf-8 -*-
"""
abmptools.gro2udf
-----------------
Converts GROMACS .gro files back to COGNAC-UDF format.

Legacy mode (GRO → UDF frames, using existing UDF as template)::

    from abmptools.gro2udf import Exporter
    Exporter().export("test.udf", "output.gro")
    # writes: test_groout.udf  in current directory

New mode (TOP + GRO → full COGNAC UDF)::

    from abmptools.gro2udf import TopExporter
    TopExporter().export("system.top", "output.gro",
                         template_path="template.udf",
                         out_path="result.udf")

CLI usage (via __main__.py or standalone script)::

    # Legacy mode
    python -m abmptools.gro2udf test.udf output.gro

    # From-TOP mode
    python -m abmptools.gro2udf --from-top system.top output.gro \\
        [--template template.udf] [--out result.udf]

    python gro2udf.py test.udf output.gro   # legacy wrapper
"""
from .exporter import Exporter
from .top_exporter import TopExporter

__all__ = ["Exporter", "TopExporter"]
