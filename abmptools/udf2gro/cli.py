# -*- coding: utf-8 -*-
"""
cli.py
------
Command-line entry point for udf2gro.

Maintains full backward compatibility with the original::

    python udf2gro.py <input.udf> <output_prefix>

Can also be called via the installed package::

    python -m abmptools.udf2gro <input.udf> <output_prefix>
"""
from __future__ import annotations
import sys
import os


def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) != 3:
        print("Usage: {} in_udf_name output_file_base".format(
            os.path.basename(argv[0])
        ))
        raise RuntimeError("")

    in_udf_name      = argv[1]
    output_file_base = argv[2]

    from .exporter import Exporter
    return Exporter().export(in_udf_name, output_file_base)


if __name__ == "__main__":
    try:
        main(sys.argv)
    except RuntimeError:
        print("ERROR: Export Gromacs failed.")
        print("may be parameter error.")
