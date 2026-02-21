# -*- coding: utf-8 -*-
"""
cli.py
------
Command-line entry point for gro2udf.

Maintains full backward compatibility with the original::

    python gro2udf.py <input.udf> <input.gro>
    python gro2udf.py <input.udf> <input.gro> <input.xvg>   (xvg ignored)

Can also be called via the installed package::

    python -m abmptools.gro2udf <input.udf> <input.gro>
"""
from __future__ import annotations

import os
import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) < 3:
        print("Usage: {} udffile grofile [xvg]".format(
            os.path.basename(argv[0])
        ))
        raise RuntimeError("Illegal arguments.")

    udf_path = argv[1]
    gro_path = argv[2]
    # argv[3] (xvg) is accepted for CLI compatibility but not processed here

    from .exporter import Exporter
    return Exporter().export(udf_path, gro_path)


if __name__ == "__main__":
    try:
        main(sys.argv)
    except RuntimeError:
        print("ERROR: gro2udf failed.")
        sys.exit(1)
