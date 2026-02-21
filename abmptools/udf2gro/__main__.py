# -*- coding: utf-8 -*-
"""Allows: python -m abmptools.udf2gro <input.udf> <output_prefix>"""
import sys
from .cli import main

if __name__ == "__main__":
    try:
        main(sys.argv)
    except RuntimeError:
        print("ERROR: Export Gromacs failed.")
        print("may be parameter error.")
