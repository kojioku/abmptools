# -*- coding: utf-8 -*-
"""
Allows running as:  python -m abmptools.gro2udf <input.udf> <input.gro>
"""
import sys
from .cli import main

if __name__ == "__main__":
    try:
        main(sys.argv)
    except RuntimeError:
        print("ERROR: gro2udf failed.")
        sys.exit(1)
