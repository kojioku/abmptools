#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Standalone CLI for the abmptools amorphous structure builder.

Usage examples::

    # SDF input
    python build_amorphous.py \\
        --mol api.sdf polymer.sdf \\
        --n_mol 200 40 \\
        --box 6.0 --temperature 300

    # SMILES input
    python build_amorphous.py \\
        --smiles "CCCCC" "c1ccccc1" \\
        --name pentane benzene \\
        --n_mol 200 50 --density 0.8

    # JSON config
    python build_amorphous.py --config mixture.json
"""
from abmptools.amorphous.cli import main

if __name__ == "__main__":
    main()
