# -*- coding: utf-8 -*-
"""
abmptools.crystal.legacy
-------------------------
Thin re-exports of the existing flat-layout crystal modules.

This namespace exists so Phase A can publish the crystal subpackage
without touching the legacy implementation. Phase C will refactor the
top-level files (``abmptools.readcif``, ``abmptools.pdb2fmo``, ...) and
this namespace will continue to expose the same callables for backward
compatibility.

Re-exported modules (no behaviour change):

- :mod:`abmptools.readcif`        — CIF → supercell PDB, hard-coded
                                     symmetry/layer dispatch
- :mod:`abmptools.pdb2fmo`        — PDB → for_abmp/*.{ajf,pdb}
- :mod:`abmptools.ajf2config`     — AJF template → segment_data.dat
- :mod:`abmptools.pdbmodify`      — PDB residue/chain renumber + sort
- :mod:`abmptools.getifiepieda`   — FMO log → IFIE/PIEDA CSV

The legacy CLI invocations remain unchanged:

    python -m abmptools.readcif    -i XXXI.cif --atomnum 32 -l 5
    python -m abmptools.pdb2fmo    -i XXXI*.pdb -p input_param
    python -m abmptools.ajf2config -i UNK.ajf
    python -m abmptools.pdbmodify  -i in.pdb -o out.pdb ...
    python -m abmptools.getifiepieda --multi 1 -dimeres ...
"""
from abmptools import (
    ajf2config,
    getifiepieda,
    pdb2fmo,
    pdbmodify,
    readcif,
)

__all__ = [
    "ajf2config",
    "getifiepieda",
    "pdb2fmo",
    "pdbmodify",
    "readcif",
]
