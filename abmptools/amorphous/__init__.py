# -*- coding: utf-8 -*-
"""
abmptools.amorphous
--------------------
Amorphous multi-component structure builder for GROMACS MD.

Public API::

    from abmptools.amorphous import AmorphousBuilder, BuildConfig, ComponentSpec

    config = BuildConfig(
        components=[
            ComponentSpec(smiles="CCCCC", name="pentane", n_mol=200),
            ComponentSpec(smiles="c1ccccc1", name="benzene", n_mol=50),
        ],
        density_g_cm3=0.8,
        temperature=300,
        output_dir="./output",
    )
    builder = AmorphousBuilder(config)
    result = builder.build()
"""
from .models import BuildConfig, ComponentSpec
from .builder import AmorphousBuilder

__all__ = [
    "AmorphousBuilder",
    "BuildConfig",
    "ComponentSpec",
    # Lazily-exposed PubChem helpers (via module __getattr__):
    "fetch_3d_sdf",
    "fetch_smiles",
    "download_3d_sdf",
    "PubChemError",
    "PubChemNotFoundError",
    "PubChemNo3DError",
]

_PUBCHEM_NAMES = {
    "fetch_3d_sdf",
    "fetch_smiles",
    "download_3d_sdf",
    "PubChemError",
    "PubChemNotFoundError",
    "PubChemNo3DError",
}


def __getattr__(name):
    """Lazy re-export of pubchem helpers.

    Keeping these out of the eager import list avoids the
    ``RuntimeWarning: ... found in sys.modules`` when the module is
    invoked as ``python -m abmptools.amorphous.pubchem``.
    """
    if name in _PUBCHEM_NAMES:
        from . import pubchem as _pc
        return getattr(_pc, name)
    raise AttributeError(f"module 'abmptools.amorphous' has no attribute {name!r}")
