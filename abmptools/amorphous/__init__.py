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

__all__ = ["AmorphousBuilder", "BuildConfig", "ComponentSpec"]
