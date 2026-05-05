# -*- coding: utf-8 -*-
"""
abmptools.genesis
-----------------
GENESIS-related sub-packages.

Currently the only occupant is :mod:`abmptools.genesis.grest`
(generalized Replica-Exchange with Solute Tempering -- gREST_SSCR
builder + analysis). Future siblings (REUS, FEP via GENESIS) may
land here without flattening the tree.

GENESIS itself (spdyn / atdyn / remd_convert) is licensed
LGPL-3.0-or-later and is **not** bundled. abmptools shells out to
the user-supplied binaries via :mod:`subprocess`, which is
"mere aggregation" under LGPL and compatible with the upcoming
abmptools Apache-2.0 relicense.

Build instructions: https://github.com/genesis-release-r-ccs/genesis
"""
