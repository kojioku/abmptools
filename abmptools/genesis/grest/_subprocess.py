# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest._subprocess
-----------------------------------
Re-exports private subprocess/file helpers from :mod:`abmptools.cg.peptide`.

The :mod:`abmptools.genesis.grest` module needs the same helpers
(``run_command``, ``ensure_dir``, ``write_text``, ``CommandError``) as
the cg/peptide and cg/membrane sub-packages. To avoid duplicating
private surface, we re-export. ``setup_logging`` is overridden so the
grest module's logger root is configured separately.

When a fourth consumer appears, these shared helpers should move to a
neutral location (``abmptools/_subprocess_common.py``); deferred until
then.
"""
from __future__ import annotations

import logging
import sys

from abmptools.cg.peptide._subprocess import (
    CommandError,
    ensure_dir,
    run_command,
    write_text,
)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for abmptools.genesis.grest CLI invocations.

    Idempotent: safe to call multiple times.
    """
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger("abmptools.genesis.grest")
    root.setLevel(level)
    if not root.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter(
                "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                datefmt="%H:%M:%S",
            )
        )
        root.addHandler(handler)


__all__ = [
    "CommandError",
    "ensure_dir",
    "run_command",
    "write_text",
    "setup_logging",
]
