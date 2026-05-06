# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa._subprocess
------------------------------------
Re-exports private subprocess/file helpers from :mod:`abmptools.cg.peptide`.

Mirrors the pattern established by :mod:`abmptools.cg.membrane._subprocess`
and :mod:`abmptools.genesis.grest._subprocess` -- avoids duplicating the
private surface, while overriding ``setup_logging`` so the mmgbsa module
configures its own logger root.

When a fifth consumer appears, these shared helpers should move to a
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
    """Configure logging for abmptools.genesis.mmgbsa CLI invocations.

    Idempotent: safe to call multiple times.
    """
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger("abmptools.genesis.mmgbsa")
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
