# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide._subprocess
---------------------------------
Internal helpers: subprocess wrapper, file/dir helpers, logging setup.

Private to abmptools.cg.peptide. The leading ``_`` is intentional --
this module is not part of the public API. Future cg modules
(polymer, smallmol) may share these via ``abmptools/cg/_common/``;
keep the private surface minimal until that refactor.
"""
from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence

logger = logging.getLogger(__name__)


class CommandError(Exception):
    """Raised when an external command fails."""

    def __init__(self, cmd: str, returncode: int, stderr: str):
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(
            f"Command failed (exit {returncode}): {cmd}\nstderr: {stderr}"
        )


def run_command(
    cmd: Sequence[str],
    cwd: Optional[Path] = None,
    check: bool = True,
    capture: bool = True,
    stdin_text: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """Run an external command with logging.

    Parameters
    ----------
    cmd
        Command arguments.
    cwd
        Working directory for the child process.
    check
        Raise CommandError on non-zero exit (default True).
    capture
        Capture stdout/stderr as text (default True).
    stdin_text
        Optional text fed via stdin (used by ``gmx genion``/``make_ndx``).
    """
    cmd_list = list(cmd)
    cmd_str = " ".join(cmd_list)
    logger.info("Running: %s", cmd_str)

    result = subprocess.run(
        cmd_list,
        cwd=cwd,
        capture_output=capture,
        text=True,
        input=stdin_text,
    )

    if check and result.returncode != 0:
        logger.error("Command failed: %s", cmd_str)
        if result.stderr:
            logger.error("stderr: %s", result.stderr)
        raise CommandError(cmd_str, result.returncode, result.stderr or "")

    return result


def ensure_dir(path: Path) -> Path:
    """Create directory (and parents) if missing."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_text(path: Path, content: str) -> None:
    """Write text content to file, creating parent dirs as needed."""
    path = Path(path)
    ensure_dir(path.parent)
    path.write_text(content)
    logger.debug("Wrote: %s", path)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for abmptools.cg.peptide CLI invocations.

    Idempotent: safe to call multiple times.
    """
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger("abmptools.cg.peptide")
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
