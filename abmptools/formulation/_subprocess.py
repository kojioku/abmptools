# -*- coding: utf-8 -*-
"""Thin subprocess helper for abmptools.formulation.

We do not depend on a heavier framework because subprocess invocation
patterns differ slightly across existing sub-packages
(cg.peptide._subprocess vs. genesis.mmgbsa._subprocess); each module
owns its own thin wrapper.
"""
from __future__ import annotations

import logging
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

logger = logging.getLogger(__name__)


@dataclass
class CommandError(RuntimeError):
    cmd: str
    returncode: int
    stderr: str = ""

    def __post_init__(self) -> None:  # pragma: no cover - trivial
        super().__init__(
            f"command failed (rc={self.returncode}): {self.cmd}\n{self.stderr}"
        )


def run_command(
    argv: Sequence[str],
    *,
    cwd: Optional[Union[str, Path]] = None,
    env: Optional[Dict[str, str]] = None,
    input_text: Optional[str] = None,
    check: bool = True,
    capture: bool = True,
    allow_returncodes: Sequence[int] = (0,),
) -> subprocess.CompletedProcess:
    """Run *argv* and return the CompletedProcess.

    Parameters
    ----------
    allow_returncodes
        Return codes that are not treated as errors. Useful e.g. for
        ``packmol`` which exits 173 on "soft" failure but still writes
        a usable mixture PDB.
    """
    cmd_str = " ".join(shlex.quote(str(a)) for a in argv)
    logger.info("running: %s", cmd_str)
    result = subprocess.run(
        list(argv),
        cwd=str(cwd) if cwd is not None else None,
        env=env,
        input=input_text,
        capture_output=capture,
        text=True,
    )
    if check and result.returncode not in allow_returncodes:
        raise CommandError(
            cmd=cmd_str,
            returncode=result.returncode,
            stderr=(result.stderr or "") + (result.stdout or ""),
        )
    return result


def ensure_dir(p: Union[str, Path]) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_text(p: Union[str, Path], content: str) -> Path:
    p = Path(p)
    ensure_dir(p.parent)
    p.write_text(content)
    return p
