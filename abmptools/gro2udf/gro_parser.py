# -*- coding: utf-8 -*-
"""
gro_parser.py
-------------
Parses GROMACS .gro files frame by frame.

Each call to parse_frames() is a generator that yields one GROFrame per
frame found in the file.

The column-reading logic faithfully replicates the original gro2udf.py so
that numeric precision is identical:
  - x, y, z are located by finding the decimal-point positions starting
    at column 20 (handles wider-than-standard fields correctly).
  - vx, vy, vz are read at fixed offsets relative to ipz.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterator, List, Optional, Tuple


@dataclass
class GROAtomLine:
    """Data for one atom in a GRO frame."""
    residue_col: str    # first 5 chars of line (used to detect mol boundary)
    x: float            # [nm]
    y: float            # [nm]
    z: float            # [nm]
    vx: float           # [nm/ps]  (0.0 if velocity columns absent)
    vy: float           # [nm/ps]
    vz: float           # [nm/ps]


@dataclass
class GROFrame:
    """One frame from a (possibly multi-frame) .gro file."""
    title: str                          # full header line (first line)
    step: int                           # from "step= N" field
    time: float                         # [ps] from "t= N" field
    atoms: List[GROAtomLine] = field(default_factory=list)
    box_vals: List[float] = field(default_factory=list)  # 3 or 9 values [nm]


class GROParser:
    """Reads a .gro file and yields GROFrame objects one by one."""

    def parse_frames(self, filepath: str) -> Iterator[GROFrame]:
        """Yield GROFrame for each frame found in *filepath*."""
        with open(filepath, "r") as f:
            while True:
                frame = self._read_one_frame(f)
                if frame is None:
                    break
                yield frame

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _read_one_frame(self, f) -> Optional[GROFrame]:
        # --- header line ---
        s = f.readline()
        if s == "":
            return None  # EOF
        title = s.rstrip("\n")
        step, time = self._parse_header(s)

        # --- atom count ---
        s = f.readline()
        if s == "":
            return None
        num_atoms = int(s.strip())

        # --- atom lines ---
        atoms: List[GROAtomLine] = []
        for _ in range(num_atoms):
            s = f.readline()
            atoms.append(self._parse_atom_line(s))

        # --- box line ---
        s = f.readline()
        box_vals = [float(v) for v in s.split()]

        return GROFrame(
            title=title,
            step=step,
            time=time,
            atoms=atoms,
            box_vals=box_vals,
        )

    @staticmethod
    def _parse_header(s: str) -> Tuple[int, float]:
        """Extract step and time from a GRO header line.

        Expected format (example)::

            test.udf t=   0.00000 step= 0

        Falls back to (0, 0.0) if the fields are absent.
        """
        temp = s.split()
        if "t=" in temp:
            npos = temp.index("t=")
            time = float(temp[npos + 1])
        else:
            time = 0.0
            print("warning : no record")

        if "step=" in temp:
            npos = temp.index("step=")
            step = int(temp[npos + 1])
        else:
            step = 0

        return step, time

    @staticmethod
    def _parse_atom_line(s: str) -> GROAtomLine:
        """Parse one atom line using the decimal-point search from original.

        The original gro2udf.py locates x/y/z by finding '.' characters
        starting at column 20.  This is more robust than fixed-width slicing
        for GRO files that have extra-wide number columns.
        """
        residue_col = s[0:5]

        # Locate decimal points for x, y, z
        ipx = s.find(".", 20)
        ipy = s.find(".", ipx + 1)
        ipz = s.find(".", ipy + 1)

        x = float(s[20: ipx + 4])
        y = float(s[ipx + 4: ipy + 4])
        z = float(s[ipy + 4: ipz + 4])

        # Velocity columns are optional (present when len > 45)
        if len(s) > 45:
            if ipz == 40:
                vx = float(s[44:52])
                vy = float(s[52:60])
                vz = float(s[60:68])
            else:
                inc = ipz - 40
                vx = float(s[44 + inc: 52 + inc])
                vy = float(s[52 + inc: 60 + inc])
                vz = float(s[60 + inc: 68 + inc])
        else:
            vx = vy = vz = 0.0

        return GROAtomLine(
            residue_col=residue_col,
            x=x, y=y, z=z,
            vx=vx, vy=vy, vz=vz,
        )
