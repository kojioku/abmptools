# -*- coding: utf-8 -*-
"""
gro2udf writes ``Statistics_Data.Pressure = Pressure + Pres. DC``.

GROMACS reports the virial ``Pressure`` WITHOUT the long-range dispersion
correction; the barostat-controlled (true) pressure is ``Pressure + Pres. DC``.
Both xvg legends fold to ``Statistics_Data.Pressure`` and are summed by
``_aggregate_statistics_per_frame`` (the same folding used for
Proper+Improper Dih. -> Torsion), matching J-OCTA's converter.

Reference check (DRO10-PVP10 NPT, ref_p = 1 bar, 100k steps):
  avg Pressure alone      = 44.1 MPa   (wrong; virial only)
  avg(Pressure + Pres.DC) = 0.22 MPa   (~ 0.1 MPa barostat target; matches J-OCTA)
"""
from types import SimpleNamespace

from abmptools.gro2udf.top_exporter import (
    _XVG_TO_UDF_STATS, _aggregate_statistics_per_frame,
)


def test_pres_dc_mapped_to_pressure():
    """Both `Pressure` and `Pres. DC` legends must target Statistics_Data.Pressure."""
    assert _XVG_TO_UDF_STATS["Pressure"][0] == "Pressure"
    assert _XVG_TO_UDF_STATS["Pres. DC"][0] == "Pressure"
    # same (class, leaf, unit) key so they get summed, not written separately
    assert _XVG_TO_UDF_STATS["Pressure"] == _XVG_TO_UDF_STATS["Pres. DC"]


def test_pressure_is_virial_plus_dc():
    """Per-frame Pressure = xvg Pressure + xvg Pres. DC (summed)."""
    frames = [SimpleNamespace(time=0.0), SimpleNamespace(time=1.0)]
    energy_times = [0.0, 1.0]
    # frame 0 mirrors the DRO reference (440.95 + (-438.80) = 2.15 bar = 0.215 MPa)
    energy_series = {
        "Pressure": [440.95, 100.0],
        "Pres. DC": [-438.80, -60.0],
    }
    per_frame = _aggregate_statistics_per_frame(frames, energy_times, energy_series)

    key = ("Pressure", "", "[bar]")
    inst0 = per_frame[0][key][0]
    inst1 = per_frame[1][key][0]
    assert abs(inst0 - (440.95 - 438.80)) < 1e-6      # 2.15 bar (0.215 MPa)
    assert abs(inst1 - (100.0 - 60.0)) < 1e-6         # 40 bar
