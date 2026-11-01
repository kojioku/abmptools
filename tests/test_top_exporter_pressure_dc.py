# -*- coding: utf-8 -*-
"""gro2udf writes GROMACS' ``Pressure`` as it is: it already includes ``Pres. DC``.

With ``DispCorr = EnerPres`` the .edr lists ``Pres. DC`` next to ``Pressure``,
but ``Pressure`` already contains it. gro2udf used to add the two, counting
the correction twice. Measured at step 0 from the same coordinates and
velocities (SPC water, GROMACS 2021.3 and 2026.3):

  DispCorr = no        Pressure = 84180.04 bar
  DispCorr = EnerPres  Pressure = 83912.38 bar, Pres. DC = -267.67 bar
  difference                    = -267.66 bar  (= Pres. DC)

Both cases are below, with those numbers.
"""
from types import SimpleNamespace

from abmptools.gro2udf.top_exporter import (
    _XVG_TO_UDF_STATS, _aggregate_statistics_per_frame,
)

KEY = ("Pressure", "", "[bar]")


def _pressure(series):
    frames = [SimpleNamespace(time=0.0)]
    per_frame = _aggregate_statistics_per_frame(frames, [0.0], series)
    return per_frame[0][KEY][0]


def test_pres_dc_is_not_mapped():
    assert _XVG_TO_UDF_STATS["Pressure"][0] == "Pressure"
    assert "Pres. DC" not in _XVG_TO_UDF_STATS


def test_without_dispcorr_pressure_is_written_as_is():
    """DispCorr = no: the .edr has no Pres. DC column."""
    assert _pressure({"Pressure": [84180.04]}) == 84180.04


def test_with_dispcorr_pressure_is_written_as_is_not_plus_dc():
    """DispCorr = EnerPres: Pres. DC is listed, and already inside Pressure."""
    p = _pressure({"Pressure": [83912.38], "Pres. DC": [-267.67]})
    assert p == 83912.38
    assert abs((p - 84180.04) - (-267.67)) < 0.05      # the correction is in it once
