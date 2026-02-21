# -*- coding: utf-8 -*-
"""
exporter.py
-----------
Orchestrator: GROParser → GROAdapter → UDFWriter → output UDF.

Usage::

    from abmptools.gro2udf import Exporter
    Exporter().export("test.udf", "output.gro")
    # writes: test_groout.udf  in the current directory

The output filename is always ``{udf_basename}_groout.udf`` placed in the
current working directory, matching the original gro2udf.py behaviour.

Step-filter logic (replicates importStructure / ConvertStructure):
  - read Output_Interval_Steps from the UDF's static section
  - erase all existing records
  - for each GRO frame: write it as a new record only when
      frame.step % output_interval == 0
  - stop after writing more than max_record records (= Total_Steps // interval)
"""
from __future__ import annotations

import os


class Exporter:
    """Thin coordinator: wires together GROParser, GROAdapter, UDFWriter."""

    def export(self, udf_path: str, gro_path: str) -> int:
        """
        Convert *gro_path* frames and write updated UDF to
        ``{udf_basename}_groout.udf`` in the current directory.

        Returns 0 on success, 1 on error.
        """
        try:
            return self._run(udf_path, gro_path)
        except Exception as exc:
            print("ERROR:", exc)
            return 1

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _run(self, udf_path: str, gro_path: str) -> int:
        from UDFManager import UDFManager
        from .gro_parser import GROParser
        from .gro_adapter import GROAdapter
        from .udf_writer import UDFWriter

        print("## gro2udf ")

        udf = UDFManager(udf_path)

        # Read simulation parameters from static data (before erasing records)
        dt            = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Time.delta_T", "[ps]"
        )
        total_steps   = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps"
        )
        output_interval = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps"
        )
        max_record = int(total_steps // output_interval)

        # Erase all existing records (same as importStructure)
        udf.eraseRecord(0, udf.totalRecord())

        parser  = GROParser()
        adapter = GROAdapter()
        writer  = UDFWriter()

        written = 0
        for frame in parser.parse_frames(gro_path):
            # Guard: stop if we already wrote more than max_record records
            # (replicates the "while j <= maxRecord" condition)
            if written > max_record:
                break

            # Determine step number (fall back to written count when dt==0)
            steps = frame.step if dt != 0 else written

            # Skip frames not on an output boundary
            if steps % output_interval != 0:
                continue

            udf.newRecord()
            print("steps = {}, record = {}".format(steps, udf.currentRecord()))

            positions, cell = adapter.to_positions_and_cell(frame)
            writer.write_frame(udf, positions, cell, steps, frame.time)

            written += 1

        print("Total number of records : ", udf.totalRecord())

        # Output file: {udf_basename}_groout.udf in current directory
        output_file = (
            os.path.basename(os.path.splitext(udf_path)[0]) + "_groout.udf"
        )
        print("output file : ", output_file)
        udf.write(output_file)

        print(" Finished!! ")
        return 0
