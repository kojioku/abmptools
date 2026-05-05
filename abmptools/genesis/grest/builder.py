# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest.builder
-------------------------------
End-to-end GENESIS gREST_SSCR system builder.

5-stage pipeline (:meth:`GrestBuilder.build`):

    1. ``_stage1_amber_parameterize``  -- ``tleap`` -> prmtop + coor
                                            + reference PDB
    2. ``_stage2_resolve_rest``        -- ``RESTSelectionSpec`` ->
                                            list of resnos +
                                            ``rno:`` selection string
    3. ``_stage3_temperature_ladder``  -- generate / passthrough ladder
    4. ``_stage4_inp_files``           -- render 4 ``.inp`` files
    5. ``_stage5_run_script``          -- emit ``run.sh`` + HPC scaffolds

Output layout::

    output_dir/
      input/config.json                          # config dump
      build/
        system.tleap                             # auto-rendered tleap script
        system.prmtop                            # AMBER topology
        system.coor                              # AMBER inpcrd
        system_ref.pdb                           # reference PDB
        system_tleap.log
        rest_residues.txt                        # mode=around only
      inp/
        step1_minimize.inp
        step2_equilibrate.inp
        step3_grest.inp
        step5_remd_convert.inp
      run.sh                                     # local mpirun
      run_pjsub.sh                               # Fugaku scaffold (not tested)
      run_sbatch.sh                              # SLURM scaffold (not tested)
      logs/                                      # populated by run.sh
"""
from __future__ import annotations

import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List

from . import (
    grest_runner,
    inp_writer,
    rest_selection,
    system_builder,
)
from ._subprocess import ensure_dir, write_text
from .inp_writer import SystemMetadata
from .models import GrestBuildConfig
from .replica_temperatures import generate_ladder

logger = logging.getLogger(__name__)


class GrestBuilder:
    """Build a GENESIS gREST_SSCR system from a :class:`GrestBuildConfig`."""

    def __init__(self, config: GrestBuildConfig):
        self.config = config
        self.output_dir = ensure_dir(Path(config.output_dir).resolve())
        self.input_dir = ensure_dir(self.output_dir / "input")
        self.build_dir = ensure_dir(self.output_dir / "build")
        self.inp_dir = ensure_dir(self.output_dir / "inp")
        ensure_dir(self.output_dir / "logs")
        # Filled in across stages.
        self._tleap_result: system_builder.TleapBuildResult | None = None
        self._rest_result: rest_selection.RESTSelectionResult | None = None
        self._ladder: List[float] = []
        self._meta: SystemMetadata | None = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(self) -> Dict[str, Any]:
        """Run the full 5-stage pipeline. Returns a result dict.

        The pipeline produces all GENESIS inputs but does not run
        ``spdyn``. The user invokes ``bash run.sh`` separately.
        """
        logger.info(
            "=== GrestBuilder.build start (project=%s) ===",
            self.config.project_name,
        )
        self._dump_config()
        self._stage1_amber_parameterize()
        self._stage2_resolve_rest()
        self._stage3_temperature_ladder()
        self._stage4_inp_files()
        self._stage5_run_script()
        result = self._collect_result()
        logger.info("=== GrestBuilder.build complete ===")
        return result

    # ------------------------------------------------------------------
    # Stage helpers
    # ------------------------------------------------------------------

    def _dump_config(self) -> None:
        """Persist the resolved config for reproducibility."""
        cfg_path = self.input_dir / "config.json"
        cfg_path.write_text(
            json.dumps(asdict(self.config), indent=2, ensure_ascii=False)
        )

    def _stage1_amber_parameterize(self) -> None:
        logger.info("=== Stage 1/5: tleap parameterisation ===")
        self._tleap_result = system_builder.build_amber_system(
            cfg=self.config,
            workdir=self.build_dir,
        )
        logger.info(
            "  prmtop=%s, box=%.2f x %.2f x %.2f A, n_protein_residues=%d",
            self._tleap_result.prmtop.name,
            *self._tleap_result.box_size_A,
            self._tleap_result.n_protein_residues,
        )

    def _stage2_resolve_rest(self) -> None:
        logger.info("=== Stage 2/5: REST residue resolution ===")
        assert self._tleap_result is not None
        self._rest_result = rest_selection.resolve_rest_selection(
            spec=self.config.rest_selection,
            prmtop=self._tleap_result.prmtop,
            workdir=self.build_dir,
            cpptraj_path=self.config.cpptraj_path,
        )
        logger.info(
            "  REST: %s (%d residues)",
            self._rest_result.selection_string,
            self._rest_result.n_residues,
        )

    def _stage3_temperature_ladder(self) -> None:
        logger.info("=== Stage 3/5: replica temperature ladder ===")
        self._ladder = generate_ladder(self.config.replica_temperatures)
        logger.info(
            "  ladder (%d replicas): %s",
            len(self._ladder),
            ", ".join(f"{t:.2f}" for t in self._ladder),
        )

    def _stage4_inp_files(self) -> None:
        logger.info("=== Stage 4/5: GENESIS .inp file rendering ===")
        assert self._tleap_result is not None
        assert self._rest_result is not None
        self._meta = SystemMetadata(
            box_size_A=self._tleap_result.box_size_A,
            n_protein_residues=self._tleap_result.n_protein_residues,
            rest_selection_string=self._rest_result.selection_string,
        )
        cfg = self.config
        meta = self._meta

        write_text(
            self.inp_dir / "step1_minimize.inp",
            inp_writer.render_minimize_inp(cfg, meta),
        )
        write_text(
            self.inp_dir / "step2_equilibrate.inp",
            inp_writer.render_equilibrate_inp(cfg, meta),
        )
        write_text(
            self.inp_dir / "step3_grest.inp",
            inp_writer.render_grest_inp(cfg, meta, self._ladder),
        )
        write_text(
            self.inp_dir / "step5_remd_convert.inp",
            inp_writer.render_remd_convert_inp(cfg, meta, self._ladder),
        )
        logger.info("  4 .inp files written -> %s", self.inp_dir)

    def _stage5_run_script(self) -> None:
        logger.info("=== Stage 5/5: run.sh + HPC scaffolds ===")
        run_path = self.output_dir / "run.sh"
        grest_runner.write_run_script(
            cfg=self.config,
            n_replicas=len(self._ladder),
            out_path=run_path,
        )
        grest_runner.write_pjsub_scaffold(self.output_dir / "run_pjsub.sh")
        grest_runner.write_sbatch_scaffold(self.output_dir / "run_sbatch.sh")
        logger.info("  run.sh -> %s", run_path)

    def _collect_result(self) -> Dict[str, Any]:
        assert self._tleap_result is not None
        assert self._rest_result is not None
        assert self._meta is not None
        return {
            "output_dir": self.output_dir,
            "build_dir": self.build_dir,
            "inp_dir": self.inp_dir,
            "config_json": self.input_dir / "config.json",
            "prmtop": self._tleap_result.prmtop,
            "coor": self._tleap_result.coor,
            "ref_pdb": self._tleap_result.ref_pdb,
            "leap_log": self._tleap_result.leap_log,
            "box_size_A": self._tleap_result.box_size_A,
            "n_protein_residues": self._tleap_result.n_protein_residues,
            "rest_residues": self._rest_result.residues,
            "rest_selection_string": self._rest_result.selection_string,
            "n_rest_residues": self._rest_result.n_residues,
            "temperature_ladder": self._ladder,
            "n_replicas": len(self._ladder),
            "inp_files": {
                "minimize": self.inp_dir / "step1_minimize.inp",
                "equilibrate": self.inp_dir / "step2_equilibrate.inp",
                "grest": self.inp_dir / "step3_grest.inp",
                "remd_convert": self.inp_dir / "step5_remd_convert.inp",
            },
            "run_script": self.output_dir / "run.sh",
            "run_pjsub": self.output_dir / "run_pjsub.sh",
            "run_sbatch": self.output_dir / "run_sbatch.sh",
        }


__all__ = ["GrestBuilder"]
