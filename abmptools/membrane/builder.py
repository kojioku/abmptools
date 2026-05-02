# -*- coding: utf-8 -*-
"""
abmptools.membrane.builder
--------------------------
Orchestrator for the membrane-US pipeline.

Stages
------
1. Bilayer construction      (bilayer.build_bilayer)
2. Force-field parameterise  (parameterize_amber / parameterize_charmm)
3. Equilibration MDPs        (mdp_us_protocol.write_equilibration_mdps)
4. Pulling MDP               (pulling.write_pulling_mdp)
5. Umbrella window MDPs      (umbrella.write_window_mdps)
6. PMF analysis (post-MD)    (pmf.run_wham)

Stages 1–5 produce all input files; the MD itself is run by an external
script generated alongside (run.sh). Stage 6 is invoked separately as
``builder.analyze_pmf()`` once windows have finished.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict

from .models import MembraneConfig

logger = logging.getLogger(__name__)


class MembraneUSBuilder:
    """Build a peptide-in-bilayer umbrella-sampling system for GROMACS.

    Usage::

        from abmptools.membrane import (
            MembraneUSBuilder, MembraneConfig,
            LipidSpec, PeptideSpec,
        )

        cfg = MembraneConfig(
            backend="amber",
            lipids=[LipidSpec(resname="POPC", n_per_leaflet=64)],
            peptide=PeptideSpec(name="poly_ala", sequence="AAAAA"),
            output_dir="./run01",
        )
        builder = MembraneUSBuilder(cfg)
        result = builder.build()
        # ... run MD via result['run_script'] ...
        builder.analyze_pmf()
    """

    def __init__(self, config: MembraneConfig) -> None:
        self.config = config
        self._build_dir: Path = Path()
        self._equil_dir: Path = Path()
        self._pull_dir: Path = Path()
        self._windows_dir: Path = Path()
        self._analysis_dir: Path = Path()
        # Stage-2 outputs cached for downstream stages.
        self._top: str = ""
        self._gro: str = ""
        self._ndx: str = ""

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(self) -> Dict[str, Any]:
        """Run stages 1–5 (everything up to but not including MD).

        Returns
        -------
        dict
            Keys: ``"output_dir"``, ``"top"``, ``"gro"``, ``"ndx"``,
            ``"equil_mdps"``, ``"pull_mdp"``, ``"window_mdps"``,
            ``"run_script"``, ``"config_json"``.
        """
        self._setup_dirs()
        self._save_config()

        logger.info("=== Stage 1: bilayer construction ===")
        bilayer_pdb = self._stage1_bilayer()

        logger.info("=== Stage 2: force-field parameterisation (%s) ===",
                    self.config.backend)
        top, gro, ndx = self._stage2_parameterize(bilayer_pdb)
        self._top, self._gro, self._ndx = top, gro, ndx

        logger.info("=== Stage 3: equilibration MDPs ===")
        equil_mdps = self._stage3_equilibration_mdps()

        logger.info("=== Stage 4: pulling MDP ===")
        pull_mdp = self._stage4_pulling_mdp()

        logger.info("=== Stage 5: umbrella window MDPs ===")
        window_mdps = self._stage5_window_mdps()

        logger.info("=== Stage 6 deferred (run MD first, then analyze_pmf) ===")
        run_script = self._write_run_script(
            top=top, gro=gro, ndx=ndx,
            equil_mdps=equil_mdps, pull_mdp=pull_mdp, window_mdps=window_mdps,
        )

        return {
            "output_dir": str(Path(self.config.output_dir).resolve()),
            "top": top,
            "gro": gro,
            "ndx": ndx,
            "equil_mdps": equil_mdps,
            "pull_mdp": pull_mdp,
            "window_mdps": window_mdps,
            "run_script": run_script,
            "config_json": str((Path(self.config.output_dir) / "input"
                                / "config.json").resolve()),
        }

    def analyze_pmf(self) -> Dict[str, Any]:
        """Run stage 6 (PMF analysis) after MD windows have completed.

        Reads ``windows/*/pullx.xvg`` and ``pullf.xvg``, runs ``gmx wham``
        and/or PyMBAR, writes ``analysis/pmf.xvg`` + plots.
        """
        from . import pmf
        return pmf.run_wham(
            windows_dir=str(self._windows_dir),
            analysis_dir=str(self._analysis_dir),
            config=self.config,
        )

    # ------------------------------------------------------------------
    # Stage delegates (kept short — heavy lifting lives in submodules)
    # ------------------------------------------------------------------

    def _setup_dirs(self) -> None:
        out = Path(self.config.output_dir)
        self._build_dir = out / "build"
        self._equil_dir = out / "equil"
        self._pull_dir = out / "pull"
        self._windows_dir = out / "windows"
        self._analysis_dir = out / "analysis"
        for d in (
            out / "input", self._build_dir, self._equil_dir,
            self._pull_dir, self._windows_dir, self._analysis_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)

    def _save_config(self) -> None:
        path = Path(self.config.output_dir) / "input" / "config.json"
        self.config.to_json(str(path))

    def _stage1_bilayer(self) -> str:
        from . import bilayer
        return bilayer.build_bilayer(
            config=self.config, build_dir=str(self._build_dir),
        )

    def _stage2_parameterize(self, bilayer_pdb: str) -> tuple[str, str, str]:
        if self.config.backend == "amber":
            from . import parameterize_amber as backend
        elif self.config.backend == "charmm36":
            from . import parameterize_charmm as backend
        else:
            raise ValueError(
                f"Unknown backend: {self.config.backend!r}"
            )
        return backend.parameterize(
            config=self.config,
            input_pdb=bilayer_pdb,
            build_dir=str(self._build_dir),
        )

    def _stage3_equilibration_mdps(self) -> Dict[str, str]:
        from . import mdp_us_protocol
        return mdp_us_protocol.write_equilibration_mdps(
            config=self.config, equil_dir=str(self._equil_dir),
        )

    def _stage4_pulling_mdp(self) -> str:
        from . import pulling
        # Estimate the initial pull-coord from the *unequilibrated* system.
        # Equilibration typically shifts the peptide by < 0.5 nm; if the
        # post-NPT shift is significant, regenerate pull.mdp before stage 4
        # via pulling.estimate_initial_pull_coord(equil/npt.gro, ...).
        init_z = pulling.estimate_initial_pull_coord(
            gro_path=self._gro,
            ndx_path=self._ndx,
            pull_group1=self.config.umbrella.pull_group1,
            pull_group2=self.config.umbrella.pull_group2,
        )
        logger.info("estimated initial pull-coord: %+.3f nm", init_z)
        # Bilayer always needs an explicit pbcatom (group spans most of
        # the box xy); pre-compute the atom closest to the bilayer COM.
        pbc_g1 = pulling.find_pbc_center_atom(
            gro_path=self._gro, ndx_path=self._ndx,
            group_name=self.config.umbrella.pull_group1,
        )
        logger.info("pbc-atom for %s: %d",
                    self.config.umbrella.pull_group1, pbc_g1)
        self._pbc_atom_g1 = pbc_g1
        return pulling.write_pulling_mdp(
            config=self.config, pull_dir=str(self._pull_dir),
            pull_init_nm=init_z, pbc_atom_g1=pbc_g1,
        )

    def _stage5_window_mdps(self) -> Dict[int, str]:
        from . import umbrella
        return umbrella.write_window_mdps(
            config=self.config, windows_dir=str(self._windows_dir),
            pbc_atom_g1=getattr(self, "_pbc_atom_g1", None),
        )

    def _write_run_script(
        self, *, top: str, gro: str, ndx: str,
        equil_mdps: Dict[str, str], pull_mdp: str,
        window_mdps: Dict[int, str],
    ) -> str:
        """Generate a top-level run script tying all stages together.

        The actual implementation lives in :mod:`umbrella` (since per-window
        invocation logic dominates); this is a thin pass-through.
        """
        from . import umbrella
        return umbrella.write_run_script(
            config=self.config,
            top=top, gro=gro, ndx=ndx,
            equil_mdps=equil_mdps, pull_mdp=pull_mdp,
            window_mdps=window_mdps,
            output_dir=self.config.output_dir,
        )
