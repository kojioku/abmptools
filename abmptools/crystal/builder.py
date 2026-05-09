# -*- coding: utf-8 -*-
"""
abmptools.crystal.builder
--------------------------
Top-level orchestrator for the crystal-FMO pipeline.

:class:`CrystalOrchestrator` runs the 5-stage workflow:

1. ``expand_cif``    -- per-input CIF -> ``cifout/layer<L>/{pdb,xyz}/``
                         (legacy or ase backend, selected by
                         :attr:`CIFEngineConfig.engine`)
2. ``generate_fmo``  -- copy XYZ + drivers next to the supercell PDB,
                         invoke ``python -m abmptools.pdb2fmo -xyz``,
                         and collect ``for_abmp/*.{ajf,pdb}``
3. ``write_jobs``    -- render per-AJF jobscripts (PJM/SLURM/PBS/local)
                         and a batch submitter
4. ``run_abinit``    -- (Phase D) optional local ``abinitmp`` smoke run
5. ``postprocess``   -- (Phase C-6) optional IFIE/PIEDA + nearest-atom
                         annotation

Stage 4 and 5 are stubbed out in Phase C-5 (``NotImplementedError``);
the orchestrator always executes stages 1-3 and they are sufficient
to drop into the historical csp7 workflow.

The orchestrator drives :mod:`abmptools.pdb2fmo` via the in-process
:func:`abmptools.pdb2fmo.run_pdb2fmo` API. Each input PDB is processed
with a fresh ``abmp.setfmo`` instance (matching the legacy CLI loop),
so per-input state isolation is preserved without paying the Python
interpreter start-up cost of ``subprocess`` per file. The CLI form
(``python -m abmptools.pdb2fmo``) remains untouched and continues to
delegate to the same ``run_pdb2fmo`` implementation.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from abmptools.pdb2fmo import run_pdb2fmo

from . import job_templates
from ._subprocess import ensure_dir
from .cif_engine_legacy import run_legacy
from .models import (
    CIFInputSpec,
    CrystalBuildConfig,
    HPCJobSpec,
)

logger = logging.getLogger("abmptools.crystal.builder")


# ---------------------------------------------------------------------------
# Driver-file emission (input_param / segment_data.dat)
# ---------------------------------------------------------------------------

def _emit_input_param(cfg: CrystalBuildConfig) -> str:
    """Render the ``input_param`` Python-dict literal consumed by
    :mod:`abmptools.pdb2fmo`.

    Mirrors the schema validated by :func:`setfmo.setrfmoparam`.
    """
    p = cfg.fragment
    f = cfg.fmo
    return (
        "param = {\n"
        f"    'cutmode': {p.cutmode!r},\n"
        f"    'solutes': {list(p.solutes)!r},\n"
        f"    'criteria': {p.criteria!r},\n"
        f"    'tgtpos': {list(p.tgtpos)!r},\n"
        f"    'molname': {list(p.molname)!r},\n"
        f"    'getmode': {p.getmode!r},\n"
        f"    'piedaflag': {p.pieda!r},\n"
        f"    'cmmflag': {p.cmm!r},\n"
        f"    'solv_flag': {p.solv_flag!r},\n"
        f"    'ajf_method': {f.method!r},\n"
        f"    'ajf_basis_set': {f.basis_set!r},\n"
        f"    'cpfflag': {f.cpfflag!r},\n"
        f"    'abinit_ver': {f.abinit_ver!r},\n"
        f"    'memory': {f.memory!r},\n"
        f"    'npro': {f.npro!r},\n"
        "}\n"
    )


def _emit_segment_data(spec: CIFInputSpec, molname: str) -> str:
    """Render a minimal ``segment_data.dat`` from one CIF input.

    Single-segment, single-fragment-per-mol convention (``connect=[]``,
    each molecule is one fragment containing all atoms). Suitable for
    the csp7 demo where each cocrystal molecule is treated as one
    fragment; users with custom fragmentation should set
    :attr:`FragmentTemplate.seg_data_path` instead.
    """
    n_atoms = spec.atoms_in_mol[0]
    seg_info = list(range(1, n_atoms + 1))
    return (
        "seg_data = [\n"
        "    {\n"
        f"    'name': {molname!r},\n"
        f"    'atom':  [{n_atoms}] ,\n"
        "    'charge':  [0] ,\n"
        "    'connect_num':  [0] ,\n"
        f"    'seg_info':  [{seg_info!r}] ,\n"
        "    'connect':  [] ,\n"
        "    'nummol_seg': [1],\n"
        "    'repeat': [1],\n"
        "    'term': [],\n"
        "    'main_seg': [1],\n"
        "    },\n"
        "]\n"
    )


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

class CrystalOrchestrator:
    """Run the crystal-FMO pipeline end to end.

    Parameters
    ----------
    config
        Validated :class:`CrystalBuildConfig`.
    config_dir
        Directory where the config file lives. Used to resolve relative
        CIF paths and ``output_dir``. Defaults to current working dir.
    """

    def __init__(
        self,
        config: CrystalBuildConfig,
        *,
        config_dir: Optional[str] = None,
    ) -> None:
        self.config = config
        self.config_dir = Path(config_dir or os.getcwd()).resolve()
        self.output_dir = self._resolve(config.output_dir)
        ensure_dir(self.output_dir)

        # Per-input subdirectory cache (populated stage by stage).
        self.expand_results: Dict[str, Dict[str, Path]] = {}
        self.fmo_results: Dict[str, Dict[str, List[Path]]] = {}
        self.job_scripts: List[Path] = []

    def _resolve(self, p: str) -> Path:
        path = Path(p)
        if not path.is_absolute():
            path = (self.config_dir / path).resolve()
        return path

    def _input_name(self, spec: CIFInputSpec) -> str:
        return spec.name or Path(spec.cif).stem

    def _per_input_dir(self, spec: CIFInputSpec) -> Path:
        return self.output_dir / self._input_name(spec)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self) -> Dict[str, Any]:
        """Run all stages and return a summary dict."""
        self._dump_config()
        self.expand_cif()
        self.generate_fmo()
        self.write_jobs()
        if self.config.run_local:
            self.run_abinit()
        if self.config.postproc.enable:
            self.postprocess()
        return self._summary()

    # ------------------------------------------------------------------
    # Stage 1 — CIF expansion
    # ------------------------------------------------------------------

    def expand_cif(self) -> Dict[str, Dict[str, Path]]:
        """Expand each CIF into ``layer<L>/{pdb,xyz}/`` outputs.

        ``CIFEngineConfig.engine='legacy'`` (default) drives the
        validated :func:`abmptools.readcif.run_legacy_cif_pipeline`
        wrapper. ``engine='ase'`` is reserved for the Phase D PDB
        emitter; in Phase C-5 it raises :class:`NotImplementedError`
        because a byte-equivalent PDB writer is not yet available.
        """
        engine = self.config.cif_engine.engine
        for spec in self.config.inputs:
            cif_path = self._resolve(spec.cif)
            in_dir = ensure_dir(self._per_input_dir(spec))

            # Stage outputs: <in_dir>/cifout/layer<L>/{pdb,xyz}
            shutil.copy(cif_path, in_dir / cif_path.name)
            if engine == "legacy":
                layer_dir = run_legacy(
                    cif=cif_path.name,
                    layer=spec.layer,
                    atoms_in_mol=spec.atoms_in_mol,
                    odir="cifout",
                    asymmetric_only=spec.asymmetric_only,
                    cwd=str(in_dir),
                )
            elif engine == "ase":
                from .cif_engine_ase import run_ase_pipeline

                layer_dir = run_ase_pipeline(
                    cif=cif_path.name,
                    layer=spec.layer,
                    atoms_in_mol=spec.atoms_in_mol,
                    bond_tolerance=self.config.cif_engine.bond_tolerance,
                    odir="cifout",
                    molname=self.config.fragment.molname[0],
                    cwd=str(in_dir),
                )
            else:
                raise ValueError(f"unknown engine: {engine!r}")

            self.expand_results[self._input_name(spec)] = {
                "layer_dir": layer_dir,
                "pdb_dir": layer_dir / "pdb",
                "xyz_dir": layer_dir / "xyz",
                "input_dir": in_dir,
            }
            logger.info(
                "expand_cif done: %s -> %s",
                self._input_name(spec), layer_dir,
            )
        return self.expand_results

    # ------------------------------------------------------------------
    # Stage 2 — FMO input generation
    # ------------------------------------------------------------------

    def generate_fmo(self) -> Dict[str, Dict[str, List[Path]]]:
        """Invoke :func:`abmptools.pdb2fmo.run_pdb2fmo` on each supercell PDB.

        Drivers (``input_param`` / ``segment_data.dat`` / ``UNK.ajf``)
        are still written into the layer<L>/pdb directory so a user can
        rerun the legacy CLI (``python -m abmptools.pdb2fmo``) on the
        same outputs interchangeably; ``run_pdb2fmo`` itself takes the
        param dict directly. The XYZ file with matching basename is
        copied alongside so :func:`pdb_io.exportardxyzfull` can read it.
        """
        if not self.expand_results:
            self.expand_cif()

        for spec in self.config.inputs:
            name = self._input_name(spec)
            paths = self.expand_results[name]
            pdb_dir: Path = paths["pdb_dir"]
            xyz_dir: Path = paths["xyz_dir"]

            self._stage_drivers(spec, pdb_dir)
            xyzs = list(xyz_dir.glob("*.xyz"))
            for xyz in xyzs:
                shutil.copy(xyz, pdb_dir / xyz.name)

            pdbs = sorted(pdb_dir.glob("*layer*.pdb"))
            if not pdbs:
                raise RuntimeError(
                    f"no layer PDB files in {pdb_dir} for input {name!r}"
                )

            param = self._build_param_dict()
            logger.info(
                "generate_fmo (in-process): %s, %d pdb(s), is_xyz=%s",
                name, len(pdbs), self.config.fmo.is_xyz,
            )
            prev_cwd = os.getcwd()
            os.chdir(pdb_dir)
            try:
                run_pdb2fmo(
                    pdb_files=[p.name for p in pdbs],
                    param=param,
                    xyz=self.config.fmo.is_xyz,
                    output_dir="for_abmp",
                    verbose=False,
                )
            finally:
                os.chdir(prev_cwd)

            for_abmp_dir = pdb_dir / "for_abmp"
            ajfs = sorted(for_abmp_dir.glob("*.ajf"))
            out_pdbs = sorted(for_abmp_dir.glob("*.pdb"))
            self.fmo_results[name] = {
                "ajf_files": ajfs,
                "pdb_files": out_pdbs,
                "for_abmp_dir": for_abmp_dir,
            }
            logger.info(
                "generate_fmo done: %s -> %d ajf(s)", name, len(ajfs),
            )
        return self.fmo_results

    def _build_param_dict(self) -> Dict[str, Any]:
        """Build the param dict consumed by :meth:`abmp.setfmo.setrfmoparam`.

        Mirrors :func:`_emit_input_param` (the file form) but skips the
        textual round-trip — we go directly from ``CrystalBuildConfig``
        to the dict that ``run_pdb2fmo`` accepts.
        """
        p = self.config.fragment
        f = self.config.fmo
        return {
            "cutmode": p.cutmode,
            "solutes": list(p.solutes),
            "criteria": p.criteria,
            "tgtpos": list(p.tgtpos),
            "molname": list(p.molname),
            "getmode": p.getmode,
            "piedaflag": p.pieda,
            "cmmflag": p.cmm,
            "solv_flag": p.solv_flag,
            "ajf_method": f.method,
            "ajf_basis_set": f.basis_set,
            "cpfflag": f.cpfflag,
            "abinit_ver": f.abinit_ver,
            "memory": f.memory,
            "npro": f.npro,
        }

    def _stage_drivers(self, spec: CIFInputSpec, pdb_dir: Path) -> None:
        """Write input_param / segment_data.dat / UNK.ajf into *pdb_dir*."""
        # input_param (always synthesised from config).
        (pdb_dir / "input_param").write_text(_emit_input_param(self.config))

        # segment_data.dat (custom path overrides the synthesised default).
        seg_path = self.config.fragment.seg_data_path
        if seg_path:
            shutil.copy(self._resolve(seg_path), pdb_dir / "segment_data.dat")
        else:
            molname = self.config.fragment.molname[0]
            (pdb_dir / "segment_data.dat").write_text(
                _emit_segment_data(spec, molname)
            )

        # UNK.ajf — required by setfmo as a template for the AJF body.
        tmpl_path = self.config.fragment.template_ajf
        if not tmpl_path:
            raise ValueError(
                "FragmentTemplate.template_ajf is required (path to a "
                "minimal ABINIT-MP AJF template, e.g. UNK.ajf). "
                "See sample/crystal/csp7_smoke/UNK.ajf for the canonical "
                "csp7 template."
            )
        shutil.copy(self._resolve(tmpl_path), pdb_dir / "UNK.ajf")

    # ------------------------------------------------------------------
    # Stage 3 — HPC job scripts
    # ------------------------------------------------------------------

    def write_jobs(self) -> List[Path]:
        """Render per-AJF jobscripts and a batch submitter."""
        if not self.fmo_results:
            self.generate_fmo()

        scripts: List[Path] = []
        for spec in self.config.inputs:
            name = self._input_name(spec)
            for_abmp_dir = self.fmo_results[name]["for_abmp_dir"]
            for ajf in self.fmo_results[name]["ajf_files"]:
                body = job_templates.render_jobscript(
                    self.config.hpc, str(ajf),
                )
                base = ajf.stem
                spec_hpc: HPCJobSpec = self.config.hpc
                suffix = (
                    f"_{spec_hpc.nodes}n-{spec_hpc.proc_per_node}p-"
                    f"{spec_hpc.omp_threads}t"
                )
                script_path = for_abmp_dir / f"{base}{suffix}.sh"
                script_path.write_text(body)
                script_path.chmod(0o755)
                scripts.append(script_path)

            submit_cmd = self._submit_command(self.config.hpc.scheduler)
            ajf_basenames = [
                ajf.name for ajf in self.fmo_results[name]["ajf_files"]
            ]
            job_templates.write_batch_runner(
                ajf_paths=ajf_basenames,
                output_dir=str(for_abmp_dir),
                submit_command=submit_cmd,
            )
        self.job_scripts = scripts
        logger.info("write_jobs done: %d script(s)", len(scripts))
        return scripts

    @staticmethod
    def _submit_command(scheduler: str) -> Optional[str]:
        return {
            "PJM": "pjsub",
            "SLURM": "sbatch",
            "PBS": "qsub",
            "local": "bash",
        }.get(scheduler)

    # ------------------------------------------------------------------
    # Stage 4 — local ABINIT-MP run
    # ------------------------------------------------------------------

    def run_abinit(self, *, timeout_s: Optional[int] = None) -> List[Path]:
        """Invoke ``abinitmp`` locally on each rendered AJF.

        For each AJF in ``fmo_results``:

        1. Resolve the ``abinitmp`` binary. Order of resolution:
           a. ``hpc.abinit_dir / hpc.binary_name`` if both fields set
              and the file exists.
           b. ``shutil.which(hpc.binary_name)`` (PATH lookup).
           c. ``shutil.which("abinitmp")`` (legacy plain name).
           Missing binary raises :class:`FileNotFoundError`.
        2. ``subprocess.run([abinitmp], stdin=ajf, stdout=log,
           stderr=err, cwd=for_abmp_dir)``.
        3. Append ``log`` to the return list; raise on non-zero exit
           when ``config.fail_fast`` is True (default), warn-and-continue
           otherwise.

        Designed for smoke / debug runs (1-handful structures). Production
        runs of the csp7 1500-structure pipeline should still go through
        the HPC scheduler via the rendered jobscripts; ``--run-local``
        is not intended for that scale.

        Parameters
        ----------
        timeout_s
            Per-invocation timeout. ``None`` means no timeout (the
            FMO calculations can take from minutes to hours depending
            on system size).

        Returns
        -------
        List[Path]
            One log file per executed AJF.
        """
        if not self.fmo_results:
            self.generate_fmo()

        abinit = self._resolve_abinit_binary()
        logger.info("run_abinit: using binary %s", abinit)

        log_files: List[Path] = []
        for spec in self.config.inputs:
            name = self._input_name(spec)
            for_abmp_dir: Path = self.fmo_results[name]["for_abmp_dir"]
            for ajf in self.fmo_results[name]["ajf_files"]:
                log_path = for_abmp_dir / f"{ajf.stem}.log"
                err_path = for_abmp_dir / f"{ajf.stem}.err"
                logger.info(
                    "run_abinit: %s < %s > %s",
                    abinit.name, ajf.name, log_path.name,
                )
                with open(ajf, "rb") as fin, \
                     open(log_path, "wb") as flog, \
                     open(err_path, "wb") as ferr:
                    result = subprocess.run(
                        [str(abinit)],
                        stdin=fin, stdout=flog, stderr=ferr,
                        cwd=str(for_abmp_dir),
                        timeout=timeout_s,
                        check=False,
                    )
                if result.returncode != 0:
                    msg = (
                        f"abinitmp exited non-zero ({result.returncode}) "
                        f"for {ajf.name}; log={log_path}, err={err_path}"
                    )
                    if self.config.fail_fast:
                        raise RuntimeError(msg)
                    logger.warning("%s — continuing (fail_fast=False)", msg)
                log_files.append(log_path)
        return log_files

    def _resolve_abinit_binary(self) -> Path:
        """Locate the ``abinitmp`` executable for ``--run-local``."""
        spec = self.config.hpc
        candidates: List[Path] = []
        if spec.abinit_dir and spec.binary_name:
            candidates.append(Path(spec.abinit_dir) / spec.binary_name)
        for name in (spec.binary_name, "abinitmp"):
            which = shutil.which(name) if name else None
            if which:
                candidates.append(Path(which))
        for c in candidates:
            if c.exists():
                return c
        raise FileNotFoundError(
            "abinitmp binary not found. Set HPCJobSpec.abinit_dir + "
            f"binary_name (got abinit_dir={spec.abinit_dir!r}, "
            f"binary_name={spec.binary_name!r}) or place `abinitmp` on PATH."
        )

    # ------------------------------------------------------------------
    # Stage 5 — postproc (Phase C-6)
    # ------------------------------------------------------------------

    def postprocess(self) -> Dict[str, Any]:  # pragma: no cover - C-6
        raise NotImplementedError(
            "IFIE/PIEDA postprocessing lands in Phase C-6; until then "
            "use `python -m abmptools.getifiepieda` directly on the run logs."
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _dump_config(self) -> None:
        """Persist the resolved config next to the outputs."""
        cfg_path = self.output_dir / "crystal_config.resolved.json"
        cfg_path.write_text(
            __import__("json").dumps(
                asdict(self.config), indent=2, ensure_ascii=False,
            )
        )
        logger.info("dumped resolved config -> %s", cfg_path)

    def _summary(self) -> Dict[str, Any]:
        return {
            "output_dir": str(self.output_dir),
            "n_inputs": len(self.config.inputs),
            "n_ajf_total": sum(
                len(r["ajf_files"]) for r in self.fmo_results.values()
            ),
            "n_jobscripts": len(self.job_scripts),
            "engine": self.config.cif_engine.engine,
            "is_xyz": self.config.fmo.is_xyz,
        }


__all__ = ["CrystalOrchestrator", "_emit_input_param", "_emit_segment_data"]
