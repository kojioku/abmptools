# -*- coding: utf-8 -*-
"""Unit tests for ``abmptools.crystal.job_templates``."""
from __future__ import annotations

from pathlib import Path

import pytest

from abmptools.crystal.job_templates import (
    PJM_TEMPLATE,
    SLURM_TEMPLATE,
    PBS_TEMPLATE,
    LOCAL_TEMPLATE,
    render_jobscript,
    write_batch_runner,
)
from abmptools.crystal.models import HPCJobSpec


# ---------------------------------------------------------------------------
# render_jobscript
# ---------------------------------------------------------------------------

def test_pjm_render_matches_runab23q2_shape():
    """PJM render must include the Fugaku-style pjsub directives."""
    spec = HPCJobSpec(
        scheduler="PJM",
        queue="small",
        group="hp190133",
        nodes=12,
        proc_per_node=2,
        omp_threads=24,
        elapse="24:00:00",
        abinit_dir="/data/hp190133/programs/ABINIT-MP",
        binary_name="abinitmp_smp",
    )
    out = render_jobscript(spec, "XXXI-MMFF-R00001layer5Zp1-around_ar6.0.ajf")

    assert '#PJM -L "rscgrp=small"' in out
    assert '#PJM -L "node=12"' in out
    assert '#PJM --mpi "proc=24,max-proc-per-node=2"' in out  # 12*2=24
    assert '#PJM -g "hp190133"' in out
    assert "OMP_NUM_THREADS=24" in out
    assert "/data/hp190133/programs/ABINIT-MP/abinitmp_smp" in out
    assert "XXXI-MMFF-R00001layer5Zp1-around_ar6.0.ajf" in out


def test_slurm_render_basics():
    spec = HPCJobSpec(
        scheduler="SLURM", queue="gpu",
        nodes=4, proc_per_node=8, omp_threads=4,
        elapse="12:00:00",
        abinit_dir="/opt/abinitmp",
    )
    out = render_jobscript(spec, "test.ajf")
    assert "#SBATCH --partition=gpu" in out
    assert "#SBATCH --nodes=4" in out
    assert "#SBATCH --ntasks-per-node=8" in out
    assert "mpirun -np 32" in out  # 4*8


def test_pbs_render_basics():
    spec = HPCJobSpec(
        scheduler="PBS", queue="long",
        nodes=2, proc_per_node=24, omp_threads=2,
        elapse="48:00:00",
        abinit_dir="/usr/local/abinitmp",
    )
    out = render_jobscript(spec, "test.ajf")
    assert "#PBS -q long" in out
    assert "select=2:ncpus=48:mpiprocs=24" in out
    assert "$PBS_O_WORKDIR" in out  # double-dollar escape verified


def test_local_render_basics():
    spec = HPCJobSpec(
        scheduler="local",
        nodes=1, proc_per_node=1, omp_threads=4,
        elapse="01:00:00",
        abinit_dir="/usr/local/abinitmp",
    )
    out = render_jobscript(spec, "test.ajf")
    assert "set -euo pipefail" in out
    assert "/usr/local/abinitmp/abinitmp_smp" in out
    assert "OMP_NUM_THREADS=4" in out


def test_template_override(tmp_path):
    """Custom template paths bypass the bundled scheduler templates."""
    custom = tmp_path / "custom.tmpl"
    custom.write_text("CUSTOM ${queue} ${nodes} ${total_proc}\n")

    spec = HPCJobSpec(
        scheduler="PJM", queue="custom_q", nodes=3, proc_per_node=4,
        template_override=str(custom),
    )
    out = render_jobscript(spec, "x.ajf")
    assert out == "CUSTOM custom_q 3 12\n"


def test_extra_lines_are_appended():
    spec = HPCJobSpec(
        scheduler="local", abinit_dir="/x",
        extra_lines=["echo done", "echo bye"],
    )
    out = render_jobscript(spec, "t.ajf")
    assert out.rstrip().endswith("echo done\necho bye")


# ---------------------------------------------------------------------------
# write_batch_runner
# ---------------------------------------------------------------------------

def test_write_batch_runner_pjsub(tmp_path):
    runner = write_batch_runner(
        ajf_paths=["a.ajf", "b.ajf"],
        output_dir=str(tmp_path),
        submit_command="pjsub",
    )
    body = Path(runner).read_text()
    assert body.startswith("#!/bin/bash")
    assert 'pjsub "a.sh"' in body
    assert 'pjsub "b.sh"' in body
    # Executable bit set.
    import stat
    mode = Path(runner).stat().st_mode
    assert mode & stat.S_IXUSR


def test_write_batch_runner_dry_run(tmp_path):
    runner = write_batch_runner(
        ajf_paths=["a.ajf"],
        output_dir=str(tmp_path),
    )
    body = Path(runner).read_text()
    assert 'echo "would submit a.sh"' in body
