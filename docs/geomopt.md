# abmptools.geomopt — Structure Optimization and Minimization

## Overview

`abmptools.geomopt` provides two parallel classes for PDB structure
relaxation, each using a different backend:

| Class | Backend | Convergence criterion | Energy unit |
|-------|---------|----------------------|-------------|
| `MacePdbOptimizer` | MACE ML potential + ASE | `fmax` [eV/Å] | eV |
| `OpenFFOpenMMMinimizer` | OpenFF force field + OpenMM | `tolerance` [kJ/mol/nm] | kJ/mol |

Both classes follow the same conventions:
- All heavy dependencies are imported lazily (method call time only)
- Importable in any environment without optional packages
- `in_pdb` is never overwritten; `out_pdb` is required
- Failed calls raise exceptions; callers decide whether to skip or abort
- Logging via the standard `logging` module

---

## Installation of optional dependencies

### MACE optimizer

```bash
pip install ase
pip install mace-torch          # pulls in torch automatically
# GPU: install the CUDA-enabled torch build first
# See https://pytorch.org for platform-specific instructions
```

### OpenFF + OpenMM minimizer

```bash
# Core (required)
conda install -c conda-forge openmm openff-toolkit

# Bond inference from PDB files (strongly recommended)
conda install -c conda-forge rdkit

# AM1-BCC partial charge assignment (required by most OpenFF force fields)
conda install -c conda-forge ambertools
# Alternative: pip install openff-nagl
```

---

## API reference

### `MacePdbOptimizer`

```python
from abmptools.geomopt import MacePdbOptimizer

opt = MacePdbOptimizer(
    model_path=None,        # path to custom .model/.pt file (optional)
    model_name="small",     # MACE-OFF pretrained: "small"|"medium"|"large"
    device="auto",          # "auto"|"cpu"|"cuda"|"mps"
    dtype="float64",        # "float32"|"float64"
    optimizer="BFGS",       # "BFGS"|"FIRE"|"LBFGS"
    fmax=0.05,              # force convergence criterion [eV/Å]
    steps=500,              # max optimization steps
    trajectory=None,        # ASE .traj output path (used when write_traj=True)
    logfile=None,           # ASE optimizer log (None = silent, "-" = stdout)
    seed=None,              # random seed for reproducibility
)

result = opt.optimize_pdb("input.pdb", "output.pdb")
# result = {
#     "energy":    -1234.56,    # eV  (final potential energy)
#     "steps":     42,          # optimization steps taken
#     "converged": True,
#     "out_pdb":   "/abs/path/to/output.pdb",
# }
```

#### Key parameters

| Parameter    | Default     | Description |
|--------------|-------------|-------------|
| `model_path` | `None`      | Path to a custom MACE model file. When set, `model_name` is ignored. |
| `model_name` | `"small"`   | MACE-OFF pre-trained model size. |
| `device`     | `"auto"`    | Auto-selects CUDA > MPS > CPU. |
| `optimizer`  | `"BFGS"`    | ASE optimizer. BFGS is recommended for most systems. |
| `fmax`       | `0.05`      | Maximum force component [eV/Å] for convergence. |
| `steps`      | `500`       | Hard step limit. Structures are written even if not converged. |

---

### `OpenFFOpenMMMinimizer`

```python
from abmptools.geomopt import OpenFFOpenMMMinimizer

mini = OpenFFOpenMMMinimizer(
    forcefield="openff_unconstrained-2.1.0.offxml",  # OpenFF force field
    platform="auto",          # "auto"|"CUDA"|"OpenCL"|"CPU"|"HIP"
    tolerance=10.0,           # RMS force convergence [kJ/mol/nm]
    max_iterations=0,         # 0 = run until convergence
    constraints=None,         # None|"HBonds"|"AllBonds"|"HAngles"
    hydrogen_mass=None,       # HMR [Da] — stored, not yet applied
    logfile=None,             # log file path (used when write_log=True)
    seed=None,                # seeds random and numpy.random
)

result = mini.minimize_pdb("input.pdb", "output.pdb")
# result = {
#     "energy_before": -100.0,  # kJ/mol (pre-minimization PE)
#     "energy_after":  -120.0,  # kJ/mol (post-minimization PE)
#     "energy":        -120.0,  # kJ/mol  (alias for energy_after)
#     "converged":     True,
#     "elapsed":       1.5,     # wall-clock seconds
#     "out_pdb":       "/abs/path/to/output.pdb",
# }
```

#### Key parameters

| Parameter        | Default                                    | Description |
|------------------|--------------------------------------------|-------------|
| `forcefield`     | `"openff_unconstrained-2.1.0.offxml"`      | OpenFF force field name or path. Use "unconstrained" for minimization so all DOF can relax. |
| `platform`       | `"auto"`                                   | Auto-selects CUDA > OpenCL > CPU. |
| `tolerance`      | `10.0`                                     | RMS force threshold [kJ/mol/nm]. Lower = tighter convergence. |
| `max_iterations` | `0`                                        | Step limit. `0` = unlimited (run until convergence). |
| `constraints`    | `None`                                     | Bond constraints. `None` recommended for minimization. |

#### PDB → OpenFF molecule loading

PDB files often lack bond/connectivity information.  The class tries two
strategies in order:

1. **`Molecule.from_file`** — works when the PDB has `CONECT` records.
2. **RDKit `MolFromPDBFile` → `Molecule.from_rdkit`** — infers bonds for
   most organic small molecules (requires `rdkit`).

If both fail, a `RuntimeError` is raised with a diagnosis and suggested fix.

**Atom count check**: after loading, the class verifies that the OpenFF
molecule and OpenMM PDB reader report the same atom count.  A mismatch
(e.g., from missing hydrogens) raises `ValueError` with a clear message.

---

## Side-by-side comparison

```python
import logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

from abmptools.geomopt import MacePdbOptimizer, OpenFFOpenMMMinimizer

# MACE: ML-level accuracy, requires GPU for large systems
r1 = MacePdbOptimizer(model_name="small", fmax=0.05).optimize_pdb(
    "mol.pdb", "mol_mace.pdb"
)
print(f"MACE : {r1['energy']:.4f} eV  steps={r1['steps']}  converged={r1['converged']}")

# OpenFF: classical force field, fast on CPU, interpretable energy terms
r2 = OpenFFOpenMMMinimizer(tolerance=10.0).minimize_pdb(
    "mol.pdb", "mol_openff.pdb"
)
print(f"OpenFF: {r2['energy']:.4f} kJ/mol  converged={r2['converged']}  elapsed={r2['elapsed']:.2f}s")
```

---

## Smoke tests

### MacePdbOptimizer

```bash
python - <<'EOF'
import logging
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s %(name)s: %(message)s")

from abmptools.geomopt import MacePdbOptimizer

opt = MacePdbOptimizer(model_name="small", device="auto", fmax=0.05, steps=100)
result = opt.optimize_pdb("water.pdb", "water_mace.pdb")

print("energy   :", result["energy"], "eV")
print("steps    :", result["steps"])
print("converged:", result["converged"])
print("output   :", result["out_pdb"])
EOF
```

### OpenFFOpenMMMinimizer

```bash
python - <<'EOF'
import logging
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s %(name)s: %(message)s")

from abmptools.geomopt import OpenFFOpenMMMinimizer

mini = OpenFFOpenMMMinimizer(
    forcefield="openff_unconstrained-2.1.0.offxml",
    platform="CPU",
    tolerance=10.0,
)
result = mini.minimize_pdb("mol.pdb", "mol_min.pdb", write_log=True)

print("energy_before:", result["energy_before"], "kJ/mol")
print("energy_after :", result["energy_after"], "kJ/mol")
print("converged    :", result["converged"])
print("elapsed      :", result["elapsed"], "s")
print("output       :", result["out_pdb"])
EOF
```

---

## Notes

- `in_pdb` and `out_pdb` must differ; passing the same path raises `ValueError`.
- Unconverged structures are still written (check `result["converged"]`).
- **MACE**: write an ASE trajectory with `optimize_pdb(..., write_traj=True)`.
- **OpenFF**: write a text log with `minimize_pdb(..., write_log=True)`.
- Logging is controlled via the standard `logging` module.
  Set `logging.basicConfig(level=logging.DEBUG)` for verbose output.
- **OpenFF multi-molecule**: currently only single-molecule PDB is supported.
  Multi-component systems (e.g., solute + solvent) require additional topology
  handling and are planned for a future release.
