# abmptools.mlopt — ML-Based Structure Optimization

## Overview

`abmptools.mlopt` provides PDB structure optimization using
[MACE](https://github.com/ACEsuit/mace) machine-learning potentials
via [ASE](https://wiki.fysik.dtu.dk/ase/).

The subpackage is importable even when the optional dependencies
(ase, mace-torch, torch) are not installed; dependency errors are
deferred to the point where optimization is actually attempted.

---

## Installation of optional dependencies

```bash
pip install ase
pip install mace-torch          # pulls in torch automatically
# GPU: install the CUDA-enabled torch build first
# See https://pytorch.org for platform-specific instructions
```

---

## API

### `MacePdbOptimizer`

```python
from abmptools.mlopt import MacePdbOptimizer

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
#     "energy":    -1234.56,    # eV
#     "steps":     42,
#     "converged": True,
#     "out_pdb":   "/abs/path/to/output.pdb",
# }
```

`optimize_pdb` raises exceptions on failure (never overwrites `in_pdb`).
Callers are responsible for catch/continue logic.

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

## Smoke test

Requires a PDB file (e.g., `water.pdb`) and optional dependencies installed.

```bash
python - <<'EOF'
import logging
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s %(name)s: %(message)s")

from abmptools.mlopt import MacePdbOptimizer

opt = MacePdbOptimizer(model_name="small", device="auto", fmax=0.05, steps=100)
result = opt.optimize_pdb("water.pdb", "water_opt.pdb")

print("energy   :", result["energy"], "eV")
print("steps    :", result["steps"])
print("converged:", result["converged"])
print("output   :", result["out_pdb"])
EOF
```

Expected output (values will vary by system):

```
INFO abmptools.mlopt.mace_optimizer: Optimizing: /abs/water.pdb -> /abs/water_opt.pdb
INFO abmptools.mlopt.mace_optimizer: Converged in 7 steps, energy = -14.123456 eV
INFO abmptools.mlopt.mace_optimizer: Written: /abs/water_opt.pdb
energy   : -14.123456 eV
steps    : 7
converged: True
output   : /abs/water_opt.pdb
```

---

## Notes

- `in_pdb` and `out_pdb` must differ; passing the same path raises `ValueError`.
- Unconverged structures are still written (check `result["converged"]`).
- To write a trajectory (for inspection with `ase gui`):
  ```python
  result = opt.optimize_pdb("in.pdb", "out.pdb", write_traj=True)
  ```
- Logging is controlled via the standard `logging` module.
  Set `logging.basicConfig(level=logging.DEBUG)` for verbose output.
