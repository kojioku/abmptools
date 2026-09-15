# pentane / benzene mixture (from a JSON config)

The same system as [`../pentane_benzene/`](../pentane_benzene/) — 200
pentane + 50 benzene at 0.8 g/cm^3, 300 K — but the composition comes
from **`mixture.json`** instead of command-line flags.

Use this form when the component list outgrows a command line, or when
the composition should sit in version control beside the results.

## Quick start

```bash
cd sample/amorphous/mixture_json
bash run_sample.sh

# After the build completes:
cd md
bash run_all.sh       # 5-stage annealing MD
python wrap_pbc.py    # produces *_pbc.xtc for VMD
```

## The config

```json
{
  "components": [
    {"smiles": "CCCCC", "name": "pentane", "n_mol": 200},
    {"smiles": "c1ccccc1", "name": "benzene", "n_mol": 50}
  ],
  "density_g_cm3": 0.8,
  "temperature": 300,
  "T_high": 600,
  "seed": 42,
  "forcefield": "openff_unconstrained-2.1.0.offxml",
  "output_dir": "."
}
```

`components` takes any number of entries. `T_high` is the melt
temperature the annealing protocol heats to before cooling to
`temperature`; `seed` fixes the packmol packing so a run is
reproducible.

## Notes

- **Linux / macOS only.** There is no `.bat`: packmol has no conda-forge
  win-64 build, so the box cannot be packed on Windows. The reasoning is
  in [the README one level up](../README.md), section
  "Windows で組めない理由".
