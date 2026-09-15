# pentane / benzene mixture (from SMILES)

A two-component amorphous box built straight from **SMILES**, with no
input file to prepare. Useful as the smallest end-to-end check of the
builder, and as the template for any mixture whose components RDKit can
generate a conformer for.

200 pentane + 50 benzene at 0.8 g/cm^3, 300 K.

## Quick start

```bash
cd sample/amorphous/pentane_benzene
bash run_sample.sh

# After the build completes:
cd md
bash run_all.sh       # 5-stage annealing MD
python wrap_pbc.py    # produces *_pbc.xtc for VMD
```

Everything the sample needs is on the command line in `run_sample.sh`
(`--smiles "CCCCC" "c1ccccc1"`), so there is no `input/` directory to
ship. The build writes `./input`, `./build` and `./md` here.

## The same system, driven by a file

[`../mixture_json/`](../mixture_json/) builds this system from a JSON
config instead of CLI flags. Use that one when the component list grows
past what is comfortable on a command line, or when the composition
should live in version control next to the results.

## Notes

- **Linux / macOS only.** There is no `.bat`: packmol has no conda-forge
  win-64 build, so the box cannot be packed on Windows. The reasoning is
  in [the README one level up](../README.md), section
  "Windows で組めない理由".
- The build itself is quick because both components are small and
  `--charge_method` defaults to AM1-BCC on a handful of unique molecules,
  not on all 250.
