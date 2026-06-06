# `_postprocess/` — DEPRECATED (機能は `abmptools.trajectory` に移行済み)

旧 bash script (`trajectory_thin_nojump.sh`) は廃止予定。 同等処理は
Windows でも動く Python module **`abmptools.trajectory`** から呼ぶ。

## 移行先

```bash
# 環境
micromamba activate abmptoolsenv     # abmptools が editable install されていれば OK
# (gmx は abmptoolsenv に含まれないので gmxcudaenv 経由が安全:)
micromamba activate gmxcudaenv && pip install -e <abmptools repo>

# 使い方 (CLI、 旧 trajectory_thin_nojump.sh の置換):
python -m abmptools.trajectory thin_nojump \
    --traj prod/prod.xtc \
    --tpr  prod/prod.tpr \
    --skip 10
# → prod/prod_nojump_skip10.xtc (100 frame、 1 ns stride、 -pbc nojump 済)
```

```python
# 使い方 (Python API):
from abmptools.trajectory import thin_and_nojump
thin_and_nojump(
    trajectory="prod/prod.xtc",
    tpr="prod/prod.tpr",
    skip=10,
)
# → prod/prod_nojump_skip10.xtc
```

## 提供 subcommand

| Subcommand | gmx trjconv オプション | 用途 |
|---|---|---|
| `thin_nojump` (default) | `-pbc nojump -skip N` | aggregation 系の基本セット (= 旧 .sh) |
| `nojump` | `-pbc nojump` | OCTA / J-OCTA Viewer 用 (frame 数そのまま) |
| `thin` | `-skip N` | 純粋に間引きだけ |
| `wrap_pbc` | `-pbc mol -ur compact [-center]` | VMD 向け wrap (= amorphous 旧 `wrap_pbc.sh`) |

## なぜ Python 化したか

- **Windows でも動く** (旧 .sh は Linux/WSL 専用)
- subprocess の error handling と stdin group 選択を統一
- **重複 source 撤廃**: amorphous の `wrap_pbc.sh` / `gen_for_udf.sh`、
  formulation の `trajectory_thin_nojump.sh` が同じ gmx trjconv パターン
  を別実装で持っていた → `abmptools.trajectory` に集約

## 参考

- `abmptools/trajectory/postprocess.py` — 関数本体
- `abmptools/trajectory/cli.py` — argparse subcommand 定義
- `tests/test_trajectory_postprocess.py` — 14 unit tests
- `python -m abmptools.trajectory --help` — 全 subcommand 一覧
