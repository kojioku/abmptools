#!/usr/bin/env python3
"""切り出した液滴 PDB が「複合体が接触した状態」で切り出せているかを検証する。

MD の後処理で周期境界 (PBC) を組み直す際、複合体が数十 A 離れた構造や
水和殻に穴が空いた構造が、見た目では気づきにくいまま生成されることがある。
FMO に投入する前にこのスクリプトで機械的にチェックする。

判定内容:
  1. 溶質フラグメント間の実座標最小原子間距離 (--frag で指定した各ペア)
     すべて --contact 未満なら接触 OK。
  2. 水和殻の最大空隙 (溶媒露出原子から最近接水までの距離の最大値)
     --hole を超えると PBC イメージング失敗を疑う。

依存: numpy, scipy

使い方 (残基番号でフラグメントを指定; cpptraj の :N-M と同じ番号):
    python3 check_contact.py \\
        --frag A:1-840 --frag B:841-1184 --frag lig:1185 \\
        CR8_AF3_MD_pr-optedpdb/*_fmo_mask.pdb

--frag を省略すると、溶質(非水)を単一フラグメント扱いにし、
水和殻チェックのみ行う。
"""
import argparse
import sys

import numpy as np
from scipy.spatial import cKDTree

WATER_RESNAMES = ('WAT', 'HOH', 'SOL', 'TIP3', 'T3P')


def parse_frag(spec):
    """'A:1-840' -> ('A', 1, 840)  /  'lig:1185' -> ('lig', 1185, 1185)"""
    name, rng = spec.split(':')
    if '-' in rng:
        a, b = rng.split('-')
    else:
        a = b = rng
    return name, int(a), int(b)


def read_pdb(path):
    """PDB を (xyz, resseq, resname, atomname) で返す。"""
    xyz, resseq, resname, atname = [], [], [], []
    with open(path) as f:
        for l in f:
            if l.startswith(('ATOM', 'HETATM')):
                xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                resseq.append(int(l[22:26]))
                resname.append(l[17:20].strip())
                atname.append(l[12:16].strip())
    return (np.array(xyz), np.array(resseq),
            np.array(resname), np.array(atname))


def mindist(a, b):
    return cKDTree(b).query(a, k=1)[0].min()


def check(path, frags, contact, hole):
    xyz, resseq, resname, atname = read_pdb(path)
    iswat = np.isin(resname, WATER_RESNAMES)
    sol_mask = ~iswat
    sol = xyz[sol_mask]
    wo = xyz[iswat & (atname == 'O')]

    # フラグメント間距離
    pair_ok = True
    pair_str = []
    if frags:
        seg = {}
        for name, r0, r1 in frags:
            m = sol_mask & (resseq >= r0) & (resseq <= r1)
            seg[name] = xyz[m]
        names = [f[0] for f in frags]
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                if len(seg[a]) == 0 or len(seg[b]) == 0:
                    continue
                d = mindist(seg[a], seg[b])
                pair_str.append(f'{a}-{b}={d:6.2f}')
                if d >= contact:
                    pair_ok = False

    # 水和殻の穴
    gap = np.nan
    if len(wo) and len(sol):
        tsol = cKDTree(sol)
        ncon = np.array([len(x) for x in tsol.query_ball_point(sol, 8.0)])
        exposed = sol[ncon < np.percentile(ncon, 25)]
        if len(exposed) == 0:          # 全原子の接触数が均一なとき (小さな溶質等)
            exposed = sol
        gap = float(cKDTree(wo).query(exposed, k=1)[0].max())
    hole_ok = np.isnan(gap) or gap < hole

    ok = pair_ok and hole_ok
    name = path.split('/')[-1]
    print(f'{name:45s} {" ".join(pair_str):40s} '
          f'water={len(wo):6d} max_gap={gap:6.2f}  '
          f'{"OK" if ok else "*** NG ***"}')
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pdb', nargs='+', help='検証する液滴 PDB')
    ap.add_argument('--frag', action='append', default=[],
                    metavar='NAME:START-END',
                    help='溶質フラグメントを残基番号で指定 (複数可)')
    ap.add_argument('--contact', type=float, default=5.0,
                    help='接触とみなす最大距離 [A] (default: 5.0)')
    ap.add_argument('--hole', type=float, default=8.0,
                    help='水和殻の穴とみなす最大空隙 [A] (default: 8.0)')
    args = ap.parse_args()

    frags = [parse_frag(s) for s in args.frag]
    results = [check(p, frags, args.contact, args.hole) for p in args.pdb]
    nbad = results.count(False)
    print(f'\n{len(results) - nbad}/{len(results)} OK')
    sys.exit(1 if nbad else 0)


if __name__ == '__main__':
    main()
