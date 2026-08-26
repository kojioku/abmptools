#!/usr/bin/env python3
"""切り出した液滴 PDB が「複合体が接触した状態」で切り出せているかを検証する。

MD の後処理で周期境界 (PBC) を組み直す際、複合体が数十 A 離れた構造や
水和殻に穴が空いた構造が、見た目では気づきにくいまま生成されることがある。
FMO に投入する前にこのスクリプトで機械的にチェックする。

判定内容:
  1. 溶質部分どうしの実座標最小原子間距離 (--residues で指定した各ペア)
     すべて --contact 未満なら接触 OK。
  2. 水和殻の最大空隙 (溶媒露出原子から最近接水までの距離の最大値)
     --hole を超えると PBC イメージング失敗を疑う。

**ここで指定するのは残基番号であって、FMO のフラグメント番号ではない。**
FMO の分割はこの後の 3_fmosetup.sh (ajf 生成) で決まるので、この時点では
まだ存在しない。両者は溶質の範囲では一致することが多いが、一致する保証は
ない。実データ (EGFR) の AUTOMATIC FRAGMENTATION 表:

    Seq. Frag. Residue
     324   324    HYZ
     335   325    WAT      <- 残基通番 335 が FMO フラグメント 325

液滴に残らなかった水の分だけ番号がずれる。cpptraj の `:N-M` と同じ番号
(= prmtop の残基番号) を渡すこと。

依存: numpy, scipy

使い方:
    python3 check_contact.py \\
        --residues A:1-840 --residues B:841-1184 --residues lig:1185 \\
        <head>-optedpdb/*_fmo_mask.pdb

--residues を省略すると、溶質(非水)をひとまとまりとして扱い、
水和殻チェックのみ行う。
"""
import argparse
import sys

import numpy as np
from scipy.spatial import cKDTree

WATER_RESNAMES = ('WAT', 'HOH', 'SOL', 'TIP3', 'T3P')


def parse_group(spec):
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


def check(path, groups, contact, hole, expect_water=True):
    """液滴 1 つを判定する。(ok, 対の文字列, 水分子数, 空隙, 理由) を返す。

    「測れなかった」を合格にしないことが要点。空のグループも、水の無い液滴も、
    以前は静かに素通りしていた。原子 50 個・水 0 に切り詰めた液滴が
    `water= 0 max_gap= nan OK` と報告された実例がある。
    """
    xyz, resseq, resname, atname = read_pdb(path)
    problems = []

    if len(xyz) == 0:
        return False, '', 0, float('nan'), ['原子が 1 つも無い']

    iswat = np.isin(resname, WATER_RESNAMES)
    sol_mask = ~iswat
    sol = xyz[sol_mask]
    wo = xyz[iswat & (atname == 'O')]

    # --- 溶質部分どうしの距離 ---
    pair_ok = True
    pair_str = []
    if groups:
        seg = {}
        for name, r0, r1 in groups:
            m = sol_mask & (resseq >= r0) & (resseq <= r1)
            seg[name] = xyz[m]
            # 空のグループは「測れなかった」であって「問題なし」ではない。
            # レンジの指定間違い (FMO フラグメント番号を渡した等) がここに出る。
            if len(seg[name]) == 0:
                problems.append(f'{name} に該当する残基が無い')
        names = [g[0] for g in groups]
        npair = 0
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                if len(seg[a]) == 0 or len(seg[b]) == 0:
                    continue
                d = mindist(seg[a], seg[b])
                pair_str.append(f'{a}-{b}={d:6.2f}')
                npair += 1
                if d >= contact:
                    pair_ok = False
        if len(names) > 1 and npair == 0:
            problems.append('比較できた組が 1 つも無い')

    # --- 水和殻の穴 ---
    gap = np.nan
    if len(wo) and len(sol):
        tsol = cKDTree(sol)
        ncon = np.array([len(x) for x in tsol.query_ball_point(sol, 8.0)])
        exposed = sol[ncon < np.percentile(ncon, 25)]
        if len(exposed) == 0:          # 全原子の接触数が均一なとき (小さな溶質等)
            exposed = sol
        gap = float(cKDTree(wo).query(exposed, k=1)[0].max())
    elif expect_water:
        # 液滴なのに水が無いのは切り出しの失敗。nan を合格にしない。
        problems.append('水が 1 分子も無い')

    hole_ok = np.isnan(gap) or gap < hole
    return ((pair_ok and hole_ok and not problems),
            ' '.join(pair_str), len(wo), gap, problems)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pdb', nargs='+', help='検証する液滴 PDB')
    ap.add_argument('--residues', '--frag', action='append', default=[],
                    dest='residues', metavar='NAME:START-END',
                    help='溶質の部分を残基番号で指定 (複数可)。'
                         'FMO のフラグメント番号ではない')
    ap.add_argument('--contact', type=float, default=5.0,
                    help='接触とみなす最大距離 [A] (default: 5.0)')
    ap.add_argument('--hole', type=float, default=8.0,
                    help='水和殻の穴とみなす最大空隙 [A] (default: 8.0)')
    ap.add_argument('--allow-dry', action='store_true',
                    help='水を含まない構造も許す (液滴化前の検証など)')
    args = ap.parse_args()

    groups = [parse_group(s) for s in args.residues]
    nbad = 0
    for path in args.pdb:
        ok, pair_str, nwat, gap, problems = check(
            path, groups, args.contact, args.hole,
            expect_water=not args.allow_dry)
        name = path.split('/')[-1]
        note = ('  <- ' + ' / '.join(problems)) if problems else ''
        print(f'{name:45s} {pair_str:40s} '
              f'water={nwat:6d} max_gap={gap:6.2f}  '
              f'{"OK" if ok else "*** NG ***"}{note}')
        if not ok:
            nbad += 1
    print(f'\n{len(args.pdb) - nbad}/{len(args.pdb)} OK')
    sys.exit(1 if nbad else 0)


if __name__ == '__main__':
    main()
