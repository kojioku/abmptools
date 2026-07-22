#!/usr/bin/env python3
"""3 フラグメント複合体の極小 toy 系を生成する。

実系 (タンパク A + タンパク B + リガンド) の再現ポイント:
  - 3 フラグメントが「1 つの GROMACS moleculetype」に同居し、
    フラグメント間には結合が無い
    (Amber -> acpype -> GROMACS で溶質が 1 分子に潰れた状況を模す)。
  - 3 フラグメントは互いに接触している (compact な複合体)。
  - 生構造では B が周期境界の向こう側のイメージに置かれている
    (MD 中に分子が境界を越えて拡散した状況を模す)。

生成物:
  toy_good.gro  ... 正解構造 (3 フラグメント接触・箱内)
  toy.gro       ... 生構造 (B を境界の向こうへラップ; FMO 前処理の入力に相当)
  toy.top       ... 最小 GROMACS トポロジ (A,B,L が単一 moleculetype)
"""
import numpy as np

L = 5.0          # nm, cubic box
BOND = 0.153     # nm
NPF = 6          # atoms per fragment


def blob(center):
    """center 周りの小さな塊 (直鎖 6 原子を折り返して compact に)。"""
    c = np.array(center)
    offs = np.array([
        [0.0, 0.0, 0.0], [BOND, 0.0, 0.0], [2 * BOND, 0.0, 0.0],
        [2 * BOND, BOND, 0.0], [BOND, BOND, 0.0], [0.0, BOND, 0.0],
    ])
    return c + offs


# 3 フラグメントの重心を三角形に配置し、互いに接触させる
cen = L / 2
A = blob([cen - 0.30, cen - 0.20, cen])       # タンパク A
B = blob([cen + 0.30, cen - 0.20, cen])       # タンパク B
Lg = blob([cen - 0.02, cen + 0.25, cen])      # リガンド

# 生構造: B を -x 方向に 1 箱ぶんラップ
B_wrap = B.copy()
B_wrap[:, 0] -= L

FRAGS = [('FRA', 1, 'A'), ('FRB', 2, 'B'), ('LIG', 3, 'L')]


def write_gro(path, title, a, b, lg):
    coords = {'A': a, 'B': b, 'L': lg}
    n = sum(len(v) for v in coords.values())
    with open(path, 'w') as f:
        f.write(title + '\n')
        f.write(f'{n}\n')
        idx = 0
        for resname, resnr, key in FRAGS:
            for xyz in coords[key]:
                idx += 1
                f.write(f'{resnr:5d}{resname:<5s}{"C":>5s}{idx:5d}'
                        f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}\n')
        f.write(f'{L:10.5f}{L:10.5f}{L:10.5f}\n')


write_gro('toy_good.gro', 'toy 3-fragment complex (contact, in box)', A, B, Lg)
write_gro('toy.gro', 'toy 3-fragment complex (B wrapped across PBC)', A, B_wrap, Lg)

# --- toy.top : A,B,L を単一 moleculetype に同居、フラグメント間結合なし ---
bonds = []
for fi in range(3):
    off = fi * NPF
    for i in range(off + 1, off + NPF):     # 各フラグメント内は直鎖結合
        bonds.append((i, i + 1))

with open('toy.top', 'w') as f:
    f.write('; toy 3-fragment complex in a SINGLE moleculetype\n')
    f.write('; (mirrors an Amber->acpype->GROMACS solute collapsed into one molecule)\n\n')
    f.write('[ defaults ]\n1 2 yes 0.5 0.8333\n\n')
    f.write('[ atomtypes ]\n')
    f.write(';name at.num mass    charge ptype sigma      epsilon\n')
    f.write(' C    6      12.011  0.0    A     3.4e-01    3.6e-01\n\n')
    f.write('[ moleculetype ]\n;name   nrexcl\nCOMPLEX  3\n\n')
    f.write('[ atoms ]\n;  nr type resnr residue atom cgnr charge  mass\n')
    idx = 0
    for resname, resnr, _ in FRAGS:
        for _i in range(NPF):
            idx += 1
            f.write(f'{idx:6d}  C {resnr:5d} {resname:>5s} {"C":>5s} '
                    f'{idx:5d} {0.0:8.4f} {12.011:8.3f}\n')
    f.write('\n[ bonds ]\n')
    for a, b in bonds:
        f.write(f'{a:6d}{b:6d}  1  {BOND:8.4f}  200000.0\n')
    f.write('\n[ system ]\ntoy complex\n\n[ molecules ]\nCOMPLEX 1\n')

print(f'wrote toy_good.gro / toy.gro / toy.top : {3 * NPF} atoms, box {L} nm')
print('  A=res1(1-6)  B=res2(7-12)  L=res3(13-18)  single moleculetype, no inter-fragment bond')
