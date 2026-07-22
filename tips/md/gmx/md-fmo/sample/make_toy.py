#!/usr/bin/env python3
"""溶媒つき 3 フラグメント複合体の toy 系を生成する。

実系 (タンパク A + タンパク B + リガンド + 水) を数百原子で模す:
  - 溶質 = 3 フラグメント (A,B,L) が単一 moleculetype に同居・相互に無結合
    (Amber -> acpype -> GROMACS で溶質が 1 分子に潰れた状況)
  - 3 フラグメントは互いに接触した compact 複合体
  - 水 (WAT: O,H1,H2 の 3 点) が箱を満たす
  - 複合体は箱の隅 (周期境界をまたぐ位置) に置く
    -> MD 中にタンパクが箱端へ拡散した現実的な状況

箱は複合体に対してやや小さめ (tight box) にしてある。これにより
gmx trjconv -pbc cluster が「接触は戻すが水和殻を置き去りにする
(溶質が脱水する)」実系の失敗モードが再現される (run_demo.sh 参照)。

生成物:
  toy_good.gro : 正解 (複合体が接触・箱内にまとまった参照)
  toy.gro      : 生構造 (複合体を箱の隅へずらして境界でラップ; FMO 前処理の入力に相当)
  toy.top      : GROMACS トポロジ (COMPLEX 1 + WAT N)
"""
import numpy as np

L = 2.8            # nm, cubic box (複合体に対しやや小さめ = tight box)
BOND = 0.153
ND = 3             # タンパク blob の 1 辺あたり原子数 (ND^3 原子)
WSPACING = 0.40    # 水グリッド間隔 [nm]
WEXCL = 0.28       # 溶質からこの距離 [nm] 内には水を置かない
OH = 0.09572
HOH = np.deg2rad(104.52)


def cube(center, n):
    pts = np.array([[i * BOND, j * BOND, k * BOND]
                    for i in range(n) for j in range(n) for k in range(n)], float)
    return pts - pts.mean(0) + center


c = L / 2
half = (ND - 1) * BOND / 2
sep = 2 * half + 0.32                      # A,B の表面が接触する重心間距離
A = cube([c - sep / 2, c, c], ND)          # タンパク A
B = cube([c + sep / 2, c, c], ND)          # タンパク B
Lg = cube([c, c + half + 0.20, c], 2)      # リガンド (界面上部に接触)
solute = np.vstack([A, B, Lg])
sol_c = {'A': A, 'B': B, 'L': Lg}


def water_at(o):
    return np.array([o, o + [OH, 0, 0],
                     o + [OH * np.cos(HOH), OH * np.sin(HOH), 0]])


waters = []
grid = np.arange(WSPACING / 2, L, WSPACING)
for x in grid:
    for y in grid:
        for z in grid:
            o = np.array([x, y, z])
            if np.min(np.linalg.norm(solute - o, axis=1)) < WEXCL:
                continue
            waters.append(water_at(o))
waters = np.array(waters)
NW = len(waters)

FRAGS = [('FRA', 1, 'A'), ('FRB', 2, 'B'), ('LIG', 3, 'L')]


def wrap(v):
    return v - L * np.floor(v / L)


def write_gro(path, title, sol, wat, shift):
    n = len(solute) + 3 * len(wat)
    with open(path, 'w') as f:
        f.write(title + '\n')
        f.write(f'{n}\n')
        idx = 0
        for resname, rn, key in FRAGS:
            coords = wrap(sol[key] + shift) if shift is not None else sol[key]
            for xyz in coords:
                idx += 1
                f.write(f'{rn:5d}{resname:<5s}{"C":>5s}{idx % 100000:5d}'
                        f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}\n')
        rr = 3
        wl = (wrap(wat.reshape(-1, 3) + shift).reshape(wat.shape)
              if shift is not None else wat)
        for w in wl:
            rr += 1
            for aname, xyz in zip(('O', 'H1', 'H2'), w):
                idx += 1
                f.write(f'{rr % 100000:5d}{"WAT":<5s}{aname:>5s}{idx % 100000:5d}'
                        f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}\n')
        f.write(f'{L:10.5f}{L:10.5f}{L:10.5f}\n')


# 正解構造
write_gro('toy_good.gro', 'toy solvated 3-frag complex (contact, in box)',
          sol_c, waters, None)
# 生構造: 複合体の界面が境界にかかるよう隅へずらしてラップ
shift = np.array([-(c - sep / 2) + 0.15, -c + 0.15, 0.0])
write_gro('toy.gro', 'toy solvated 3-frag complex (wrapped across PBC)',
          sol_c, waters, shift)

# --- toy.top ---
NA, NB, NL = len(A), len(B), len(Lg)
with open('toy.top', 'w') as f:
    f.write('; toy solvated 3-fragment complex\n')
    f.write('; solute A,B,L share ONE moleculetype (no inter-fragment bond); water = WAT\n\n')
    f.write('[ defaults ]\n1 2 yes 0.5 0.8333\n\n')
    f.write('[ atomtypes ]\n')
    f.write(';name at.num mass    charge ptype sigma      epsilon\n')
    f.write(' C    6      12.011  0.0    A     3.4e-01    3.6e-01\n')
    f.write(' OW   8      16.000  0.0    A     3.15e-01   6.4e-01\n')
    f.write(' HW   1       1.008  0.0    A     0.0e+00    0.0e+00\n\n')
    f.write('[ moleculetype ]\n;name   nrexcl\nCOMPLEX  3\n\n[ atoms ]\n')
    idx = 0
    for resname, rn, _ in FRAGS:
        n = NA if rn == 1 else NB if rn == 2 else NL
        for _i in range(n):
            idx += 1
            f.write(f'{idx:6d}  C {rn:5d} {resname:>5s} {"C":>5s} '
                    f'{idx:5d} {0.0:8.4f} {12.011:8.3f}\n')
    f.write('\n[ bonds ]\n')
    off = 0
    for fr in (A, B, Lg):
        for i in range(off + 1, off + len(fr)):
            f.write(f'{i:6d}{i + 1:6d}  1  {BOND:8.4f}  200000.0\n')
        off += len(fr)
    f.write('\n[ moleculetype ]\n;name   nrexcl\nWAT  2\n\n[ atoms ]\n')
    f.write('     1  OW 1 WAT  O 1  0.0 16.000\n')
    f.write('     2  HW 1 WAT H1 1  0.0  1.008\n')
    f.write('     3  HW 1 WAT H2 1  0.0  1.008\n')
    f.write('\n[ bonds ]\n     1     2  1  0.09572  400000.0\n'
            '     1     3  1  0.09572  400000.0\n')
    f.write('\n[ angles ]\n     2     1     3  1  104.52  400.0\n')
    f.write('\n[ system ]\ntoy solvated complex\n\n[ molecules ]\n')
    f.write('COMPLEX 1\n')
    f.write(f'WAT {NW}\n')

print('wrote toy_good.gro / toy.gro / toy.top')
print(f'  solute {len(solute)} atoms (A,B,L single moleculetype) + WAT x{NW}'
      f' = {len(solute) + 3 * NW} atoms, box {L} nm')
