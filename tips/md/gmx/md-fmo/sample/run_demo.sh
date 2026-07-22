#!/bin/bash
#
# md-fmo PBC 再構成デモ (極小 toy 系)。
#
# 単一 moleculetype に同居した 3 フラグメント複合体 (タンパク A + B + リガンドを模す) で、
# B が周期境界を越えた「生構造」を各手法で組み直し、結果を check_contact.py で判定する。
#
# 期待される結果:
#   raw        : NG  (B が境界の向こうのイメージにある)
#   -pbc whole : NG  (A,B,L は単一 moleculetype・相互に無結合なので再結合できない)
#   -pbc cluster : OK (この小さな系では修復できる。ただし実系では初期配置依存で
#                      破綻することがある -- ../README.md の実データ表を参照)
#   autoimage  : OK  (prmtop の分子情報から決定論的に組み直す)
#
# 必要: gmx, cpptraj (AmberTools), python3 + numpy + scipy + parmed
#
set -e

PY=${PYTHON:-python3}
CHECK=../check_contact.py
FR="--frag A:1-1 --frag B:2-2 --frag L:3-3"

command -v gmx     >/dev/null || { echo "gmx が見つかりません"; exit 1; }
command -v cpptraj >/dev/null || { echo "cpptraj が見つかりません"; exit 1; }

# --- 0. toy 系を生成 (toy_good.gro, toy.gro, toy.top) ---
$PY make_toy.py

# --- 1. prmtop を生成 (0_parmed.sh の gromber 相当; cpptraj autoimage 用) ---
$PY - <<'PYEOF'
import parmed as pmd
g = pmd.load_file('toy.top', xyz='toy.gro')
g.save('toy.prmtop', format='amber', overwrite=True)
print('toy.prmtop written (%d molecules)' % len(list(g.split())))
PYEOF

# --- 2. 参照 tpr (gmx の -pbc 用トポロジ) ---
printf 'integrator=steep\nnsteps=0\nrvdw=0.9\nrlist=0.9\nrcoulomb=0.9\npbc=xyz\n' > em.mdp
gmx grompp -f em.mdp -c toy.gro -p toy.top -o toy_ref.tpr -maxwarn 2 >/dev/null 2>&1

run_all () {   # $1 = 入力 gro, $2 = タグ
    local inp=$1 tag=$2
    gmx trjconv -f "$inp" -s toy_ref.tpr -o ${tag}_raw.pdb           <<< "0" >/dev/null 2>&1
    gmx trjconv -f "$inp" -s toy_ref.tpr -pbc whole -o ${tag}_whole.pdb <<< "0" >/dev/null 2>&1
    printf "0\n0\n" | gmx trjconv -f ${tag}_whole.pdb -s toy_ref.tpr \
        -pbc cluster -o ${tag}_cluster.pdb >/dev/null 2>&1
    cat > ${tag}_ai.in <<EOF
parm toy.prmtop
trajin ${tag}_whole.pdb
autoimage anchor :1 fixed :2,3
trajout ${tag}_autoimage.pdb pdb
run
quit
EOF
    cpptraj -i ${tag}_ai.in >${tag}_ai.log 2>&1
}

# --- 3. 生構造 (B ラップ; FMO 前処理の入力に相当) を 4 手法で組み直し ---
run_all toy.gro W
echo
echo "###### 生構造 (B が周期境界の向こう) を組み直した結果 ######"
$PY $CHECK $FR W_raw.pdb W_whole.pdb W_cluster.pdb W_autoimage.pdb || true

# --- 4. 対照: 正解構造はどの手法でも壊れないこと ---
run_all toy_good.gro G
echo
echo "###### 対照: 正解構造 (最初から接触) ######"
$PY $CHECK $FR G_raw.pdb G_cluster.pdb G_autoimage.pdb || true

echo
echo "ポイント: -pbc whole は NG のまま (単一 moleculetype なので再結合不可)。"
echo "         autoimage は prmtop の分子情報で決定論的に修復する。"
