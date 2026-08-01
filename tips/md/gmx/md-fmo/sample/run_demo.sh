#!/bin/bash
#
# md-fmo PBC 再構成デモ (溶媒つき toy 系)。
#
# 単一 moleculetype に同居した 3 フラグメント複合体 (タンパク A + B + リガンド) を
# 水で溶媒和し、B 等が周期境界を越えた「生構造」を各手法で組み直して
# check_contact.py で判定する。実系の 2_optmask-frame.sh と同じ 4 段
# (whole -> cluster -> mol/compact) を再現し、cpptraj autoimage と比較する。
#
# 期待される結果 (接触距離 [A] と 水和殻の最大空隙 max_gap [A]):
#   good input   : OK           正解構造 (接触・水和済み)
#   -pbc whole   : NG           複合体が四散 (単一 moleculetype で再結合できない)
#   -pbc cluster : NG (接触は戻るが max_gap 大) <-- 溶質が脱水! 水和殻を置き去りにする
#   -pbc mol/compact : OK       後段でようやく水を戻す (全段が成功すればの話)
#   autoimage    : OK           接触も水和殻も 1 手で決定論的に回復
#
# 要点:
#   - -pbc whole/mol だけでは複合体を組み直せない。
#   - -pbc cluster は接触を戻せても水和殻に穴を開ける (溶媒系特有の失敗)。
#     実系ではさらに接触自体を壊すフレームがある (../README.md の 11/20)。
#   - autoimage は prmtop の分子情報で 1 手・決定論的に組み直す。
#
# 必要: gmx, cpptraj (AmberTools), python3 + numpy + scipy + parmed
#
set -e

PY=${PYTHON:-python3}
CHECK=../check_contact.py
FR="--frag A:1-1 --frag B:2-2 --frag L:3-3"

command -v gmx     >/dev/null || { echo "gmx が見つかりません"; exit 1; }
command -v cpptraj >/dev/null || { echo "cpptraj が見つかりません"; exit 1; }

NSOL=56   # 溶質原子数 (make_toy.py の 3^3*2 + 2^3 = 62 ... ND/リガンドを変えたら要調整)

# --- 0. toy 系を生成 ---
$PY make_toy.py
NSOL=$($PY - <<'EOF'
# toy.gro から溶質(非WAT)原子数を数える
n=0
for l in open('toy.gro').read().splitlines()[2:-1]:
    if l[5:10].strip()!='WAT': n+=1
print(n)
EOF
)
echo "solute atoms = $NSOL"

# --- 1. prmtop (0_parmed.sh 相当; cpptraj autoimage 用) ---
$PY - <<'PYEOF'
import parmed as pmd
g = pmd.load_file('toy.top', xyz='toy.gro')
g.save('toy.prmtop', format='amber', overwrite=True)
PYEOF

# --- 2. 参照 tpr と solute index ---
printf 'integrator=steep\nnsteps=0\nrvdw=0.9\nrlist=0.9\nrcoulomb=0.9\npbc=xyz\n' > em.mdp
gmx grompp -f em.mdp -c toy.gro -p toy.top -o toy_ref.tpr -maxwarn 5 >/dev/null 2>&1
printf "keep 0\na 1-${NSOL}\nname 1 solute\nq\n" | gmx make_ndx -f toy.gro -o toy.ndx >/dev/null 2>&1

topdb () { gmx trjconv -f "$1" -s toy_ref.tpr -o "$2" <<< "System" >/dev/null 2>&1; }

# --- 3. 実パイプライン (whole -> cluster -> mol/compact; .gro 中間) ---
gmx trjconv -f toy.gro -s toy_ref.tpr -pbc whole -o c1.gro <<< "System" >/dev/null 2>&1
printf "solute\nSystem\n" | gmx trjconv -f c1.gro -s toy_ref.tpr -n toy.ndx \
    -pbc cluster -o c2.gro >/dev/null 2>&1
printf "solute\nSystem\n" | gmx trjconv -f c2.gro -s toy_ref.tpr -n toy.ndx \
    -pbc mol -ur compact -center -o c3.gro >/dev/null 2>&1
topdb c1.gro c1.pdb ; topdb c2.gro c2.pdb ; topdb c3.gro c3.pdb
topdb toy_good.gro good.pdb

# --- 4. autoimage パイプライン (whole -> autoimage) ---
gmx trjconv -f toy.gro -s toy_ref.tpr -pbc whole -o aiw.gro <<< "System" >/dev/null 2>&1
topdb aiw.gro aiw.pdb
cat > ai.in <<EOF
parm toy.prmtop
trajin aiw.pdb
autoimage anchor :1 fixed :2,3 mobile :WAT
trajout ai.pdb pdb
run
quit
EOF
cpptraj -i ai.in >ai.log 2>&1

echo
echo "###### 生構造 (複合体が周期境界をまたぐ) を組み直した結果 ######"
printf "%-14s" "good input :"; $PY $CHECK $FR good.pdb 2>/dev/null | head -1 | sed 's#good.pdb##'
printf "%-14s" "whole      :"; $PY $CHECK $FR c1.pdb   2>/dev/null | head -1 | sed 's#c1.pdb##'
printf "%-14s" "cluster    :"; $PY $CHECK $FR c2.pdb   2>/dev/null | head -1 | sed 's#c2.pdb##'
printf "%-14s" "mol/compact:"; $PY $CHECK $FR c3.pdb   2>/dev/null | head -1 | sed 's#c3.pdb##'
printf "%-14s" "autoimage  :"; $PY $CHECK $FR ai.pdb   2>/dev/null | head -1 | sed 's#ai.pdb##'
echo
echo "ポイント: whole では複合体を再結合できず (NG)、cluster は接触を戻しても"
echo "         水和殻に穴を開ける (max_gap 大)。autoimage は 1 手で接触も水和も回復する。"
