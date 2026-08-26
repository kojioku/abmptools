#!/bin/bash

### setting ###
templatepdb='egfr-HYZ_pr_xxx_fmo_mask.pdb'
method='MP2'
basis='6-31G*'

ndigits=1
memory=7000
np=1

fgkflag=true
runsh='setupbindsv1-l39-bulk.sh' #runabfmogw.sh

starttime=0
endtime=4
interval=1

sortflag=false
sortlig='1H8'
sortsolv='SOL'

Ligchargeflag=false
ligand='1H8'
chgval=1
################

namehead=${templatepdb%.pdb}

# ジョブ投入スクリプトを探して手元へ持ってくる。
#
# 以前はカレントに置いてあることを要求し、無ければ即座に落ちていた。
# 実際には 3_fmosetup.sh を流すのは -optedpdb の中で、投入スクリプトは
# 配布物の一式と一緒に上の階層にある。手で cp させる手順が 1 つ増えるだけで、
# 忘れると理由の分かりにくいエラーになる。
#
# 探す順: カレント -> このスクリプトと同じ場所 -> 上へ 3 階層。
# 見つかったらカレントへ複製する (元は動かさない)。
if [ ! -f "${runsh}" ]; then
  _here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  for _d in "${_here}" .. ../.. ../../..; do
    if [ -f "${_d}/${runsh}" ]; then
      cp "${_d}/${runsh}" .
      echo "[copy] ${_d}/${runsh} を持ってきた"
      break
    fi
  done
fi

# 名前が違うだけのことがある (setupv1dd2024-bulk.sh / setupbindsv1-*.sh 等)。
# 候補が 1 つに決まるなら、それを使って先へ進む。
if [ ! -f "${runsh}" ]; then
  _here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  _cand=$(ls -1 *bulk.sh runabfmogw.sh "${_here}"/*bulk.sh \
                "${_here}"/runabfmogw.sh ../*bulk.sh ../runabfmogw.sh \
          2>/dev/null | sort -u)
  if [ "$(printf '%s\n' "${_cand}" | grep -c .)" = 1 ] && [ -n "${_cand}" ]; then
    cp "${_cand}" .
    runsh=$(basename "${_cand}")
    echo "[copy] ${_cand} を使う (runsh を ${runsh} に読み替えた)"
  fi
fi

if [ ! -f "${runsh}" ]; then
  echo "Error: ${runsh} が見つからない。" >&2
  echo "       カレント・スクリプトと同じ場所・上位 3 階層を探した。" >&2
  echo "       投入スクリプトを置くか、冒頭の runsh= を実名に直すこと。" >&2
  exit 1
fi
echo "runsh = ${runsh}"

# step1 rename CYX to CYS
python -m abmptools.pdbmodify -mode rename -str 'CYX' 'CYS' -i *.pdb
mkdir orig
mv *mask.pdb orig

# step2 sort Res and trivial modify
if "${sortflag}"; then
    python -m abmptools.pdbmodify -i *renamed.pdb -s ${sortlig} ${sortsolv}
else
    python -m abmptools.pdbmodify -i *renamed.pdb
fi
mkdir renamed
mv *renamed.pdb renamed

# step3 generate ajftemplate
if "${Ligchargeflag}"; then
    python -m abmptools.generateajf -i ${namehead}-renamed-mod.pdb -cmm -mem ${memory} -np ${np} -lc ${ligand} ${chgval} -basis ${basis} --method ${method} -cpfv 23
else
    python -m abmptools.generateajf -i ${namehead}-renamed-mod.pdb -cmm -mem ${memory} -np ${np} -basis ${basis} --method ${method} -cpfv 23
fi

# step4 generate serial ajfs
if [ "${basis}" = '6-31G*' ]; then
    basisstr='6-31Gd'
else
    basisstr=${basis}
fi

python -m abmptools.ajfserial -i ${namehead}-renamed-mod-${method}-${basisstr}-nbo.ajf -t ${starttime} ${endtime} ${interval} -str xxx

mkdir ${namehead}-${method}-${basisstr}-fmoset
mv *pdb *ajf ${namehead}-${method}-${basisstr}-fmoset

# step5 generate createbachrun.sh
cd ${namehead}-${method}-${basisstr}-fmoset

if "$fgkflag"; then
    bash ../${runsh} ${namehead}-renamed-mod-${method}-${basisstr}-nbo.ajf

else
    echo "#!/bin/bash
    ########
    script=\$1

    stime=\$2
    etime=\$3
    interval=${interval}
    zeropad=${ndigits}
    temp='${namehead}-renamed-mod-${method}-${basisstr}-nbo.ajf'
    ########

    fhead=\${temp%xxx*}
    ftail=\${temp##*xxx}

    for i in \`seq \${stime} \${interval} \${etime}\`
    do
        num=\`printf \"%0\${zeropad}d\" \${i}\`
        bash \${script} \${fhead}\${num}\${ftail}
    done
    " > 4_batchrun.sh
    echo "4_batchrun.sh is generated"
    echo "Next: start calc e.g) bash 4_batchrun.sh ../runabfmogw.sh 0 100"
fi

