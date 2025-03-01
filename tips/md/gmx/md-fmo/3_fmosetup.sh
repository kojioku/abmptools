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

if [ ! -f ${runsh} ]; then
  echo "Error: ${runsh} not found in the current directory."
  exit 1
fi

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

