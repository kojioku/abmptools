#!/bin/bash

### setting ###
templatepdb='4j0p-1H8_md-samp_int10ps_xxxps-optmasked-4.0ang-noion.pdb'

sortflag=true
sortlig='1H8'
sortsolv='SOL'

method='MP2'
basis='6-31G*'

memory=5000
np=4

Ligchargeflag=true
ligand='1H8'
chgval=1

starttime=10000
endtime=20000
interval=200
################

namehead=${templatepdb%.pdb}

# step1 rename CYX to CYS
python -m abmptools.pdbmodify -mode rename -str 'CYX' 'CYS' -i *-noion.pdb
mkdir orig
mv *-noion.pdb orig

# step2 sort Res and trivial modify
if [ "${sortflag}" ]; then
    python -m abmptools.pdbmodify -i *renamed.pdb -s ${sortlig} ${sortsolv}
else
    python -m abmptools.pdbmodify -i *renamed.pdb
fi
mkdir renamed
mv *renamed.pdb renamed

# step3 generate ajftemplate
if [ "${Ligchargeflag}" ]; then
    python -m abmptools.generateajf -i ${namehead}-renamed-mod.pdb -cmm -mem ${memory} -np ${np} -lc ${ligand} ${chgval} -basis ${basis} --method ${method}
else
    python -m abmptools.generateajf -i ${namehead}-renamed-mod.pdb -cmm -mem ${memory} -np ${np} -basis ${basis} --method ${method}
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

echo "#!/bin/bash
########
script=''

stime=${starttime}
etime=${endtime}
interval=${interval}
zeropad=5
temp='${namehead}-renamed-mod-${method}-${basisstr}-nbo.ajf'
########

fhead=\${temp%xxx*}
ftail=\${temp##*xxx}

for i in \`seq \${stime} \${interval} \${etime}\`
do
    num=\`printf \"%0\${zeropad}d\" \${i}\`
    bash \${script} \${fhead}\${num}\${ftail}
done
" > ${namehead}-${method}-${basisstr}-fmoset/createbatchrun.sh

