#!/bin/sh
# WARNING! CURRENT DIRECTORY MUST BE THAT OF THIS FILE
# Date: 2020/3/24
# Author: Koji Okuwaki
# ver2.1: Calculations up to equilibration is applied to the MOE.
# ver2.0: fix for setting cutoff, density, temp, and mask
set -e

# -- user setting --
fbase=$1
np=128

# prepare ps
heatps=50
equilps=1000
anealps=100

# run-step
# prodps=100000
# interval=1000

# md param
cutval=12
pval=1.01
dtval=0.001
tempval=50
restmask="@N,CA,C,O3',C3',C4',C5',O5',P"
# -------------------------

mpi="mpirun -np $np"
ambermin="sander.MPI"
amber="sander.MPI"
coor=$fbase.a.0.coor
prmtop=$fbase.z.prmtop
cfg=$fbase.z.cfg

t0=.a.0
rstrt=rstrt
n=1

#$ -l f_node=5                 ## * [資源タイプ名]=[使用ノード数個数]
#$ -l h_rt=12:00:00             ## * 経過時間 -> Max 24時間までです。
#$ -p -5                       ## * ジョブの優先度 *デフォルトは-5です。

# カレントディレクトリに移動
#$ -cwd

# Moduleコマンドの初期化
. /etc/profile.d/modules.sh

# CUDA 環境の読み込み
module load amber/16up10_cuda

# Intel Compiler環境の読み込み
module load intel/18.0.1.163
# module load intel/19.0.0.117

# Intel MPI 環境の読み込み * コンパイル時の設定環境を指定
module load intel-mpi/18.1.163
# module load intel-mpi/19.0.117


# AMBER configuration.
# Set the amber executable in $amber
# - subshells will _MOE_AMBER_EXE to override $amber
# - Default to $AMBERHOME which we assume has been set up properly.
# - CUDA_VISIBLE_DEVICES is required for GPU usage.
# - insert mpirun -np # if needed for mpi use


function gethlab(){
if [ $tmpst -lt 10 ]; then
    hlab='a'
elif [ $tmpst -lt 100 ]; then
    hlab='b'
elif [ $tmpst -lt 1000 ]; then
    hlab='c'
elif [ $tmpst -lt 10000 ]; then
    hlab='d'
elif [ $tmpst -lt 100000 ]; then
    hlab='e'
elif [ $tmpst -lt 1000000 ]; then
    hlab='f'
elif [ $tmpst -lt 10000000 ]; then
    hlab='g'
elif [ $tmpst -lt 100000000 ]; then
    hlab='h'
fi
}

function gethlab_e(){
if [ $tmpet -lt 10 ]; then
    hlab_e='a'
elif [ $tmpet -lt 100 ]; then
    hlab_e='b'
elif [ $tmpet -lt 1000 ]; then
    hlab_e='c'
elif [ $tmpet -lt 10000 ]; then
    hlab_e='d'
elif [ $tmpet -lt 100000 ]; then
    hlab_e='e'
elif [ $tmpet -lt 1000000 ]; then
    hlab_e='f'
elif [ $tmpet -lt 10000000 ]; then
    hlab_e='g'
elif [ $tmpet -lt 100000000 ]; then
    hlab_e='h'
fi
}

function settimelab(){
tmpst=$tmpet
echo tmpst = $tmpst
gethlab
tmpet=$(($tmpst + $1))
echo tmpet = $tmpet
gethlab_e
}

function gencfg(){

echo "# Name;Protocol;Init;Stage" > $cfg
# heat

tmpet=0
settimelab $heatps
echo "heat;z.1.in;${hlab}.${tmpst}.min2.rstrt;${hlab_e}.$tmpet" >> $cfg
# den
settimelab $equilps
echo "equil;z.2.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg

settimelab $anealps
echo "aneal;z.3.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg

# prod
# echo "prod;z.1.in;Equilibration.restrt;c.$interval" >> $cfg
# num=8
# # prodstart=$(($interval * 2))
# for i in `seq $interval $interval $prodps`
# do
#     gen_prodin
#     num=$(($num+1))
#     settimelab $interval
# 
#     echo "prod;z.${num}.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg
# done
}

function run_cp_md(){
\cp -f $3 inpcrd
$mpi $2 -O -ref inpcrd -r restrt_tmp
\cp -f restrt_tmp restrt
ambpdb -c restrt > pdb
mv -f mdin ${1}.mdin
mv -f mdout ${1}.mdout
mv -f restrt ${1}.restrt
mv -f mdcrd ${1}.mdcrd
mv -f pdb ${1}.pdb
mv -f restrt_tmp restrt
}

function runmin1() {
echo FALSE > FINISH
echo start min1
cat << EOF > mdin
Minimization1
&cntrl
  imin=1, maxcyc=100000, ncyc=3000, drms=0.1,
  ntr=1, restraint_wt=3.0, restraintmask='!@H=',

ntb=0,
cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000, ntwf=1000,
ig=0,
vlimit=-1,
/

EOF

run_cp_md Minimization1 $ambermin $coor
}


function runmin2(){
echo start min2

cat << EOF > mdin
Minimization2
&cntrl
  imin=1, maxcyc=200000, ncyc=3000,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

ntb=0,
cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000, ntwf=1000,
ig=0,
vlimit=-1,

/

EOF

run_cp_md Minimization2 $ambermin restrt
}


function runminlast(){
echo start min_last

cat << EOF > mdin
Minimization3
&cntrl
  imin=1, maxcyc=3000, ncyc=1500,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

ntb=0,
cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000, ntwf=1000,
ig=0,
vlimit=-1,

/

EOF

run_cp_md Minimize_${2} $ambermin ${1}
}


# ntwf option was removed for all md mdin because of error (GPU)

# heat (default 100ps)
function getin_heat(){
echo set Heat

cat << EOF > $1
Heating
&cntrl
  nstlim=`echo "$heatps/$dtval" | bc`,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

ntb=0,
cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=0,
vlimit=-1,


nmropt=1,
/
&wt type='TEMP0', istep1=0, istep2=`echo "$heatps/$dtval" | bc`, value1=0., value2=${tempval}, /
&wt type='END', /

EOF

# run_cp_md Heating $amber restrt
}

function getin_aneal(){
echo set aneal

cat << EOF > $1
aneal
&cntrl
  nstlim=`echo "$heatps/$dtval" | bc`,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

ntb=0,
cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=10,
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=0,
vlimit=-1,


/
EOF

# run_cp_md Heating $amber restrt
}


function getin_equil(){
echo set equil

cat << EOF > $1
Equilibration
&cntrl
  irest=1, ntx=5,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",
  nstlim=`echo "$equilps/$dtval" | bc`,
  ntb=0,
cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=0,
vlimit=-1,

/

EOF

# run_cp_md Equilibration $amber restrt
}


# -- main --
# -- setup --
\cp $prmtop prmtop
gencfg
# -- minimization --

getin_heat $fbase.z.1.in
getin_equil $fbase.z.2.in
getin_aneal $fbase.z.3.in


# -----------------------------------------
# production run (copy from MOE)
# -----------------------------------

protocol=${fbase}.z.cfg
extensions="${t0}.coor .z.prmtop .z.cfg"
# cp Equilibration.restrt $fbase.Equilibration.restrt

if [ ! -e Minimization1.mdout ]; then
    runmin1
else
    echo "already min1"
fi
if [ ! -e Minimization2.mdout ]; then
    runmin2
else
    echo " already min2"
fi

cp Minimization2.restrt $fbase.a.0.min2.rstrt
# cp  $fbase.a.0.coor $fbase.a.0.min2.rstrt


while read ext; do
    case $ext in
    [!#]*)
    ext=${ext#*;}
    ;;
    esac
done < $protocol

stdoutfiles=""

while read line; do
    case $line in
    [#]*) continue ;;
    esac
    name=${line%%;*}
    line=${line#*;}
    a=${line%%;*}
    b=${line%;*}
    b=${b##*;}
    c=${line##*;}

    if [ -f ${fbase}.${c}.${rstrt} ] && [ -f ${fbase}.${c}.valid ]; then
        echo Stage $n Complete: $name
    else
        rm -f ${fbase}.${c}.*
        echo Starting Stage $n: $name
            $mpi $amber -O -i ${fbase}.${a} -o ${fbase}.${c}.stdout \
            -p ${fbase}.z.prmtop -c ${fbase}.${b} \
            -ref ${fbase}.${b} -r ${fbase}.${c}.${rstrt} \
            -inf ${fbase}.${c}.mdinfo -x ${fbase}.${c}.mdcrd  \
            < /dev/null
        if [ -f ${fbase}.${c}.${rstrt} ]; then
            touch ${fbase}.${c}.valid
        else
            echo Failure Stage $n: $name
            exit 1
        fi
    fi
    stdoutfiles="${stdoutfiles} ${fbase}.${c}.stdout"
    n=$((n+1))
done < $protocol

if [ -e ${fbase}.${c}.stdout ]; then
    runminlast ${fbase}.${c}.rstrt ${fbase}.${c}
fi

# Scan all of of the .stdout files and extract the simulation
# measurements and store in sim.z.out; the format is a titled csv

awkscript=$(cat <<"EOT"
    BEGIN {
        _ti = 0; _avg = 0; _res = 0;
        stage = "?"; n = -1; Eele = 0; dVdL = 0;
        printf "stage,t,T,P,V,U,K,Estr,Eang,Edih,Evdw,Eele,Eres,dVdL\n"
    }
    /[|] TI region  1/   { _ti = 1 }
    /[|] TI region  2/   { _ti = 2 }
    /RESULTS/            { _res = 1; _avg = 0 }
    /A V E R A G E S/    { _avg = 1 }
    /^ Etot   =/         { U = $3; K = $6 }
    /^ EKCMT  =/         { V = $9 }
    /^ BOND   =/         { Estr = $3; Eang = $6; Edih = $9 }
    /^ 1-4 NB =/         { Evdw = $3 + $9; Eele = Eele + $6 }
    /^ EELEC  =/         { Eele = Eele + $3; Eres = $9 }
    /^ DV\/DL  =/        { dVdL = $3 }
    /^#MOE PROTOCOL:/    { stage = $3 }
    /^ NSTEP =/ {
        if (n >= 0) {
        printf "%s,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
            stage,t,T,P,V,U,K,Estr,Eang,Edih,Evdw,Eele,Eres,dVdL;
        n = -1;
        } else if (_res && ! _avg) {
        n = $3; t = $6; T = $9; P = $12 * 100
        Eele = 0;
        }
    }
EOT
); awk "${awkscript}" ${stdoutfiles} > ${fbase}.z.out

exit 0;

