#!/bin/sh
# Date: 2020/4/30
# Author: Koji Okuwaki
# ver2.6: add anneal and 3step min
# ver2.4: fix .in file number
# ver2.3: fix prodction run ps bug
# ver2.2: modify heat interval
# ver2.1: Calculations up to equilibration is applied to the MOE.
# ver2.0: fix for setting cutoff, density, temp, and mask
set -e
set -x

# -- user setting --
fbase=$1
np=4 #np for md: default=4
npmin=56   # np for minimization: default=num core

# prepare ps
heatps=50
den1ps=5
den2ps=5
den3ps=10
den4ps=10
den5ps=10
den6ps=10
equilps=1000

# run-step
prodps=130000
interval=1000

# md param
cutval=12.0
pval=1.01
dtval=0.001
tempval=310
restmask="@N,CA,C,O3',C3',C4',C5',O5',P"
# -------------------------

mpi="mpirun -np $np"
mpimin="mpirun -np $npmin"
ambermin="pmemd.MPI"
amber="pmemd.cuda.MPI"
coor=$fbase.a.0.coor
prmtop=$fbase.z.prmtop
cfg=$fbase.z.cfg

t0=.a.0
rstrt=rstrt
n=1

#$ -l f_node=2                 ## * [資源タイプ名]=[使用ノード数個数]
#$ -l h_rt=24:00:00             ## * 経過時間 -> Max 24時間までです。
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
echo "heat;z.1.in;${hlab}.${tmpst}.min3.rstrt;${hlab_e}.$tmpet" >> $cfg
# den
settimelab $den1ps
echo "heat;z.2.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg
settimelab $den2ps
echo "heat;z.3.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg
settimelab $den3ps
echo "den;z.4.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg
settimelab $den4ps
echo "den;z.5.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg
settimelab $den5ps
echo "den;z.6.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg
settimelab $den6ps
echo "den;z.7.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg
settimelab $equilps
echo "equil;z.8.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg

# prod
# echo "prod;z.1.in;Equilibration.restrt;c.$interval" >> $cfg
num=9
# prodstart=$(($interval * 2))
for i in `seq $interval $interval $prodps`
do
    gen_prodin
    settimelab $interval
    echo "prod;z.${num}.in;${hlab}.${tmpst}.rstrt;${hlab_e}.$tmpet" >> $cfg
    num=$(($num+1))
done
}

function gen_prodin(){

cat << EOF > ${fbase}.z.${num}.in
Production
&cntrl
  irest=1, ntx=5,
  nstlim=`echo "$interval/$dtval" | bc`,
  ntp=1,

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF
}


function run_cp_md(){
# run_cp_md name(fbase,out), bin, coor(in)
\cp -f $3 inpcrd
$mpimin $2 -O -ref inpcrd -r restrt_tmp
\cp -f restrt_tmp restrt
ambpdb -c restrt > pdb
mv -f mdin ${1}.mdin
mv -f mdout ${1}.mdout
mv -f restrt ${1}.restrt
if [ -e mdcrd ]; then
    mv -f mdcrd ${1}.mdcrd
fi
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
  ntr=1, restraint_wt=100, restraintmask='!@H=',

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000, ntwf=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

run_cp_md Minimization1 $ambermin $coor
}


function runmin2(){
echo start min2

cat << EOF > mdin
Minimization2
&cntrl
  imin=1, maxcyc=3000, ncyc=1500,
  ntr=1, restraint_wt=100, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000, ntwf=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

run_cp_md Minimization2 $ambermin restrt
}


function runmin3(){
echo start min3

cat << EOF > mdin
Minimization3
&cntrl
  imin=1, maxcyc=6000, ncyc=3000,

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000, ntwf=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

run_cp_md Minimization3 $ambermin restrt
}



# ntwf option was removed for all md mdin because of error (GPU)

# heat (default 100ps)
function getin_heat(){
echo start Heat

cat << EOF > $1
Heating
&cntrl
  nstlim=`echo "$heatps/$dtval" | bc`,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,


nmropt=1,
/
&wt type='TEMP0', istep1=0, istep2=`echo "$heatps/$dtval" | bc`, value1=0., value2=${tempval}, /
&wt type='END', /

EOF

# run_cp_md Heating $amber restrt
}


function getin_density1(){
echo start density1

cat << EOF > $1
Density1
&cntrl
  irest=1, ntx=5,
  nstlim=`echo "$den1ps/$dtval" | bc`,
  ntp=1,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

# run_cp_md Density1 $amber restrt
}


function getin_density2(){
echo start density2

cat << EOF > $1
Density2
&cntrl
  irest=1, ntx=5,
  nstlim=`echo "$den2ps/$dtval" | bc`,
  ntp=1,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

# run_cp_md Density2 $ambermin restrt
}


function getin_density3(){
echo start density3

cat << EOF > $1
Density3
&cntrl
  irest=1, ntx=5,
  nstlim=`echo "$den3ps/$dtval" | bc`,
  ntp=1,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

# run_cp_md Density3 $amber restrt
}


function getin_density4(){
echo start density4

cat << EOF > $1
Density4
&cntrl
  irest=1, ntx=5,
  nstlim=`echo "$den4ps/$dtval" | bc`,
  ntp=1,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

# run_cp_md Density4 $amber restrt
}


function getin_density5(){
echo start density5

cat << EOF > $1
Density5
&cntrl
  irest=1, ntx=5,
  nstlim=`echo "$den5ps/$dtval" | bc`,
  ntp=1,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

# run_cp_md Density5 $amber restrt
}


function getin_density6(){
echo start density6

cat << EOF > $1
Density6
&cntrl
  irest=1, ntx=5,
  nstlim=`echo "$den6ps/$dtval" | bc`,
  ntp=1,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

# run_cp_md Density6 $amber restrt
}


function getin_equil(){
echo start equil

cat << EOF > $1
Equilibration
&cntrl
  irest=1, ntx=5,
  nstlim=`echo "$equilps/$dtval" | bc`,
  ntp=1,
cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

# run_cp_md Equilibration $amber restrt
}

# run_cp_md name(fbase,out), bin, coor(in)
function runanneal(){
infile=$1
outfile=$2
echo start anneal
cat << EOF > mdin
Heating
&cntrl
  nstlim=`echo "$heatps/$dtval" | bc`,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

nmropt=1,
/
&wt type='TEMP0', istep1=0, istep2=`echo "$heatps/$dtval" | bc`, value1=${tempval}, value2=0., /
&wt type='END', /

EOF

run_cp_md $outfile $ambermin $infile
}


function runminlast1() {
infile=$1
outfile=$2
echo "step1 Hopt"
cat << EOF > mdin
Minimization1
&cntrl
  imin=1, maxcyc=100000, ncyc=3000, drms=0.1,
  ntr=1, restraint_wt=10000, restraintmask='!@H=',
  ig=-1,
  vlimit=-1,
/
EOF

run_cp_md $outfile $ambermin $infile
}


function runminlast2() {
infile=$1
outfile=$2
echo "step2(backbone:constrain SideChain:restraint)"
cat << EOF > mdin
Minimization2
&cntrl
  imin=1, maxcyc=100000, ncyc=3000, drms=0.3,
  !ntr=1, restraint_wt=10000.0, restraintmask="${restmask}",
  ig=-1,
  vlimit=-1,
/
EOF

run_cp_md $outfile $ambermin $infile
}


function runminlast3() {
infile=$1
outfile=$2
echo "step3 Hopt"
cat << EOF > mdin
Minimization3
&cntrl
  imin=1, maxcyc=100000, ncyc=3000, drms=0.1,
  ntr=1, restraint_wt=10000, restraintmask='!@H=',
  ig=-1,
  vlimit=-1,
/
EOF

run_cp_md $outfile $ambermin $infile
}


# -- main --
# -- setup --
\cp $prmtop prmtop
gencfg
# -- minimization --

getin_heat $fbase.z.1.in
getin_density1 $fbase.z.2.in
getin_density2 $fbase.z.3.in
getin_density3 $fbase.z.4.in
getin_density4 $fbase.z.5.in
getin_density5 $fbase.z.6.in
getin_density6 $fbase.z.7.in
getin_equil $fbase.z.8.in

# -----------------------------------------
# production run
# -----------------------------------

protocol=${fbase}.z.cfg
extensions="${t0}.coor .z.prmtop .z.cfg"
# cp Equilibration.restrt $fbase.Equilibration.restrt

if [ ! -e Minimization1.mdout ]; then
    runmin1
fi
if [ ! -e Minimization2.mdout ]; then
    runmin2
fi
if [ ! -e Minimization3.mdout ]; then
    runmin3
fi

cp Minimization3.restrt $fbase.a.0.min3.rstrt

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

#fbase: name
#c: number
# run_cp_md: name(fbase), bin, coor
echo "start lastmin(anneal, minimize)"
if [ -e ${fbase}.${c}.stdout -a ! -e ${fbase}.${c}.anneal.restrt ]; then
    runanneal ${fbase}.${c}.rstrt ${fbase}.${c}.anneal
    # in: fbase.c.rstrt out: fbase.c_anneal.restrt
fi
if [ -e ${fbase}.${c}.anneal.restrt -a ! -e ${fbase}.${c}.minlast1.restrt ]; then
    runminlast1 ${fbase}.${c}.anneal.restrt ${fbase}.${c}.minlast1
    # in: fbase.c.anneal out: fbase.c.minlast1.restrt
fi
if [ -e ${fbase}.${c}.minlast1.restrt -a ! -e ${fbase}.${c}.minlast2.restrt ]; then
    runminlast2 ${fbase}.${c}.minlast1.restrt ${fbase}.${c}.minlast2
    # in: fbase.c.minlast1 out: fbase.c.minlast2.restrt
fi
if [ -e ${fbase}.${c}.minlast2.restrt -a ! -e ${fbase}.${c}.minlast3.restrt ]; then
    runminlast3 ${fbase}.${c}.minlast2.restrt ${fbase}.${c}.minlast3
    # out: fbase.c.minlast2 out: fbase.c.minlast3.restrt
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

