#!/bin/bash

# 2020/6/12 v3.3
# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)
module load amber 2> /dev/null

# Warning: cpptraj trajout(startframe) have a bug. please check a MODEL number in tgt.pdb
# v3.3: add error handling for cpptraj bug
# v3.2: add method to exclude minimization.mdcrd
# v3.0: add tempolary solution for cpp traj trajout(start) bug

#--user setting--
# caputure time info(ps)
stime=1100  # please specify the start time of production run
etime=51100
interval=10000
sampletime=1.0

centerinfo=":1-306"
maskinfo=":1-312"
stripdist=4.0
maskflag=true
#--user setting end--


# tempolary solution for cpp traj trajout(start) bug
stime_inprod=0
etime_inprod=`echo "$etime - $stime" | bc`
sframe_inprod=0
eframe_inprod=`echo "$etime_inprod / $sampletime" | bc`
# -----

sframe=`echo "$stime / $sampletime" | bc`
eframe=`echo "$etime / $sampletime" | bc`
ivframe=`echo "$interval / $sampletime" | bc`

name=`ls *.a.0.coor`
head=${name%%.*}

dir='mdout_centerd_foropt'
mkdir $dir 2>/dev/null

prmtop=`ls *.prmtop`
init=`ls *.a.0.coor`
prodinit=`ls ${head}*.${stime}.rstrt`
if [ "$prodinit" == "" ]; then
    prodinit="xxxxx"
    echo "Error!!: check start time setting"
    exit 0
fi

trajs=`ls *.mdcrd`
# echo $trajs

# input filter
trajs_filt=''
for traj in $trajs
do
    if [ "$traj" == "Minimization1.mdcrd" ] || [ "$traj" == "Minimization2.mdcrd" ] || [ "$traj" == "Minimization3.mdcrd" ]; then
        echo "$traj skip"
        continue
    fi
    buf=${traj%.*}
    number=${buf##*.}
    if [ $number -le $stime ]; then
        echo $traj skip
    else
        trajs_filt="$trajs_filt $traj"
    fi

    if [ $number -ge $etime ]; then
        echo $traj: last
        break
    fi
done
echo "input trajs: $trajs_filt"



#-- print section --
echo stime: $stime, etime: $etime, getinterval: $getinterval, sampletime: $sampletime
echo centerinfo, $centerinfo
echo stripdist, $stripdist

echo sframe, $sframe, eframe, $eframe, ivframe, $ivframe

sleep 3

# input section
if [ $stime_inprod -eq 0 ]; then
    echo 'add init file'
    ls $prodinit > inbuf.txt
    ls $trajs_filt >> inbuf.txt
    sframe_inprod=`echo "$sframe_inprod + 1" | bc`
    eframe_inprod=`echo "$eframe_inprod + 1" | bc`
else
    ls $trajs_filt > inbuf.txt
fi

echo "parm $prmtop" > cpptraj.in

# write section
while read line
do
    echo $line
    echo "trajin $line" >> cpptraj.in
done < inbuf.txt

echo "parminfo $prmtop" >> cpptraj.in
echo "autoimage anchor $centerinfo origin" >> cpptraj.in

for i in `seq $sframe_inprod $ivframe $sframe_inprod`
do
    label=`echo "($i + $sframe -1) * $sampletime" | bc | awk '{printf("%d\n",$1)}'`
    newtraj=tgt${label}.rst
    echo "trajout $newtraj onlyframes $i restart" >> cpptraj.in
done
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in


mv tgt*.rst $dir
cp $prmtop $dir
echo "$dir/tgt*.rst ($stime,$etime) was generated."

