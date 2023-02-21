#!/bin/bash

# updated on 2023/02/20
# generated on 2020/06
# Author Koji Okuwaki

# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)

# Warning: cpptraj trajout(startframe) have a bug. please check a MODEL number in tgt.pdb
# v3.3: add error handling for cpptraj bug
# v3.2: add method to exclude minimization.mdcrd
# v3.0: add tempolary solution for cpp traj trajout(start) bug

# --- user setting ---
## caputure time info(ps)
stime=0  # please specify the start time of production run
etime=2000
interval=200
sampletime=1.0
# --- user setting end ---

# load module (for tsubame)
module load amber 2> /dev/null

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

dirhead='mdout_foropt'
nowtime=`date +"%Y%m%d%I%M"`
dir=${dirhead}.${nowtime}
mkdir ${dir} 2>/dev/null

parm=`ls *.psf`
xsc=`ls *.a.0.xsc`
init=`ls *.a.0.coor`
trajs=`ls *.coor.dcd`

# prodinit=`ls ${head}*.${stime}.coor`
# if [ "$prodinit" == "" ]; then
#     prodinit="xxxxx"
#     echo "Error!!: check start time setting"
#     exit 0
# fi

# --- input filter ---
trajs_filt=''
for traj in $trajs
do
    buf=${traj%.*}
    buf2=${buf%.*}
    number=${buf2##*.}
    # echo $number, $buf2
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

# --- parm print section ---
echo stime: $stime, etime: $etime, interval: $interval, sampletime: $sampletime
echo centerinfo, $centerinfo

echo sframe, $sframe, eframe, $eframe, ivframe, $ivframe
sleep 3

# --- input section ---
if [ $stime_inprod -eq 0 ]; then
    initbox=`tail -n 1 $xsc | awk '{printf ("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%16i\n", $2, $6, $10, 90, 90, 90, 1)}'`
    sed -e "9i $initbox" $init > init.coor
    ls init.coor > inbuf.txt
    ls $trajs_filt >> inbuf.txt
    sframe_inprod=`echo "$sframe_inprod + 1" | bc`
    eframe_inprod=`echo "$eframe_inprod + 1" | bc`
else
    ls $trajs_filt > inbuf.txt
fi

# --- write section ---
echo "parm $parm" > cpptraj.in
while read line
do
    # echo $line
    echo "trajin $line" >> cpptraj.in
done < inbuf.txt

echo "parminfo $parm" >> cpptraj.in
# echo "autoimage anchor $centerinfo origin" >> cpptraj.in
#echo "center $centerinfo" >> cpptraj.in
#echo "image" >> cpptraj.in

for i in `seq $sframe_inprod $ivframe $eframe_inprod`
do
    label=`echo "($i + $sframe -1) * $sampletime" | bc | awk '{printf("%d\n",$1)}'`
    newtraj=tgt${label}ps.dcd
    echo "trajout $newtraj dcd onlyframes $i restart" >> cpptraj.in
done
echo 'run' >> cpptraj.in

today=$(date "+%Y%m%d%I%M")
ofile="cuttraj-foropt.${today}.log"

cpptraj < cpptraj.in | tee $ofile

# move generated structure (each frame)
mv tgt*.dcd $dir 
cp $parm $dir
echo "$dir/tgt*.dcd ($stime,$etime) was generated."


echo `date` >> $ofile
echo $parm >> $ofile
echo $traj >> $ofile
echo $starttime, $endtime, $interval >> $ofile
