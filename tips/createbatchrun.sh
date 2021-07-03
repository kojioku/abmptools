#!/bin/bash
stime=1
etime=20
interval=1
script='runab24_fgk202107.sh'
temp='1eo8-gmxtrr-wat4ang-xxx-renamed-MP2-6-31Gd-nbo.ajf'
fhead=${temp%xxx*}
ftail=${temp##*xxx}

for i in `seq $stime $interval $etime`
do
    bash $script $fhead$i$ftail
done
