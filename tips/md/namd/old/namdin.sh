#!/bin/sh
# WARNING! CURRENT DIRECTORY MUST BE THAT OF THIS FILE
# Automatically generated on Tue Aug 27 14:51:04 2019.
# Written by MOE 2019.01 (Chemical Computing Group ULC)

fb="namd_bex"
t0=".a.0"
namd="namd2"
charmrun=""

dir=`/usr/bin/dirname "$0"`
cd "$dir"
if [ ! -f "$fb.z.cfg" ]; then
    echo configuration file "$fb.z.cfg" does not exist
    exit 1;
fi
for x in coor psf par fix
do
    if [ ! -f "$fb$t0.$x" ]; then
	echo configuration file "$fb$t0.$x" does not exist
	exit 1;
    fi
done
while :
do
    touch "$fb.z.err"
        $namd  +p$1   $fb.z.cfg >> "$fb.z.stdout"
    if [ -f "$fb.z.err" ]; then exit 1; fi
    if [ -f "$fb.z.coor" -a -f "$fb.z.vel" -a -f "$fb.z.xsc" ]; then
	break;
    fi
done
exit 0;
