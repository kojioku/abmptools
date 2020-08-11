import sys
import os
import pandas as pd
import itertools
import copy
import csv
import abmptools as ampt

# get_ifiepieda.py
# Author: Koji Okuwaki
# v3.0(2020.05.10): reform system
# v2.2(2020.03.25): add multi-molname mode
# v2.1(2020.03.19): only gas phase ifie and pieda

# Warning:
# note: start from label 1(frag, mol)

aobj = ampt.anlfmo()
# --- user setting ---
aobj.anlmode= 'multi' #frag, 'mol', 'fraginmol', 'ff-multi'
aobj.fragmode = 'manual'  #'hybrid', 'auto', 'manual'
aobj.dist = 8.0
aobj.abinit_ver='rev20'

aobj.start = 1700
aobj.end = 4200
aobj.interval = 500

aobj.ilog_head = 'sbecd7_50nsdynamics_namd'
aobj.ilog_tail = '-moved-sed-around-8.0-for_abmp.log'

aobj.pynp = 4
aobj.tgt2type = 'molname'
# aobj.tgt2molname = 'WAT'

print(aobj.tgt2type)
# ---- user setting end ---

tgt1 = sys.argv[1]
tgt2 = sys.argv[2]

aobj = aobj.readifiewrap(tgt1, tgt2)

aobj.writecsvwrap()



