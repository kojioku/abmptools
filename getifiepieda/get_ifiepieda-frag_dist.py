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
aobj.anlmode= 'frag' #frag, 'mol', 'fraginmol', 'multi'
aobj.fragmode = 'manual'  #'hybrid', 'auto', 'manual'
aobj.dist = 8.0
aobj.abinit_ver='rev15'
aobj.tgt2type = 'dist'
# aobj.tgt2molname = 'WAT'

print(aobj.tgt2type)
# ---- user setting end ---

logname = sys.argv[1]
tgtfragid = sys.argv[2]

aobj = aobj.readifiewrap(logname, int(tgtfragid))
aobj = aobj.filterifiewrap()

aobj.writecsvwrap()


