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
aobj.anlmode= 'mol' #frag, 'mol', 'fraginmol', 'ff-multi'
aobj.fragmode = 'manual'  #'hybrid', 'auto', 'manual'
aobj.dist = 2.5
aobj.abinit_ver='rev15'

aobj.selecttype = 'molid'
# aobj.tgt2molname = 'WAT'

print('tgtselecttype', aobj.selecttype)
# ---- user setting end ---

logname = sys.argv[1]
aobj.tgtmolid = sys.argv[2]

aobj = aobj.readifiewrap(logname)
aobj = aobj.filterifiewrap()

aobj.writecsvwrap()

