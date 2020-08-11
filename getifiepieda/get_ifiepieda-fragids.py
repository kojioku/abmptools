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
aobj.tgt2type = 'frag'
aobj.abinit_ver='rev15'
# aobj.tgt2molname = 'WAT'

# ---- user setting end ---

logname = sys.argv[1]
tgtfrag1 = sys.argv[2]
tgtfrag2 = sys.argv[3]

aobj = aobj.readifiewrap(logname, tgtfrag1, tgtfrag2)
aobj = aobj.filterifiewrap()

# print('ifdf\n', aobj.ifdfs)
# print('ifdf_filter\n', aobj.ifdf_filters)
# print('pidf\n', aobj.pidfs)
# # print('pitgtdf\n', aobj.pitgtdfs)
# print('pitgtdf\n', aobj.pidf_filters)

aobj.writecsvwrap()




