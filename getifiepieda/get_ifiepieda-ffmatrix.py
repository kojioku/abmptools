import sys
import os
import pandas as pd
import itertools
import copy
import csv
import abmptools as ampt

# get_ifiepieda.py
# Author: Koji Okuwaki

aobj = ampt.anlfmo()
# --- user setting ---
aobj.anlmode= 'frag' #frag, 'mol', 'fraginmol', 'multi'
aobj.fragmode = 'manual'  #'hybrid', 'auto', 'manual'
aobj.abinit_ver='rev15'

aobj.tgt2type = 'frag'
aobj.matrixtype='frags-frags'
aobj.f90soflag = True
aobj.logMethod = 'MP3'
aobj.addresinfo = True

print('tgt2type', aobj.tgt2type)
# ---- user setting end ---

logname = sys.argv[1]
tgtfrag1 = sys.argv[2]
tgtfrag2 = sys.argv[3]

aobj = aobj.readifiewrap(logname, tgtfrag1, tgtfrag2)
aobj = aobj.filterifiewrap()

aobj.writecsvwrap()

