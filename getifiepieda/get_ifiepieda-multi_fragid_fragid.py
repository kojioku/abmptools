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
aobj.anlmode= 'multi' #frag, 'mol', 'fraginmol', 'multi'
aobj.fragmode = 'manual'  #'hybrid', 'auto', 'manual'
aobj.dist = 8.0
aobj.abinit_ver='rev20'

aobj.start = 100
aobj.end = 19100
aobj.interval = 1000

# ./6lu7_minHSCG_163hip100-hopt--mod_forabmp_192n-2p-24t.log

aobj.ilog_head = '6lu7orig_md040j8_163neu-'
aobj.ilog_tail = '-hopt-ps-mod_forabmp_192n-2p-24t.log'

aobj.tgt2type = 'frag'
aobj.pynp = 4
aobj.f90soflag = True
aobj.addresinfo = True

print(aobj.tgt2type)
# ---- user setting end ---

tgt1 = sys.argv[1]
tgt2 = sys.argv[2]

aobj = aobj.readifiewrap(tgt1, tgt2)
#out
aobj.writecsvwrap()


