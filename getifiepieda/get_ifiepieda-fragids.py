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

print(aobj.tgt2type)
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




'''
# parameters

anlmode:
fragmode
dist
tgt2type

(multi mode)
    argv[1]: tgt1frag id,  [2]: tgtmolname or fragid
(others)
    argv[1]: logname, [2]:tgtid(mol or frag)


# readpdb = False
# pdbname='sbecd7_50nsdynamics_namd2200-moved-sed-around-8.0-for_abmp.pdb'   # 'iss2-spg2-ok20200130opt-for_abmp.pdb'
# abinit_ver = 16
#
# # -- for mol mode or multi mode--
# tgt2type = 'frag' #frag: mol-frag, mol: mol-mol
#
# # -- fraginmol mode --
# tgt1_lofrag = 2
# tgt2molname = '000'
# tgt2_lofrag = 4
#
# # ------ multi mode ------
# # if tgt2type == 'frag':
# tgt1frag = eval(sys.argv[1])
# tgt2frag = int(sys.argv[2])
#
# # if tgt2type == 'molname':
# tgt2molname = sys.argv[2]
# tgt2dist = 3.0
#
# # multi file setting
# name_head = 'sbecd7_50nsdynamics_namd'
# name_tail = '-moved-sed-around-8.0-for_abmp.log'
# pdb_head = 'sbecd7_50nsdynamics_namd'
# pdb_tail = '-moved-sed-around-8.0-for_abmp.pdb'
#
# start = 1700
# end = 4200
# interval = 500
#
# # --hybrid mode
# hyfrag = 317 #320
# hynum = 3

'''

