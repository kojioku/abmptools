import sys
import os
import pandas as pd
import itertools
import copy
import csv
import ampt.abinit_io as abio
import ampt.cutset_fmo as cutf

# get_ifiepieda.py
# Author: Koji Okuwaki
# v2.2(2020.03.25): add multi-molname mode
# v2.1(2020.03.19): only gas phase ifie and pieda

# --- user setting ---
# note: start from label 1(frag, mol)

molmode= 'ff-multi' #frag, 'mol', 'fraginmol', 'ff-multi'
fragmode = 'manual'  #'hybrid', 'auto', 'manual'
dist = 8.0

# for mol molmode
molselect = 'mol' # mol or frag

# for ff-multi
if molmode == 'ff-multi':
    frag1 = eval(sys.argv[1])
    print(type(frag1))

    tgt2type = 'molname' # frag or molname
    if tgt2type == 'frag':
        frag2 = int(sys.argv[2])

    if tgt2type == 'molname':
        tgt2molname = sys.argv[2]
        tgt2dist = 3.0
    name_head = 'sbecd7_50nsdynamics_namd'
    name_tail = '-moved-sed-around-8.0-for_abmp.log'

    pdb_head = 'sbecd7_50nsdynamics_namd'
    pdb_tail = '-moved-sed-around-8.0-for_abmp.pdb'

    start = 1700
    end = 4200
    interval = 500

readpdb = False
abinit_ver = 15

# fraginmol mode
if molmode == 'fraginmol':
    tgt1_lofrag = 2
    tgt2molname = '000'
    tgt2_lofrag = 4

pdbname='sbecd7_50nsdynamics_namd2200-moved-sed-around-8.0-for_abmp.pdb'   # 'iss2-spg2-ok20200130opt-for_abmp.pdb'

# --hybrid setting
hyfrag = 317 #320
hynum = 3
# --------------------

# load module
obj = cutf.cut_fmo()

## step1 setup param
if molmode == 'ff-multi':

    lognames = []
    pdbnames = []
    times = []
    for i in range(start, end+1, interval):
        times.append(str(i))
        lognames.append(name_head + str(i) + name_tail)
        pdbnames.append(pdb_head + str(i) + pdb_tail)
    print('lognames', lognames)
    print('molmode:' , molmode)
    print('fragmode:', fragmode)
    if tgt2type == 'frag':
        print('frag1, frag2', frag1, frag2)
    if tgt2type == 'molname':
        print('frag1, frag2mol', frag1, tgt2molname)

else:
    argvs = sys.argv
    logname = argvs[1]
    tgtid = int(sys.argv[2])

    # --- print setting ---
    print('logname:', logname)
    print('dist:', dist)
    print('tgtid:', tgtid)
    print('molmode:' , molmode)
    print('fragmode:', fragmode)
    print('molselect:', molselect)

## read pdb
if molmode == 'fraginmol':
    readpdb = True
if  readpdb == True:
    obj.assignmolname = False
    totalMol, atomnameMol, molnames, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfo(pdbname)
    print(molnames)

if molmode == 'ff-multi' and tgt2type == 'molname':
    nfs = []
    molnames_perrec = []
    obj.assignmolname = False
    for i in range(len(pdbnames)):
        totalMol, atomnameMol, molnames, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfo(pdbnames[i])
        nf = obj.getnf(lognames[i], fragmode)
        nfs.append(nf)
        molnames_perrec.append(molnames)

# ------------
# --- ifie ---
# ------------

print('--- ifie ---')
column = ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE', 'MP2-IFIE', 'PR-TYPE1', 'GRIMME', 'JUNG', 'HILL']

### read ifie (from log) section
if molmode == 'ff-multi':
    dfs = []
    skips = []
    for i in range(len(lognames)):
        ifie = obj.read_ifie(lognames[i])
        dfs.append(obj.getifiedf(ifie, column))
        if len(ifie) == 0:
            print('Warning: time', times[i], 'not converged: skip data')
            skips.append(i)

    dellist = lambda items, indexes: [item for index, item in enumerate(items) if index not in indexes]
    times = dellist(times,skips)
    # print(dfs)

    if tgt2type == 'molname':
        print('### read frag info ###')
        molfragss_perrec = []
        for i in range(len(pdbnames)):
            molfragss = obj.getallmolfrags(lognames[i], dfs[i], nfs[i], molmode)
            # print('frags_permol\n', molfragss)
            molfragss_perrec.append(molfragss)
        # print(molfragss_perrec)
        print(len(molnames), len(molfragss))

else:
    ifie = obj.read_ifie(logname)
    df = obj.getifiedf(ifie, column)
    # print(df.head())

### filter section
# frag mode
if molmode == 'frag':
    tgtdf, tgtdf_filter = obj.gettgtdf(df, dist)
    print(tgtdf_filter)

# multi mode
if molmode == 'ff-multi':
    tgtdf_filters = pd.DataFrame()
    print('times', times)

    count = 0
    if tgt2type == 'frag':
        for df in dfs:
            # print('read', lognames[count])
            tgtdf_filters = tgtdf_filters.append(obj.gettgtdf_ffmulti(df, frag1, frag2))
            count += 1
        tgtdf_filters['TIME'] = times
        print(tgtdf_filters)


    if tgt2type == 'molname':
        # get tgt frag id
        tgtmolfrags_perrec = []
        HF_IFIE_sums = []
        MP2_IFIE_sums = []
        PR_TYPE1_sums = []
        GRIMME_sums = []
        JUNG_sums = []
        HILL_sums = []

        for i in range(len(molnames_perrec)):
            tgtmolfrags = []
            for j in range(len(molnames_perrec[i])):
                try:
                    if molnames_perrec[i][j] == tgt2molname:
                        tgtmolfrags.append(molfragss_perrec[i][j])
                except:
                    continue

            tgtmolfrags_perrec.append(tgtmolfrags)
        print(len(tgtmolfrags_perrec))

        if type(frag1) == int:
            frag1s = [frag1]
            print(frag1s)
        else:
            frag1s = copy.deepcopy(frag1)

        for i in range(len(tgtmolfrags_perrec)):
            df = dfs[i]
            if df.empty:
                # print(i, 'empty!!')
                continue
            tgtdf_filters = pd.DataFrame()
            for frag1p in frag1s:
                print('frag1p', frag1p)
                for tgtmolfrags in tgtmolfrags_perrec[i]:
                    for frag2 in tgtmolfrags:
                        tgtdf_filters = tgtdf_filters.append(obj.gettgtdf_ffmulti(df, frag1p, frag2))
                        count += 1
                # print(tgtdf_filters.tail())

            HF_IFIE_sums.append(tgtdf_filters['HF-IFIE'].sum())
            MP2_IFIE_sums.append(tgtdf_filters['MP2-IFIE'].sum())
            PR_TYPE1_sums.append(tgtdf_filters['PR-TYPE1'].sum())
            GRIMME_sums.append(tgtdf_filters['GRIMME'].sum())
            JUNG_sums.append(tgtdf_filters['JUNG'].sum())
            HILL_sums.append(tgtdf_filters['HILL'].sum())

        tgtdfsum = pd.DataFrame()
        tgtdfsum['HF-IFIE'] = HF_IFIE_sums
        tgtdfsum['MP2-IFIE'] = MP2_IFIE_sums
        tgtdfsum['PR-TYPE1'] = PR_TYPE1_sums
        tgtdfsum['GRIMME'] = GRIMME_sums
        tgtdfsum['JUNG'] = JUNG_sums
        tgtdfsum['HILL'] = HILL_sums
        tgtdfsum['TIME'] = times

        print(tgtdfsum)


# mol-mol mode
if molmode == 'mol':
    if molselect == 'frag':
        molfrags = obj.getmolfrags(tgtid, df, molmode)
        print('target-frags:', molfrags)
    elif molselect == 'mol':
        nf = obj.getnf(logname, fragmode)
        molfragss = obj.getallmolfrags(logname, df, nf, molmode)
        print('frags_permol\n', molfragss)
        molfrags = molfragss[tgtid-1]

        contactmolfrags, ifie_permols, ifiesum = obj.getifiesummol(df, column)

# fraginmol mode
if molmode == 'fraginmol':
    nf = obj.getnf(logname, fragmode)
    print(nf)
    molfragss = obj.getallmolfrags(logname, df, nf, molmode)
    print(molfragss)
    tgtmol = tgtid - 1

    tgt2_glofrags = []
    tgt1_glofrag =molfragss[tgtmol][tgt1_lofrag - 1]
    print('centermolfrag:', tgt1_glofrag)
    print('tgt2molname', tgt2molname)
    for i in range(len(molnames)):
        if molnames[i] == tgt2molname:
            tgt2frag = molfragss[i][tgt2_lofrag - 1]
            tgt2_glofrags.append(tgt2frag)
    print('tgt2_glofrags', tgt2_glofrags)
    tgtdf = df[df['I'] == tgt1_glofrag]
    tgtdf = tgtdf.append(df[df['J'] == tgt1_glofrag])
    tgtdf = tgtdf[tgtdf['DIST'] < dist]

    tgt_new2 = pd.DataFrame()
    for tgt2_glofrag in tgt2_glofrags:
        tgt_new = tgtdf[(tgtdf['I'] == tgt2_glofrag) |(tgtdf['J'] == tgt2_glofrag)]
        tgt_new2 = tgt_new2.append(tgt_new)

    print('tgt_new2\n', tgt_new2)


# -------------
# --- PIEDA ---
# -------------
if abinit_ver >= 17:
    pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'Solv(ES)', 'DI(MP2)', 'q(I=>J)']
else:
    pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']

print('--- pieda ---')

### read pieda (from log) section
if molmode == 'ff-multi':
    pidfs = []
    for logname in lognames:
        pieda = obj.read_pieda(logname)
        pidfs.append(obj.getpiedadf(pieda, pcolumn))
else:
    pieda = obj.read_pieda(logname)
    pidf = obj.getpiedadf(pieda, pcolumn)

### read fraginfo section
frags = []
if fragmode != 'manual':
    frags = obj.read_fraginfo(logname)
    print('frags', frags)

if fragmode == 'hybrid':
    getf = frags.pop(hyfrag-1)
    for i in range(hynum):
        frags.append(getf)
    print('frags', frags)

### filter
# frag-frag mode
if molmode == 'frag':
    pitgtdf, pitgtdf_filter = obj.getpitgtdf(pidf, tgtdf_filter)
    if fragmode != 'manual':
        print('len_frags', len(frags))
        #assign resname(e.g. Gly6)
        for i in range(1, len(frags) + 1):
            pitgtdf.I = pitgtdf.I.replace(i, frags[i-1])
            pitgtdf.J = pitgtdf.J.replace(i, frags[i-1])

# ff-multi mode
if molmode == 'ff-multi':
    if tgt2type == 'frag':
        pitgtdf_filters = pd.DataFrame()
        count = 0
        for pidf in pidfs:
            # print('read', lognames[count])
            pitgtdf_filters = pitgtdf_filters.append(obj.gettgtdf_ffmulti(pidf, frag1, frag2))
            count += 1
        pitgtdf_filters['TIME'] = times
        print(pitgtdf_filters)

    if tgt2type == 'molname':
        ES_sums = []
        EX_sums = []
        CT_sums = []
        DI_sums = []
        q_sums = []
        for i in range(len(tgtmolfrags_perrec)):
            pidf = pidfs[i]
            if pidf.empty:
                continue
            pitgtdf_filters = pd.DataFrame()
            for frag1p in frag1s:
                for tgtmolfrags in tgtmolfrags_perrec[i]:
                    for frag2 in tgtmolfrags:
                        pitgtdf_filters = pitgtdf_filters.append(obj.gettgtdf_ffmulti(pidf, frag1p, frag2))
                        count += 1
            ES_sums.append(pitgtdf_filters['ES'].sum())
            EX_sums.append(pitgtdf_filters['EX'].sum())
            CT_sums.append(pitgtdf_filters['CT-mix'].sum())
            DI_sums.append(pitgtdf_filters['DI(MP2)'].sum())
            q_sums.append(pitgtdf_filters['q(I=>J)'].sum())

        pitgtdfsum = pd.DataFrame()
        pitgtdfsum['ES'] = ES_sums
        pitgtdfsum['EX'] = EX_sums
        pitgtdfsum['CT-mix'] = CT_sums
        pitgtdfsum['DI(MP2)'] = DI_sums
        pitgtdfsum['q(I=>J)'] = q_sums
        pitgtdfsum['TIME'] = times

        print(pitgtdfsum)


# mol-mol mode
if molmode == 'mol':
    pieda_permols = []
    for contactmolfrag in contactmolfrags:
        pieda_permol = pd.DataFrame(columns=pcolumn)
        for contact in contactmolfrag:
            for tgtfrag in molfrags:
                # print(contact, tgtfrag)
                pieda_permol = pieda_permol.append(pidf[((pidf['I'] == contact) & (pidf['J'] == tgtfrag)) | ((pidf['I'] == tgtfrag) & (pidf['J'] == contact))])
        pieda_permols.append(pieda_permol)
        # print(pieda_permol, file=plogdt)

    count = 0
    if abinit_ver >= 17:
        piedasums = [['contactmolfrag', 'tgtmolfrags', 'ES', 'EX', 'Solv(ES)', 'CT-mix', 'DI(MP2)', 'a(I=>J)']]
    else:
        piedasums = [['contactmolfrag', 'tgtmolfrags', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'a(I=>J)']]
    for datadf in pieda_permols:
        piedasums.append([contactmolfrags[count], molfrags, datadf['ES'].sum(), datadf['EX'].sum(), datadf['CT-mix'].sum(), datadf['DI(MP2)'].sum(),  datadf['q(I=>J)'].sum()])
        count += 1

# fraginmol mode
if molmode == 'fraginmol':

    print('--- pieda ----')
    pitgtdf = pidf[pidf['I'] == tgt1_glofrag]
    pitgtdf = pitgtdf.append(pidf[pidf['J'] == tgt1_glofrag])

    pitgt_new2 = pd.DataFrame()
    for tgt2_glofrag in tgt2_glofrags:
        pitgt_new = pitgtdf[(pitgtdf['I'] == tgt2_glofrag) |(pitgtdf['J'] == tgt2_glofrag)]
        pitgt_new2 = pitgt_new2.append(pitgt_new)

    if fragmode != 'manual':
        print('len_frags', len(frags))

        for i in range(1, len(frags) + 1):
            pitgt_new2.I = pitgt_new2.I.replace(i, frags[i-1])
            pitgt_new2.J = pitgt_new2.J.replace(i, frags[i-1])

    path = 'csv'
    if os.path.exists('csv') == False:
        os.mkdir('csv')
    head, ext = os.path.splitext(logname)

    ohead = head + '-' 'frag' + str(tgt1_glofrag) + '-mol' + str(tgt2molname) + 'frag' + str(tgt2_lofrag)
    print(pitgt_new2.head())

# -------------
# --- write ---
# -------------
# define out name
path = 'csv'
if os.path.exists('csv') == False:
    os.mkdir('csv')

head, ext = os.path.splitext(logname)

if molmode == 'frag':
    try:
        ohead = head + '-' + str(tgtid) + '-' + frags[tgtid - 1]
    except:
        ohead = head + '-' + str(tgtid)

    print(pitgtdf.head())

    tgtdf.to_csv(path + '/' + ohead + '-ifie.csv')
    tgtdf_filter.to_csv(path + '/' + ohead + '-ifie_' + 'dist' + str(dist) + '.csv')
    pitgtdf.to_csv(path + '/' + ohead + '-pieda.csv')

    print(ohead + '-ifie.csv', ohead + '-pieda.csv generated.')

if molmode == 'ff-multi':
    if tgt2type == 'frag':
        oifie = 'frag' + str(frag1) + '-frag' + str(frag2) + '-ifie.csv'
        opieda = 'frag' + str(frag1) + '-frag' + str(frag2) + '-pieda.csv'
        tgtdf_filters.to_csv(path + '/' + oifie)
        pitgtdf_filters.to_csv(path + '/' + opieda)

    if tgt2type == 'molname':
        oifie = 'frag' + str(frag1) + '-' + str(tgt2molname) + '-ifie.csv'
        opieda = 'frag' + str(frag1) + '-' + str(tgt2molname) + '-pieda.csv'
        tgtdfsum.to_csv(path + '/' + oifie)
        pitgtdfsum.to_csv(path + '/' + opieda)

if molmode == 'mol':
    ilogdtname = path + '/' + head + '_ifie-mol-' + molselect + str(tgtid) + 'dist' + str(dist) + '.txt'
    isumname = path + '/' + head + '_ifiesum-mol-' + molselect + str(tgtid) + 'dist' + str(dist) + '.csv'
    plogdtname = path + '/' + head + '_pieda-mol-' + molselect + str(tgtid) + 'dist' + str(dist) + '.txt'
    psumname = path + '/' + head + '_piedasum-mol-' + molselect + str(tgtid) + 'dist' + str(dist) + '.csv'

    # write section
    ilogdt = open(ilogdtname, 'w')
    for ifie_permol in ifie_permols:
        print(ifie_permol, file=ilogdt)

    plogdt = open(plogdtname, 'w')
    for pieda_permol in pieda_permols:
        print(pieda_permol, file=plogdt)

    with open(isumname, 'w') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerows(ifiesums)

    with open(psumname, 'w') as f:
        writer = csv.writer(f, lineterminator='\n')
        # writer.writerow(list)
        writer.writerows(piedasums)

    print('---out---')
    print(ilogdtname)
    print(isumname)
    print(plogdtname)
    print(psumname)

if molmode == 'fraginmol':
    tgtdf.to_csv(path + '/' + head + '-ifie.csv')
    tgt_new2.to_csv(path + '/' + ohead + '-ifie_'  + 'dist' + str(dist) + '.csv')
    pitgt_new2.to_csv(path + '/' + ohead + '-pieda.csv')
    print(ohead + '-ifie.csv', ohead + '-pieda.csv generated.')

