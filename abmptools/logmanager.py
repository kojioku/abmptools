import sys
import os
import gzip
import pprint
scrdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scrdir)

from multiprocessing import Pool
import copy
import random
import math
import re
import subprocess
import csv
import pdb_io as pdio
import time
import abmptools

try:
    import numpy as np
except:
    pass
try:
    import collections
except:
    pass
try:
    import pandas as pd
except:
    pass
try:
    import itertools
except:
    pass


class LOGManager():
    def __init__(self):
        self.cpfver = 23
        self.atominfo = {}
        self.fraginfo = []
        self.condition = []
        self.static_data = {}
        self.mominfo = {}
        self.diminfo = {}
        self.labels = {}
        self.tgtfrag = 0
        self.is_gz = False
        return

        '''
        atominfo = {
            'alabels': [],  # 原子の番号(i10)
            'elems': [],  # 元素記号(a2)
            'elemtypes': [],  # 原子タイプ(a4)
            'resnames': [],  # 残基名(a3)
            'resids': [],  # 残基番号(i10)
            'fragids': [],  # フラグメント番号(i10)
            'xcoords': [],  # x座標(f20.10)
            'ycoords': [],  # y座標(f20.10)
            'zcoords': [],  # z座標(f20.10)
            'chainids': [],  # Chain ID(a3)
            'optflags': [],  # 構造最適化オプション(i1)
            }
        '''

        '''
            charge_label = ["MUL-HF", "MUL-MP2", "NPA-HF", "NPA-MP2", "ESP-HF", "ESP-MP2"]
            DPM_label = ["DPM-HF-X", "DPM-HF-Y", "DPM-HF-Z", "DPM-MP2-X",
                         "DPM-MP2-Y", "DPM-MP2-Z"]
            monomer_label = ["NR", "HF", "MP2", "MP3"]
            dimer_label = ["NR", "HF", "ES", "MP2", "PR-MP2", "SCS-MP2(Grimme)",
                           "MP3", "SCS-MP3(MP2.5)", "HF-BSSE", "MP2-BSSE",
                           "SCS-MP2-BSSE", "MP3-BSSE", "SCS-MP3-BSSE",
                           "SOLV-ES", "SOLV-NP", "EX", "CT", "DQ"]
        '''

        '''
        labels = {
            'charge': charge_label,
            'DPM': DPM_label,
            'monomer': monomer_label,
            'dimer': dimer_label,
            }
        '''

        '''
        fraginfo = {
            'natoms': fnatoms,
            'baas': fbaas,
            'connects': fconnects
        }
        '''

        '''
        'static_data': {'natom': 38,  #総原子数
                        'ndimer': 10, # <- fraginfoから
                        'nfrag': 5, #AUTOMATIC FRAGMENTATION もしくは NF <- fraginfoからとれる
                        'ntetramer': '0',
                        'ntrimer': '0',
                        'nuclear_repulsion_energy': 1889.49334720085,
                        'total_electronic_energy': -2986.14838661631,
                        'total_energy': -1096.65503941546},

            ## FRAGMENT NUCLEAR REPULSION

               NUCLEAR REPULSION ENERGY =          1889.4933472009

         =========================
            ## FMO TOTAL ENERGY
         =========================

               FMO2-HF
               Nuclear repulsion =          1889.4933472009
               Electronic energy =         -2985.0596593788
               Total      energy =         -1095.5663121780

               E(FMO2-MP2)       =            -1.0887272375
               Electronic energy =         -2986.1483866163
               Total      energy =         -1096.6550394155
        '''

    def parse(self, filepath):
        # Initialize the data structures
        # natom = 0
        # nfrag = 0
        # atom_data = []
        charge_label = ["MUL-HF", "MUL-MP2", "NPA-HF", "ESP-HF"]
        DPM_label = ["DPM-HF-X", "DPM-HF-Y", "DPM-HF-Z"]
        monomer_label = ["NR", "HF", "MP2",]
        dimer_label = ["NR", "HF", "ES", "MP2", "PR-MP2", "SCS-MP2(Grimme)",
                       "EX", "CT", "DQ"]

        labels = {
            'charge': charge_label,
            'DPM': DPM_label,
            'monomer': monomer_label,
            'dimer': dimer_label,
            }

        self.labels = labels

        print('start read', filepath)
        # filepath の拡張子が .gz なら gzip で読み込む
        if os.path.splitext(filepath)[-1] == '.gz':
            self.is_gz = True

        if self.is_gz:
            file = gzip.open(filepath, 'rt')
        else:
            file = open(filepath, 'rt')

        # condition
        Method, ElecState, BasisSet, Laoc, Lptc, Ldimer, ReadGeom, fragmode \
            = self.getcondition(file)

        condition = {
            'aoc': Laoc,
            'basis_set': BasisSet,
            'calculation_method': Method,
            'electronic_state': ElecState,
            'ldimer': Ldimer,
            'ptc': Lptc,
            'readgeom': ReadGeom,
        }

        elems = self.readelemslog(file)

        # fraginfo
        nf, fnatoms, fbaas, fconnects, fragdatas, natom\
            = self.getfraginfo(file, fragmode)

        '''
        'static_data': {'natom': 38,  #総原子数
                        'ndimer': 10, # <- fraginfoから
                        'nfrag': 5, #AUTOMATIC FRAGMENTATION もしくは NF <- fraginfoからとれる
                        'ntetramer': '0',
                        'ntrimer': '0',
                        'nuclear_repulsion_energy': 1889.49334720085,
                        'total_electronic_energy': -2986.14838661631,
                        'total_energy': -1096.65503941546},
        '''

        fraginfo = {
            'natoms': fnatoms,
            'baas': fbaas,
            'connects': fconnects,
            'atomlabel': fragdatas,
        }

        # Create a dictionary to store the index of each number
        number_to_index = {}
        # Iterate through the data list to find the index of each number
        for index, sublist in enumerate(fragdatas):
            # print(index, sublist)
            for number in sublist:
                number_to_index[number] = index + 1
        # キーの数字が若い順に辞書のキーをソート
        sorted_keys = sorted(number_to_index.keys())
        # ソートされた順番でリストを作成
        fragids = [number_to_index[key] for key in sorted_keys]

        # number_to_index
        NR = self.getnr(file)

        static_data = {
            'natom': natom,
            'ndimer': int(nf * (nf - 1) / 2),
            'nfrag': nf,
            'ntetramer': 0,
            'ntrimer': 0,
            'nuclear_repulsion_energy': NR,
        }

        nums, heads, molnames, atypenames, labs, chains, \
            resnums, codes, xcoords, ycoords, zcoords, occs, temps, \
            amarks, charges, totalatom \
            = self.readpdb(ReadGeom)

        atominfo = {
            'alabels': nums,  # 原子の番号(i10)
            'elems': elems,  # 元素記号(a2)
            'elemtypes': atypenames,  # 原子タイプ(a4)
            'resnames': molnames,  # 残基名(a3)
            'resids': resnums,  # 残基番号(i10)
            'fragids': fragids,  # フラグメント番号(i10)
            'xcoords': xcoords,  # x座標(f20.10)
            'ycoords': ycoords,  # y座標(f20.10)
            'zcoords': zcoords,  # z座標(f20.10)
            'chainids': ["" for i in range(natom)],  # Chain ID(a3)
            'optflags': ["" for i in range(natom)],  # 構造最適化オプション(i1)
            }

        for chg in charge_label:
            atominfo[chg] = [0.0 for i in range(natom)]

        # pprint.pprint(vars(pdbobj))
        print('atominfo\n', self.atominfo)

        hf, mp2 = self.getmominfo(file, nf)
        dpm_hf_x, dpm_hf_y, dpm_hf_z = self.getdipoleinfo(file, nf)

        print('dipole moment section')
        mominfo = {
            'fragi': [i+1 for i in range(nf)],
            }

        # Initialize the data structures
        for dpm in DPM_label:
            mominfo[dpm] = []

        print('monomer section')
        # Initialize the data structures
        for mom in monomer_label:
            mominfo[mom] = []

        mominfo['NR'] = [0 for i in range(nf)]
        mominfo['HF'] = hf
        mominfo['MP2'] = mp2
        mominfo['DPM-HF-X'] = dpm_hf_x
        mominfo['DPM-HF-Y'] = dpm_hf_y
        mominfo['DPM-HF-Z'] = dpm_hf_z

        self.mominfo = mominfo
        print('mominfo\n', self.mominfo)

        if Method == 'MP2':
            ifrags, jfrags, dists, hfs, mp2s, prs, grimmes, \
                pifrags, pjfrags, ess, exs, cts, dis, dqs \
                = self.read_mp2ifiepieda(file)
        else:
            print("Methods other than MP2 are not supported.")
            sys.exit()

        ifieinfo = {
            "fragi": ifrags,
            "fragj": jfrags,
            "min-dist": dists,
            "NR": [0 for i in range(len(ifrags))],
            "HF": hfs,
            "MP2": mp2s,
            "PR-MP2": prs,
            "SCS-MP2(Grimme)": grimmes,
        }

        piedainfo = {
            "fragi": pifrags,
            "fragj": pjfrags,
            "ES": ess,
            "EX": exs,
            "CT": cts,
            "DQ": dqs,
        }


        # for dim in dimer_label:
        #     diminfo[dim] = []

        ifiedf = pd.DataFrame(ifieinfo)
        piedadf = pd.DataFrame(piedainfo)
        ifpidf = pd.merge(ifiedf, piedadf, on=['fragi', 'fragj'])
        # print(ifpidf)

        total_electronic_energy, total_energy = self.gettotalenergy(file)
        static_data['total_electronic_energy'] = total_electronic_energy
        static_data['total_energy'] = total_energy

        file.close()

        self.atominfo = pd.DataFrame(atominfo)
        self.fraginfo = fraginfo
        self.condition = condition
        self.static_data = static_data
        self.mominfo = pd.DataFrame(mominfo)
        self.diminfo = ifpidf
        self.labels = labels

        print('atominfo\n', self.atominfo)
        print('fraginfo\n', self.fraginfo)
        print('condition\n', self.condition)
        print('static_data\n', self.static_data)
        print('mominfo\n', self.mominfo)
        print('diminfo\n', self.diminfo)
        print('labels\n', self.labels)

        '''
        self.cpfver = cpfver
        self.atominfo = pd.DataFrame(atominfo)
        self.fraginfo = fraginfo
        self.condition = condition
        self.static_data = static_data
        self.mominfo = pd.DataFrame(mominfo)
        self.diminfo = pd.DataFrame(diminfo)
        self.labels = labels
        '''

    @staticmethod
    def gettotalenergy(file):
        flag = False
        eneflag = False
        count = 0
        gcount = 0
        data = []
        print('--- get totalenergy from log ---')
        while True:
            Items = file.readline().strip().split()
            if eneflag:
                data.append(Items[3])
                gcount += 1
                if gcount == 2:
                    break
                continue
            if flag:
                count += 1
                if flag and count == 8:
                    flag = False
                    eneflag = True
                    continue
                continue
            if len(Items) <= 2:
                continue
            if Items[0:3] == ['##', 'FMO', 'TOTAL']:
                flag = True
                continue
        return data[0], data[1]

    @staticmethod
    def getmominfo(file, nmom):
        flag = False
        eneflag = False
        count = 0
        gcount = 0
        hf = []
        mp2 = []
        print('--- get mominfo from log ---')
        while True:
            Items = file.readline().strip().split()
            if eneflag:
                hf.append(float(Items[1]))
                mp2.append(float(Items[2]))
                gcount += 1
                if gcount == nmom:
                    break
                continue
            if flag:
                count += 1
                if flag and count == 3:
                    flag = False
                    eneflag = True
                    continue
                continue
            if len(Items) <= 2:
                continue
            if Items[0:3] == ['##', 'MONOMER', 'ENERGY']:
                flag = True
                continue
        return hf, mp2

    @staticmethod
    def getdipoleinfo(file, nmom):
        flag = False
        eneflag = False
        count = 0
        gcount = 0
        hf_x = []
        hf_y = []
        hf_z = []
        print('--- get dipoleinfo from log ---')
        while True:
            Items = file.readline().strip().split()
            if eneflag:
                hf_x.append(float(Items[1]))
                hf_y.append(float(Items[2]))
                hf_z.append(float(Items[3]))
                gcount += 1
                if gcount == nmom:
                    break
                continue
            if flag:
                count += 1
                if flag and count == 6:
                    flag = False
                    eneflag = True
                    continue
                continue
            if len(Items) <= 2:
                continue
                # MONOMER ELECTRIC DIPOLE MOMENT
            if Items[0:3] == ['##', 'MONOMER', 'ELECTRIC']:
                flag = True
                continue
                break
        return hf_x, hf_y, hf_z

    @staticmethod
    def readelemslog(file):
        flag = False
        elems = []
        print('--- get eleminfo from log ---')
        while True:
            Items = file.readline().strip().split()
            if len(Items) <= 2:
                continue
            # READ MOLECULAR STRUCTURE FROM PDB-FILE
            if Items[0:3] == ['##', 'READ', 'MOLECULAR']:
                flag = True
                continue
            # Molecular formula
            if Items[0:3] == ['##', 'Molecular', 'formula']:
                print('End Read Elems')
                break
            if flag:
                elems.append(Items[1])
        return elems

    @staticmethod
    def getcondition(file):

        while True:
            Items = file.readline().strip().split()
            if len(Items) <= 2:
                continue
            if Items[0:2] == ['Method', '=']:
                # print('logMethod =', Items[2])
                Method = Items[2]
            if Items[0:2] == ['ElecState', '=']:
                # print('ElecState =', Items[2])
                ElecState = Items[2]
            if Items[0:2] == ['BasisSet', '=']:
                # print('BasisSet =', Items[2])
                BasisSet = Items[2]
            if Items[0:2] == ['Laoc', '=']:
                # print('Laoc =', Items[2])
                Laoc = float(Items[2])
            if Items[0:2] == ['Lptc', '=']:
                # print('Lptc =', Items[2])
                Lptc = float(Items[2])
            if Items[0] == 'ReadGeom':
                ReadGeom = Items[2]
            if Items[0] == 'AutoFrag':
                if Items[2] == 'ON':
                    fragmode = 'auto'
                else:
                    fragmode = 'manual'
            if Items[0:2] == ['Ldimer', '=']:
                # print('Ldimer =', Items[2])
                Ldimer = float(Items[2])
            if Items[0:3] == ['##', 'CHECK', 'AVAILABLE']:
                print('End Read Condition')
                break

            '''
                ElecState               = S1
                BasisSet                = STO-3G
                Method                  = MP2
                Laoc                    =     0.000
                Lptc                    =     2.000
                Ldimer                  =     2.000

            '''
            '''
              'condition': {
              'aoc': 0.0,
              'basis_set': 'STO-3G',
              'calculation_method': 'MP2',
              'electronic_state': 'S1',
              'ldimer': 2.0,
              'ptc': 2.0},
            '''

        return Method, ElecState, BasisSet, Laoc, Lptc, Ldimer, ReadGeom, fragmode

#     @staticmethod
#     def getpiedadf(self, pieda):
#         # print('l669', pieda[1])
#         pidf = pd.DataFrame(pieda, columns=self.pcolumn)
#         pidf['I'] = pidf['I'].astype(int)
#         pidf['J'] = pidf['J'].astype(int)
#         pidf['ES'] = pidf['ES'].astype(float)
#         pidf['EX'] = pidf['EX'].astype(float)
#         pidf['CT-mix'] = pidf['CT-mix'].astype(float)
#         if self.abinit_ver == 'rev16' or self.abinit_ver == 'rev17':
#             pidf['Solv(ES)'] = pidf['Solv(ES)'].astype(float)
#         pidf['DI(MP2)'] = pidf['DI(MP2)'].astype(float)
#         pidf['q(I=>J)'] = pidf['q(I=>J)'].astype(float)
#         return pidf
#
#     @staticmethod
#     def getifiedf(self, ifie, solv=[]):
#         '''get ifie data frame from ifie file
#
#         get ifie data frame from ifie file
#
#         Args:
#             ifie (str): ifie file name
#             solv (list): solvent name list
#
#         Returns:
#             df (pandas.DataFrame): ifie data frame
#         '''
#
#         df = pd.DataFrame(ifie, columns=self.icolumn)
#         df['I'] = df['I'].astype(int)
#         df['J'] = df['J'].astype(int)
#         df['DIST'] = df['DIST'].astype(float)
#         df['HF-IFIE'] = df['HF-IFIE'].astype(float) * 627.5095
#
#         if self.logMethod == 'MP2':
#             df['MP2-IFIE'] = df['MP2-IFIE'].astype(float) * 627.5095
#             df['PR-TYPE1'] = df['PR-TYPE1'].astype(float) * 627.5095
#             df['GRIMME'] = df['GRIMME'].astype(float) * 627.5095
#             df['JUNG'] = df['JUNG'].astype(float) * 627.5095
#             df['HILL'] = df['HILL'].astype(float) * 627.5095
#
#         # print('solv', solv)
#         if len(solv) != 0:
#             solvdf = pd.DataFrame(solv, columns=['I', 'J', 'Solv(ES)'])
#             solvdf['I'] = solvdf['I'].astype(int)
#             solvdf['J'] = solvdf['J'].astype(int)
#             solvdf['Solv(ES)'] = solvdf['Solv(ES)'].astype(float)
#
#             # print(df.head())
#             # print(solvdf.head())
#             df = pd.merge(df, solvdf, on=['I', 'J'], how='left')
#         print(df.head())
#         return df
#
#             df = self.getifiedf(ifie)
#             self.ifdf = df
#
#             pidf = self.getpiedadf(pieda)
#             self.pidf = pidf

    @staticmethod
    def read_mp2ifiepieda(file):
        ''' read ifie and pieda data from file

        read ifie and pieda data from file

        Args:
            fname (str): file name

        Returns:
            ifie (list): ifie data
            pieda (list): pieda data
        '''
        # print('start getifiepieda')

        ifie = []
        count = 0
        pieda = []
        pcount = 0
        pflag = False
        bsseflag = False
        bssecount = 0
        bsse = []

        flag = False
        momflag = False
        dimflag = False

        ifrags = []
        jfrags = []
        dists = []
        hfs = []
        mp2s = []
        prs = []
        grimmes = []
        pifrags = []
        pjfrags = []
        ess = []
        exs = []
        cts = []
        dis = []
        dqs = []

        # print text
        print('--- get mominfo from log ---')
        while True:
            Items = file.readline().strip().split()

            if len(Items) < 2:
                continue
            if Items[1] == 'MP2-IFIE':
                flag = True
                continue
            if Items[1] == 'PIEDA':
                flag = False
                pflag = True
                continue
            if flag:
                count += 1
            if flag and count > 2:
                # ifie.append(Items)
                ifrags.append(int(Items[0]))
                jfrags.append(int(Items[1]))
                dists.append(float(Items[2]))
                if float(Items[4]) < -2:
                    Items[4] = 0.0
                hfs.append(float(Items[4]))
                mp2s.append(float(Items[5]))
                prs.append(float(Items[6]))
                grimmes.append(float(Items[7]))

            if len(Items) < 2:
                continue
            # after pieda or BSSE (break)
            if Items[1] == 'Mulliken':
                # flag = False
                break

            # for BSSE
            # if pflag and Items[:5] == ['##', 'BSSE', 'for','non-bonding','MP2-IFIE']:
            #     pflag = False
            #     # print('pieda end! next is BSSE')
            #     continue
            # if Items[:4] == ['##', 'BSSE', 'for', 'MP2-IFIE']:
            #     bsseflag = True
            #     # print('BSSE start!')
            #     continue
            # if bsseflag:
            #     bssecount += 1
            # if bsseflag and bssecount > 2:
            #     bsse.append(Items)

            # for pieda
            if pflag:
                pcount += 1
            if pflag and pcount > 2:
                # pieda.append(Items)
                pifrags.append(int(Items[0]))
                pjfrags.append(int(Items[1]))
                ess.append(float(Items[2]))
                exs.append(float(Items[3]))
                cts.append(float(Items[4]))
                dis.append(float(Items[5]))
                dqs.append(float(Items[6]))

#         for i in range(len(ifie)):
#             if float(ifie[i][4]) < -2:
#                 ifie[i][4] = 0.0
        # print('bsse', bsse[0])
        # file.close()

        return ifrags, jfrags, dists, hfs, mp2s, prs, grimmes, \
            pifrags, pjfrags, ess, exs, cts, dis, dqs
        # print ifie

    @staticmethod
    def getfraginfo(file, fragmode):
        # print('fragmode', fragmode)
        '''
        fraginfo = {
            'natoms': fnatoms,
            'baas': fbaas,
            'connects': fconnects
        }

        Frag.   Elec.   ATOM
            1      16       1     2     5     6     7     8
            2      30       3     4     9    10    13    14    15
            3      30      11    12    16    17    20    21    22
            4      30      18    19    23    24    27    28    29
            5      54      25    26    30    31    32    33    34    35    36    37
                           38

        ALL ELECTRON =    160
        ALL ATOM     =     38

        Frag.   Bonded Atom  Proj.
            2       2     3     3
            3      10    11     3
            4      17    18     3
            5      24    25     3
        '''

        if fragmode == 'auto':
            readflag = False
            autoreadflag = False
            baaflag = True
            fragdata = []
            fragdatas = []
            elecs = []
            seqnos = []
            fragnos = []
            residuestr = []
            baareadstart = False
            baafrags = []
            fconnects = []
            fprojs = []
            elec = 0

            while True:
                line = file.readline()
                items = line.strip().split()

                # chains = line[0]
                if len(items) == 0:
                    continue

                if items[0:3] == ['Frag.', 'Elec.', 'ATOM']:
                    readflag = True
                    continue

                if items[0:2] == ["ALL", "ELECTRON"]:
                    fragdatas.append(list(map(int, fragdata)))
                    readflag = False
                    baaflag = True

                if baaflag is True and items[0:2] == ["Frag.", "Bonded"]:
                    baareadstart = True
                    continue

                if baareadstart:
                    if items[0:2] == ["##", "READ"]:
                        break
                    baafrags.append(int(items[0]))
                    fconnects.append(list(map(int, items[1:3])))
                    fprojs.append(int(items[3]))

                if items[0:2] == ["ALL", "ATOM"]:
                    natom = int(items[3])

                if readflag:
                    if line[0:21] == "                     ":
                        fragdata = fragdata + items
                    else:
                        if len(fragdata) != 0:
                            fragdatas.append(list(map(int, fragdata)))
                            elecs.append(int(elec))
                        fragdata = []
                        elec = items[1]
                        fragdata = fragdata + items[2:]

                # if items [0:2] == ['START', 'FRAGMENT']:
                #     break

                # AUTOMATIC FRAGMENTATION
                if items[0:3] == ['Seq.', 'Frag.', 'Residue']:
                    autoreadflag = True
                    continue

                if autoreadflag and items[0:2] == ['The', 'system']:
                    autoreadflag = False
                    continue

                if autoreadflag and items[0] == 'Ions':
                    autoreadflag = False
                    continue

                if autoreadflag:
                    # print(items)
                    seqnos.append(items[0])
                    fragnos.append(items[1])
                    residuestr.append(items[2])

            nf = len(fragdatas)
            fbaas = [0] * nf
            fnatoms = []
            # get natoms
            for fdata in fragdatas:
                fnatoms.append(len(fdata))

            # get fbaas
            for baafrag in baafrags:
                fbaas[baafrag - 1] += 1
        return nf, fnatoms, fbaas, fconnects, fragdatas, natom

    @staticmethod
    def readpdb(fname):

        print('--- get pdbinfo ---')
        print('infile:',  fname)
        file = open(fname, 'r')
        molnames = []
        xcoords = []
        ycoords = []
        zcoords = []
        atypenames = []
        heads = []
        molnames = []
        atypenames = []
        labs = []
        chains = []
        resnums = []
        codes = []
        occs = []
        temps = []
        amarks = []
        charges = []
        nums = []
        atomcount = 0

        while True:
            line = file.readline()
            if not line:
                break
            itemlist = line.strip().split()
            if len(itemlist) < 3:
                continue

            if line[0:6] == 'HETATM' or line[0:4] == 'ATOM':
                atomcount += 1
                try:
                    num = int(line[6:12])
                except:
                    head = line[0:6]
                    num = line[6:11]
                    atypename = line[12:16]
                    lab = line[16]
                    res = line[17:20].strip()
                    chain = line[21]
                    resnum = line[22:26]
                    code = line[26]
                    xcoord = float(line[30:38].strip())
                    ycoord = float(line[38:46].strip())
                    zcoord = float(line[46:54].strip())
                    occ = line[54:60]
                    temp = line[60:66]
                    amark = line[76:78]
                    charge = line[78:80]

                if num < 100000:
                    head = line[0:6]
                    num = line[6:11]
                    atypename = line[12:16]
                    lab = line[16]
                    res = line[17:20].strip()
                    chain = line[21]
                    resnum = line[22:26]
                    code = line[26]
                    xcoord = float(line[30:38].strip())
                    ycoord = float(line[38:46].strip())
                    zcoord = float(line[46:54].strip())
                    occ = line[54:60]
                    temp = line[60:66]
                    amark = line[76:78]
                    charge = line[78:80]

                else:
                    head = line[0:6]
                    num = line[6:12]
                    atypename = line[13:17]
                    lab = line[17]
                    res = line[18:21].strip()
                    chain = line[22]
                    resnum = line[23:27]
                    code = line[27]
                    xcoord = float(line[31:39].strip())
                    ycoord = float(line[39:47].strip())
                    zcoord = float(line[47:55].strip())
                    occ = line[55:61]
                    temp = line[61:67]
                    amark = line[77:79]
                    charge = line[79:81]

                nums.append(num)
                heads.append(head)
                molnames.append(res)
                atypenames.append(atypename)
                labs.append(lab)
                chains.append(chain)
                resnums.append(resnum)
                codes.append(code)
                xcoords.append(xcoord)
                ycoords.append(ycoord)
                zcoords.append(zcoord)
                occs.append(occ)
                temps.append(temp)
                amarks.append(amark)
                charges.append(charge.replace('\n', ''))

            totalatom = atomcount

        return nums, heads, molnames, atypenames, labs, chains, \
            resnums, codes, xcoords, ycoords, zcoords, occs, temps, amarks, charges, totalatom

    @staticmethod
    def getnr(file):
        '''

        ## FRAGMENT NUCLEAR REPULSION

           NUCLEAR REPULSION ENERGY =          1889.4933472009


     =========================
        ## FMO TOTAL ENERGY
     =========================

           FMO2-HF
           Nuclear repulsion =          1889.4933472009
           Electronic energy =         -2985.0596593788
           Total      energy =         -1095.5663121780

           E(FMO2-MP2)       =            -1.0887272375
           Electronic energy =         -2986.1483866163
           Total      energy =         -1096.6550394155
        '''

        while True:
            Items = file.readline().strip().split()
            if len(Items) <= 2:
                continue
            if Items[0:2] == ['NUCLEAR', 'REPULSION']:
                # print('logMethod =', Items[2])
                NR = float(Items[4])
                return NR

    def readlog(self, logname, fragmode):
        # print('fragmode', fragmode)
        if fragmode == 'manual':
            text = open(logname, "r").readlines()
            for i in range(len(text)):
                itemList = text[i][:-1].split()
                # print(itemList)
                if len(itemList) < 2:
                    continue
                if itemList[:2]== ['NF', '=']:
                   nf = int(itemList[2])
                   break
        if fragmode == 'auto':
            f = open(logname, 'r')
            readflag = False
            autoreadflag = False
            fragdata = []
            fragdatas = []
            elecs = []
            seqnos = []
            fragnos = []
            residuestr = []
            logreadGeom = []
            pdbabs = ""

            for line in f:
                items = line[1:].split()
                chains = line[0]
                if len(items) == 0:
                    continue

                if items[0] == 'ReadGeom':
                     logreadGeom = items[2]

                if items[0] == 'AutoFrag':
                    if items[2] == 'ON':
                        self.fragmode = 'auto'
                    else:
                        self.fragmode = 'manual'

                # read frag table
                # if len(items) == 3:
                    # print (items)
                if items[0:3] == ['Frag.', 'Elec.', 'ATOM']:
                    readflag = True
                    # print('# readflag ON #')
                    continue
                if items[0:2] == ["ALL", "ELECTRON"]:
                    fragdatas.append(fragdata)
                    readflag = False
                if items[0:2] == ["ALL", "ATOM"]:
                    natom = int(items[3])
                if readflag == True:
                    if line[0:21] == "                     ":
                        # print(line)
                        fragdata = fragdata + items
                    else:
                        if len(fragdata) != 0:
                            fragdatas.append(fragdata)
                            elecs.append(int(elec))
                        fragdata = []
                        elec = items[1]
                        fragdata = fragdata + items[2:]

                if items [0:2] == ['START', 'FRAGMENT']:
                    break

                ## AUTOMATIC FRAGMENTATION
                if items[0:3] == ['Seq.', 'Frag.', 'Residue']:
                    autoreadflag = True
                    continue

                if autoreadflag == True and items[0:2] == ['The', 'system']:
                    autoreadflag = False
                    continue

                if autoreadflag == True and items[0] == 'Ions':
                    autoreadflag = False
                    continue

                if autoreadflag == True:
                   # print(items)
                   seqnos.append(items[0])
                   fragnos.append(items[1])
                   residuestr.append(items[2])

            nf = len(fragdatas)
            print('-------nf--------', nf)

        return natom

    def getlognatom(self, fname):
        acount = 0
        flag = False
        f =open(fname, "r")
        text = f.readlines()
        for i in range(len(text)):
            itemList = text[i].split()
            if itemList[0:4] == ['##', 'READ', 'MOLECULAR', 'STRUCTURE']:
                flag = True
                continue
            if flag:
                if itemList[0:3] == ['##', 'Molecular', 'formula']:
                    break
                elif len(itemList) <= 1:
                    continue
                else:
                    acount += 1
        return acount


    def getlogchg(self, fname, natom):
        chgs = []
        f =open(fname, "r")
        text = f.readlines()
        for i in range(len(text)):
            itemList = text[i].split()
            # MO
            if itemList == ['TWO-STAGE', 'RESP', 'FITTING:', 'SECOND', 'STAGE']:
                for j in range(int(natom)):
                    chgval = text[i+20+j].split()
                    chgs.append(float(chgval[2]))

            # FMO
            if itemList == ['##', 'ESP-FITTING', 'TYPE:', 'RESP']:
                for j in range(int(natom)):
                    chgval = text[i+6+j].split()
                    chgs.append(float(chgval[4]))
        return chgs


    def getlogchgall(self, fname, natom, chgtype):
        alabs = []
        elems = []
        ress = []
        frags = []
        chgs = []
        pops = []
        f =open(fname, "r")
        text = f.readlines()
        if chgtype == 'nbo':
            for i in range(len(text)):
                itemList = text[i].split()
                # FMO
                if itemList == ['##', 'NATURAL', 'ATOMIC', 'POPULATIONS']:
                    for j in range(int(natom)):
                        chgval = text[i+6+j].split()
                        alabs.append(int(chgval[0]))
                        elems.append(str(chgval[1]))
                        ress.append(int(chgval[2]))
                        frags.append(int(chgval[3]))
                        chgs.append(float(chgval[4]))
                        pops.append(float(chgval[5]))
                    chgdf=pd.DataFrame({'AtomLabel':alabs,
                                       'Elem': elems,
                                       'Res': ress,
                                       'Frag': frags,
                                       'Charge': chgs,
                                       'Pop': pops})
        else:
            print('Options except A are not supported yet.')
            sys.exit()
        return chgdf

        '''
         ========================================================
           ## NATURAL POPULATION ANALYSIS -- Ver.2.73(20131003)
         ========================================================


           ## NATURAL ATOMIC POPULATIONS

         -----------------------------------------
             Atom  Res Frag     Charge      Pop
                               FMO2-HF    FMO2-HF
         -----------------------------------------
             1 N     1    1  -0.801541   7.801541
             2 C     1    1  -0.108078   6.108078
             3 C     1    2   0.847514   5.152486
             4 O     1    2  -0.776786   8.776786
             5 C     1    1  -0.669410   6.669410
             6 H     1    1   0.461111   0.538889
        '''

