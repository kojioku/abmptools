import sys
import os
import gzip
scrdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scrdir)
try:
    import pandas as pd
except ImportError:
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
            charge_label = ["MUL-HF", "MUL-MP2", "NPA-HF", "NPA-MP2", "ESP-HF", "ESP-MP2"]
            DPM_label = ["DPM-HF-X", "DPM-HF-Y", "DPM-HF-Z", "DPM-MP2-X",
                         "DPM-MP2-Y", "DPM-MP2-Z"]
            monomer_label = ["NR", "HF", "MP2", "MP3"]
            dimer_label = ["NR", "HF", "ES", "MP2", "PR-MP2", "SCS-MP2(Grimme)",
                           "MP3", "SCS-MP3(MP2.5)", "HF-BSSE", "MP2-BSSE",
                           "SCS-MP2-BSSE", "MP3-BSSE", "SCS-MP3-BSSE",
                           "SOLV-ES", "SOLV-NP", "EX", "CT", "DQ"]
        '''

    def parse(self, filepath):
        # Initialize the data structures
        # natom = 0
        # nfrag = 0
        # atom_data = []

        print('start read', filepath)
        # filepath の拡張子が .gz なら gzip で読み込む
        if os.path.splitext(filepath)[-1] == '.gz':
            self.is_gz = True

        if self.is_gz:
            file = gzip.open(filepath, 'rt')
        else:
            file = open(filepath, 'rt')

        # get condition
        Method, ElecState, BasisSet, Laoc, Lptc, Ldimer, ReadGeom, fragmode,\
            is_npa, is_resp = self.getcondition(file)

        condition = {
            'aoc': Laoc,
            'basis_set': BasisSet,
            'calculation_method': Method,
            'electronic_state': ElecState,
            'ldimer': Ldimer,
            'ptc': Lptc,
            'readgeom': ReadGeom,
        }

        # setup label
        charge_label = ["MUL-HF", "MUL-MP2", "NPA-HF", "ESP-HF"]
        DPM_label = ["DPM-HF-X", "DPM-HF-Y", "DPM-HF-Z"]
        monomer_label = ["NR", "HF"]
        if Method == "MP2":
            monomer_label.append("MP2")
        dimer_label = ["NR", "HF", "ES"]
        if Method == "MP2":
            dimer_label += ["MP2", "PR-MP2", "SCS-MP2(Grimme)"]
        dimer_label += ["EX", "CT", "DQ"]

        labels = {
            'charge': charge_label,
            'DPM': DPM_label,
            'monomer': monomer_label,
            'dimer': dimer_label,
            }

        self.labels = labels

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
        # print('atominfo\n', atominfo)

        if Method == 'HF':
            hf = self.getmominfo(file, nf, Method)
        if Method == 'MP2':
            hf, mp2 = self.getmominfo(file, nf, Method)
        if Method == 'MP2D':
            print('MP2D is not supported yet')
            sys.exit()

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
        if Method == 'MP2':
            mominfo['MP2'] = mp2
        mominfo['DPM-HF-X'] = dpm_hf_x
        mominfo['DPM-HF-Y'] = dpm_hf_y
        mominfo['DPM-HF-Z'] = dpm_hf_z

        # print('mominfo\n', self.mominfo)

        if Method == 'MP2':
            ifrags, jfrags, dists, hfs, mp2s, prs, grimmes, \
                pifrags, pjfrags, ess, exs, cts, dis, dqs \
                = self.readifiepieda(file, Method)
        elif Method == 'HF':
            ifrags, jfrags, dists, hfs, \
                pifrags, pjfrags, ess, exs, cts, dis, dqs \
                = self.readifiepieda(file, Method)
        else:
            print("Methods other than MP2 or HF are not supported.")
            sys.exit()

        hartree = 627.5095
        ifieinfo = {
            "fragi": ifrags,
            "fragj": jfrags,
            "min-dist": dists,
            "NR": [0 for i in range(len(ifrags))],
            "HF": [x / hartree for x in hfs],
        }
        if Method == 'MP2':
            ifieinfo["MP2"] = [x / hartree for x in mp2s]
            ifieinfo["PR-MP2"] = [x / hartree for x in prs]
            ifieinfo["SCS-MP2(Grimme)"] = [x / hartree for x in grimmes]

        piedainfo = {
            "fragi": pifrags,
            "fragj": pjfrags,
            "ES": [x / hartree for x in ess],
            "EX": [x / hartree for x in exs],
            "CT": [x / hartree for x in cts],
            "DQ": dqs,
        }

        # for dim in dimer_label:
        #     diminfo[dim] = []

        ifiedf = pd.DataFrame(ifieinfo)
        piedadf = pd.DataFrame(piedainfo)
        print('ifiedf\n', ifiedf)
        print('piedadf\n', piedadf)

        # merge ifie and pieda
        ifpidf = pd.merge(ifiedf, piedadf, on=['fragi', 'fragj'], how='left') \
            .fillna(0.0)
        ifpidf = ifpidf[['fragi', 'fragj', 'min-dist'] + dimer_label]

        # get Mulliken
        mulalabs, mulelems, mulpops, mulchgs = \
            self.getlogmul(file, natom)
        atominfo['MUL-HF'] = mulchgs
        # print('Mulliken charges\n', atominfo['MUL-HF'])

        # get NPA val
        if is_npa:
            npaalabs, npaelems, nparess, npafrags, npachgs, npapops = \
                self.getlognpa(file, natom)
            atominfo['NPA-HF'] = npachgs
            # print('NPA charges\n', atominfo['NPA-HF'])
        # get RESP val
        if is_resp:
            espalabs, espelems, espress, espfrags, espchgs = \
                self.getlogresp(file, natom)
            atominfo['ESP-HF'] = espchgs

        total_electronic_energy, total_energy = \
            self.gettotalenergy(file, Method)
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

        return

    @staticmethod
    def gettotalenergy(file, Method):
        flag = False
        eneflag = False
        count = 0
        gcount = 0
        data = []
        if Method == 'MP2':
            tgtline = 8
        elif Method == 'HF':
            tgtline = 4
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
                if flag and count == tgtline:
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

        '''
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

    @staticmethod
    def getmominfo(file, nmom, Method):
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
                if Method == 'MP2':
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
        if Method == 'MP2':
            return hf, mp2
        if Method == 'HF':
            return hf
        '''
         =======================
            ## MONOMER ENERGY
         =======================

               Frag. No.  Monomer HF energy (E'(I))  Monomer MP2 energy
                      1      -79.0445787441             -0.0790837991
                      2     -189.3725896993             -0.1822649539
                      3     -189.3632919644             -0.1827496481
                      4     -189.3625086375             -0.1828280415
                      5     -389.4033075305             -0.3520136478
        '''

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
        is_npa = False
        is_resp = False
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
            if Items[0:2] == ['NBOANL', '=']:
                if Items[2] == 'ON':
                    is_npa = True
            if Items[0:2] == ['ESPTYP', '=']:
                if Items[2] == 'RESP':
                    is_resp = True
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

        return Method, ElecState, BasisSet, Laoc, Lptc, Ldimer, \
            ReadGeom, fragmode, is_npa, is_resp

    @staticmethod
    def readifiepieda(file, Method):
        ''' read ifie and pieda data from file

        read ifie and pieda data from file

        Args:
            fname (str): file name

        Returns:
            ifie (list): ifie data
            pieda (list): pieda data
        '''

        count = 0
        pcount = 0
        flag = False
        pflag = False
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
        print('--- get IFIE and PIEDA from log ---')
        while True:
            Items = file.readline().strip().split()
            if (pflag or flag) and len(Items) != 0:
                if Items[0] == '=================================':
                    break
            if len(Items) < 2:
                continue
            if Items[1] == 'MP2-IFIE' or Items[1] == 'HF-IFIE':
                flag = True
                continue
            if Items[1] == 'PIEDA':
                flag = False
                pflag = True
                continue
            if flag:
                count += 1
            if flag and count >= 3:
                ifrags.append(int(Items[0]))
                jfrags.append(int(Items[1]))
                dists.append(float(Items[2]))
                if float(Items[4]) < -2:
                    Items[4] = 0.0
                hfs.append(float(Items[4]))
                if Method == 'MP2':
                    mp2s.append(float(Items[5]))
                    prs.append(float(Items[6]))
                    grimmes.append(float(Items[7]))

            # for pieda
            if pflag:
                pcount += 1
            if pflag and pcount >= 3:
                # pieda.append(Items)
                pifrags.append(int(Items[0]))
                pjfrags.append(int(Items[1]))
                ess.append(float(Items[2]))
                exs.append(float(Items[3]))
                cts.append(float(Items[4]))
                dis.append(float(Items[5]))
                dqs.append(float(Items[6]))

        if Method == 'MP2':
            return ifrags, jfrags, dists, hfs, mp2s, prs, grimmes, \
                pifrags, pjfrags, ess, exs, cts, dis, dqs
        elif Method == 'HF':
            return ifrags, jfrags, dists, hfs, \
                pifrags, pjfrags, ess, exs, cts, dis, dqs

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


        '''

        while True:
            Items = file.readline().strip().split()
            if len(Items) <= 2:
                continue
            if Items[0:2] == ['NUCLEAR', 'REPULSION']:
                # print('logMethod =', Items[2])
                NR = float(Items[4])
                return NR

    @staticmethod
    def getlogresp(file, natom):
        alabs = []
        elems = []
        chgs = []
        ress = []
        frags = []
        flag = False
        count = 0
        gcount = 0

        while True:
            Items = file.readline().strip().split()
            # FMO
            if Items == ['##', 'ESP-FITTING', 'TYPE:', 'RESP']:
                flag = True
                continue
            if count == 5:
                alabs.append(int(Items[0]))
                elems.append(str(Items[1]))
                ress.append(int(Items[2]))
                frags.append(int(Items[3]))
                chgs.append(float(Items[4]))
                gcount += 1
                if gcount == natom:
                    break
                continue
            if flag:
                count += 1
                continue
        return alabs, elems, ress, frags, chgs

        # MO
        # if Items == ['TWO-STAGE', 'RESP', 'FITTING:', 'SECOND', 'STAGE']:
        #     for j in range(int(natom)):
        #         chgval = text[i+20+j].split()
        #         chgs.append(float(chgval[2]))

        '''
        ============================================================
          ## ELECTROSTATIC POTENTIAL FITTING -- ver.1.1 (20120905)
        ============================================================

          ## ESP-FITTING TYPE: RESP

        ------------------------------
            Atom  Res Frag     Charge
                              FMO2-HF
        ------------------------------
            1 N     1    1  -0.840813
            2 C     1    1  -0.007067
            3 C     1    2   0.624572
            4 O     1    2  -0.416961
            5 H     1    1   0.312765
            6 H     1    1   0.335686
            7 H     1    1   0.060326
        '''

    @staticmethod
    def getlogmul(file, natom):
        print('get mulliken charge from log')
        alabs = []
        elems = []
        chgs = []
        pops = []
        flag = False
        count = 0
        gcount = 0

        while True:
            Items = file.readline().strip().split()
            # FMO
            if Items == ['##', 'Mulliken', 'atomic', 'population']:
                flag = True
                continue
            if count == 4:
                alabs.append(int(Items[0]))
                elems.append(str(Items[1]))
                pops.append(float(Items[2]))
                chgs.append(float(Items[3]))
                gcount += 1
                if gcount == natom:
                    break
                continue
            if flag:
                count += 1
                continue
        return alabs, elems, pops, chgs

        '''
        =================================
          ## Mulliken atomic population
        =================================

         No. Atom   Atomic pop.  Net charge
                        FMO2         FMO2
           1   N      7.397224    -0.397224
           2   C      6.038224    -0.038224
           3   C      5.697646     0.302354
           4   O      8.310874    -0.310874
           5   H      0.822331     0.177669
        '''

    @staticmethod
    def getlognpa(file, natom):
        print('get npa charges from log')
        alabs = []
        elems = []
        ress = []
        frags = []
        chgs = []
        pops = []
        flag = False
        count = 0
        gcount = 0

        while True:
            Items = file.readline().strip().split()
            # FMO
            if Items == ['##', 'NATURAL', 'ATOMIC', 'POPULATIONS']:
                flag = True
                continue
            if count == 5:
                alabs.append(int(Items[0]))
                elems.append(str(Items[1]))
                ress.append(int(Items[2]))
                frags.append(int(Items[3]))
                chgs.append(float(Items[4]))
                pops.append(float(Items[5]))
                gcount += 1
                if gcount == natom:
                    break
                continue
            if flag:
                count += 1
                continue
        return alabs, elems, ress, frags, chgs, pops

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
