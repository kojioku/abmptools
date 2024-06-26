import math
import copy
import os
import gzip
import sys


class CPFManager:
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
        self.is_difie = False
        return

    def parse(self, filepath):
        # Initialize the data structures
        natom = 0
        nfrag = 0
        tgtfrag = self.tgtfrag

        print('start read', filepath)
        # filepath の拡張子が .gz なら gzip で読み込む
        if os.path.splitext(filepath)[-1] == '.gz':
            self.is_gz = True

        if self.is_gz:
            file = gzip.open(filepath, 'rt')
        else:
            file = open(filepath, 'rt')

        # Read header
        header, natom, nfrag, labels, cpfver = self.read_header(file)
        # Read atom data
        atominfo = self.read_atominfo(file, labels, natom, cpfver, tgtfrag)
        # Read fragment data
        fraginfo = self.read_fraginfo(file, nfrag, atominfo, cpfver, tgtfrag)
        # Read fragment distance data
        diminfo, dimaccept = self.read_dimdist(file, nfrag, tgtfrag)
        # Read dipole data
        mominfo = self.read_dipole(file, labels, nfrag, tgtfrag)
        # Read condition data
        condition = \
            self.read_condition(file)
        # Read static data
        static_data = \
            self.read_static(file, natom, nfrag, atominfo, tgtfrag)
        # Read monomer data
        mominfo = \
            self.read_monomer(file, mominfo, labels, natom,
                              tgtfrag, cpfver, nfrag)
        # Read dimer data
        diminfo, static_data = \
            self.read_dimer(file, diminfo, labels, nfrag,
                            cpfver, dimaccept, static_data, tgtfrag)
        import pandas as pd
        self.cpfver = cpfver
        self.atominfo = pd.DataFrame(atominfo)
        self.fraginfo = fraginfo
        self.condition = condition
        self.static_data = static_data
        self.mominfo = pd.DataFrame(mominfo)
        self.diminfo = pd.DataFrame(diminfo)
        self.labels = labels

        file.close()
        return self

    def write(self, header, filename, cpfver=23):
        '''Write the CPF file.

        Args:
            header (str): Header of the CPF file.
            filename (str): Filename of the CPF file.

        Returns:
            None
        '''

        # header section
        header = self.set_header(header, cpfver, self.static_data, self.labels)
        # atom section
        atomstr = self.set_atom(self.labels, self.atominfo)
        # fragment section
        fragstr = self.setupfragstr(self.fraginfo)
        # static section
        staticstr = self.set_static(self.static_data)
        # condition section
        conditionstr = self.set_condition(self.condition)
        # mom section
        dpmstr, momstr = self.set_mom(
            self.labels, self.mominfo, self.static_data)
        # dimer section
        diststr, dimstr = self.set_dimer(
            self.diminfo, self.labels, self.static_data, self.is_difie)

        trmstr = "{:>10}".format(self.static_data['ntrimer']) + '\n'
        tetrstr = "{:>10}".format(self.static_data['ntetramer']) + '\n'

        out = header + '\n' + atomstr + fragstr + diststr + dpmstr + \
            conditionstr + staticstr + momstr + \
            dimstr + trmstr + tetrstr + 'END'
        with open(filename, 'w') as f:
            f.write(out)
        print(filename + ' is created.')

        return None

    def info(self):
        import pprint

        pprint.pprint(dir(self))
        pprint.pprint(vars(self))

        return

    @staticmethod
    def read_header(file):
        header = file.readline().strip()
        ''' \
        CPF Open1.0 rev23
        CPF Open1.0 rev10
        CPF Ver.7.0 (MIZUHO 4.0)
        CPF Ver.4.201
        '''
        # Read the number of atoms and fragments
        natom, nfrag = map(int, file.readline().split())

        ''' \
        MUL-HF, MUL-MP2, NPA-HF, NPA-MP2, ESP-HF, ESP-MP2
        DPM-HF-X DPM-HF-Y DPM-HF-Z DPM-MP2-X DPM-MP2-Y DPM-MP2-Z
        NR, HF, MP2, MP3
        NR, HF, ES, MP2, PR-MP2, SCS-MP2(Grimme), MP3, SCS-MP3(MP2.5), \
            HF-BSSE, MP2-BSSE, SCS-MP2-BSSE, MP3-BSSE, SCS-MP3-BSSE, \
                SOLV-ES, SOLV-NP, EX, CT, DQ
        '''

        if header[0:17] == 'CPF Open1.0 rev23':
            cpfver = 23
        elif header[0:17] == 'CPF Open1.0 rev10':
            cpfver = 10
        elif header[0:24] == 'CPF Ver.7.0 (MIZUHO 4.0)':
            cpfver = 7.001
        elif header[0:13] == 'CPF Ver.4.201':
            cpfver = 4.201
        else:
            print('Unknown CPF version.')
            sys.exit()

        print('cpfver', cpfver)
        if cpfver == 23:
            charge_label = file.readline().strip().split()
            DPM_label = file.readline().strip().split()
            monomer_label = file.readline().strip().split()
            dimer_label = file.readline().strip().split()
        else:
            charge_label = ["MUL-HF", "MUL-MP2", "NPA-HF", "NPA-MP2", "ESP-HF", "ESP-MP2"]
            DPM_label = ["DPM-HF-X", "DPM-HF-Y", "DPM-HF-Z", "DPM-MP2-X",
                         "DPM-MP2-Y", "DPM-MP2-Z"]
            monomer_label = ["NR", "HF", "MP2", "MP3"]

        if cpfver == 10:
            dimer_label = ["NR", "HF", "ES", "MP2", "PR-MP2", "SCS-MP2(Grimme)",
                           "MP3", "SCS-MP3(MP2.5)", "HF-BSSE", "MP2-BSSE",
                           "SCS-MP2-BSSE", "MP3-BSSE", "SCS-MP3-BSSE",
                           "SOLV-ES", "SOLV-NP", "EX", "CT", "DQ"]

        if cpfver == 7.001:
            dimer_label = ["NR", "HF", "ES", "MP2", "SCS-MP2(Grimme)",
                           "MP3", "SCS-MP3(MP2.5)", "HF-BSSE", "MP2-BSSE",
                           "SCS-MP2-BSSE", "MP3-BSSE", "SCS-MP3-BSSE",
                           "SOLV-ES", "SOLV-NP", "EX", "CT", "DQ"]

        if cpfver == 4.201:
            dimer_label = ["NR", "HF", "ES", "MP2", "SCS-MP2(Grimme)",
                           "MP3", "SCS-MP3(MP2.5)", "HF-BSSE", "MP2-BSSE",
                           "SCS-MP2-BSSE", "MP3-BSSE", "SCS-MP3-BSSE",
                           "EX", "CT", "DQ"]

        labels = {
            'charge': charge_label,
            'DPM': DPM_label,
            'monomer': monomer_label,
            'dimer': dimer_label,
            }

        return header, natom, nfrag, labels, cpfver

    @staticmethod
    def read_atominfo(file, labels, natom, cpfver, tgtfrag):

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

        # Initialize the data structures
        for chg in labels['charge']:
            atominfo[chg] = []

        # Read the atom data
        print('atom section')
        for _ in range(natom):
            atom_data = file.readline().rstrip()
            if cpfver == 23:
                if tgtfrag != 0 and \
                        int(atom_data[34:44].strip()) not in tgtfrag:
                    continue
                atominfo['alabels'].append(int(atom_data[0:10].strip()))
                atominfo['elems'].append(atom_data[11:13].strip())
                atominfo['elemtypes'].append(atom_data[14:18])
                atominfo['resnames'].append(atom_data[19:22].strip())
                atominfo['resids'].append(int(atom_data[23:33].strip()))
                atominfo['fragids'].append(int(atom_data[34:44].strip()))
                atominfo['xcoords'].append(float(atom_data[45:65].strip()))
                atominfo['ycoords'].append(float(atom_data[65:85].strip()))
                atominfo['zcoords'].append(float(atom_data[85:105].strip()))
                atominfo['chainids'].append(atom_data[106:109].strip())
                atominfo['optflags'].append(atom_data[110:111].strip())
                for chg in labels['charge']:
                    lstart = labels['charge'].index(chg)*20 + 111
                    lend = lstart + 20
                    atominfo[chg].append(float(atom_data[lstart:lend].strip()))

            elif cpfver in [10, 7.001, 4.201]:
                if tgtfrag != 0 and \
                        int(atom_data[23:27].strip()) not in tgtfrag:
                    continue
                atominfo['alabels'].append(int(atom_data[0:5].strip()))
                atominfo['elems'].append(atom_data[6:8].strip())
                atominfo['elemtypes'].append(atom_data[9:13])
                atominfo['resnames'].append(atom_data[14:17].strip())
                atominfo['resids'].append(int(atom_data[18:22].strip()))
                atominfo['fragids'].append(int(atom_data[23:27].strip()))
                atominfo['xcoords'].append(float(atom_data[28:40].strip()))
                atominfo['ycoords'].append(float(atom_data[40:52].strip()))
                atominfo['zcoords'].append(float(atom_data[52:64].strip()))
                for chg in labels['charge']:
                    lstart = labels['charge'].index(chg)*12 + 64
                    lend = lstart + 12
                    atominfo[chg].append(float(atom_data[lstart:lend].strip()))
                atominfo['chainids'].append(atom_data[137:138].strip())
                atominfo['optflags'].append("")

        return atominfo

        ''' \
         1 N   N   GLY          1          1         0.1620000000\
            -0.2020000000        0.0000000000     1       -0.3972241421
        '''
        '''"原子数"分ループ{
                原子の番号(i10), 元素記号(a2), x, \
                原子タイプ(a4), x, 残基名(a3), x, 残基番号(i10), x,フラグメント番号(i10), x,\
                x座標(f20.10), y座標(f20.10) z座標(f20.10), x, Chain ID(a3), x, \
                構造最適化オプション(i1), !Mulliken原子電荷(HF密度)(f20.10), \
                !Mulliken原子電荷(MP2密度)(f20.10), !NPA原子電荷(HF密度)(f20.10), \
                !NPA原子電荷(MP2密度)(f20.10), !ESP原子電荷(HF密度)(f20.10), \
                !ESP原子電荷(MP2密度)(f20.10)
                }
        '''

    @staticmethod
    def read_fraginfo(file, nfrag, atominfo, cpfver, tgtfrag):
        '''
        16      30      30      30      54
         0       1       1       1       1
         2         3
        10        11
        17        18
        24        25
        フラグメント数/10行ループ {フラグメントごとの電子数(10i8)}
        フラグメント数/10行ループ {
            結合電子を自分のフラグメントに割り当てたフラグメント間の結合の本数(フラグメントごとのBDA数)(10i8)
            }
        BDA総数行ループ{
            フラグメントが結合している相手の原子番号(i10), 相手フラグメントと結合している原子の番号(i10)
            }
        '''

        print('fragment section')
        fnatoms, fbaas, fconnects = CPFManager.readfragcpf(
            file, nfrag, tgtfrag, atominfo['alabels'], cpfver)

        # cut the fragment data if tgtfrag is specified
        if tgtfrag != 0:
            fnatoms = copy.deepcopy(fnatoms[:len(tgtfrag)])
            fbaas = copy.deepcopy(fbaas[:len(tgtfrag)])

        # print(fconnects)
        fraginfo = {
            'natoms': fnatoms,
            'baas': fbaas,
            'connects': fconnects
        }

        return fraginfo

    @staticmethod
    def readfragcpf(file, nf, tgtfrag=0, alabels=None, cpfver=23):
        '''Read the fragment data from the CPF file.

        Args:
            file (file): The CPF file to read from
            nf (int): The number of fragments in the CPF file

        Returns:
            fnatoms (list): The number of atoms in each fragment
            fbaas (list): The BAA of each fragment
            fconnects (list): The connect (bda baa) atom IDs of each fragment

        Note:
            The CPF file format is described in the CPF manual.

            fnatoms = [natom1, natom2, ...]
            fbaas = [baa1, baa2, ...]
            fconnects = [[atomid1, atomid2, ...], [atomid1, atomid2, ...], ...]
        '''

        # Initialize the data structures
        flag = True
        nline = 1  # line number
        fnatoms = []  # fragment natom section
        fbaas = []  # bda section
        fconnects = []  # bda-baa atom section
        typcount = 0   # type (0: fnatoms, 1: fbaas, 2: fconnects)
        lcount = 0  # line count

        if cpfver == 23 and nf > 10:
            nline = math.ceil(nf/10)
            # print('n_line', nline)
        if cpfver in [4.201, 7.001, 10] and nf > 16:
            nline = math.ceil(nf/16)

        count = 0
        while True:
            itemlist = file.readline().strip().split()
            if len(itemlist) == 0:
                continue
            # fragment natom section
            if flag is True and typcount == 0:
                fnatoms.append(itemlist)
                lcount += 1
                if lcount == nline:
                    typcount += 1
                    lcount = 0
                    fnatoms = CPFManager.flatten(fnatoms)
                    # print(fnatoms)
                    continue

            # bda section
            if flag is True and typcount == 1:
                fbaas.append(itemlist)
                lcount += 1
                if lcount == nline:
                    typcount += 1
                    lcount = 0
                    fbaas = CPFManager.flatten(fbaas)
                    # print(fbaas)
                    continue

            # bda-baa atom section
            if flag is True and typcount == 2:
                count += 1

                def split_integer(itemlist):
                    # itemlistが整数であることを確認
                    if len(itemlist) == 1:
                        if 100000000 <= int(itemlist[0]) <= 999999999:
                            # print('length = 9')
                            length = 4
                        if 1000000000 <= int(itemlist[0]) <= 9999999999:
                            # print('length = 10')
                            length = 5
                        itemlist_str = str(itemlist[0])
                        first_part = int(itemlist_str[:length])  # 最初の4桁
                        second_part = int(itemlist_str[length:])  # 残りの5桁
                        return [first_part, second_part]
                    else:
                        return itemlist

                itemlist = split_integer(itemlist)
                # print(itemlist)
                if count == sum(fbaas):
                    if tgtfrag != 0:
                        if int(itemlist[0]) not in alabels or int(itemlist[1]) not in alabels:
                            print('break!')
                            break
                    fconnects.append(itemlist)
                    fconnects = CPFManager.functor(int, fconnects)
                    break
                if tgtfrag != 0:
                    # print(alabels)
                    if int(itemlist[0]) not in alabels or int(itemlist[1]) not in alabels:
                        # print(itemlist, 'skip!!!!')
                        continue
                fconnects.append(itemlist)
                fconnects = CPFManager.functor(int, fconnects)

                # datas.append(itemlist)

        return fnatoms, fbaas, fconnects

    @staticmethod
    def read_dimdist(file, nfrag, tgtfrag):
        '''
            2           1  0.000000000000000E+000
            3           1   6.43498504089706
            3           2  0.000000000000000E+000
            4           1   6.42318412936426
            4           2   5.20966306910324
            4           3  0.000000000000000E+000
            5           1   6.69353320455502
            5           2   3.97902962070384
            5           3   5.20963256562360
            5           4  0.000000000000000E+000

            ダイマー数行ループ{
                フラグメントIの番号,フラグメントJの番号,フラグメントIJ間の最短原子間距離
            }
        '''

        print('dimer distance section')
        # read the dimer info
        diminfo = {
            'fragi': [],
            'fragj': [],
            'min-dist': []
            }

        dimaccept = []
        for _ in range(nfrag*(nfrag-1)//2):
            dim_data = file.readline().strip().split()
            if tgtfrag != 0:
                if int(dim_data[0]) not in tgtfrag or int(dim_data[1]) not in tgtfrag:
                    continue
            diminfo['fragi'].append(int(dim_data[0]))
            diminfo['fragj'].append(int(dim_data[1]))
            diminfo['min-dist'].append(float(dim_data[2]))
            dimaccept.append(_)

        return diminfo, dimaccept

    @staticmethod
    def read_dipole(file, labels, nfrag, tgtfrag):
        # read the dipole data

        print('dipole moment section')
        mominfo = {
            'fragi': []
            }

        # Initialize the data structures
        for dpm in labels['DPM']:
            mominfo[dpm] = []

        for _ in range(nfrag):
            dipole_data = file.readline().strip().split()
            if tgtfrag != 0:
                if (_+1) not in tgtfrag:
                    continue
            count = 0
            for dpm in labels['DPM']:
                mominfo[dpm].append(float(dipole_data[count]))
                count += 1

        return mominfo

        '''
        0.785970969662152E+00    0.281296712643063E+01   -0.132604581524068E+01
        0.145048892092187E+01   -0.334817286212412E+01   -0.161464626780787E+01
        0.524704704075969E+00   -0.357367905308988E+01    0.156523382974318E+01
        0.245666859488253E+01   -0.160468359833147E+01    0.280152739221378E+01
        0.368656733092432E+01   -0.478966423463728E+01   -0.248675894502097E+00

            モノマー数行ループ {
                !双極子モーメント(HF密度)のx成分,!y成分,!z成分,!双極子モーメント(MP2密度)のx成分,!y成分,!z成分
            }
        '''

    @staticmethod
    def read_condition(file):
        '''
        STO-3G
        S1
        HF
          0.000000000000000E+000   2.00000000000000        2.00000000000000
           1889.49334720085
          -2985.05965937836
          -1095.56631217751

            基底関数名(a20)
            電子状態(a4)
            計算方法(a20)
            AOpopulation近似のパラメータ,点電荷近似のパラメータ,dimer-es近似のパラメータ
            核間反発エネルギー
            全電子エネルギー
            全エネルギー
        '''

        # Read the basis set, electronic state, and calculation method
        print('condition and static data section')
        condition = {}
        condition['basis_set'] = file.readline().strip()
        condition['electronic_state'] = file.readline().strip()
        condition['calculation_method'] = file.readline().strip()

        # Read the parameters for the AOpopulation, ptc, and ldimer
        data_approxy = file.readline().strip().split()
        condition['aoc'] = float(data_approxy[0])
        condition['ptc'] = float(data_approxy[1])
        condition['ldimer'] = float(data_approxy[2])

        return condition

    @staticmethod
    def read_static(file, natom, nfrag, atominfo, tgtfrag):
        static_data = {}
        static_data['nuclear_repulsion_energy'] = float(file.readline().strip())
        static_data['total_electronic_energy'] = float(file.readline().strip())
        static_data['total_energy'] = float(file.readline().strip())
        if tgtfrag != 0:
            static_data['natom'] = len(atominfo['alabels'])
            static_data['nfrag'] = len(tgtfrag)
        else:
            static_data['natom'] = natom
            static_data['nfrag'] = nfrag

        return static_data

    @staticmethod
    def read_monomer(file, mominfo, labels, natom, tgtfrag, cpfver, nfrag):
        '''
         5
         1   0.321733218676314E+02  -0.111217900611777E+03
         2   0.105581623892680E+03  -0.294954213591969E+03
         3   0.105580760608816E+03  -0.294944052573276E+03
         4   0.105592376341217E+03  -0.294954884978768E+03
         5   0.302972088661468E+03  -0.692375396192077E+03
        '''

        print('monomer section')
        # Initialize the data structures
        for mom in labels['monomer']:
            mominfo[mom] = []

        if cpfver == 23:
            nmonomer = int(file.readline().strip())
            # if self.tgtfrag != 0:
            #     nmonomer = static_data['nfrag']
            for _ in range(nmonomer):
                monomer_data = file.readline().strip().split()
                if tgtfrag != 0:
                    if (_ + 1) not in tgtfrag:
                        continue
                count = 1
                # print(monomer_data)
                mominfo['fragi'].append(int(monomer_data[0]))
                for mom in labels['monomer']:
                    mominfo[mom].append(float(monomer_data[count]))
                    count += 1

        else:
            nmonomer = nfrag
            icount = 1
            for _ in range(nmonomer):
                monomer_data = file.readline().strip().split()
                if tgtfrag != 0:
                    if (_+1) not in tgtfrag:
                        continue
                mominfo['fragi'].append(icount)
                icount += 1

                count = 0
                for mom in labels['monomer']:
                    mominfo[mom].append(float(monomer_data[count]))
                    count += 1
        # print(mominfo)
        return mominfo

    @staticmethod
    def read_dimer(file, diminfo, labels, nfrag,
                   cpfver, dimaccept, static_data, tgtfrag):

        print('dimer section')
        # Initialize the data structures
        for dim in labels['dimer']:
            diminfo[dim] = []

        if cpfver == 23:
            ndimer = int(file.readline().strip())
        else:
            ndimer = nfrag*(nfrag-1)//2
        static_data['ndimer'] = ndimer

        dimaccept_set = set(dimaccept)
        for _ in range(ndimer):
            if cpfver == 23:
                count = 2
            else:
                count = 0
            dimer_data = file.readline().strip().split()
            if tgtfrag != 0:
                if _ not in dimaccept_set:
                    continue
            for dim in labels['dimer']:
                diminfo[dim].append(float(dimer_data[count]))
                count += 1
        static_data['ntrimer'] = "0"
        static_data['ntetramer'] = "0"

        if tgtfrag != 0:
            static_data['ndimer'] = len(dimaccept)

        return diminfo, static_data

        '''
        10
         2         1
            0.898358210380181E+02  -0.104597690503434E+03  -0.144643244417453E+02
            -0.170246083060746E+00  -0.127298940609776E+00   0.100810192882720E+02
         3         1
            0.481151423108853E+02  -0.481122012096812E+02   0.292560697485555E-02
            0.199376860194889E-04  -0.444345681671621E-05   0.257150081637292E-04
         3         2
            0.151596219157517E+03  -0.166354359041009E+03  -0.144607605291407E+02
            -0.171095817875511E+00  -0.126283536475938E+00   0.100219172897770E+02
         0
         0

        モノマーの数(i10)
        モノマー数分ループ {
        !モノマーの番号(i10), !モノマーの核間反発エネルギー(e24.15), !モノマーの電子エネルギー(e24.15)
        !モノマーのMP2相関エネルギー(e24.15), !モノマーのMP3相関エネルギー(e24.15)
        }
        ダイマーの数(i10)
        ダイマー数分ループ {I(I10) J(I10), !IFIEの核間反発エネルギー項(e24.15),
            !HF-IFIEの電子エネルギー項(e24.15), !HF-IFIEの静電エネルギー項(e24.15),
            !MP2-IFIE (e24.15) PR-MP2-IFIE (e24.15),
            !SCS-MP2(Grimme)-IFIE (e24.15),
            !MP3-IFIE(e24.15), !SCS-MP3(MP2.5)-IFIE(e24.15),
            !HF-IFIE-BSSE(e24.15),
            !MP2-IFIE-BSSE (e24.15) SCS-MP2(Grimme)-IFIE-BSSE (e24.15),
            !MP3-IFIE-BSSE (e24.15)SCS-MP3(MP2.5)-IFIE-BSSE (e24.15),
            !溶媒による遮蔽効果のES項 (e24.15) 溶媒による遮蔽効果のNP項 (e24.15),
            !PIEDAのEX項 (e24.15), !PIEDAのCT項 (e24.15), !PIEDAのΔq(e24.15)
        }
        トリマーの数(I10)
        トリマー数分ループ {I(I10), J(I10), K(I10), HF(e24.15), MP2(e24.15),
            SCS-MP2(e24.15), MP3(e24.15), SCS-MP3(e24.15)
        }
        テトラマーの数(I10),
        テトラマー数分ループ {I(I10), J(I10), K(I10), L(I10),
            HF(e24.15), MP2(e24.15), SCS-MP2(e24.15),
            MP3(e24.15), SCS-MP3(e24.15)
        }
        END(a3)
        '''

    @staticmethod
    def set_header(header, cpfver, static_data, labels):
        if cpfver == 23:
            header = 'CPF Open1.0 rev23 ' + header
        if cpfver == 10:
            header = 'CPF Open1.0 rev10 ' + header
        if cpfver == 4:
            header = 'CPF Ver.4.201 ' + header

        header += '\n' + '{:>10}'.format(static_data['natom']) \
            + '{:>10}'.format(static_data['nfrag']) + '\n'

        if cpfver == 23:
            header += ' '.join(labels['charge']) + '\n'
            header += ' '.join(labels['DPM']) + '\n'
            header += ' '.join(labels['monomer']) + '\n'
            header += ' '.join(labels['dimer'])
        return header

    @staticmethod
    def set_atom(labels, atominfo):
        # setup format dict
        formats = {
            'alabels': "{:>10}",  # 原子の番号(i10)
            'elems': "{:>2}",  # 元素記号(a2)
            'elemtypes': "{:<4}",  # 原子タイプ(a4)
            'resnames': "{:>3}",  # 残基名(a3)
            'resids': "{:>10}",  # 残基番号(i10)
            'fragids': "{:>10}",  # フラグメント番号(i10)
            'xcoords': "{:20.10f}",  # x座標(f20.10)
            'ycoords': "{:19.10f}",  # y座標(f20.10)
            'zcoords': "{:19.10f}",  # z座標(f20.10)
            'chainids': "{:>3}",  # Chain ID(a3)
            'optflags': "{:>1}",  # 構造最適化オプション(i1)
        }

        for chglabel in labels['charge']:
            formats[chglabel] = "{:19.10f}"

        # update atominfodf using fmt
        atominfo = copy.deepcopy(atominfo)
        for col, fmt in formats.items():
            atominfo[col] = atominfo[col].apply(lambda x: fmt.format(x))

        atomstr = ''
        for index, row in atominfo.iterrows():
            atomstr += ' '.join(row.astype(str).tolist()) + '\n'

        return atomstr

    @staticmethod
    def set_static(static_data):
        # static section
        staticstr = ''
        staticstr += str(static_data['nuclear_repulsion_energy']) + '\n'
        staticstr += str(static_data['total_electronic_energy']) + '\n'
        staticstr += str(static_data['total_energy']) + '\n'

        '''
            static_data = {}
            static_data['nuclear_repulsion_energy'] = float(file.readline().strip()) -
            static_data['total_electronic_energy'] = float(file.readline().strip()) -
            static_data['total_energy'] = float(file.readline().strip()) -
            static_data['natom'] = natom
            static_data['nfrag'] = nfrag
        '''

        return staticstr

    @staticmethod
    def set_condition(condition):
        # condition
        conditionstr = ''
        conditionstr += "{:<20}".format(condition['basis_set']) + '\n'
        conditionstr += "{:<4}".format(condition['electronic_state']) + '\n'
        conditionstr += "{:<20}".format(condition['calculation_method']) + '\n'
        conditionstr += "{:24.15e}".format(condition['aoc']) + \
            "{:24.15e}".format(condition['ptc']) + \
            "{:24.15e}".format(condition['ldimer']) + '\n'

        '''
        self.condition
            condition = {}
            condition['basis_set'] = file.readline().strip() a20
            condition['electronic_state'] = file.readline().strip() a4
            condition['calculation_method'] = file.readline().strip() a20
            condition['aoc'] = float(data_approxy[0]) -
            condition['ptc'] = float(data_approxy[1]) -
            condition['ldimer'] = float(data_approxy[2]) -
        '''

        return conditionstr

    @staticmethod
    def set_mom(labels, mominfo, static_data):
        formats = {
            'fragi': "{:>10}",  # フラグメント番号(i10)
        }

        for dpmlabel in labels['DPM']:
            formats[dpmlabel] = "{:24.15e}"

        for momlabel in labels['monomer']:
            formats[momlabel] = "{:24.15e}"

        # update mominfodf using fmt
        mominfo = copy.deepcopy(mominfo)
        for col, fmt in formats.items():
            mominfo[col] = mominfo[col].apply(lambda x: fmt.format(x))

        # dpmlist に含まれるカラムのみを選択
        selected_columns = [col for col in mominfo.columns if col in labels['DPM']]

        dpmstr = ''
        # 選択したカラムのデータをスペース区切りで出力
        for index, row in mominfo[selected_columns].iterrows():
            dpmstr += ''.join(row.astype(str).tolist()) + '\n'

        selected_columns = [col for col in mominfo.columns if col in labels['monomer']]
        selected_columns.insert(0, 'fragi')

        momstr = "{:>10}".format(static_data['nfrag']) + '\n'
        for index, row in mominfo[selected_columns].iterrows():
            # print ('index', index, 'row', row)
            momstr += ''.join(row.astype(str).tolist()) + '\n'

        return dpmstr, momstr
        '''
                labels = {
                    'charge': charge_label,
                    'DPM': DPM_label,
                    'monomer': monomer_label,
                    'dimer': dimer_label,
                    }

        self.mominfo
            # read the dipole data
            0.785970969662152E+00    0.281296712643063E+01   -0.132604581524068E+01
            0.145048892092187E+01   -0.334817286212412E+01   -0.161464626780787E+01
           -0.524704704075969E+00   -0.357367905308988E+01    0.156523382974318E+01
            0.245666859488253E+01   -0.160468359833147E+01    0.280152739221378E+01
            0.368656733092432E+01   -0.478966423463728E+01   -0.248675894502097E+00
        '''

    @staticmethod
    def set_dimer(diminfo, labels, static_data, is_difie):
        formats = {
            'fragi': "{:>10}",  # フラグメント番号(i10)
            'fragj': "{:>10}",  # フラグメント番号(i10)
            'min-dist': "{:24.15e}",  # 最小距離(24.15e)
        }

        if is_difie:
            formats['M-min-dist'] = "{:24.15e}"  # 平均最小距離(24.15e)
            formats['S-min-dist'] = "{:24.15e}"  # 平均最小距離(24.15e)

        diminfo = copy.deepcopy(diminfo)
        for dimlabel in labels['dimer']:
            formats[dimlabel] = "{:24.15e}"

        # update diminfodf using fmt
        for col, fmt in formats.items():
            diminfo[col] = diminfo[col].apply(lambda x: fmt.format(x))

        diststr = ''
        # standard ifie case
        if not is_difie:
            print(static_data['ndimer'])
            for i in range(static_data['ndimer']):
                diststr += str(diminfo['fragi'][i]) + str(diminfo['fragj'][i]) \
                    + str(diminfo['min-dist'][i]) + '\n'
                # print(diststr)

        # difie case
        else:
            for i in range(static_data['ndimer']):
                diststr += str(diminfo['fragi'][i]) + str(diminfo['fragj'][i]) \
                    + str(diminfo['min-dist'][i]) + str(diminfo['M-min-dist'][i]) \
                    + str(diminfo['S-min-dist'][i]) + '\n'

        # dimer に含まれるカラムのみを選択
        selected_columns = [col for col in diminfo.columns if col in labels['dimer']]
        selected_columns.insert(0, 'fragj')
        selected_columns.insert(0, 'fragi')

        dimstr = "{:>10}".format(static_data['ndimer']) + '\n'
        # 選択したカラムのデータをスペース区切りで出力
        for index, row in diminfo[selected_columns].iterrows():
            dimstr += ''.join(row.astype(str).tolist()) + '\n'

        return diststr, dimstr

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
            MUL-HF, MUL-MP2, NPA-HF, NPA-MP2, ESP-HF, ESP-MP2


        diminfo
            fragi, fragj, min-dist
            NR, HF, ES, MP2, PR-MP2, SCS-MP2(Grimme), MP3, SCS-MP3(MP2.5), \
                HF-BSSE, MP2-BSSE, SCS-MP2-BSSE, MP3-BSSE, SCS-MP3-BSSE, \
                    SOLV-ES, SOLV-NP, EX, CT, DQ

        mominfo
            fragi, DPM-HF-X DPM-HF-Y DPM-HF-Z DPM-MP2-X DPM-MP2-Y DPM-MP2-Z
            NR, HF, MP2, MP3

                モノマー数行ループ {
                    !双極子モーメント(HF密度)のx成分,!y成分,!z成分,!双極子モーメント(MP2密度)のx成分,!y成分,!z成分
                }

        # self.diminfoの内容を'fragi', 'fragj', 'min-dist'の順に出力
            diminfo = {
                'fragi': [],        space
                'fragj': [],        space
                'min-dist': []
                }
    '''

    @staticmethod
    def setupfragstr(fraginfo):
        fragstr = ''
        count = 0
        for fnatom in fraginfo['natoms']:
            fragstr += '{:>8}'.format(fnatom)
            count += 1
            if count == 10:
                fragstr += '\n'
                count = 0
        if count != 0:
            fragstr += '\n'

        count = 0
        for fnatom in fraginfo['baas']:
            fragstr += '{:>8}'.format(fnatom)
            count += 1
            if count == 10:
                fragstr += '\n'
                count = 0
        if count != 0:
            fragstr += '\n'

        for fnatom in fraginfo['connects']:
            # print(fnatom)
            fragstr += '{:>10}'.format(fnatom[0]) + '{:>10}'.format(fnatom[1]) + '\n'

        return fragstr

        '''
            fraginfo = {
                'natoms': fnatoms,      10i8
                'baas': fbaas,          10i8
                'connects': fconnects   i10
            }
        '''

    @staticmethod
    def flatten(nested_list):
        '''2重のリストをフラットにする関数

        Args:
            nested_list (list): 2重のリスト

        Returns:
            list: フラットにしたリスト
        '''

        return [int(e) for inner_list in nested_list for e in inner_list]

    @staticmethod
    def functor(f, lname):
        '''リストの中身を関数で変換する関数

        Args:
            f (function): 適用する関数
            lname (list): リスト

        Returns:
            list: 適用後のリスト
        '''

        if isinstance(lname, list):
            return [CPFManager.functor(f, i) for i in lname]
        else:
            return f(lname)

    # def listtoint(nested_list):
    #     """リストの中をintegerに"""
    #     return [int(e) for inner_list in nested_list for e in inner_list]

    @staticmethod
    def selectfrag(item1):
        if type(item1) == str:
            if '-' in item1:
                tgt = item1.split('-')
                # print('tgt', tgt)
                tgtlist = [i for i in range(int(tgt[0]), int(tgt[1]) + 1)]
            else:
                tgtlist = list(map(int, item1.split(',')))
            return tgtlist

        elif type(item1) == int:
            return int(item1)
