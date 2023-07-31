import math


def flatten(nested_list):
    '''2重のリストをフラットにする関数

    Args:
        nested_list (list): 2重のリスト

    Returns:
        list: フラットにしたリスト
    '''

    return [int(e) for inner_list in nested_list for e in inner_list]


def functor(f, lname):
    '''リストの中身を関数で変換する関数

    Args:
        f (function): 適用する関数
        lname (list): リスト

    Returns:
        list: 適用後のリスト
    '''

    if isinstance(lname, list):
        return [functor(f, i) for i in lname]
    else:
        return f(lname)

# def listtoint(nested_list):
#     """リストの中をintegerに"""
#     return [int(e) for inner_list in nested_list for e in inner_list]


def readfragcpf(file, nf):
    '''Read the fragment data from the CPF file.

    Args:
        file (file): The CPF file to read from
        nf (int): The number of fragments in the CPF file

    Returns:
        fnatoms (list): The number of atoms in each fragment
        fbaas (list): The BAA of each fragment
        connects (list): The connect (bda baa) atom IDs of each fragment

    Note:
        The CPF file format is described in the CPF manual.

        fnatoms = [natom1, natom2, ...]
        fbaas = [baa1, baa2, ...]
        connects = [[atomid1, atomid2, ...], [atomid1, atomid2, ...], ...]
    '''

    # Initialize the data structures
    flag = True
    nline = 1  # line number
    fnatoms = []  # fragment natom section
    fbaas = []  # bda section
    connects = []  # bda-baa atom section
    typcount = 0   # type (0: fnatoms, 1: fbaas, 2: connects)
    lcount = 0  # line count

    if nf > 10:
        nline = math.ceil(nf/10)
        print('n_line', nline)

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
                fnatoms = flatten(fnatoms)
                continue

        # bda section
        if flag is True and typcount == 1:
            fbaas.append(itemlist)
            lcount += 1
            if lcount == nline:
                typcount += 1
                lcount = 0
                fbaas = flatten(fbaas)
                continue

        # bda-baa atom section
        if flag is True and typcount == 2:
            count += 1
            connects.append(itemlist)
            connects = functor(int, connects)
            if count == sum(fbaas):
                break

            # datas.append(itemlist)

    return fnatoms, fbaas, connects


def read_data(filepath):
    # Initialize the data structures
    natom = 0
    nfrag = 0
    # atom_data = []
    atominfo = {
        'alabels': [],
        'elems': [],
        'elemtype': [],
        'resnames': [],
        'resids': [],
        'fragids': [],
        'xcoords': [],
        'ycoords': [],
        'zcoords': [],
        'chainids': [],
        'optflags': [],
        }

    with open(filepath, 'r') as file:
        # Read the header
        header = file.readline().strip()
        ''' \
        CPF Open1.0 rev23
        CPF Open1.0 rev10
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

        if header[0:18] == 'CPF Open1.0 rev23':
            cpfver = 23
        if header[0:18] == 'CPF Open1.0 rev10':
            cpfver = 10

        if cpfver == 23:
            charge_label = file.readline().strip().split()
            DPM_label = file.readline().strip().split()
            monomer_label = file.readline().strip().split()
            dimer_label = file.readline().strip().split()

        # Read the atom data
        ''' \
         1 N   N   GLY          1          1         0.1620000000\
            -0.2020000000        0.0000000000     1       -0.3972241421
        '''
        '''"原子数"分ループ{
                原子の番号(i10), x, 元素記号(a2), x, \
                原子タイプ(a4), x, 残基名(a3), x, 残基番号(i10), x,フラグメント番号(i10), x,\
                x座標(f20.10), y座標(f20.10) z座標(f20.10), x, Chain ID(a3), x, \
                構造最適化オプション(i1), !Mulliken原子電荷(HF密度)(f20.10), \
                !Mulliken原子電荷(MP2密度)(f20.10), !NPA原子電荷(HF密度)(f20.10), \
                !NPA原子電荷(MP2密度)(f20.10), !ESP原子電荷(HF密度)(f20.10), \
                !ESP原子電荷(MP2密度)(f20.10)
                }
        '''

        if cpfver == 23:
            # Initialize the data structures
            for chg in charge_label:
                atominfo[chg] = []

        # Read the atom data
        for _ in range(natom):
            atom_data = file.readline().rstrip()
            atominfo['alabels'].append(atom_data[0:10].strip())
            atominfo['elems'].append(atom_data[11:13].strip())
            atominfo['elemtype'].append(atom_data[14:18].strip())
            atominfo['resnames'].append(atom_data[19:22].strip())
            atominfo['resids'].append(atom_data[23:33].strip())
            atominfo['fragids'].append(atom_data[34:44].strip())
            atominfo['xcoords'].append(atom_data[45:65].strip())
            atominfo['ycoords'].append(atom_data[65:85].strip())
            atominfo['zcoords'].append(atom_data[85:105].strip())
            atominfo['chainids'].append(atom_data[106:109].strip())
            atominfo['optflags'].append(atom_data[110:111].strip())
            if cpfver == 23:
                for chg in charge_label:
                    lstart = charge_label.index(chg)*20 + 111
                    lend = lstart + 20
                    atominfo[chg].append(atom_data[lstart:lend].strip())

        # Read the fragment data
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
            フラグメントが結合している相手の原子の番号(i10), 相手のフラグメントと結合している原子の番号(i10)
            }
        '''

        fnatoms, fbaas, fatminfos = readfragcpf(file, nfrag)
        pass

        # Read the fragment distance data
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

        # read the dimer info
        diminfo = {
            'fragi': [],
            'fragj': [],
            'min-dist': []
            }

        for _ in range(nfrag*(nfrag-1)//2):
            dim_data = file.readline().strip().split()
            diminfo['fragi'].append(dim_data[0])
            diminfo['fragj'].append(dim_data[1])
            diminfo['min-dist'].append(dim_data[2])

        # read the dipole data
        '''
        0.785970969662152E+00    0.281296712643063E+01   -0.132604581524068E+01
        0.145048892092187E+01   -0.334817286212412E+01   -0.161464626780787E+01
       -0.524704704075969E+00   -0.357367905308988E+01    0.156523382974318E+01
        0.245666859488253E+01   -0.160468359833147E+01    0.280152739221378E+01
        0.368656733092432E+01   -0.478966423463728E+01   -0.248675894502097E+00

            モノマー数行ループ {
                !双極子モーメント(HF密度)のx成分,!y成分,!z成分,!双極子モーメント(MP2密度)のx成分,!y成分,!z成分
            }
        '''

        mominfo = {
            'fragi': []
            }

        if cpfver == 23:
            # Initialize the data structures
            for dpm in DPM_label:
                mominfo[dpm] = []

        for _ in range(nfrag):
            dipole_data = file.readline().strip().split()
            count = 0
            for dpm in DPM_label:
                mominfo[dpm].append(dipole_data[count])
                count += 1

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
        condition = {}
        condition['basis_set'] = file.readline().strip()
        condition['electronic_state'] = file.readline().strip()
        condition['calculation_method'] = file.readline().strip()

        # Read the parameters for the AOpopulation, ptc, and ldimer
        data_approxy = file.readline().strip().split()
        condition['aoc'] = data_approxy[0]
        condition['ptc'] = data_approxy[1]
        condition['ldimer'] = data_approxy[2]

        static_data = {}
        static_data['nuclear_repulsion_energy'] = file.readline().strip()
        static_data['total_electronic_energy'] = file.readline().strip()
        static_data['total_energy'] = file.readline().strip()

        # Read the nuclear repulsion energy,
        # total electronic energy, and total energy
        '''
         5
         1   0.321733218676314E+02  -0.111217900611777E+03
         2   0.105581623892680E+03  -0.294954213591969E+03
         3   0.105580760608816E+03  -0.294944052573276E+03
         4   0.105592376341217E+03  -0.294954884978768E+03
         5   0.302972088661468E+03  -0.692375396192077E+03
        '''

        if cpfver == 23:
            # Initialize the data structures
            for mom in monomer_label:
                mominfo[mom] = []

        nmonomer = int(file.readline().strip())
        for _ in range(nmonomer):
            count = 0
            monomer_data = file.readline().strip().split()
            if cpfver == 23:
                for mom in monomer_label:
                    mominfo['fragi'].append(monomer_data[0])
                    mominfo[mom].append(monomer_data[count])
                    count += 1

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
         4         1
         0.472340649687743E+02  -0.472327218765366E+02   0.135816158277180E-02
         0.380952857454986E-04  -0.531646308132849E-04   0.292747285648431E-04
         4         2
         0.103751129938383E+03  -0.103749009778853E+03   0.201925943345316E-02
         0.734678685546442E-03  -0.633778588650102E-03   0.186391389006602E-02
         4         3
         0.151609836272825E+03  -0.166363652509399E+03  -0.144566293381012E+02
         -0.170776291331123E+00  -0.126410607141253E+00   0.100185717607355E+02
         5         1
         0.693589957512211E+02  -0.693562420378788E+02   0.278835582868453E-02
         0.343382369294432E-05  -0.380763100480408E-04  -0.695127646110905E-05
         5         2
         0.159493576274165E+03  -0.159501166271000E+03  -0.716229202939189E-02
         0.448447060864510E-02  -0.491217541426181E-02   0.144539031120985E-01
         5         3
         0.181521409736929E+03  -0.181515319038899E+03   0.598906538814958E-02
         0.747041979622054E-03  -0.645409338090985E-03   0.184540822159818E-02
         5         4
         0.235076980380320E+03  -0.249830849163806E+03  -0.144568950278076E+02
         -0.171820856781494E+00  -0.125152898897113E+00   0.100103224179901E+02
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

        if cpfver == 23:
            # Initialize the data structures
            for dim in dimer_label:
                diminfo[dim] = []

        ndimer = int(file.readline().strip())
        for _ in range(ndimer):
            count = 0
            dimer_data = file.readline().strip().split()
            if cpfver == 23:
                for dim in dimer_label:
                    diminfo[dim].append(dimer_data[count])
                    count += 1

        print(cpfver)
        print(atominfo)
        print(fnatoms, fbaas, fatminfos)
        print(condition)
        print(static_data)
        print(mominfo)
        print(diminfo)

    return True


if __name__ == '__main__':
    read_data('gly5.cpf')
