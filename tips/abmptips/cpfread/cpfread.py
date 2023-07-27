def read_data(filepath):
    # Initialize the data structures
    cpfver = ''
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
        cpfver = file.readline().strip()
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
        NR, HF, ES, MP2, PR-MP2, SCS-MP2(Grimme), MP3, SCS-MP3(MP2.5), HF-BSSE, MP2-BSSE, SCS-MP2-BSSE, MP3-BSSE, SCS-MP3-BSSE, SOLV-ES, SOLV-NP, EX, CT, DQ
        '''

        charge_label = file.readline().strip().split()
        DPM_label = file.readline().strip().split()
        monomer_label = file.readline().strip().split()
        dimer_label = file.readline().strip().split()


        # Read the atom data
        ''' \
         1 N   N   GLY          1          1         0.1620000000       -0.2020000000        0.0000000000     1       -0.3972241421
        '''
        '''L7 ("原子数"分ループ) {原子の番号(i10), x, 元素記号(a2), x, 原子タイプ(a4), x, 残基名(a3), x, 残基番号(i10), x,フラグメント番号(i10), x, x座標(f20.10), y座標(f20.10) z座標(f20.10), x, Chain ID(a3), x, 構造最適化オプション(i1), !Mulliken原子電荷(HF密度)(f20.10), !Mulliken原子電荷(MP2密度)(f20.10), !NPA原子電荷(HF密度)(f20.10), !NPA原子電荷(MP2密度)(f20.10), !ESP原子電荷(HF密度)(f20.10), !ESP原子電荷(MP2密度)(f20.10)}
        '''

        for chg in charge_label:
            atominfo[chg] = []

        for _ in range(natom):
            atom_data = file.readline().rstrip()
            atominfo['alabels'].append(atom_data[0:10])
            atominfo['elems'].append(atom_data[11:13])
            atominfo['elemtype'].append(atom_data[14:18])
            atominfo['resnames'].append(atom_data[19:22])
            atominfo['resids'].append(atom_data[23:33])
            atominfo['fragids'].append(atom_data[34:44])
            atominfo['xcoords'].append(atom_data[45:65])
            atominfo['ycoords'].append(atom_data[65:85])
            atominfo['zcoords'].append(atom_data[85:105])
            atominfo['chainids'].append(atom_data[106:109])
            atominfo['optflags'].append(atom_data[110:111])
            for chg in charge_label:
                lstart = charge_label.index(chg)*20 + 111 
                lend = lstart + 20
                atominfo[chg].append(atom_data[lstart:lend])


          # Read the fragment data
#         for _ in range(fragment_count):
#             fragment_data = list(map(int, file.readline().split()))
#             data['fragments'].append(fragment_data)
# 
#         # Add similar loops for monomers, dimers, trimers, tetramers data reading as per your data structure
# 
#         # Read the basis function name, electronic state, and calculation method
#         data['basis_function_name'] = file.readline().strip()
#         data['electronic_state'] = file.readline().strip()
#         data['calculation_method'] = file.readline().strip()

    return data

read_data('gly5.cpf')


