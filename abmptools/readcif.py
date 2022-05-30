import sys
import os
import numpy as np
import math
import statistics
import copy
import argparse
import glob
import shutil

def getcartesiancellvec(angle, length):

    a_vec = [0] * 3
    b_vec = [0] * 3
    c_vec = [0] * 3

    a = length[0]
    b = length[1]
    c = length[2]

    alpha = angle[0]/180.0*math.pi
    beta =  angle[1]/180.0*math.pi
    gamma = angle[2]/180.0*math.pi

#     print('atom', atoms)
#     print('coords', coords)
#     print('length', length)
    print('a, b, c', a, b, c)
    print('angle', angle)
    # print('alpha, beta, gamma', alpha, beta, gamma)
    # print(math.sin(alpha))

    a_vec[0] = a
    a_vec[1] = 0
    a_vec[2] = 0
    # print('a_vec0', a_vec[0])

    b_vec[0] = b*math.cos(gamma)
    b_vec[1] = b*math.sin(gamma)
    b_vec[2] = 0

    tan_theta = (c*math.cos(alpha)-c*math.cos(gamma)*math.cos(beta))/(c*math.cos(beta)*math.sin(gamma))
    theta = math.atan(tan_theta)
    t = c*math.cos(beta) / math.cos(theta)

    print('theta', theta)
    print('theta(degree)', math.degrees(theta))
    print('t', t)
    c_vec[0] = c*math.cos(beta)
    print('c_x', c_vec[0])
    c_vec[1] = t*math.sin(theta)
    print('c_y', c_vec[1])
    c_vec[2] = math.sqrt(c**2 - (c_vec[0]**2 + c_vec[1]**2))

    # c_vec[0] = c*math.cos(beta) + (c*math.cos(alpha)*math.cos(gamma))
    # c_vec[1] = c*math.cos(alpha)*math.sin(gamma)

    # print(c_vec[0], c_vec[1], c_vec[2], c)
    # print(c**2 - (c_vec[0]**2 + c_vec[1]**2))

    print('a_vec', a_vec)
    print('b_vec', b_vec)
    print('c_vec', c_vec)

    return a_vec, b_vec, c_vec

def getcartesianmol(a, b, c, acell, bcell, ccell):

    a = list(np.array(acell) * a)
    b = list(np.array(bcell) * b)
    c = list(np.array(ccell) * c)

    return a, b, c

# def intocell_xyz(incoords, anum_mol):
#     coordss = []
#     a1 = anum_mol[0]
#     coordss.append(incoords[:a1])
#     coordss.append(incoords[a1:])
#     # print('coordss', coordss)
#
#     shiftnums = []
#
#     for coords in coordss:
#         print('mean',statistics.mean(coords))
#         if statistics.mean(coords) < 0.0:
#             shiftnum = -(math.floor(statistics.mean(coords)))
#             # print('mean', statistics.mean(coords))
#         elif statistics.mean(coords) > 1.0:
#             shiftnum = -(math.floor(statistics.mean(coords)))
#             # print('mean', statistics.mean(coords))
#         else:
#             shiftnum = 0
#         # print('shiftnum', shiftnum)
#         shiftnums.append(shiftnum)
#     # print('newcoordss', newcoordss)
#     print('shiftnums', shiftnums)
#     return shiftnums

def intocell(incoords, anum_mol):
    if len(anum_mol) == 2:
        coordss = []
        a1 = anum_mol[0]
        coordss.append(incoords[:a1])
        coordss.append(incoords[a1:])
        # print('coordss', coordss)

        newcoordss = []
        for coords in coordss:
            # print('mean',statistics.mean(coords))
            while statistics.mean(coords) < 0.0:
                coords =  list(np.array(coords) + 1.0)
                # print('mean', statistics.mean(coords))
            while statistics.mean(coords) > 1.0:
                coords =  list(np.array(coords) - 1.0)
                # print('mean', statistics.mean(coords))
            newcoordss.extend(coords)
        # print('newcoordss', newcoordss)

    elif len(anum_mol) == 1:
        coordss = []
        a1 = anum_mol[0]
        coordss.append(incoords[:a1])
        # print('coordss', coordss)

        newcoordss = []
        for coords in coordss:
            # print('mean',statistics.mean(coords))
            while statistics.mean(coords) < 0.0:
                coords =  list(np.array(coords) + 1.0)
                # print('mean', statistics.mean(coords))
            while statistics.mean(coords) > 1.0:
                coords =  list(np.array(coords) - 1.0)
                # print('mean', statistics.mean(coords))
            newcoordss.extend(coords)
        # print('newcoordss', newcoordss)

    else:
        print('Not supported yet.')
        sys.exit()

    return newcoordss

def getsymcoord(sympos, incoord, label):

    if sympos in ['x']:
        coord = incoord[0]
    if sympos in ['y']:
        coord = incoord[1]
    if sympos in ['z']:
        coord = incoord[2]
    if sympos in ['-x']:
        coord = -incoord[0]
    if sympos in ['-y']:
        coord = -incoord[1]
    if sympos in ['-z']:
        coord = -incoord[2]

    if sympos in ['1/2+z', 'z+1/2']:
        coord = 1/2 + incoord[2]
    if sympos in ['1/2+y', 'y+1/2']:
        coord = 1/2 + incoord[1]
    if sympos in ['1/2+x', 'x+1/2']:
        coord = 1/2 + incoord[0]
    if sympos in ['1/2-z', '-z+1/2']:
        coord = 1/2 - incoord[2]
    if sympos in ['1/2-y', '-y+1/2']:
        coord = 1/2 - incoord[1]
    if sympos in ['1/2-x', '-x+1/2']:
        coord = 1/2 - incoord[0]

    if sympos in ['-1/2+z', 'z-1/2']:
        coord = -1/2 + incoord[2]
    if sympos in ['-1/2+y', 'y-1/2']:
        coord = -1/2 + incoord[1]
    if sympos in ['-1/2+x', 'x-1/2']:
        coord = -1/2 + incoord[0]
    if sympos in ['-1/2-z', '-z-1/2']:
        coord = -1/2 - incoord[2]
    if sympos in ['-1/2-y', '-y-1/2']:
        coord = -1/2 - incoord[1]
    if sympos in ['-1/2-x', '-x-1/2']:
        coord = -1/2 - incoord[0]



    if sympos in ['1/3+z']:
        coord = 1/3 + incoord[2]
    if sympos in ['1/3+y']:
        coord = 1/3 + incoord[1]
    if sympos in ['1/3+x']:
        coord = 1/3 + incoord[0]
    if sympos in ['2/3+z']:
        coord = 2/3 + incoord[2]
    if sympos in ['2/3+y']:
        coord = 2/3 + incoord[1]
    if sympos in ['2/3+x']:
        coord = 2/3 + incoord[0]

    if sympos in ['1/3-z']:
        coord = 1/3 - incoord[2]
    if sympos in ['1/3-y']:
        coord = 1/3 - incoord[1]
    if sympos in ['1/3-x']:
        coord = 1/3 - incoord[0]
    if sympos in ['2/3-z']:
        coord = 2/3 - incoord[2]
    if sympos in ['2/3-y']:
        coord = 2/3 - incoord[1]
    if sympos in ['2/3-x']:
        coord = 2/3 - incoord[0]

    if sympos in ['1/6+z']:
        coord = 1/6 + incoord[2]
    if sympos in ['1/6+y']:
        coord = 1/6 + incoord[1]
    if sympos in ['1/6+x']:
        coord = 1/6 + incoord[0]
    if sympos in ['1/6-z']:
        coord = 1/6 - incoord[2]
    if sympos in ['1/6-y']:
        coord = 1/6 - incoord[1]
    if sympos in ['1/6-x']:
        coord = 1/6 - incoord[0]

    if sympos in ['1/4+z']:
        coord = 1/4 + incoord[2]
    if sympos in ['1/4+y']:
        coord = 1/4 + incoord[1]
    if sympos in ['1/4+x']:
        coord = 1/4 + incoord[0]
    if sympos in ['1/4-z']:
        coord = 1/4 - incoord[2]
    if sympos in ['1/4-y']:
        coord = 1/4 - incoord[1]
    if sympos in ['1/4-x']:
        coord = 1/4 - incoord[0]

    if sympos in ['3/4+z']:
        coord = 3/4 + incoord[2]
    if sympos in ['3/4+y']:
        coord = 3/4 + incoord[1]
    if sympos in ['3/4+x']:
        coord = 3/4 + incoord[0]
    if sympos in ['3/4-z']:
        coord = 3/4 - incoord[2]
    if sympos in ['3/4-y']:
        coord = 3/4 - incoord[1]
    if sympos in ['3/4-x']:
        coord = 3/4 - incoord[0]


    if sympos in ['5/6+z']:
        coord = 5/6 + incoord[2]
    if sympos in ['5/6+y']:
        coord = 5/6 + incoord[1]
    if sympos in ['5/6+x']:
        coord = 5/6 + incoord[0]
    if sympos in ['5/6-z']:
        coord = 5/6 - incoord[2]
    if sympos in ['5/6-y']:
        coord = 5/6 - incoord[1]
    if sympos in ['5/6-x']:
        coord = 5/6 - incoord[0]

    if sympos in ['x-y']:
        coord = incoord[0] - incoord[1]

    if sympos in ['-x+y']:
        coord = - incoord[0] + incoord[1]

    if sympos in ['1/3+x-y']:
        coord = 1/3 + incoord[0] - incoord[1]

    if sympos in ['1/3-x+y']:
        coord = 1/3 - incoord[0] + incoord[1]

    if sympos in ['2/3-x+y']:
        coord = 2/3 - incoord[0] + incoord[1]

    if sympos in ['2/3+x-y']:
        coord = 2/3 + incoord[0] - incoord[1]

    return coord


# def shiftatom_xyz(a_xyz, b_xyz, c_xyz, a_num, b_num, c_num, cella, cellb, cellc, atomid):
#     print('abc_shiftnum', a_num, b_num, c_num)
#     if atomid < anum_inmol[0]:
#         # print('beforea', a_xyz)
#         a_xyz = np.array(a_xyz) + np.array(cella) * a_num[0]
#         b_xyz = np.array(b_xyz) + np.array(cellb) * b_num[0]
#         c_xyz = np.array(c_xyz) + np.array(cellc) * c_num[0]
#         # print('aftera', a_xyz)
#     if atomid >= anum_inmol[0]:
#         a_xyz = np.array(a_xyz) + np.array(cella) * a_num[1]
#         b_xyz = np.array(b_xyz) + np.array(cellb) * b_num[1]
#         c_xyz = np.array(c_xyz) + np.array(cellc) * c_num[1]
#     return list(a_xyz), list(b_xyz), list(c_xyz)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                prog='readcif.py', # program name
                usage='python readcif.py -i xxx.cif --atomnum xx -l x', # program usage
                description='readcif script',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--input',
                        help='input coordinate file (pdb)',
                        nargs='*',
                        action='append',
                        required=True)

    parser.add_argument('-calcd', '--calcdist',
                        help='calcflag',
                        action='store_true',
                        )

    parser.add_argument('-intra', '--intra',
                        help='nointra',
                        action='store_true',
                        )

    parser.add_argument('-d', '--dist',
                        help='dist',
                        type=float,
                        default=2.0,
                        )

    parser.add_argument('-a', '--tgtatom',
                        help='tgtatm',
                        nargs='*',
                        action='append',
                        default=[])

    parser.add_argument('-od', '--odir',
                        help='outdir',
                        default='cifout'
                        )

    parser.add_argument('-an', '--atomnum',
                        help='atom num',
                        nargs='*',
                        type=int,
                        action='append',
                        required=True)


    parser.add_argument('-l', '--layer',
                        help='layer',
                        type=int,
                        default=1
                        )

    parser.add_argument('-nopdb', '--nopdb',
                        help='flag nopdb',
                        action='store_false')

    parser.add_argument('-noout', '--noout',
                        help='flag nopdb',
                        action='store_true')

    parser.add_argument('-min', '--min',
                        help='min',
                        action='store_true')

    parser.add_argument('-au', '--asymonly',
                        help='assymmetric unit only',
                        action='store_true')


    # get args
    args = parser.parse_args()

    print('coord(cif) =', args.input)
    print('odir = ', args.odir)
    print('atomnum = ', args.atomnum)
    print('layer =', args.layer)
    print('pdbflag =', args.nopdb)
    print('calcdist', args.calcdist)
    print('intra', args.intra)
    print('dist', args.dist)
    print('tgtatom', args.tgtatom)
    print('out', args.noout)
    print('min', args.min)
    print('asymmetric only', args.asymonly)

    calc_dist = args.calcdist
    tgtdist = args.dist

    if args.tgtatom:
        tgtatoms = args.tgtatom[0]

    if args.intra:
        nointra = False
    else:
        nointra = True
    maxnum = 10

    ## -- user setting end
    # argvs = sys.argv
    # print(argvs)

    infiles = args.input
    odir = args.odir
    xyzdir = odir + '/layer' + str(args.layer) + '/xyz'
    pdbdir = odir + '/layer' + str(args.layer) + '/pdb'

    anum_inmol = args.atomnum[0]
    image = args.layer
    pdbflag = args.nopdb

    # symlist = ['P21/N', 'P21/M', 'P21', 'PNA21', 'P212121', 'P21212', 'PCA21', 'P21/C', 'C2/C', 'C2', 'CC', 'P-1', 'P1', 'PC', 'P2/N', 'IC', 'I2', 'IA',  'I2/A', 'I2/C', 'P2/C',  'PBCN', 'PCCN', 'PBCA', 'PNMA', 'PN']

    if pdbflag:
        import abmptools as ampt
        aobj = ampt.mol_io()

    print('atom num mol:', anum_inmol)
    # print('tgtatom:', tgtatoms)
    # print('tgtdist:', tgtdist)
    print('########## Read Start #########')

    if os.path.exists(odir) is False:
        os.makedirs(odir)
    if os.path.exists(xyzdir) is False:
        os.makedirs(xyzdir)
    if not os.path.exists(pdbdir) and pdbflag:
        os.makedirs(pdbdir)

    minvals = []
    mindatas = []
    for infile in infiles[0]:
        out, ext = os.path.splitext(infile)
        out = out.split('/')[-1]

        print('\n##  Start Read', infile)
        lines = open(infile, 'r').readlines()
        ## initialize
        num = 1
        coordflag = False
        coordend = False

        cellend = False
        atoms = []
        coords = []
        atomsmol = []
        coordsmol = []

        symflag = False
        symxyzs = []
        symxyz = []

        angle = []
        length = []
        lengthmol = []
        anglemol = []
        paxnums = []
        znum = []
        acolumns = []
        acolumnsflag = False
        scolumns = []
        scolumnsflag = False
        is_eq_pax = False
        paxnum = 0

        ## get line
        for i in range(len(lines)):
            line = lines[i].split()
            if len(line) == 0:
                continue
            # CONFLEX case
            if line[0:2] == ['#', 'CONFLEX8']:
                if num > maxnum:
                    break
                print('start data', num)
                num += 1
            if '_space_group_name_H-M' in line[0]:
                if len(line) == 3:
                    line[1] = line[1] + line[2]
                if len(line) == 4:
                    line[1] = line[1] + line[2] + line[3]
                if len(line) == 5:
                    line[1] = line[1] + line[2] + line[3] + line[4]
                sym = line[1].replace("'", '').upper()
                print('space group:', sym)

            ## get _symmetry
            if '_symmetry_equiv_pos' in line[0]:
                print(line[0])
                scolumns.append(line[0])
                scolumnsflag = True
                if line[0] == '_symmetry_equiv_pos_as_xyz':
                    print('start!!!!')
                    is_eq_pax = True
                    paxnum = 0
                continue
            if scolumnsflag:
                if '_symmetry_equiv_pos' not in line[0]:
                    scolumnsflag = False
                    symflag = True
                    xyzsym_idx = scolumns.index('_symmetry_equiv_pos_as_xyz')
                    print('xyzsym_idx', xyzsym_idx)

            if symflag == True:
                print(line[0:1])
                if 'loop_' in line[0] or '_cell' in line[0]:
                    symflag = False
                    symxyzs.append(symxyz)
                    symxyz = []
                else:
                    symxyz.append(line[xyzsym_idx].split(','))
                    print(symxyz)

            # stop count equiv pos as xyz
            if is_eq_pax:
                if line[0:1] == ['_cell_length_a']:
                    print('stop!')
                    is_eq_pax = False
                    paxnums.append(paxnum)
                else:
                    paxnum += 1

            ## get cell length and line
            if line[0:1] == ['_cell_length_a']:
                length.append(float(line[1]))
                continue
            if line[0:1] == ['_cell_length_b']:
                length.append(float(line[1]))
                continue
            if line[0:1] == ['_cell_length_c']:
                length.append(float(line[1]))
                continue
            if line[0:1] == ['_cell_angle_alpha']:
                angle.append(float(line[1]))
                continue
            if line[0:1] == ['_cell_angle_beta']:
                angle.append(float(line[1]))
                continue
            if line[0:1] == ['_cell_angle_gamma']:
                angle.append(float(line[1]))
                continue
            if line[0:1] == ['_cell_formula_units_Z']:
                znum.append(int(line[1]))
                cellend = True
                continue

            ## save length and angle for mol list
            if cellend == True:
                cellend = False
                # coords = np.array(coords)
                lengthmol.append(length)
                anglemol.append(angle)
                # print('aaaaa', lengthmol)
                length = []
                angle = []

            ## get coord
            if '_atom_site' in line[0] and '_geom' not in line[0]:
                print(line[0])
                acolumns.append(line[0])
                acolumnsflag = True
                continue
            if acolumnsflag == True:
                if '_atom_site' not in line[0]:
                    acolumnsflag = False
                    coordflag = True
                    symbol_idx = acolumns.index('_atom_site_type_symbol')
                    x_idx = acolumns.index('_atom_site_fract_x')
                    y_idx = acolumns.index('_atom_site_fract_y')
                    z_idx = acolumns.index('_atom_site_fract_z')
                    print(symbol_idx, x_idx, y_idx, z_idx)
            if coordflag == True:
                if line[0:1] == ['loop_']:
                    coordflag = False
                    coordend = True
                    continue

                atoms.append(line[symbol_idx])
                coords.append([float(line[x_idx]), float(line[y_idx]), float(line[z_idx])])

            # end read coord xyz and save for list
            if coordend == True:
                coordend = False
                atomsmol.append(atoms)
                coordsmol.append(coords)
                atoms = []
                coords = []
    #     print(atomsmol)
    #     print(coordsmol)
    #     print(len(atomsmol[0]))
    #     print(len(coordsmol[0]))
    #     print(lengthmol)
    #     print(anglemol)
    #     print(len(lengthmol))
    #     print(len(anglemol))

        print('znum', znum)
        print('paxnums', paxnums)
        # initialize
        a_xyz = [0] * 3
        b_xyz = [0] * 3
        c_xyz = [0] * 3
        xyz = [0] * 3
        a_xyzs = []
        b_xyzs = []
        c_xyzs = []
        xyzs = []
        xyzs_sym = []
        xyzsmol = []
        xyzs27mols = []
        truemolid = []
        zprime = []
        molnum = len(atomsmol)
        atomnum = len(atomsmol[0])
        print('molnum', molnum)
        print('atomnum', atomnum)

        ## get cartesian  coordinate
        # i: datanum in each cif file
        for i in range(molnum):
            print('\nrun mol', i)
            tgtflag = False
            length = lengthmol[i]
            angle = anglemol[i]
            # print(length)
            # print(angle)
            cella, cellb, cellc = copy.deepcopy(getcartesiancellvec(angle, length))

            as_list = []
            bs_list = []
            cs_list = []

            # compare znum and num(equiv_pos_as_xyz)
            # zprime = znum / eq_pax
            zprime.append(int(znum[i] / paxnums[i]))

            # f = open('zprime.log', 'a')
            # print("Z'=", zprime[i], file=f)

            if args.asymonly:
                paxnums[i] = 1

            for z_id in range(paxnums[i]):
                a_s = []
                b_s = []
                c_s = []
                # j: atomnum
                for j in range(atomnum):
                    atoms = atomsmol[i][j]
                    coords = coordsmol[i][j]

                    # get symmetry_equiv_pos_as_xyz from cif file
                    acoord = getsymcoord(symxyzs[i][z_id][0], coords, 0)
                    bcoord = getsymcoord(symxyzs[i][z_id][1], coords, 1)
                    ccoord = getsymcoord(symxyzs[i][z_id][2], coords, 2)

                    a_s.append(copy.deepcopy(acoord))
                    b_s.append(copy.deepcopy(bcoord))
                    c_s.append(copy.deepcopy(ccoord))

                if len(anum_inmol) == 1:
                    anum_inmol *= zprime[i]
                elif len(anum_inmol) == 2 and zprime[i] == 1:
                    print("Z' != 1 and multi-type mol is not supported.")

                # print(len(a_s))
                a_s = copy.deepcopy(intocell(a_s, anum_inmol))
                b_s = copy.deepcopy(intocell(b_s, anum_inmol))
                c_s = copy.deepcopy(intocell(c_s, anum_inmol))
                # print(len(a_s))

                # print('check shift (z_id:', z_id, ')')
                # a_num = copy.deepcopy(intocell_xyz(a_s, anum_inmol))
                # b_num= copy.deepcopy(intocell_xyz(b_s, anum_inmol))
                # c_num = copy.deepcopy(intocell_xyz(c_s, anum_inmol))

                # print('shiftnum', a_num, b_num, c_num)
                for k in range(len(a_s)):
                    a_xyz, b_xyz, c_xyz = copy.deepcopy(getcartesianmol(a_s[k], b_s[k], c_s[k], cella, cellb, cellc))

                    # a_xyz, b_xyz, c_xyz = copy.deepcopy(shiftatom(a_xyz, b_xyz, c_xyz, a_num, b_num, c_num, cella, cellb, cellc, k))
                    a_xyzs.append(copy.deepcopy(a_xyz))
                    b_xyzs.append(copy.deepcopy(b_xyz))
                    c_xyzs.append(copy.deepcopy(c_xyz))

    #             print('cartesian section')
    #             print('a_xyzs', a_xyzs)
    #             print('b_xyzs', b_xyzs)
    #             print('c_xyzs', c_xyzs)
    #             print('lenaxyzs', len(a_xyzs))

                for k in range(len(a_xyzs)):
                    xyz[0] = a_xyzs[k][0] + b_xyzs[k][0] + c_xyzs[k][0]
                    xyz[1] = a_xyzs[k][1] + b_xyzs[k][1] + c_xyzs[k][1]
                    xyz[2] = a_xyzs[k][2] + b_xyzs[k][2] + c_xyzs[k][2]

                    xyzs.append(copy.deepcopy(xyz))
    #             print('xyzs', xyzs)
    #             print(len(xyzs))
    #             print(atomsmol[i])

                xyzs_sym.append(copy.deepcopy(xyzs))
                xyzs = []
                a_xyzs = []
                b_xyzs = []
                c_xyzs = []

                as_list = as_list + (copy.deepcopy(a_s))
                bs_list = bs_list + (copy.deepcopy(b_s))
                cs_list = cs_list + (copy.deepcopy(c_s))

            xyzsmol.append(xyzs_sym)
            xyzs_sym = []

            # print('as_list', as_list)

            def getpermol(coords, anum_inmol):
                permols = []
                permol = []
                count = 0
                i = 0
                for coord in coords:
                    permol.append(coord)
                    count += 1
                    if count == anum_inmol[i % len(anum_inmol)]:
                        permols.append(permol)
                        permol = []
                        count = 0
                        i += 1
                return permols

            # calc dist
            as_27 = []
            bs_27 = []
            cs_27 = []
            if image == 5:
                for xshift in [0, 1, -1, 2, -2]:
                    for yshift in [0, 1, -1, 2, -2]:
                        for zshift in [0, 1, -1, 2, -2]:
                            as_27 = as_27 + list(np.array(as_list) + xshift)
                            bs_27 = bs_27 + list(np.array(bs_list) + yshift)
                            cs_27 = cs_27 + list(np.array(cs_list) + zshift)

            if image == 4:
                for xshift in [0, 1, -1, 2]:
                    for yshift in [0, 1, -1, 2]:
                        for zshift in [0, 1, -1, 2]:
                            as_27 = as_27 + list(np.array(as_list) + xshift)
                            bs_27 = bs_27 + list(np.array(bs_list) + yshift)
                            cs_27 = cs_27 + list(np.array(cs_list) + zshift)

            if image == 3:
                for xshift in [0, 1, -1]:
                    for yshift in [0, 1, -1]:
                        for zshift in [0, 1, -1]:
                            as_27 = as_27 + list(np.array(as_list) + xshift)
                            bs_27 = bs_27 + list(np.array(bs_list) + yshift)
                            cs_27 = cs_27 + list(np.array(cs_list) + zshift)

            if image == 2:
                for xshift in [0, 1]:
                    for yshift in [0, 1]:
                        for zshift in [0, 1]:
                            as_27 = as_27 + list(np.array(as_list) + xshift)
                            bs_27 = bs_27 + list(np.array(bs_list) + yshift)
                            cs_27 = cs_27 + list(np.array(cs_list) + zshift)

            if image == 1:
                for xshift in [0]:
                    for yshift in [0]:
                        for zshift in [0]:
                            as_27 = as_27 + list(np.array(as_list) + xshift)
                            bs_27 = bs_27 + list(np.array(bs_list) + yshift)
                            cs_27 = cs_27 + list(np.array(cs_list) + zshift)

            if image >= 6:
                print('Error!! layer over 5 is not supported yet.')
                sys.exit()
            # print(as_27, bs_27, cs_27)
            # sys.exit()

            # --- calc_dist case ---
            # print(a_s)
            distdatas = []
            distvals = []
            if calc_dist == True:
                # print(as_27)
                a_permol = getpermol(as_27, anum_inmol)
                b_permol = getpermol(bs_27, anum_inmol)
                c_permol = getpermol(cs_27, anum_inmol)
                atoms_permol = getpermol(atomsmol[0], anum_inmol)
                # print('a_permol', len(a_permol))

                for moli in range(len(a_permol)):
                    for molj in range(len(a_permol)):
                        if nointra == True:
                            if moli == molj:
                                continue
                        for atomi in range(len(a_permol[moli])):
                            for atomj in range(len(a_permol[molj])):
                                atmname_i = atoms_permol[moli % len(anum_inmol)][atomi]
                                atmname_j = atoms_permol[molj % len(anum_inmol)][atomj]
                                if (atmname_i in tgtatoms) and (atmname_j in tgtatoms):

                                    adist = abs(a_permol[moli][atomi] - a_permol[molj][atomj])
#                                     if adist > 0.5:
#                                         adist = 1.0 - adist
                                    bdist = abs(b_permol[moli][atomi] - b_permol[molj][atomj])
#                                     if bdist > 0.5:
#                                         bdist = 1.0 - bdist
                                    cdist = abs(c_permol[moli][atomi] - c_permol[molj][atomj])
#                                     if cdist > 0.5:
#                                         cdist = 1.0 - cdist
                                    adistxyz, bdistxyz, cdistxyz = copy.deepcopy(getcartesianmol(adist, bdist, cdist, cella, cellb, cellc))

                                    xvec = adistxyz[0] + bdistxyz[0] + cdistxyz[0]
                                    yvec = adistxyz[1] + bdistxyz[1] + cdistxyz[1]
                                    zvec = adistxyz[2] + bdistxyz[2] + cdistxyz[2]
                                    dist = math.sqrt(xvec**2 + yvec**2 + zvec**2)

                                    if dist < tgtdist and dist > 0.1:
                                        print('mol', moli+1, 'atom', atomi+1, atmname_i, '- mol', molj+1, 'atom', atomj+1, atmname_j, "{:6.3f}".format(dist))
                                        tgtflag = True
                                        distdatas.append([moli+1, atomi+1, atmname_i, molj+1, atomj+1, atmname_j])
                                        distvals.append(dist)
                if tgtflag == True:
                    truemolid.append(i)

            if args.min:
                minidx = distvals.index(min(distvals))
                minval = min(distvals)
                mindata = distdatas[minidx]
                print('minval', minval)
                print('distdata', mindata)
                minvals.append(minval)
                mindatas.append(['file', infile, 'molecules', i, mindata])
            # print('as_list', as_list)
            # print('len as_list', len(as_list))

            a27_xyzs = []
            b27_xyzs = []
            c27_xyzs = []

            for k in range(len(as_27)):
                a27_xyz, b27_xyz, c27_xyz = copy.deepcopy(getcartesianmol(as_27[k], bs_27[k], cs_27[k], cella, cellb, cellc))
                a27_xyzs.append(a27_xyz)
                b27_xyzs.append(b27_xyz)
                c27_xyzs.append(c27_xyz)
                # print(a27_xyz, b27_xyz, c27_xyz)

            m27_xyz = [0, 0, 0]
            # print(len(as_27))
            # print('a27_xyzs[0]', a27_xyzs[0])
            # print('a27_xyzs[1]', a27_xyzs[1])
            # print('a27_xyzs[2]', a27_xyzs[2])

            xyzs27mol = []
            for k in range(len(as_27)):
                m27_xyz[0] = a27_xyzs[k][0] + b27_xyzs[k][0] + c27_xyzs[k][0]
                m27_xyz[1] = a27_xyzs[k][1] + b27_xyzs[k][1] + c27_xyzs[k][1]
                m27_xyz[2] = a27_xyzs[k][2] + b27_xyzs[k][2] + c27_xyzs[k][2]
                xyzs27mol.append(copy.deepcopy(m27_xyz))

            print(xyzs27mol[0], xyzs27mol[1])
            xyzs27mols.append(copy.deepcopy(xyzs27mol))
            print('len xyzmol', len(xyzs27mols[0]))

        ## write for file
        if calc_dist:
            print('satisfy criteria mol:', truemolid)

        if args.noout:
            continue

        print('## Write ', args.layer, ' layer sell ##')
        for i in range(molnum):
            print('mol -- ', i, ' --')
            if molnum ==1:
                number = ''
            else:
                number = '%05d' % i
            if args.asymonly:
                austr = 'asymonly'
            else:
                austr = ''
            oname = xyzdir + '/' + out + number + 'layer' + str(args.layer) + 'Zp' + str(zprime[i]) + austr + '.xyz'
            fo = open(oname, 'w')
            print(len(xyzs27mols[i]), file=fo)
            print('', file=fo)
            for j in range(len(xyzs27mols[i])):
                # for k in range(len(xyzs27mols[i][j])):
                    # print(i)
                print((atomsmol[0][j % len(atomsmol[0])]), xyzs27mols[i][j][0], xyzs27mols[i][j][1], xyzs27mols[i][j][2], file=fo)
            fo.close()

            if pdbflag:
                aobj.convert_xyz_pdb(oname)

    if calc_dist and args.min:
        fomin = open('mindist.log', 'w')
        for i in range(len(minvals)):
            print(mindatas[i], minvals[i], file=fomin)
        print('generated mindist.log')

    if pdbflag and not args.noout:
        pdbs = glob.glob(xyzdir + '/*pdb')
        for pdb in pdbs:
            shutil.move(pdb, pdbdir + '/' + pdb.split('/')[-1])
