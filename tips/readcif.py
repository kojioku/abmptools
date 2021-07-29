import sys
import os
import numpy as np
import math
import statistics
import copy

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

    print(c_vec)


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
    ## -- user setiing
    # anum_inmol = [24]
    anum_inmol =[7, 55]
    # anum_inmol =[47, 14]
    # anum_inmol =[24, 1] clonidine_hydrochloride2-PCS.cif
    # anum_inmol =[21, 1] hydralazine_hydrochloride200-PCS.cif
    # anum_inmol =[24, 17] CA_GA_2csp-PCS.cif
    # anum_inmol =[24, 0]

    calc_dist = True
    tgtdist = 1.5
    tgtatoms = ['O', 'H']
    image = 2
    nointra = True

    maxnum = 10
    ## -- user setting end

    argvs = sys.argv
    print(argvs)
    print('atom num mol:', anum_inmol)
    print('tgtatom:', tgtatoms)
    print('tgtdist:', tgtdist)

    infile = argvs[1]
    odir = argvs[2]
    if os.path.exists(odir) is False:
        os.makedirs(odir)
    out, ext = os.path.splitext(infile)
    out = out.split('/')[-1]

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
    angle = []
    length = []
    lengthmol = []
    anglemol = []
    znum = []


    ## get line
    for i in range(len(lines)):
        line = lines[i].split()
        # print(line)
        if line[0:2] == ['#', 'CONFLEX8']:
            if num > maxnum:
                break
            print('start data', num)
            num += 1
        if line[0:1] == ['_symmetry_space_group_name_H-M']:
            sym = line[1].replace("'", '')
            print('space group:', sym)
            symlist = ['P21/N', 'PNA21', 'P212121', 'P21/C', 'C2/C', 'P-1']
            if sym in symlist:
                pass
            else:
                print(sym, 'Not supported yet.')

        ## get length and line
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
        if line[0:1] == ['_atom_site_occupancy']:
            coordflag = True
            continue
        if coordflag == True:
            if line[0:1] == ['loop_']:
                coordflag = False
                coordend = True

                continue
            atoms.append(line[1])
            coords.append([float(line[2]), float(line[3]), float(line[4])])
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
    molnum = len(atomsmol)
    atomnum = len(atomsmol[0])
    print('molnum', molnum)
    print('atomnum', atomnum)

    ## get cartesian  coordinate
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
        for z_id in range(znum[i]):
            a_s = []
            b_s = []
            c_s = []
            for j in range(atomnum):
                # print('atom', j)
                atoms = atomsmol[i][j]
                coords = coordsmol[i][j]
                if sym == 'P21/N':
                    if z_id == 0:
                        acoord = coords[0]
                        bcoord = coords[1]
                        ccoord = coords[2]

                    if z_id == 1:
                        acoord = ((-coords[0]) + 0.5)
                        bcoord = (coords[1] + 0.5)
                        ccoord = ((-coords[2]) + 0.5)

                    if z_id == 2:
                        acoord = (-coords[0])
                        bcoord = (-coords[1])
                        ccoord = (-coords[2])

                    if z_id == 3:
                        acoord = (coords[0] + 0.5)
                        bcoord = ((-coords[1]) + 0.5)
                        ccoord = (coords[2] + 0.5)

                if sym == 'PNA21':
                    if z_id == 0:
                        acoord = coords[0]
                        bcoord = coords[1]
                        ccoord = coords[2]

                    if z_id == 1:
                        acoord = (-coords[0])
                        bcoord = (-coords[1])
                        ccoord = (coords[2] + 0.5)

                    if z_id == 2:
                        acoord = (coords[0] + 1/2)
                        bcoord = (-coords[1]) + 1/2
                        ccoord = coords[2]

                    if z_id == 3:
                        acoord = (-coords[0] + 0.5)
                        bcoord = (coords[1] + 0.5)
                        ccoord = (coords[2] + 0.5)

                if sym == 'P212121':
                    if z_id == 0:
                        acoord = coords[0]
                        bcoord = coords[1]
                        ccoord = coords[2]

                    if z_id == 1:
                        acoord = (-coords[0]) + 0.5
                        bcoord = (-coords[1])
                        ccoord = coords[2] + 0.5

                    if z_id == 2:
                        acoord = (-coords[0])
                        bcoord = coords[1] + 0.5
                        ccoord = (-coords[2]) + 0.5

                    if z_id == 3:
                        acoord = coords[0] + 0.5
                        bcoord = (-coords[1]) + 0.5
                        ccoord = -coords[2]

                if sym == 'P21/C':
                    if z_id == 0:
                        acoord = coords[0]
                        bcoord = coords[1]
                        ccoord = coords[2]

                    if z_id == 1:
                        acoord = -coords[0]
                        bcoord = coords[1] + 0.5
                        ccoord = -coords[2] + 0.5

                    if z_id == 2:
                        acoord = -coords[0]
                        bcoord = -coords[1]
                        ccoord = -coords[2]

                    if z_id == 3:
                        acoord = coords[0]
                        bcoord = -coords[1] + 0.5
                        ccoord = coords[2] + 0.5

                if sym == 'C2/C':
                    if z_id == 0:  #x,y,z
                        acoord = coords[0]
                        bcoord = coords[1]
                        ccoord = coords[2]

                    if z_id == 1:  #-x,y,-z+1/2
                        acoord = (-coords[0])
                        bcoord = coords[1]
                        ccoord = -coords[2] + 0.5

                    if z_id == 2:  #-x, -y, -z
                        acoord = -coords[0]
                        bcoord = -coords[1]
                        ccoord = -coords[2]

                    if z_id == 3:  #x, -y, z+1/2
                        acoord = coords[0]
                        bcoord = -coords[1]
                        ccoord = coords[2]+ 0.5

                    if z_id == 4:  #x, -y, z+1/2
                        acoord = coords[0] + 0.5
                        bcoord = coords[1] + 0.5
                        ccoord = coords[2]

                    if z_id == 5:  #-x+1/2,y+1/2,-z+1/2
                        acoord = -coords[0] + 0.5
                        bcoord = coords[1] + 0.5
                        ccoord = -coords[2] + 0.5

                    if z_id == 6:  #-x+1/2,-y+1/2,-z
                        acoord = -coords[0] + 0.5
                        bcoord = -coords[1] + 0.5
                        ccoord = -coords[2]

                    if z_id == 7:  #x+1/2,-y+1/2,z+1/2
                        acoord = coords[0] + 0.5
                        bcoord = -coords[1] + 0.5
                        ccoord = coords[2] + 0.5

                if sym == 'P-1':
                    if z_id == 0:  #x,y,z
                        acoord = coords[0]
                        bcoord = coords[1]
                        ccoord = coords[2]

                    if z_id == 1:  #-x,y,-z+1/2
                        acoord = -coords[0]
                        bcoord = -coords[1]
                        ccoord = -coords[2]

                a_s.append(copy.deepcopy(acoord))
                b_s.append(copy.deepcopy(bcoord))
                c_s.append(copy.deepcopy(ccoord))

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

        # print(a_s)
        if calc_dist == True:
            a_permol = getpermol(as_list, anum_inmol)
            b_permol = getpermol(bs_list, anum_inmol)
            c_permol = getpermol(cs_list, anum_inmol)
            atoms_permol = getpermol(atomsmol[0], anum_inmol)
        # print('a_permol', a_permol)

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
                            if adist > 0.5:
                                adist = 1.0 - adist
                            bdist = abs(b_permol[moli][atomi] - b_permol[molj][atomj])
                            if bdist > 0.5:
                                bdist = 1.0 - bdist
                            cdist = abs(c_permol[moli][atomi] - c_permol[molj][atomj])
                            if cdist > 0.5:
                                cdist = 1.0 - cdist
                            adistxyz, bdistxyz, cdistxyz = copy.deepcopy(getcartesianmol(adist, bdist, cdist, cella, cellb, cellc))

                            xvec = adistxyz[0] + bdistxyz[0] + cdistxyz[0]
                            yvec = adistxyz[1] + bdistxyz[1] + cdistxyz[1]
                            zvec = adistxyz[2] + bdistxyz[2] + cdistxyz[2]
                            dist = math.sqrt(xvec**2 + yvec**2 + zvec**2)

                            if dist < tgtdist:
                                print('mol', moli+1, 'atom', atomi+1, atmname_i, '- mol', molj+1, 'atom', atomj+1, atmname_j, "{:6.3f}".format(dist))
                                tgtflag = True
        if tgtflag == True:
            truemolid.append(i)

        # print('as_list', as_list)
        # print('len as_list', len(as_list))

        as_27 = []
        bs_27 = []
        cs_27 = []
        if image == 3:
            for xshift in [1, 0, -1]:
                for yshift in [1, 0, -1]:
                    for zshift in [1, 0, -1]:
                        as_27 = as_27 + list(np.array(as_list) + xshift)
                        bs_27 = bs_27 + list(np.array(bs_list) + yshift)
                        cs_27 = cs_27 + list(np.array(cs_list) + zshift)

        if image == 2:
            for xshift in [1, 0]:
                for yshift in [1, 0]:
                    for zshift in [1, 0]:
                        as_27 = as_27 + list(np.array(as_list) + xshift)
                        bs_27 = bs_27 + list(np.array(bs_list) + yshift)
                        cs_27 = cs_27 + list(np.array(cs_list) + zshift)

        # print(as_27, bs_27, cs_27)
        # sys.exit()

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
    for i in range(molnum):
        print('mol -- ', i, ' --')
        number = '%05d' % i
        fo = open(odir + '/' + out + number + '.xyz', 'w')
        print(len(xyzsmol[i][0]) * znum[i], file=fo)
        print('', file=fo)
        for j in range(len(xyzsmol[i])):
            for k in range(len(xyzsmol[i][j])):
                # print(i)
                print(atomsmol[0][k], xyzsmol[i][j][k][0], xyzsmol[i][j][k][1], xyzsmol[i][j][k][2], file=fo)
        fo.close()

    for i in range(molnum):
        print('mol -- ', i, ' --')
        number = '%05d' % i
        fo = open(odir + '/' + out + number + '_big.xyz', 'w')
        print(len(xyzs27mols[i]), file=fo)
        print('', file=fo)
        for j in range(len(xyzs27mols[i])):
            # for k in range(len(xyzs27mols[i][j])):
                # print(i)
            print((atomsmol[0][j % len(atomsmol[0])]), xyzs27mols[i][j][0], xyzs27mols[i][j][1], xyzs27mols[i][j][2], file=fo)
        fo.close()

    print('satisfy criteria mol:', truemolid)

#         if line[0] == '_symmetry_equiv_pos_as_xyz':
#         if line[0] =='x,y,z'
#         if line[0] =='-x+1/2,y+1/2,-z+1/2'
#         if line[0] =='-x,-y,-z'
#         if line[0] =='x+1/2,-y+1/2,z+1/2'
