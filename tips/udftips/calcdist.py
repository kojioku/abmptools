import UDFManager
import argparse
import os
import numpy as np


def getargs():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('-i', '--input',
                           type=str, required=True,
                           help='input dpm')

    argparser.add_argument('-t', '--target',
                           type=int, required=True,
                           nargs=4,
                           help='input dpm')

    args = argparser.parse_args()
    return args


def getudfinfo(uobj):
    uobj.jump(0)
    totalMol = uobj.size("Set_of_Molecules.molecule[]")
    totalRec = uobj.totalRecord()
    cell = uobj.get("Structure.Unit_Cell.Cell_Size")

    return totalMol, totalRec, cell


if __name__ == "__main__":
    args = getargs()
    infile = args.input
    target = args.target
    print('Input:', infile)
    print('Target mol{0} atom{1}, mol{2} atom{3}'.format(
        target[0], target[1], target[2], target[3]))
    uobj = UDFManager.UDFManager(infile)

    totalmol, totalrec, cell = getudfinfo(uobj)

    print('Total Mol:', totalmol)
    print('Total Rec:', totalrec)
    print('Cell:', cell)

    # Calculate distance
    dist = []
    for i in range(totalrec):
        uobj.jump(i)
        pos1 = uobj.get("Structure.Position.mol[].atom[]",
                        [target[0], target[1]])
        pos2 = uobj.get("Structure.Position.mol[].atom[]",
                        [target[2], target[3]])

        dist_scalar = [0, 0, 0]
        for j in range(3):
            tmpdist = abs(pos1[j] - pos2[j])
            if tmpdist > cell[j]:
                # print('Warning: distance is larger than cell size')
                while tmpdist > cell[j]:
                    tmpdist -= cell[j]
            if tmpdist < cell[j]/2.:
                dist_scalar[j] = tmpdist
            else:
                # print('Warning: distance is larger than half of cell size')
                dist_scalar[j] = cell[j] - tmpdist

        dist.append(np.linalg.norm(np.array(dist_scalar)))

    # write result
    outname = os.path.splitext(infile)[0] + '-dist' + \
        str(target[0]) + '-' + str(target[1]) + '_' + \
        str(target[2]) + '-' + str(target[3]) + '.csv'
    f = open(outname, 'w')
    print('Rec, Distance(sigma)', file=f)
    for i, dist in enumerate(dist):
        print('{:>4d}, {:>12.6f}'.format(i, dist), file=f)
    f.close()
    # ファイル中身を表示
    with open(outname, 'r') as f:
        print(f.read())

    print('------------------\nOutput:', outname)
