from UDFManager import *
import sys
import abmptools
import numpy as np


def gettotalmol_rec(uobj):
    totalMol = uobj.size("Set_of_Molecules.molecule[]")
    totalRec = uobj.totalRecord()
    return totalMol, totalRec


def gettotalAtm(uobj):
    totalMol = uobj.size("Set_of_Molecules.molecule[]")
    numlist = []
    for i in range(totalMol):
        numAtm = uobj.size("Set_of_Molecules.molecule[].atom[]", [i])
        numlist.append(numAtm)
    totalAtm = sum(numlist)
    return totalAtm


def getcellsize(uobj, rec):
    uobj.jump(rec)
    cell = uobj.get("Structure.Unit_Cell.Cell_Size")
    return cell


def getudfinfowrap(uobj):
    totalMol, totalRec = gettotalmol_rec(uobj)
    totalAtm = gettotalAtm(uobj)
    cell = getcellsize(uobj, totalRec-1)

    print('totalMol:', totalMol)
    print('totalAtm:', totalAtm)
    print('cellinfo:', cell)

    return totalMol, totalRec, totalAtm, cell


def gettotalAtm(uobj):
    totalMol = uobj.size("Set_of_Molecules.molecule[]")
    numlist = []
    for i in range(totalMol):
        numAtm = uobj.size("Set_of_Molecules.molecule[].atom[]", [i])
        numlist.append(numAtm)
    totalAtm = sum(numlist)
    return totalAtm


def getmolname(i, uobj):
    # --get used mol infomation--
    molname = uobj.get("Set_of_Molecules.molecule[" +
                       str(i) + "].Mol_Name")
    return molname


def getposmolrec(uobj, indexMol, record):
    # UDF operation
    # get the position of a molecule
    uobj.jump(record)
    i = indexMol
    aaa = uobj.get("Structure.Position.mol[" + str(i) + "].atom[]")
    position = np.array(aaa)
    return position

def getvelmolrec(uobj, indexMol, record):
    # UDF operation
    # get the position of a molecule
    uobj.jump(record)
    i = indexMol
    aaa = uobj.get("Structure.Velocity.mol[" + str(i) + "].atom[]")
    velocity = np.array(aaa)
    return velocity


def getCenter(posVec):
    # get center coordinates
    center = np.average(posVec, 0)
    return center


def getdist_list(p1, p2):
    dist = np.linalg.norm(np.array(p1) - np.array(p2))
    return dist


if __name__ == "__main__":

    fname = sys.argv[1]
    # tgtname = "Cholesterol"
    tgtname = "POPC"
    tgtmolid = 10

    print('target:', tgtname, "mol", tgtmolid)
    uobj = UDFManager(fname)

    totalmol, totalrec, totalatm, cell = getudfinfowrap(uobj)

    # get molname list
    molnames = []
    for i in range(totalmol):
        molnames.append(getmolname(i, uobj))
    print(totalmol, totalrec)
    # print(molnames)

    # get target mol index
    targets = []
    for i in range(totalmol):
        if molnames[i] == tgtname:
            targets.append(i)
    print('target id', targets)


    # get target posmol for all rec
    poss = []
    mol = targets[tgtmolid]
    print("--- target mol id:", mol, "---")

    for record in range(totalrec):
        poss.append(getposmolrec(uobj, mol, record))
    # print('poss0', poss[0])

    vels = []
    for record in range(totalrec):
        vels.append(getvelmolrec(uobj, mol, record))

    # centerOfMol: com of each molecules
    cocs = []
    for i in range(len(poss)):
        cocs.append(getCenter(poss[i]).tolist())
    # print('cocs', cocs)
    velmols = []
    for i in range(len(vels)):
        velmols.append(getCenter(vels[i]).tolist())


    fn1 = 'coc_mol' + str(tgtmolid) + '.csv'
    fn2 = 'vel_mol' + str(tgtmolid) + '.csv'
    fn3 = 'dist_mol' + str(tgtmolid) + '.csv'
    f1 = open(fn1, 'w')
    f2 = open(fn2, 'w')
    f3 = open(fn3, 'w')

    print('rec x y z', file=f1)
    for i in range(len(cocs)):
        print(str(i), '{0[0]:.3f} {0[1]:.3f} {0[2]:.3f}'.format(cocs[i]), file=f1)

    print('rec x y z', file=f2)
    for i in range(len(vels)):
        print(str(i), '{0[0]:.3f} {0[1]:.3f} {0[2]:.3f}'.format(velmols[i]), file=f2)

    # getdist
    dists = []
    for i in range(totalrec-1):
        dists.append(getdist_list(cocs[i], cocs[i+1]))
    for i in range(len(dists)):
        print(str(i) + '-' + str(i+1), '{:.3f}'.format(dists[i]), file=f3)

    print('out:', fn1, fn2, fn3)
    f1.close()
    f2.close()
    f3.close()
    # get cholesterol mol pos at one record

#     countss = []
#     for record in range(totalrec):
#         if record % 10 == 0:
#             print('rec', record)
#         tgtposs = []
#         for mol in targets:
#             tgtposs.append(cobj.getposmolrec(uobj, mol, record))
#         # print('tgtposs0', tgtposs[0])
# 
#         # centerOfMol: com of each molecules
#         tgtcocs = []
#         for i in range(len(tgtposs)):
#             tgtcocs.append(cobj.getCenter(tgtposs[i]).tolist())
#         # print('tgtcocs', tgtcocs)
# 
#         for k in range(3):
#             # print('take', k)
#             for i in range(len(tgtcocs)):
#                 for j in range(3):
#                     if tgtcocs[i][j] > cell[0]:
#                         # print("big", tgtcocs[i][j])
#                         tgtcocs[i][j] -= cell[0]
#                     if tgtcocs[i][j] < 0:
#                         # print("small", tgtcocs[i][j])
#                         tgtcocs[i][j] += cell[0]
# 
#         # # check contact
#         counts = []
#         for i in range(len(tgtposs)):
#             count = 0
#             for j in range(len(tgtposs)):
#                 if i == j:
#                     continue
#                 dist = cobj.getdist_list(tgtcocs[i], tgtcocs[j])
#                 # print(i, j, dist)
#                 if dist <= contactdist:
#                     # print(i, j, "contact")
#                     count += 1
#             counts.append(count)
#         countss.append(counts)
# 
#     onecount = []
#     for i in range(len(countss)):
#         onecount.append(countss[i][tgtmolid])
#     print("number_mol_around", onecount)
