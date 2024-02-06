# -*- coding: utf-8 -*-
import os
import sys
import copy
import re
import subprocess
import math
from multiprocessing import Pool
try:
    import numpy as np
except:
    pass


class molcalc():

    def getoriginpos(self, pos,cell):
        x = y = z = 0
        for i in range(len(pos)):
            x = x + float(pos[i][0])
            y = y + float(pos[i][1])
            z = z + float(pos[i][2])
        x = x / len(pos)
        y = y / len(pos)
        z = z / len(pos)
        #print x,y,z
        #print pos

        # move to origin
        for i in range(len(pos)):
            pos[i][0] = float(pos[i][0]) - x # + cell[0]/2
            pos[i][1] = float(pos[i][1]) - y # + cell[1]/2
            pos[i][2] = float(pos[i][2]) - z # + cell[2]/2
        #print pos
        return pos

    def getdist(self, p1, p2):
        dist = math.sqrt(sum((p1 - p2)**2))
        return dist

    def getdist_list(self, p1, p2):
        dist = np.linalg.norm(np.array(p1) - np.array(p2))
        return dist

    def getpos(self, pos,molnum,cell):
    #    vol=getvolume(pos)
    #    vol = vol * 1.2
        poslist=[[],[]]
        #print vol
        vec=[]
        raw=math.ceil(math.pow(sum(molnum),(1.0/3.0)))
        dist=cell[0]/raw

        #print "raw",raw,"dist",dist
        #get center pos of each pos for put
        count=0
        flag = False
        centerflag = False
        raw=int(raw)
        for i in range(raw):
            for j in range(raw):
                for k  in range(raw):
                    if centerflag == True:
                        centcount=count
                    x=i*dist
                    y=j*dist
                    z=k*dist
                    vec.append([x,y,z])
                    count=count+1
                    if count==sum(molnum):
                        flag = True
                        break
                if flag == True:
                    break
            if flag ==True:
                break

        #print "pos",pos
        count = 0
        for i in range(2):
            for j in range(molnum[i]):
                tpos=self.gettranspos(pos[i],vec[count])
                poslist[i].append(tpos)
                count += 1
        #print len(poslist[1][2])
        #for i in range(len(pos[0])):
        #    poslist[0][0][i] = [pos[0][i][0] + cell[0]/2 ,pos[0][i][1] + cell[0]/2 ,pos[0][i][2] + cell[0]/2]

        #print poslist
        return poslist

    def gettranspos(self, pos,vec):
            data=[]
            for i in range(len(pos)):
                 aaa=[pos[i][0] + vec[0] , pos[i][1] + vec[1] , pos[i][2] + vec[2] ]
                 data.append(aaa)
            return data

    def getmolradius(self, pos):
        aaa = 0
        data = []
        for i in range(len(pos)):
            aaa = math.sqrt(pow(pos[i][0], 2) + pow(pos[i][1], 2)
                            + pow(pos[i][2], 2))
            data.append(aaa)
        #print max(data)
        return max(data)

    def getvolume(self, pos):
        aaa = 0
        data = []
        for i in range(len(pos)):
            aaa = math.sqrt(pow(pos[i][0], 2) + pow(pos[i][1], 2)
                            + pow(pos[i][2], 2))
            data.append(aaa)
        #print max(data)
        return max(data)

    def getchg(self, fname,natom,molchg):
        chg=[]
    #    print fname,natom,molname
        f =open(fname, "r")
        text = f.readlines()
        for i in range(len(text)):
            itemList = text[i].split()
            if itemList == ['TWO-STAGE', 'RESP', 'FITTING:', 'SECOND', 'STAGE']:
                chg = []
                for j in range(int(natom)):
                    aaa=text[i+20+j].split()
                    chg.append(float(aaa[2]) * 18.2208933832)

            # if fmo
            if itemList == ['##', 'ESP-FITTING', 'TYPE:', 'RESP']:
                chg = []
                for j in range(int(natom)):
                    aaa=text[i+6+j].split()
                    chg.append(float(aaa[4]) * 18.2208933832)

        # print ("sumchg:", sum(chg))

        if sum(chg) != float(molchg):
            value=(sum(chg) - float(molchg)) /len(chg)
            #print value
            for i in range(len(chg)):
                chg[i] = chg[i] - value
            # print ("modchg:",sum(chg))
        return chg

    def getmolnum(self, dirname, fname, nummol_seg, repeats):
        molnum=[]

        f =open(dirname + "/" + "molnum.dat","r")
        text = f.readlines()
        for i in range(len(text)):
            itemList = text[i][:-1].split()
            if fname[0] == fname[1]:
                mol0=int(math.ceil(int(itemList[0])/float(repeats[0])))
                mol1=int(math.ceil(int(itemList[1])/float(repeats[0])))
            else:
                mol0=int(math.ceil(1*nummol_seg[0]/float(repeats[0])))
                mol1=int(math.ceil(int(itemList[1])*2/float(repeats[0])))

            molnum.append(mol0)
            molnum.append(mol1)
            print ("molnum;",molnum)
        return molnum

    def getmolmass(self, ffname,atom_list):
        data =[]
        for i in range(len(ffname)):
            for j in range(len(atom_list)):
                if ffname[i][1] == atom_list[j][0]:
                    data.append(atom_list[j][1])
        #print data
        return data

    def gettotalmass(self, mass,molnum):
        data=0
        print ("molnum",molnum)
        for i in range(len(mass)):
            data += sum(mass[i])*molnum[i]
        #print "mass=",data
        return data

    def Exportardpos(self, path, iname,  molindex, posMol, nameAtom):
        # export molindex mol data from whole posMol

        if os.path.exists(path) is False:
            print(path)
            subprocess.call(["mkdir", path])

            # print molindex
        # # Export position of mol
        out_file = path + "/" + str(iname) + ".xyz"
        # print out_file

        numlist = []
        # get total atom
        for i in molindex:
            numAtm = len(posMol[i])
            numlist.append(numAtm)
        totalAtm = sum(numlist)
        # print "totalAtom", totalAtm

        f = open(out_file, "w")
        print(totalAtm, file=f)
        print(str(iname), file=f)
        f.close()

        f = open(out_file, "a+")

        for i in molindex:
            for j in range(len(posMol[i])):
                print(nameAtom[i][j], \
                    posMol[i][j][0], posMol[i][j][1], posMol[i][j][2], file=f)
        f.close()

        # subprocess.call(["babel", "-ixyz", out_file, "-opdb",
                         # path + "/" + str(iname) + ".pdb"])

    def exportplus1pos(self, path, pos, name, atom, molindex):
        # # Export position of mol
        # head, ext = os.path.splitext(str(iname))
        if os.path.exists(path + "/pdb") is False:
            subprocess.call(["mkdir", path + "/pdb"])

        out_head = path + "/pdb/" + name
        out_file = out_head + ".xyz"
        print(out_file)

        totalnum = 0
        for i in molindex:
            totalnum += len(pos[i])
        f = open(out_file, "w")
        print(totalnum, file=f)
        print(out_file, file=f)
        f.close()

        # print pos
        f = open(out_file, "a+")

        for i in molindex:
            for j in range(len(pos[i])):
                print(atom[i][j][0:1].upper(), \
                    pos[i][j][0], pos[i][j][1], pos[i][j][2], file=f)
        f.close()

        subprocess.call(["babel", "-ixyz", out_file,
                         "-opdb", out_head + ".pdb"])

    def getatomisite(self, isitelist, typenameMol):
        # get the position of a molecule
        isitemol = []
        # print "isitelist",isitelist
        # print "typenameMol",typenameMol
        for i in typenameMol:
            for j in range(len(isitelist[0])):
                if i == isitelist[0][j]:
                    # print "match",i,isitelist[0][j],isitelist[1][j]
                    isitemol.append(isitelist[1][j])
                    break
        # print "isitemol",isitemol
        return isitemol

    def exportxyz(self, path, pos, num):
        # # Export position of mol
        # head, ext = os.path.splitext(str(iname))
        out_head = path + "/tmp" + str(num)
        out_file = out_head + ".xyz"
        print(out_file)

        f = open(out_file, "w")
        print(len(pos), file=f)
        print(out_file, file=f)
        f.close()

        # print pos
        f = open(out_file, "a+")

        for i in range(len(pos)):
            print("H", \
                pos[i][0], pos[i][1], pos[i][2], file=f)
        f.close()

        subprocess.call(["babel", "-ixyz", out_file,
                         "-opdb", out_head + ".pdb"])

    def exportdata(self, path, oname, data):
        if os.path.exists(path) is False:
            print(path)
            subprocess.call(["mkdir", path])

        out_file = path + "/" + str(oname) + ".dat"

        f = open(out_file, "w")

        for i in range(len(data)):
            print(data[i], file=f)
        f.close()



    def read_xyz(self, filename):
        lines = [l.rstrip() for l in open(filename, "U")]
        atom = []
        coord = []
        for l in range(2, len(lines)):
            temp = []
            temp.append(re.sub(' {1,}', ' ', lines[l]).split())
            atom.append(temp[0][0])
            temp[0].pop(0)
            temp2 = []
            for ll in range(len(temp[0])):
                temp2.append(float(temp[0][ll]))
            coord.append(temp2)
        return [atom, coord]


    def calcLJPairInteraction(self, pos1, pos2, param):

        # return LJ r1 r2
        sigma = param[0]
        epsilon = param[1]
    #    print "pos1",pos1
    #    print "pos2",pos2
    #    print "sigma:", sigma
    #    print "epsilon:", epsilon
        r12 = getdist(pos1, pos2)
    #    print "r12", r12

        ratio = sigma / r12
        r6 = ratio*ratio*ratio
        r6 *= r6

        LJ = 4.0 * epsilon * (r6*r6 - r6)
    #    print "LJ=",LJ
        return LJ


    def calcCoulombInteraction(self, pos1, pos2, q1, q2, epsilon):

        VALENCE = 18.220893
        r12 = getdist(pos1, pos2)

        if (epsilon != 0.0):
            # energy = 322.0637 * (q1/VALENCE) * (q2/VALENCE) / epsilon / r12_len;
            # //DREIDING paper(37)
            # energy = q1*q2/epsilon/r12_len;

            # COGNAC energy equation
            # energy = q1*q2/VALENCE/VALENCE/epsilon/r12
            energy = q1*q2/epsilon/r12
    #        print  "ENERGY", energy, "charge", q1/VALENCE ,q2/VALENCE,"length",r12
        else:
            print("dielectric constant is zero !!!")

        return energy

    def getcontactlist(self, inmol, posMol, site, neighborMol):
        contactlist = []
        #print neighborMol
        print(len(site))
        for i in range(len(neighborMol)):
            flag = False
            s1 = site[neighborMol[i][0]]
            s2 = site[neighborMol[i][1]]
            for j in range(len(s1)):
                if flag is True:
                    break
                r1 = s1[j]
                p1 = posMol[neighborMol[i][0]][j]
                for k in range(len(s2)):
                    r2 = s2[k]
                    p2 = posMol[neighborMol[i][1]][k]
                    dist = self.getdist(p1, p2)
                    if dist < (r1 + r2)/1.5:
                        contactlist.append([neighborMol[i][0], neighborMol[i][1]])
                        flag = True
                        break

        # --arrange contact list per mol--
        clistall = []
        for i in range(inmol):
            clistmol = []
            for j in contactlist:
                if j[1] == i:
                    clistmol.append(j[0])
                if j[0] == i:
                    clistmol.append(j[1])
            clistall.append(clistmol)
        # print "clistall",clistall

        # --export clistall--
        # exportdata(path[0] + "/" + path[1],clistall)

        # add center mol info for clistall
        for i in range(len(clistall)):
            clistall[i].insert(0, i)
        # print("clistall", clistall)

        return clistall


    def getrenumindex(self, index, clistall):
        # index_renum:molid arranged
        index_renum = []
        for i in range(len(index)):
            index_renum.append(i)
        #print "index", index
        # print("index_renum", index_renum)

        for i in range(len(clistall)):
            for j in range(len(clistall[i])):
                for k in range(len(index)):
                    if clistall[i][j] == index[k]:
                        clistall[i][j] = index_renum[k]

        # print("contact_list_renum", clistall)
        return index_renum, clistall


    def getrenumfrag(self, index, clistall,fragids):
        flag = False
        count = 0
        aaa = []
        for k in range(len(fragids)):
            flag = False
            for i in range(len(fragids[k])):
                if flag == True:
                    break
                for j in range(len(index)):
                    if index[j] == fragids[k][i]:
                        #print j, fragids[k][i]
                        # print "mol",str(k), "exist"
                        aaa.append(k)
                        flag = True
                        # count += 1
                        break
        #print "index", index
        # print "aaa", aaa # used mol
        # print "len", len(aaa)

        aaa_frag = []
        for i in aaa:
            aaa_frag.append(fragids[i])
        # print aaa_frag

        count = 0
        frags_renum = []
        for i in range(len(aaa_frag)):
            aaa_frags = []
            for j in range(len(aaa_frag[i])):
                count += 1
                aaa_frags.append(count)
            frags_renum.append(aaa_frags)
        # print "frags_renum", frags_renum

        convlist = []
        for i in index:
                for k in range(len(aaa_frag)):
                    for l in range(len(aaa_frag[k])):
                        if aaa_frag[k][l] == i:
                            convlist.append([aaa_frag[k][l], frags_renum[k][l]])

        clist_renum = copy.deepcopy(clistall)
        for i in range(len(clistall)):
            for j in range(len(clistall[i])):
                for k in range(len(convlist)):
                        if clistall[i][j] == convlist[k][0]:
                            clist_renum[i][j] = convlist[k][1] - 1

        # print("contact_listfrag_renum", clist_renum)

        # index_renum:molid arranged
        index_renum = copy.deepcopy(index)

        for i in range(len(index)):
            for k in range(len(convlist)):
                    if index[i] == convlist[k][0]:
                        index_renum[i] = convlist[k][1] - 1

        # print("index_renum", index_renum)

        return index_renum, clist_renum

    def getindex(self, clistall):
        # --get used mol index --
        index = []
        for i in range(len(clistall)):
            for j in clistall[i]:
                flag = False
                for k in index:
                    if j == k:
                        flag = True
                if flag is False:
                    index.append(j)

        index.sort()
        # print("index = ", index)

        return index

    def getcontactfrag(self, clist, posMol, site, fragids, infrag):
    # clist i:mol j: contact to clist[i][0]
    # posfrag_mols i:molid j:fragid k:atomid l:x,y,z(3)
    # site i:molid j:fragid k:atomid
        contactlists = []
        for i in range(len(clist)):
            for j in range(1,len(clist[i])):
                flag = False
                s1 = site[clist[i][0]]
                s2 = site[clist[i][j]]
                # print clist[i][0], clist[i][j]
                for k in range(len(s1)):
                    flag = False
                    for m in range(len(s2)):
                        flag = False
                        # print "m =", m
                        # print "1"
                        for l in range(len(s1[k])): # k: fragid l:atomid
                            if flag == True:
                                break
                            r1 = s1[k][l]
                            p1 = posMol[clist[i][0]][k][l]
                            for n in range(len(s2[m])): # m; fragid n:atomid
                                r2 = s2[m][n]
                                p2 = posMol[clist[i][j]][m][n]
                                dist = self.getdist(p1, p2)
                                # print r1, r2
                                # print "2"
                                if dist < (r1 + r2)/1.5:
                                    contactlists.append([clist[i][0], k, clist[i][j], m])
                                    # print [clist[i][0], k, clist[i][j], m]
                                    # print "break"
                                    flag = True
                                    break
        #return contactlist

        contactfrag = []
        for alist in contactlists:
            # print alist
            # print [fragids[alist[0]][alist[1]], fragids[alist[2]][alist[3]]]
            contactfrag.append([fragids[alist[0]][alist[1]], fragids[alist[2]][alist[3]]])
        # print("contactfrag", contactfrag)

        # --arrange contact list per mol--
        clistall = []
        for i in range(1,infrag+1):
            clistmol = []
            for j in contactfrag:
                flag = False
                if j[1] == i:
                    for k in clistmol:
                        if k == j[0]:
                            flag = True
                    if flag is False:
                        clistmol.append(j[0])
                if j[0] == i:
                    for k in clistmol:
                        if k == j[1]:
                            flag = True
                    if flag is False:
                        clistmol.append(j[1])
            clistall.append(clistmol)
        # print("clistfrag_all",clistall)

        # --export clistall--
        # exportdata(path[0] + "/" + path[1],clistall)

        # add center mol info for clistall
        for i in range(len(clistall)):
            clistall[i].insert(0, i+1)
        # print("clistfrag_all", clistall)

        return clistall

    def getCenter(self, posVec):
        # get center coordinates
        # print "posVec",posVec
        center = np.average(posVec, 0)
        # print "center",center
        return center

    def moveMolTrans(self, posVec, transVec):
        # Parallel shift
        posVec = posVec + transVec
        return posVec

    def moveMolEuler(self, posVec, rotatRdn):
        # ## Rotation (specify the quaternion, Euler angle ZXZ '[rad])
        a = np.cos(rotatRdn[1]/2.)*np.cos((rotatRdn[2]+rotatRdn[0])/2.)
        b = np.sin(rotatRdn[1]/2.)*np.cos((rotatRdn[2]-rotatRdn[0])/2.)
        c = np.sin(rotatRdn[1]/2.)*np.sin((rotatRdn[2]-rotatRdn[0])/2.)
        d = np.cos(rotatRdn[1]/2.)*np.sin((rotatRdn[2]+rotatRdn[0])/2.)
        A = np.array([[-c*c+b*b-d*d+a*a, 2.0*(d*a-c*b), 2.0*(b*d+c*a)],
                      [-2.0*(c*b+d*a), c*c-b*b-d*d+a*a, 2.0*(b*a-c*d)],
                      [2.0*(b*d-c*a), -2.0*(c*d+b*a), -c*c-b*b+d*d+a*a]])
        numPos = len(posVec)
        for i in range(numPos):
            posVec[i] = np.dot(A, posVec[i])
        return posVec

    def moveMolRotat(self, posVec, rotatDeg):
        # Rotation (specify the Cardin angle XYZ [deg])
        if rotatDeg[1] == 90 or rotatDeg[1] == -90:
            a = 0
        rotatRad = np.array([2*pi/360*rotatDeg[0],
                             2*pi/360*rotatDeg[1],
                             2*pi/360*rotatDeg[2]])
        a = rotatRad[0]
        b = rotatRad[1]
        c = rotatRad[2]
        A = np.array([[np.cos(c)*np.cos(b),
                       np.cos(c)*np.sin(b)*np.sin(a)-np.sin(c)*np.cos(a),
                       np.cos(c)*np.sin(b)*np.cos(a)+np.sin(c)*np.sin(a)],
                      [np.sin(c)*np.cos(b),
                       np.sin(c)*np.sin(b)*np.sin(a)+np.cos(c)*np.cos(a),
                       np.sin(c)*np.sin(b)*np.cos(a)-np.cos(c)*np.sin(a)],
                      [-np.sin(b), np.cos(b)*np.sin(a), np.cos(b)*np.cos(a)]])
        numPos = len(posVec)
        for i in range(numPos):
            posVec[i] = np.dot(A, posVec[i])
        return posVec

    def getdummyatom(self, connect, end):
        flag = False
        for i in range(len(connect)):
            if connect[i][1] == end:
                for j in range(2, len(connect[i])):
                    for k in range(len(connect)):
                        if connect[i][j] == connect[k][1]:
                            # print connect[i][j], connect[k][1]
                            flag = True
                        if flag is False:
                            return connect[i][j]

        return

    def getangle(self, p1, p2, p3, length):
        # p1 - p2 - p3 angle
        p1 = p1 - p2
        p3 = p3 - p2
        # print p1, p2, p3
        angle = (p1[0]*p3[0]+p1[1]*p3[1]+p1[2]*p3[2])/(length[0]*length[1])
        # print ('angle',angle)
        # print (angle)
        try:
            theta = math.acos(angle) * 180 / math.pi
        except:
            theta = 0.0
        return theta

    def getnppos(self, pos):
        for i in range(len(pos)):
            for j in range(len(pos[i])):
                pos[i][j] = float(pos[i][j])
        pos = np.array(pos)
        return pos


    def get2vec(self, pos1, pos2, p2l1, p1l2, p1dum2, plus):
        vec = pos1[p1l2] - pos2[p2l1]
        vec2 = pos1[p1dum2] - pos1[p1l2]
        dist = self.getdist(pos1[p1dum2], pos1[p1l2])
        vec2 = vec2 / dist * (dist + plus)
        return vec, vec2


    def xymatch(self, pos, pos2, l1, l2, dum1):
        p1 = copy.deepcopy(pos[l2])
        p2 = copy.deepcopy(pos2[l1])
        p3 = copy.deepcopy(pos2[dum1])
        p1[2] = 0
        p2[2] = 0
        p3[2] = 0
        d1 = self.getdist(p1, p2)
        d2 = self.getdist(p3, p2)
        # print d1, d2
        theta = self.getangle(p1, p2, p3, [d1, d2])
        # print (theta)

        if theta == 0.0:
            theta = 0.00001

        aaa = self.rotate_ardz(theta, pos2)
        bbb = self.rotate_ardz(-theta, pos2)
        value1 = abs(p1[0]/p1[1] - aaa[dum1, 0]/aaa[dum1, 1])
        value2 = abs(p1[0]/p1[1] - bbb[dum1, 0]/bbb[dum1, 1])
        # print (value1, value2)
        if (value2 > value1):
            rotated = copy.deepcopy(aaa)
        else:
            rotated = copy.deepcopy(bbb)

        return rotated


    def rotate_ardz(self, theta, pos):
        # ------回転角の定義------
        rot = math.pi * theta/180  # ラジアンに変換
        # ------回転行列の成分の定義---------
        rotate_z = np.arange(9.0).reshape((3, 3))

        rotate_z[0, 0] = math.cos(rot)           # z軸周りに回転させる回転行列
        rotate_z[1, 0] = math.sin(rot)
        rotate_z[2, 0] = 0.0
        rotate_z[0, 1] = -(math.sin(rot))
        rotate_z[1, 1] = math.cos(rot)
        rotate_z[2, 1] = 0.0
        rotate_z[0, 2] = 0.0
        rotate_z[1, 2] = 0.0
        rotate_z[2, 2] = 1.0
        # ------回転処理----------------

        print ("theta", theta)
        # print rotate_z
        rotpos = np.arange(len(pos) * 3.0).reshape(-1, 3)
        rotpos[:] = 0.0
    #    print rotpos
        for i in range(len(pos)):
            for j in range(3):
                for k in range(3):
                    rotpos[i][j] = rotpos[i][j] + rotate_z[j, k] * pos[i, k]

        # print rotpos
        return rotpos


    def label1match(self, pos2in, pos1):

        normal = np.arange(3.0)
        normal[0] = pos2in[dum1, 1] / math.sqrt(pos2in[dum1, 0]**2 +
                                                pos2in[dum1, 1]**2)
        normal[1] = -(pos2in[dum1, 0] / math.sqrt(pos2in[dum1, 0]**2 +
                                                  pos2in[dum1, 1]**2))
        normal[2] = 0
        # print normal
        # refと、回転させてるmolの原子1のラベルを一致させる。xy平面上に投影した二点と原点を通る直線の法線ベクトルまわりに回転
        # 今回はラベル1の原子の座標で見る
        d1 = self.getdist(pos1[l2], pos2in[l1])
        d2 = self.getdist(pos2in[dum1], pos2in[l1])
        theta = getangle(pos1[l2], pos2in[l1], pos2in[dum1], [d1, d2])
        # print ("step2 angle=", theta)

        pos2out1 = rotate_ardvec(theta, normal, pos2in)
        d1 = self.getdist(pos1[l2], pos2out1[l1])
        d2 = self.getdist(pos2out1[dum1], pos2out1[l1])
        v1 = (pos1[l2] - pos2out1[l1]) / d1
        v2 = (pos2out1[dum1] - pos2out1[l1]) / d2
        value1 = self.getdist(v1, v2)

        pos2out2 = rotate_ardvec(-theta, normal, pos2in)
        d1 = self.getdist(pos1[l2], pos2out2[l1])
        d2 = self.getdist(pos2out2[dum1], pos2out2[l1])
        v1 = (pos1[l2] - pos2out2[l1]) / d1
        v2 = (pos2out2[dum1] - pos2out2[l1]) / d2
        value2 = self.getdist(v1, v2)

        # print (value1, value2)

        if (value2 > value1):
            rotated = copy.deepcopy(pos2out1)
        else:
            rotated = copy.deepcopy(pos2out2)

        return rotated

    def rotate_ardvec(self, theta, vec, pos):
        # print ("vec", vec)
        # ------回転角の定義------
        rot = math.pi * theta/180  # ラジアンに変換
        # ------回転行列の成分の定義---------
        matrix = np.arange(9.0).reshape((3, 3))

        # ------回転行列の成分の定義---------
        matrix[0, 0] = vec[0]**2*(1-math.cos(rot))+math.cos(rot)
        matrix[1, 0] = vec[0]*vec[1]*(1-math.cos(rot))+vec[2]*math.sin(rot)
        matrix[2, 0] = vec[2]*vec[0]*(1-math.cos(rot))-vec[1]*math.sin(rot)
        matrix[0, 1] = vec[0]*vec[1]*(1-math.cos(rot))-vec[2]*math.sin(rot)
        matrix[1, 1] = vec[1]**2*(1-math.cos(rot))+math.cos(rot)
        matrix[2, 1] = vec[1]*vec[2]*(1-math.cos(rot))+vec[0]*math.sin(rot)
        matrix[0, 2] = vec[2]*vec[0]*(1-math.cos(rot))+vec[1]*math.sin(rot)
        matrix[1, 2] = vec[1]*vec[2]*(1-math.cos(rot))-vec[0]*math.sin(rot)
        matrix[2, 2] = vec[2]**2*(1-math.cos(rot))+math.cos(rot)

        # print ("theta", theta)
        # print rotate_z
        rotpos = np.arange(len(pos) * 3.0).reshape(-1, 3)
        rotpos[:] = 0.0
    #    print rotpos
        for i in range(len(pos)):
            for j in range(3):
                for k in range(3):
                    rotpos[i][j] = rotpos[i][j] + matrix[j, k] * pos[i, k]

        # print rotpos
        return rotpos


    def dihed_rotate(self, pos2in, pos1, l1, l2):
        vec = np.arange(3.0)

    #   normalization vec
        lvec = self.getdist(pos1[l1[1]], pos2in[l2[0]])
        vec = pos1[l1[1]] / lvec

        # get crossproduct
        cross1 = self.getcrossprod(pos1[l1[1]], pos2in[l2[0]], pos1[l1[0]])
        cross2 = self.getcrossprod(pos1[l1[1]], pos2in[l2[0]], pos2in[l2[1]])

        # get angle cross1-origin-cross2
        d1 = self.getdist(cross1, pos2in[l2[0]])
        d2 = self.getdist(cross2, pos2in[l2[0]])
        theta = self.getangle(cross1, pos2in[l2[0]], cross2, [d1, d2])
        # print (theta)

        theta = 180 - theta
        # out = rotate_ardvec(180 - theta, vec, pos2in)

        pos2out1 = self.rotate_ardvec(theta, vec, pos2in)
        value1 = self.getdist(pos1[l1[0]], pos2out1[l2[1]])

        pos2out2 = self.rotate_ardvec(-theta, vec, pos2in)
        value2 = self.getdist(pos1[l1[0]], pos2out2[l2[1]])

        # print ("*********", value1, value2)
        if (value2 < value1):
            rotated = copy.deepcopy(pos2out1)
        else:
            rotated = copy.deepcopy(pos2out2)

        return rotated


    def dum1_rotate(self, pos2in, pos1, l1, l2, dum1):
        # rotate around crossproduct of pos1[l2]-pos2in[l1]-pos2in[dum1]
        vec = np.arange(3.0)
        # get cross product
        cross1 = self.getcrossprod(pos1[l2], pos2in[l1], pos2in[dum1])

        # normalize
        dd = self.getdist(cross1, pos2in[l1])
        vec = cross1 / dd

        # get rotate angle
        d1 = self.getdist(pos1[l2], pos2in[l1])
        d2 = self.getdist(pos2in[dum1], pos2in[l1])
        theta = self.getangle(pos1[l2], pos2in[l1], pos2in[dum1], [d1, d2])
        # print (theta)

        # rotate
        pos2out1 = self.rotate_ardvec(theta, vec, pos2in)
        d1 = self.getdist(pos1[l2], pos2out1[l1])
        d2 = self.getdist(pos2out1[dum1], pos2out1[l1])
        v1 = (pos1[l2] - pos2out1[l1]) / d1
        v2 = (pos2out1[dum1] - pos2out1[l1]) / d2
        value1 = self.getdist(v1, v2)

        pos2out2 = self.rotate_ardvec(-theta, vec, pos2in)
        d1 = self.getdist(pos1[l2], pos2out2[l1])
        d2 = self.getdist(pos2out2[dum1], pos2out2[l1])
        v1 = (pos1[l2] - pos2out2[l1]) / d1
        v2 = (pos2out2[dum1] - pos2out2[l1]) / d2
        value2 = self.getdist(v1, v2)

        # print (value1, value2)

        # check
        if (value2 > value1):
            rotated = copy.deepcopy(pos2out1)
        else:
            rotated = copy.deepcopy(pos2out2)

        return rotated


    def babelxyzpdb(self, head):
        cmd = "babel -ixyz " + head + ".xyz -opdb " + head + ".pdb"
        subprocess.call(cmd.split(" "))


    def writexyzpoly(self, head, natomsum, atoms, pos, dums, seq, cnct_count, cncted_count):
        # print atoms
        print (dums)
        print (natomsum)
        f = open(head + ".xyz", "w")
        print (natomsum, file=f)
        print (" @@ " + head + " @@", file=f)
        for i in range(len(pos)):
            molid = seq[i]
            flag = cnct_count[i]
            flag2 = cncted_count[i]
            print (molid, flag, flag2)
            for j in range(len(pos[i])):
                if j == dums[molid][0] and flag != 0:
                    # print "CONT!"
                    continue
                skipflag = False
                for k in range(1, len(dums[molid])):
                    if j == dums[molid][k] and flag2 != 0:
                        # print "CONTi!"
                        skipflag = True
                if skipflag is True:
                    continue
                print (atoms[molid][j][1],  pos[i][j][0], \
                    pos[i][j][1], pos[i][j][2], file=f)

        f.close()

    def getcrossprod(self, p1, p2, p3):
        v1 = p1 - p2
        v2 = p3 - p2
        cross = np.arange(3.0)
        cross[0] = v1[1] * v2[2] - v1[2] * v2[1]
        cross[1] = v1[2] * v2[0] - v1[0] * v2[2]
        cross[2] = v1[0] * v2[1] - v1[1] * v2[0]

        return cross

    def getrepeatpos(self, pos_orig, l1, l2, dum1, dum2, d_plus):
        # target l1 and dum2 + 1.5 vec
        vec, vec2 = self.get2vec(pos_orig, l1, l2, dum2, d_plus)
        # print "vec2", vec2
        pos2_orig = pos_orig + vec + vec2
        pos1 = pos_orig - pos2_orig[l1]
        pos2 = pos2_orig - pos2_orig[l1]

        # step1: rotate pos2[dum1] to pos1[l2]
        pos2_dum1_rot = self.dum1_rotate(pos2, pos1, l1, l2)

        # step2: dihedral rotate
        pos2_dihed_rot = self.dihed_rotate(pos2_dum1_rot, pos1, l1, l2)

        return pos2_dihed_rot + pos2_orig[l1]

