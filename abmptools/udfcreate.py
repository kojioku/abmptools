# -*- coding: utf-8 -*-
import os
import sys
import copy
import re
import subprocess
import math
from multiprocessing import Pool
try:
    from UDFManager import *
except:
    pass
import time

class udfcreate():
    def __init__(self):
        self.algo = ['NPT_Andersen_Kremer_Grest']
        self.cellsize = [30, 30, 30]
        self.totalstep = 20000
        self.outstep = 200
        self.tempes = [300]


        # self.fname=""
        # self.atom=""
        # self.ffname=""
        # self.molname=""
        # self.bondff=""
        # self.bond=""
        # self.angleff=""
        # self.angle=""
        # self.torsionff=""
        # self.torsion=""
        # self.chg=""

        # self.cellsize=""
        # self.ljparam=""
        # self.bondparam=""
        # self.angleparam=""
        # self.torsionparam=""
        # self.atom_list=""
        # self.totalstep=""
        # self.outstep=""
        # self.totalmass=""
        # self.algo=""
        # self.poslist=""

    def setudfparam(self, param_udf):
        # -- read for paramdata --
        self.algo = param_udf['algo']

        try:
            self.nvtalgo = param_udf['nvtalgo']
        except KeyError:
            self.nvtalgo = param_udf['NVT_Nose_Hoover']

        self.cellsize = param_udf['cellsize']
        self.totalstep = param_udf['totalstep']
        self.outstep = param_udf['outstep']

        try:
            self.nvtstep = param_udf['nvtstep']
        except KeyError:
            self.nvtstep = 1000

        try:
            self.nvtoutstep = param_udf['nvtoutstep']
        except KeyError:
            self.nvtoutstep = 10

        try:
            self.timestep = param_udf['timestep']
        except KeyError:
            self.timestep = 0.001 # 1fs

        try:
            self.pressure = param_udf['pressure']
        except KeyError:
            self.pressure = 1.0 # 1atm

        self.tempes = param_udf['temperature']
        self.octahome = param_udf['octahome']
        self.cognacpath = param_udf['cognacpath']
        try:
            self.mdlocalnp = param_udf['mdlocalnp']
        except KeyError:
            self.mdlocalnp = 24

        try:
            self.mdnp = param_udf['mdnp']
        except KeyError:
            self.mdnp = 48



    def getconnectdata(self, fname):
        data=[]
        atomdata=[]
        for line in open(fname, 'r'):
                itemList = line[:-1].split()
                #print itemList[0]
                if itemList[0]=="CONECT" and len(itemList)>=4:
                    #print (len(itemList))
                    data.append(itemList)
                #if itemList[0] in ["HETATM","ATOM"]:
                   #atomdata.append([itemList[1],itemList[-1]])
        for i in range(len(data)):
            for j in range(1,len(data[i])):
                data[i][j] = int(data[i][j]) -1
        #print "connect data", data
        return data

    def getbatdata(self, data):
        #get bond info
        #print "connectdata",data
        bond=[]
        for i in range(len(data)):
            #print (i)
            for j in range(2,len(data[i])):
                #print j,len(i)
                #print "bond:",data[i][1] + "-" + data[i][j]
                bond.append([data[i][1],data[i][j]])
        #print "bond",bond

        #get angle info
        angle=[]
        for i in range(len(data)):
            for j in range(2,len(data[i])):
                for k in range(j+1,len(data[i])):
                    #print "angle:",data[i][j],data[i][1],data[i][k]
                    angle.append([data[i][j],data[i][1],data[i][k]])
        #print "angle",angle


        torsion=[]
        for i in range(len(data)-1):
            flag = False
            for j in range(i + 1, len(data)):
                for k in range (2,len(data[i])):
                    if data[i][k] in data[j]:
                        flag = True
            if flag == True:

                for j in range(2, len(data[i])):
                    for k in range (2,len(data[i+1])):
                        if data[i][j] in data[i+1] or data[i+1][k] in data[i]:
                            #print "(Duplication:",data[i][j],"-",data[i][1],"-",data[i+1][1],"-",data[i+1][k],")"
                            continue
                        else:
                            #print "torsion:",data[i][j],"-",data[i][1],"-",data[i+1][1],"-",data[i+1][k]
                            torsion.append([data[i][j],data[i][1],data[i+1][1],data[i+1][k]])

        # print ("torsionmap",torsion)
        return bond, angle, torsion


    def getffname(self, atom, xyzfile, fffile):
        # xyzfile="monomer/pdb/" + molname + ".xyz"
        # fffile="monomer/" + molname + ".ff"

        if os.path.exists(fffile) is False:
            outff = open(fffile, 'w')
            cmd = "obminimize -n 1 -ff GAFF -onul " + xyzfile
            ps = subprocess.Popen(cmd, shell=True, stderr=outff)
            ps.wait()
            outff.close()

        ffname = []
        i = 0
        # print('get', fffile)
        # time.sleep(2)
        for line in open(fffile, 'r'):
            itemList = line[:-1].split()
            i += 1
            if i >= 5:
                ffname.append([int(itemList[0])-1, itemList[1]])
            if i == len(atom) + 4:
                break
        return ffname

#         if os.path.exists(fffile) is False:
#             out = open(fffile, "w")
#             cmd = "obminimize -n 1 -ff GAFF -onul " + xyzfile
#             subprocess.call(cmd.split(" "),stderr=out)
#             out.close()
#             time.sleep(5)
# 
#         ffname=[]
#         i=0
#         for line in open(fffile, 'r'):
#             itemList = line[:-1].split()
#             #print itemList
#             i=i+1
#             if i >= 5:
#                 ffname.append([int(itemList[0])-1,itemList[1]])
#             if i == len(atom) + 4:
#                 break
#         return ffname

    def getbatff(self, ffname,bond):
        fff=[]
        #print "bond",bond
        #print "ffname",ffname
        for i in range(len(bond)):
            dbuf=[]
            for j in range(len(bond[i])):
                for k in range(len(ffname)):
                    if bond[i][j] == ffname[k][0]:
                        #print ffname[k][1]
                        dbuf.append(ffname[k][1])
            fff.append(dbuf)
        #print "fff",fff
        return fff

    def getfflist(self, infile):
        data=[]
        for line in open(infile, 'r'):
            if line[0] != "#":
                data.append(eval(line))
        return data


    def gettorsionfflist(self, infile):
        data=[]
        data1=[]
        data2=[]
        data3=[]
        data4=[]
        flag=0
        for line in open(infile, 'r'):
            aaa= line.split()
            if aaa[0] == '#"X-@-@-@"':
                flag = 1
            if aaa[0] == '#"X-X-@-@"':
                flag = 2
            if aaa[0] == '#"@-@-@-@"':
                flag = 3
            if line[0] == "#":
                continue
            if flag==0:
                data1.append(eval(line))
            if flag==1:
                data2.append(eval(line))
            if flag==2:
                data3.append(eval(line))
            if flag==3:
                data4.append(eval(line))
        data.append(data1)
        data.append(data2)
        data.append(data3)
        data.append(data4)
        return data


    def getffid(self, batff,atom_list):
        fff=[]
        for i in range(len(batff)):
            dbuf=[]
            for j in range(len(batff[i])):
                for k in range(len(atom_list)):
                    if batff[i][j] == atom_list[k][0]:
                        #print ffname[k][1]
                        dbuf.append(atom_list[k][2])
            fff.append(dbuf)
        return fff


    def getljparam(self, ffname,atom_list):
        data=[]
        for i in range(len(ffname)):
            for j in range(len(ffname[i])):
                for k in range(len(atom_list)):
                    #print ffname[i][1], atom_list[j][0]
                    if ffname[i][j][1] == atom_list[k][0]:
                        data.append(atom_list[k])
    #    print "data",data

        dellist=[]
    #    print len(data)
        for i in range(len(data)-1):
            for j in range(i+1,len(data)):
                if data[i][0] == data[j][0]:
                    dellist.append(i)
                    break

    #    print "dellist",dellist
        for i in range(len(dellist)):
    #        print dellist[-(i+1)]
            data.pop(dellist[-(i+1)])
    #    print data
    #    print ffname
    #    print atom_list
        return data


    def getbondffparam(self, bondffid,ff_list,bond,bondff):
    #    print bondffid
        iddata=[]
        fff=[]
        for i in range(len(bondffid)):
            for k in range(len(ff_list)):
                if bondffid[i][0] == ff_list[k][0] and bondffid[i][1] == ff_list[k][1]:
                    fff.append(ff_list[k])
                    break
                if bondffid[i][0] == ff_list[k][1] and bondffid[i][1] == ff_list[k][0]:
                        #print ffname[k][1]
                    fff.append(ff_list[k])
                    bondffid[i][0],bondffid[i][1]=bondffid[i][1],bondffid[i][0]
                    bond[i][0],bond[i][1]=bond[i][1],bond[i][0]
                    bondff[i][0],bondff[i][1]=bondff[i][1],bondff[i][0]

                    break

        #print "fff",fff
        for i in range(len(fff)):
            flag=0
            for k in range(i+1,len(fff)):
                if fff[i][0] == fff[k][0] and fff[i][1] == fff[k][1]:
                    flag=1
                if fff[i][0] == fff[k][1] and fff[i][1] == fff[k][0]:
                    flag=1
            if flag==0:
                iddata.append(fff[i])
       # print "iddata",iddata

       # print "bondffparam",iddata
        return iddata,bondffid,bond,bondff


    def getangleffparam(self, angleffid,ff_list,angle,angleff):
        iddata=[]
        fff=[]
        for i in range(len(angleffid)):
            for k in range(len(ff_list)):
                if angleffid[i][1] == ff_list[k][1]:
                    if angleffid[i][0] == ff_list[k][0] and angleffid[i][2] == ff_list[k][2]:
                        fff.append(ff_list[k])
                        break
                    if angleffid[i][0] == ff_list[k][2] and angleffid[i][2] == ff_list[k][0]:
                            #print ffname[k][1]
                        fff.append(ff_list[k])
                        angleffid[i][0],angleffid[i][2]=angleffid[i][2],angleffid[i][0]
                        angle[i][0],angle[i][2]=angle[i][2],angle[i][0]
                        angleff[i][0],angleff[i][2]=angleff[i][2],angleff[i][0]
                        break

    #    print "anglefff", fff
        for i in range(len(fff)):
            flag=0
            for k in range(i+1,len(fff)):
                if fff[i][1] == fff[k][1]:
                    if fff[i][0] == fff[k][0] and fff[i][2] == fff[k][2]:
                        flag=1
                    if fff[i][0] == fff[k][2] and fff[i][2] == fff[k][0]:
                        flag=1
            if flag==0:
                iddata.append(fff[i])
    #    print "angleid",angleffid
    #    print "angleffid",iddata
    #    print len(angleffid),len(fff)
        return iddata,angleffid,angle,angleff


    def gettorsionffparam(self, torsionffid,ff_list,torsion,torsionff):
        iddata=[]
        fff=[]
        # print ff_list[1]
        # print torsionffid
        for i in range(len(torsionffid)):
            flag=0
            #"@-@-@-@"
            for k in ff_list[3]:
                if torsionffid[i][0] == k[0] and torsionffid[i][1] == k[1] and torsionffid[i][2] == k[2] and torsionffid[i][3] == k[3]:
                    fff.append(k)
                    flag=1
                    break
                if torsionffid[i][0] == k[3] and torsionffid[i][1] == k[2] and torsionffid[i][2] == k[1] and torsionffid[i][3] == k[0]:
                    fff.append(k)
                    flag=1
                    torsionffid[i],torsion[i],torsionff[i]=self.inversetorsion(torsionffid[i],torsion[i],torsionff[i])
                    break
                if flag == 1: break

            #"X-@-@-@"
            for k in ff_list[1]:
                if torsionffid[i][1] == k[1] and torsionffid[i][2] == k[2] and torsionffid[i][3] == k[3]:
                    fff.append(k)
                    flag=1
                    torsionffid[i][0]=1000
                    torsionff[i][0]="X"
                    # print (torsionff[i])
                    break
                if torsionffid[i][2] == k[1] and torsionffid[i][1] == k[2] and torsionffid[i][0] == k[3]:
                    fff.append(k)
                    flag=1
                    torsionffid[i][3] = 1000
                    torsionff[i][3] = "X"
                    torsionffid[i],torsion[i],torsionff[i]=self.inversetorsion(torsionffid[i],torsion[i],torsionff[i])
                    # print (torsionffid[i])
                    break
                if flag == 1: break

            #"X-X-@-@"
            for k in ff_list[2]:
                if torsionffid[i][2] == k[2] and torsionffid[i][3] == k[3]:
                    fff.append(k)
                    flag=1
                    torsionffid[i][0] = torsionffid[i][1] = 1000
                    torsionff[i][0] = torsionff[i][1] = "X"
                    # print (torsionff[i])
                    break
                if torsionffid[i][1] == k[2] and torsionffid[i][0] == k[3]:
                    fff.append(k)
                    flag=1
                    torsionffid[i][2] = torsionffid[i][3] = 1000
                    torsionff[i][2] = torsionff[i][3] = "X"
                    torsionffid[i],torsion[i],torsionff[i]=self.inversetorsion(torsionffid[i],torsion[i],torsionff[i])
                    # print (torsionff[i])
                    break
                if flag == 1: break

            #"X-@-@-X"
            for k in ff_list[0]:
                if torsionffid[i][1] == k[1] and torsionffid[i][2] == k[2]:
                    fff.append(k)
                    flag=1
                    torsionffid[i][0] = torsionffid[i][3] = 1000
                    torsionff[i][0] = torsionff[i][3] = "X"
                    break
                if torsionffid[i][2] == k[1] and torsionffid[i][1] == k[2]:
                    fff.append(k)
                    flag=1
                    torsionffid[i][0] = torsionffid[i][3] = 1000
                    torsionff[i][0] = torsionff[i][3] = "X"
                    torsionffid[i],torsion[i],torsionff[i] = self.inversetorsion(torsionffid[i],torsion[i],torsionff[i])
                    break
                if flag == 1: break

        for i in range(len(fff)):
            flag=0
            for k in range(i+1,len(fff)):
                if fff[i][1] == fff[k][1] and fff[i][2] == fff[k][2]:
                    if fff[i][0] == fff[k][0] and fff[i][3] == fff[k][3]:
                        flag=1
                if fff[i][1] == fff[k][2] and fff[i][2] == fff[k][1]:
                    if fff[i][0] == fff[k][3] and fff[i][3] == fff[k][0]:
                        flag=1
            if flag==0:
                iddata.append(fff[i])

    #    print "torsionffid",torsionffid
    #    print "iddata",iddata
    #    print "ff_list",ff_list

    #    print "torsionffparam",iddata
        return iddata,torsionffid,torsion,torsionff


    def inversetorsion(self, aaa,bbb,ccc):
        aaa[0],aaa[1],aaa[2],aaa[3] = aaa[3],aaa[2],aaa[1],aaa[0]
        bbb[0],bbb[1],bbb[2],bbb[3] = bbb[3],bbb[2],bbb[1],bbb[0]
        ccc[0],ccc[1],ccc[2],ccc[3] = ccc[3],ccc[2],ccc[1],ccc[0]
        return aaa,bbb,ccc


    def putatomparam(self, paramfile):
        atomtable ="    ["
        for i in paramfile:
            atomtable += '{"' + i[0] + '",' + str(i[1]) + "}"
        atomtable +="]\n"
        #print atomtable
        return atomtable


    def putljparam(self, paramfile):
        ljtable="    [\n"
        for i in range(len(paramfile)):
            for j in range(i,len(paramfile)):
                sigma=(paramfile[i][4] + paramfile[j][4]) / 2.0 ** (1.0/6.0)
                epsilon=math.sqrt(paramfile[i][3] * paramfile[j][3])
                if sigma == 0 or epsilon == 0:
                    continue
                            #print "sigma=", sigma
                ljtable+='        {\n'
                ljtable+='            "' + paramfile[i][0] + '-' + paramfile[j][0] + '"' + ',"Lennard_Jones","' + paramfile[i][0] + '","' +paramfile[j][0] + '",' + str(sigma*2.5) + ',0.5,{' + str(sigma) + "," + str(epsilon) + '}\n'
                ljtable+="            {0.0,0.0,0.0}\n"
                ljtable+="            {0.0,0.0,0.0,0.0,0,0}\n"
                ljtable+="            {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}\n"
                ljtable+="            {0.0,0.0,0.0,0.0,0.0,0.0,0.0}\n"
                ljtable+="            {0.0,0.0}\n"
                ljtable+="            {0.0,0.0,0.0}\n"
                ljtable+="            {0.0,0.0,0.0}\n"
                ljtable+='            {"","",{""}{0}}{0, []}\n'
                ljtable+="        }\n"
        ljtable+= "    ]\n"
        #print ljtable
        return ljtable

    #        {
    #            "c3-c3","Lennard_Jones","c3","c3",8.49917377105882,0.50000000000000,{3.39966950842353,0.10940000000000}
    #            {0.0,0.0,0.0}
    #            {0.0,0.0,0.0,0.0,0,0}
    #            {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}
    #            {0.0,0.0,0.0,0.0,0.0,0.0,0.0}
    #            {0.0,0.0}
    #            {0.0,0.0,0.0}
    #            {0.0,0.0,0.0}
    #            {"","",{""}{0}}{0, []}
    #        }


    def putbondparam(self, paramfile, atom_list):
    #    print "paramfile",paramfile
    #    print "atomlist",atom_list
        mol=[0,0]
        bondtable="    [\n"
        for i in range(len(paramfile)):
            for j in range(len(atom_list)):
                if paramfile[i][0] == atom_list[j][2]:
                    mol[0]=j
                if paramfile[i][1] == atom_list[j][2]:
                    mol[1]=j
            k_bond=paramfile[i][2]*2
            bondtable += '        {\n'
            bondtable += '            "' + atom_list[mol[0]][0]+ '-' + atom_list[mol[1]][0] +'"'+ ',"Harmonic",' + str(paramfile[i][3]) + ',{' + str(k_bond) + '}\n'
            bondtable += '            {0.0,0.0,0.0,0.0}{0.0}{0.0,0.0,""}{0, []}{0.0}{"","",{""}{0}}{0, []}\n'
            bondtable += '        }\n'
        bondtable += "    ]\n"
        #print bondtable
        return bondtable
    #    [
    #        {
    #            "c3-ca","Harmonic",1.51300000000000,{647.000000000000}
    #            {0.0,0.0,0.0,0.0}{0.0}{0.0,0.0,""}{0, []}{0.0}{"","",{""}{0}}{0, []}
    #        }
    #    ]


    def putangleparam(self, paramfile, atom_list):
        #print paramfile
        #print atom_list
        mol=[0,0,0]
        angletable= "    [\n"
        for i in range(len(paramfile)):
            for j in range(len(atom_list)):
                if paramfile[i][0] == atom_list[j][2]:
                    mol[0]=j
                if paramfile[i][1] == atom_list[j][2]:
                    mol[1]=j
                if paramfile[i][2] == atom_list[j][2]:
                    mol[2]=j
            k_bond=paramfile[i][3]*2
            e_angle=180-paramfile[i][4]
            angletable += '        {\n'
            angletable += '            "' + atom_list[mol[0]][0]+ '-' + atom_list[mol[1]][0] + '-' + atom_list[mol[2]][0] + '"' + ',"Theta",' + str(e_angle) + ',{' + str(k_bond) + '}\n'
            angletable += '            {0.0}{0.0}{0, []}{"","",{""}{0}}{0, []}\n'
            angletable += '        }\n'
        angletable += "    ]\n"
        return angletable
    #    [
    #        {"ca-ca-ha","Theta",59.9900000000000,{97.0000000000000}{0.0}{0.0}{0, []}{"","",{""}{0}}{0, []}}
    #    ]


    def puttorsionparam(self, paramfile, atom_list):
        #print paimport numpy as npramfile
        #print atom_list
        mol=[0,0,0,0]
        torsiontable= "    [\n"
        for i in range(len(paramfile)):
            for j in range(len(atom_list)):
                if paramfile[i][0] == atom_list[j][2]:
                    mol[0]=j
                if paramfile[i][1] == atom_list[j][2]:
                    mol[1]=j
                if paramfile[i][2] == atom_list[j][2]:
                    mol[2]=j
                if paramfile[i][3] == atom_list[j][2]:
                    mol[3]=j
            torsiontable += '        {\n'
            torsiontable += '            "' + atom_list[mol[0]][0]+ '-' + atom_list[mol[1]][0] + '-' + atom_list[mol[2]][0] + '-' + atom_list[mol[3]][0] + '"' + ',"Cosine_Polynomial",\n'
            torsiontable += '            {' + str(paramfile[i][4]) + ',' + str(paramfile[i][5]) + ',' + '[' + str(paramfile[i][6]) + ',' + str(paramfile[i][7]) + ',' + str(paramfile[i][8]) + ',' + str(paramfile[i][9]) + ']}\n'
            torsiontable += '            {0.0,0,0,0.0,0}{0.0,0.0,0,0}{"","",{""}{0}}{0, []}\n'
            torsiontable += '        }\n'
        torsiontable += "    ]\n"
        return torsiontable

    #    [
    #        {
    #            "X-c3-c3-X","Cosine_Polynomial",
    #            {1.00000000000000, 4, [0.15555555555556,0.46666666666667,0.0,-0.62222222222222]}
    #            {0.0,0,0,0.0,0}{0.0,0.0,0,0}{"","",{""}{0}}{0, []}
    #        }
    #    ]


    def putinteractions(self, ljparam):

        ljtable=self.putljparam(ljparam)
        interactions='Interactions:{\n' + str(ljtable) +'''
        []
        [
            {
                "POINT_CHARGE","Ewald",0.83333333333333,
                {0.0,0.0}{0.0,0.0,0.0}{0.0,0.0}{0.0,11.3038833052088,"Auto",{0.0,{0,0,0}}}
                {0.0,0.0,0.0,0.0,0.0,0}{0.0,0.0,"",{0.0,{0,0,0}}0.0,0}
            }
        ]
    }
    '''
        #print interactions
        return interactions

    def putinteractionsitetype(self, paramfile):
        #print paramfile
        isitetable="    ["
        for i in paramfile:
            s_range=i[4] * 2 / 2**(1.0/6.0) * 1.3
            isitetable+= '{"' + i[0] + '",1,' + str(s_range) + "}"
        isitetable+= "]\n"
        return isitetable
    #(Interaction_Site_Type.Range=atomsigma*2*1.3)
    #    [{"c3",1,4.41957044601440}{"ha",1,3.37953519821167}{"hc",1,3.44439268112183}{"ca",1,4.41957044601440}]


    def putatom(self, atom, ffname, molname, tnum):
        atomtable=""
        atomtable+="            [\n"
        for i in range(len(atom)):
            atomtable+='                {'+str(tnum) + ',"' +atom[i][1] + '","' +str(ffname[i][1])+'",0,1,[{"1","' + molname + ':0"}]}\n'
            tnum=tnum + 1
        atomtable+="            ]\n"
        return tnum,atomtable
    #            [
    #                {0,  "C",  "c3",  0,  1,  [{"1","Polyethylene:0"}]}
    #                {1,  "H",  "hc",  0,  1,  [{"1","Polyethylene:0"}]}
    #                {2,  "H",  "hc",  0,  1,  [{"1","Polyethylene:0"}]}
    #                {3,  "C",  "c3",  0,  1,  [{"1","Polyethylene:0"}]}
    #            ]

    def putbond(self, bondff, bond):
    #    print bondff
    #    print bond
        bondtable="            [\n"
        for i in range(len(bondff)):
            bondtable+='                {"' + bondff[i][0] + "-" + bondff[i][1] + '",' + str(bond[i][0]) + ',' + str(bond[i][1]) + ',' + "1.0}\n"
        bondtable+="            ]\n"
        return bondtable

    #            [
    #                {"c3-hc",0,1,1.00000000000000}
    #                {"c3-hc",0,2,1.00000000000000}
    #                {"c3-hc",3,4,1.00000000000000}
    #            ]


    def putangle(self, angleff, angle):
        angletable="            [\n"
        for i in range(len(angleff)):
            angletable+='                {"' + angleff[i][0] + "-" + angleff[i][1] + "-" + angleff[i][2]+ '",' + str(angle[i][0]) + ',' + str(angle[i][1]) + ',' + str(angle[i][2]) + '}\n'
        angletable+="            ]\n"
        return angletable


    #            [
    #                {"hc-c3-hc",1,0,2}
    #                {"c3-c3-hc",1,0,3}
    #                {"hc-c3-hc",1,0,6}
    #            ]

    def puttorsion(self, torsionff, torsion):
        torsiontable="            [\n"
        for i in range(len(torsionff)):
            torsiontable+='                {"' + torsionff[i][0] + "-" + torsionff[i][1] + "-" + torsionff[i][2] + "-" + torsionff[i][3] + '",' + str(torsion[i][0]) + ',' + str(torsion[i][1]) + ',' + str(torsion[i][2]) + ',' + str(torsion[i][3]) + '}\n'
        torsiontable+="            ]\n"
        return torsiontable

    #            [
    #                {"hc-c3-c3-hc",1,0,3,4}
    #                {"hc-c3-c3-hc",1,0,3,5}
    #                {"hc-c3-c3-hc",1,0,3,7}
    #            ]

    def putinteractionsite(self, ffname):
        #print paramfile
        isitetable= "            [\n"
        for i in range(len(ffname)):
            isitetable+= '                {"' + ffname[i][1] + '",[' + str(i) + "]}\n"
        isitetable+= "            ]\n"
        return isitetable

    #            [
    #                {"c3", [0]}
    #                {"hc", [1]}
    #                {"hc", [2]}
    #            ]

    def putpointcharge(self, chg):
        #f =open(chgfile, "r")
        #text = f.readlines()
        chgtable= "            [\n"
        for i in range(len(chg)):
        #    itemList=text[i].split()
            chgtable+= '                {"POINT_CHARGE",[' + str(i) + '] ' + str(chg[i]) + "}\n"
        chgtable+= "            ]\n"
        return chgtable

    #            [
    #                {"POINT_CHARGE",  [0] -0.44582881929900}
    #                {"POINT_CHARGE",  [1]  0.14860960643300}
    #                {"POINT_CHARGE",  [2]  0.14860960643300}
    #            ]

    def putmolecularattributes(self, ljparam, bondparam, angleparam, torsionparam, atom_list):
        atomtable=self.putatomparam(ljparam)
        bondtable=self.putbondparam(bondparam,atom_list)
        angletable=self.putangleparam(angleparam,atom_list)
        torsiontable=self.puttorsionparam(torsionparam,atom_list)
        isitetable=self.putinteractionsitetype(ljparam)

        molattr='Molecular_Attributes:{\n' + str(atomtable) + str(bondtable) + str(angletable) +str(torsiontable) + str(isitetable) + '    []\n' + '}\n'
        #print molattr
        return molattr

    def putsetofmolecules(self, poslist, param):
        fname=param[0]
        atom=param[1]
        ffname=param[2]
        molname=param[3]
        bondff=param[4]
        bond=param[5]
        angleff=param[6]
        angle=param[7]
        torsionff=param[8]
        torsion=param[9]
        chg=param[10]

        somolecules=""
        tnum=0
        somolecules+='Set_of_Molecules:{\n'
        somolecules+='    [\n'
        for i in range(len(poslist[0])):
            somolecules+='        {\n'
            somolecules+='            "'+fname[0]+'"\n'

            tnum,atomtable=self.putatom(atom[0],ffname[0],molname[0],tnum)
            bondtable=self.putbond(bondff[0],bond[0])
            angletable=self.putangle(angleff[0],angle[0])
            torsiontable=self.puttorsion(torsionff[0],torsion[0])
            isitetable=self.putinteractionsite(ffname[0])
            chgtable=self.putpointcharge(chg[0])

            somolecules+=str(atomtable)+str(bondtable)+str(angletable)+str(torsiontable)+str(isitetable)+str(chgtable)
            somolecules+='        }\n'
        for i in range(len(poslist[1])):
            somolecules+='        {\n'
            somolecules+='            "'+fname[1]+'"\n'
            tnum,atomtable=self.putatom(atom[1],ffname[1],molname[1],tnum)
            bondtable=self.putbond(bondff[1],bond[1])
            angletable=self.putangle(angleff[1],angle[1])
            torsiontable=self.puttorsion(torsionff[1],torsion[1])
            isitetable=self.putinteractionsite(ffname[1])
            chgtable=self.putpointcharge(chg[1])
            somolecules+=str(atomtable)+str(bondtable)+str(angletable)+str(torsiontable)+str(isitetable)+str(chgtable)
            somolecules+='        }\n'

        somolecules+='    ]\n'
        somolecules+='}\n'
        #print somolecules
        return somolecules

    #Set_of_Molecules:{
    #    [
    #        {
    #            "molA",

    def putpos(self, pos):
        #print pos
        #print len(pos)
        #print len(pos[0])
        postable=""
        for j in range(len(pos)):
            postable+='            {\n'
            postable+='                [\n'
            for i in range(len(pos[j])):
                postable+='                    {' + str(pos[j][i][0]) + ',' + str(pos[j][i][1]) + ',' + str(pos[j][i][2]) + '}\n'
            postable+='                ]\n'
            postable+='            }\n'
        return postable

    def putstructure(self, poslist,cell):
        postable=[]
        structure='Structure:{\n'
        structure+='    {\n'
        structure+='        [\n'
        for i in range(2):
            postable.append(self.putpos(poslist[i]))
        structure+=postable[0] + postable[1]
        structure+='        ]\n'
        structure+='    }\n'
        structure+='    {[]}{[]}\n'
        structure+='    {0.0,{' + str(cell[0]) + ',' +  str(cell[1]) + ',' + str(cell[2]) + ',90.0000000000000,90.0000000000000,90.0000000000000}0.0}\n'
        structure+='}\n'
        structure+='Unit_Parameter:{"","",1.00000000000000,4.18550000000000,1.00000000000000e-001}\n'
        structure+='\end{data}\n'
        #print structure
        return structure

    def putheader(self):
        header=\
        '''COGNAC INPUT UDF DATA.

\include{"cognac90.udf"}

\\begin{header}
\\begin{data}
EngineType:"COGNAC"
EngineVersion:"V90"
IOType:"INOUT"
ProjectName:"TEST"
Comment:"FF_TYPE[GAFF]"
Action:"cognac_draw.act;cognac_info.act;cognac_plot.act;cognac_anal.act;cognac_edit.act;cognac_lammps.act"
\end{data}
\end{header}

\\begin{data}
'''

        return header

    def putsimulationcondition(self, totalstep, outstep, totalmass, algo):
        simucondition=\
            '''Simulation_Conditions:{
    {
        1000000.00000000,
        {2.04584957182253e-002,''' + str(totalstep) + ',' + str(outstep) + '''}
        {0.59595103713955,0}
        {1.43880938955919e-005,{0.0,0.0,0.0,0.0,0.0,0.0}}
        {
            "",
            {"",{0.0}{0.0,0.0}}
            {"",{0.0,0.0,0.0,0.0,0.0,0.0}{0.0,0.50000000000000,""}{0.0,0.0,0.0}100,0}
        }
        {100,0,0,0,0}
        {0,0,1.00000000000000e-005}
    }

    {
        "Dynamics",
        {
            "''' + algo[0] + '''",
            {20.0000000000000}
            {10.0000000000000}
            {0.50000000000000}
            {''' + str(totalmass) + '''}
            {''' + str(totalmass) + ''',0,""}
            {80.0000000000000,0,""}
            {''' + str(totalmass) + ''',20.0000000000000}
            {''' + str(totalmass) + ''',0.50000000000000}
            {''' + str(totalmass) + ''',0,"",25.0000000000000}
            {''' + str(totalmass) + ''',0,"",0.50000000000000}
            {20.0000000000000,1.00000000000000}
            {20.0000000000000,0,"",1.00000000000000}
            {0.0}{0.0,0.0}{0.0}
        }
        {"",0,0.0,0,{0,0.0}[]}
    }

    {"PERIODIC","PERIODIC","PERIODIC",0}
    {1,1,1,1,0,1,1,0,1,1}
    {{1,1,1,1,1,1,1,1,1}{1,1,1}{0,0,0}{0}}
    {0,0,{0.0,0.0,0.0}0.0,0}
    {"NO",[]}

}
'''
        return simucondition

    def putinitialstructure(self, cell):
        initialstructure=\
        '''Initial_Structure:{
    {0.0,{'''+ str(cell[0]) + ',' + str(cell[1]) + ',' + str(cell[2]) + ''',90.0000000000000,90.0000000000000,90.0000000000000}0.0}
    {"",-1}
    {
        "Restart",
        {"",-1,1,0}
        {0,0,0.0,0.0,[][]0.0,[]}
        {{"",{0.0,0.0}}0,0,0,[]0.0,0.0,[]0.0,""}
        {0.0,0,0.0,{0.0,0,0.0,0.0,0.0}{0.0, [] 0.0}}
        {"",0,0,0,0.0}
    }
    {1,"DYNAMICS",300.000000000000,10000}
}
'''
        return initialstructure


    def getbondparampair(self, data):
    #    print "data",data
        fff=[]
        iddata=[]
        for i in range(len(data)):
            for j in range(len(data[i])):
                fff.append(data[i][j])

    #    print fff
        for i in range(len(fff)):
            flag=0
            for k in range(i+1,len(fff)):
                if fff[i][0] == fff[k][0] and fff[i][1] == fff[k][1]:
                    flag=1
                if fff[i][0] == fff[k][1] and fff[i][1] == fff[k][0]:
                    flag=1
            if flag==0:
                iddata.append(fff[i])
    #    print "bondparampair",iddata
        return iddata

    def getangleparampair(self, data):
    #    print "data",data
        fff=[]
        iddata=[]
        for i in range(len(data)):
            for j in range(len(data[i])):
                fff.append(data[i][j])

        for i in range(len(fff)):
            flag=0
            for k in range(i+1,len(fff)):
                if fff[i][1] == fff[k][1]:
                    if fff[i][0] == fff[k][0] and fff[i][2] == fff[k][2]:
                        flag=1
                    if fff[i][0] == fff[k][2] and fff[i][2] == fff[k][0]:
                        flag=1
            if flag==0:
                iddata.append(fff[i])

    #    print "angleparampair",iddata
        return iddata

    def gettorsionparampair(self, data):
    #    print "data",data
        fff=[]
        iddata=[]
        for i in range(len(data)):
            for j in range(len(data[i])):
                fff.append(data[i][j])

        for i in range(len(fff)):
            flag=0
            for k in range(i+1,len(fff)):
                if fff[i][1] == fff[k][1] and fff[i][2] == fff[k][2]:
                    if fff[i][0] == fff[k][0] and fff[i][3] == fff[k][3]:
                        flag=1
                if fff[i][1] == fff[k][2] and fff[i][2] == fff[k][1]:
                    if fff[i][0] == fff[k][3] and fff[i][3] == fff[k][0]:
                        flag=1
            if flag==0:
                iddata.append(fff[i])

    #    print "torsionparampair",iddata
        return iddata


    def gen_udf(self, udf_param, out_name, som_param):
        cellsize = udf_param[0]
        ljparam = udf_param[1]
        bondparam = udf_param[2]
        angleparam = udf_param[3]
        torsionparam = udf_param[4]
        atom_list = udf_param[5]
        totalstep = udf_param[6]
        outstep = udf_param[7]
        totalmass = udf_param[8]
        algo = udf_param[9]
        poslist = udf_param[10]

        header = self.putheader()
        simucondition = self.putsimulationcondition(totalstep, outstep, totalmass, algo)
        initialstructure = self.putinitialstructure(cellsize)
        molattr = self.putmolecularattributes(ljparam, bondparam,
                                              angleparam, torsionparam, atom_list)
        interactions = self.putinteractions(ljparam)
        somolecules = self.putsetofmolecules(poslist,som_param)
        structure = self.putstructure(poslist,cellsize)

        udf_body = str(header) + str(simucondition) + str(initialstructure) + \
            str(molattr) + str(interactions) + str(somolecules) + str(structure)
        out_file = open(out_name, "w")
        print(udf_body, file=out_file)
        out_file.close()


    def putclusterpos(self, uobj, poslist, molnum):
        #molnum= uobj.size('Set_of_Molecules.molecule[]')
        count = 0
        #print poslist
        for i in range(molnum):
            atom = uobj.size('Set_of_Molecules.molecule[].atom[]',[i])
            for j in range(atom):
                #print count
                uobj.put(float(poslist[count][0]) + 3.0, 'Structure.Position.mol[' + str(i) + '].atom[' + str(j) + '].x')
                uobj.put(float(poslist[count][1]) + 3.0, 'Structure.Position.mol[' + str(i) + '].atom[' + str(j) + '].y')
                uobj.put(float(poslist[count][2]) + 3.0, 'Structure.Position.mol[' + str(i) + '].atom[' + str(j) + '].z')
                count +=1

    def clusterfix(self, uobj,atomlen,nummol_seg,fix_label):

        ca_id = uobj.size("Simulation_Conditions.Constraint_Conditions.Constraint_Atom[]")
        for i in range (nummol_seg):
            #print fix_label[i],atomlen
            target=fix_label[i]- atomlen*(i)
            #print "target=",target
            #uobj.put('NO',"Simulation_Conditions.Constraint_Conditions.Read_from_Restart")
            uobj.put(i,"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Index.Mol_Index")
            uobj.put(target,"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Index.Atom_Index")
            uobj.put('YES',"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Constraint_Axis.x")
            uobj.put('YES',"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Constraint_Axis.y")
            uobj.put('YES',"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Constraint_Axis.z")
            uobj.put('Steady',"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Method")
            uobj.put(0,"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Steady.Velocity.x")
            uobj.put(0,"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Steady.Velocity.y")
            uobj.put(0,"Simulation_Conditions.Constraint_Conditions.Constraint_Atom[" + str(ca_id) + "].Steady.Velocity.z")
            ca_id += 1
        print ("cluster fixed.")


