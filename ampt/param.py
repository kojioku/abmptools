# -*- coding: utf-8 -*-
import sys
import os
scrdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scrdir)

import copy
import subprocess

class setparam:

    def __init__(self):
        # -- assign param to variable --
        self.solv_flag = False

        # submit param
        self.submit_system = 'none'
        self.npro = 4
        self.nnode = 0
        self.para_job = 2
        self.memory = "1800"
        self.queue = None

        # abinitmp param
        self.ajf_method = "MP2"
        self.ajf_basis_set = "6-31G*"
        self.abinit_ver = 'rev10'
        self.abinitmp_path = 'abinitmp'
        self.mpipath = 'mpirun'
        self.pbmolrad = 'vdw'
        self.pbcnv = 1.0
        self.piedaflag = False
        # distlitname

        # others
        self.python_path = 'python'
        self.pynp = 4
        self.babelflag = False
        self.targets = []
        self.mainpath = '.'
        self.nofzc_flag = False
        self.fzcnum = 0
        self.ksubcmd = 'pjsub'

        # param for result

        return self

    def setupparam(self, param_global, argc, argvs):

        # assign param to variable
        try:
            self.ajf_method = param_global['method']
        except KeyError:
            self.ajf_method = 'MP2'

        try:
            self.ajf_basis_set = param_global['basis']
        except KeyError:
            self.ajf_basis_set = '6-31G*'

        try:
            self.solv_flag = param_global['solv_flag']
        except KeyError:
            self.solv_flag = 'False'

        try:
            self.piedaflag = param_global['piedaflag']
        except KeyError:
            self.piedaflag = False

        try:
            self.pbmolrad = param_global['pbmolrad']
        except KeyError:
            self.pbmolrad = 'vdw'

        try:
            self.pbcnv = param_global['pbcnv']
        except KeyError:
            self.pbcnv = 1.0
            # self.pbcnv = 1.0E-2

        try:
            self.ksubcmd = param_global['ksubcmd']
        except KeyError:
            self.ksubcmd = 'pjsub'

        try:
            self.fzclists = param_global['nofzclists']
        except KeyError:
            self.nofzc_flag = False
            self.fzclists = []

        if self.fzclists:
            self.nofzc_flag = True

        try:
            self.abinitmp_path = param_global['abinitmp_path']
        except KeyError:
            self.abinitmp_path = 'abinitmp'

        try:
            self.submit_system = param_global['submit_system']
        except KeyError:
            self.submit_system = 'none'

        try:
            self.npro = param_global['num_processor']
        except KeyError:
            self.npro = 8

        try:
            self.nnode = param_global['num_node']
        except KeyError:
            self.nnode = 1

        try:
            self.para_job = param_global['para_job']
        except KeyError:
            self.para_job = 1

        try:
            self.queue = param_global['queue_name']
        except KeyError:
            self.queue = None

        try:
            self.dlisttype = param_global['dist_list_type']
        except KeyError:
            self.dlisttype = 'vdW'

        try:
            self.abinit_ver = param_global['abinitmp_ver']
        except KeyError:
            self.abinit_ver = 'rev10'

        try:
            self.mpipath = param_global['mpipath']
        except KeyError:
            self.mpipath = 'mpirun'


        try:
            self.gname = param_global['group_name']
        except KeyError:
            self.gname = ''


        try:
            self.memory = param_global['memory']
        except:
            print("Memory is set to 1800MB")
            self.memory = "1800"


        # check abinitmp ver
        abinit_vers = ['rev5', 'rev10', 'rev11', 'rev15', 'rev17']
        if self.abinit_ver not in abinit_vers:
            print("abinitmp ver error!")
            sys.exit()

        return self


    def getbinaryflag(self, cmd):
        flag = True
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.kill()
        except:
            flag = False

        return flag

    def config_read(self, seg_name, seg_atom):
            fdata = {}
            exec(open(self.mainpath + "/segment_data.dat", "r").read(), fdata)
            seg_data = fdata['seg_data']
            # print (seg_data)
            # default data
            seg_conf = {'name': seg_name,
                        'atom': [seg_atom],
                        'charge': [0],
                        'connect_num': [0],
                        'connect': [],
                        'seg_info': [[i + 1 for i in range(seg_atom)]],
                        'nummol_seg': [1],
                        'repeat': [1],
                        'pair_file': [],
                        'multi_xyz': 'none'}

            # print seg_data[0]
            for i in range(len(seg_data)):
                if seg_name == seg_data[i]['name']:
                    for key, data in seg_data[i].items():
                        seg_conf[key] = data

            # print "seg_conf =", seg_conf
            return seg_conf

