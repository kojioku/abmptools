import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import abmptools as abmp
import argparse


if __name__ == "__main__":
    # main
    # create parser
    parser = argparse.ArgumentParser(
                prog='pdb2fmo.py', # program name
                usage='python pdb2fmo.py -i xxx.pdb', # program usage 
                description='generate ABNITMP input (ajf,pdb set) from orig pdb and segment_data file',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--incoord',
                        help='coordinate file (pdb)',
                        nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-p', '--parameter',
                        help='parameter file',
                        default='input_param')

    parser.add_argument('-noreid', '--norefreshresid',
                        help='refreshresid',
                        action='store_false'
                        )

    parser.add_argument('-noreatm', '--norefreshatmtype',
                        help='refreshatmtype',
                        action='store_false'
                        )

    # get args
    args = parser.parse_args()

    print('coord(pdb) =', args.incoord)
    print('parameter = ', args.parameter)

    for i in range(len(args.incoord)):
        fname = args.incoord[i]
        oname, ext = os.path.splitext(fname)

        aobj = abmp.setfmo()
        path = 'for_abmp'

        param_read = {}
        exec(open(args.parameter, 'r').read(), param_read)
        param_rfmo = param_read['param']
        aobj.setrfmoparam(param_rfmo)
        if args.norefreshresid:
            aobj.refreshresid = True
        if args.norefreshatmtype:
            aobj.refreshatmtype = True

        print('--- info ---')
        print('setup mode', aobj.cutmode)
        print('piedaflag', aobj.piedaflag)
        print('cmmflag', aobj.cmmflag)
        print('refreshresid', aobj.refreshresid)
        print('refreshatmtype', aobj.refreshatmtype)


        if aobj.cutmode == 'sphere':
            oname = oname + '-' + aobj.cutmode + 'p' + str(aobj.tgtpos[0]) + '_' + str(aobj.tgtpos[1]) + '_' + str(aobj.tgtpos[2]) + '_ar' + str(aobj.criteria)

        elif aobj.cutmode == 'around':
            oname = oname + '-' + aobj.cutmode + '_ar' + str(aobj.criteria)

        elif aobj.cutmode == 'cube':
            oname = oname + '-' + aobj.cutmode + 'p' + str(aobj.tgtpos[0]) + '_' + str(aobj.tgtpos[1]) + '_' + str(aobj.tgtpos[2]) + '_x' + str(aobj.criteria[0]) + '_y' + str(aobj.criteria[1]) + '_z' + str(aobj.criteria[2])

        elif aobj.cutmode == 'neutral':
            oname = oname + '-' + aobj.cutmode + '_ar' + str(aobj.criteria)

        if aobj.cutmode == 'none':
            oname = oname +  '-for_abmp'

        if aobj.cutmode == 'around':
            print('molset', aobj.molname)
            print('solute', aobj.solutes)

        if len(aobj.ionname) != 0:
            print('ion mode: ', aobj.ionmode)
            print('ion name: ', aobj.ionname)

        aobj.getcontact_rmapfmopdb(path, fname, oname)
