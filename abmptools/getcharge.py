import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import abmptools
import collections
import argparse
# Matrix operation

if __name__ == "__main__":
    ## -- user setting --
    parser = argparse.ArgumentParser(
                prog='getcharge.py', # program name
                usage='''e.g.)
                python -m abmptools.getcharge -i 6lu7-covneu-nowat-hinagata0514.log -f 125 126 128 -t nbo
                ''', # program usage
                description='description',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--incoord',
                        help='input frag info',
                        # nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-t', '--type',
                        help='template file',
                        required=True)

    parser.add_argument('-f', '--frag',
                        help='output config file',
                        nargs='*',
                        required=True
                        # default=['HOH', 'WAT', 'NA']
                        )
    ###################################

    # get args
    args = parser.parse_args()

    print('incoord(log) =', args.incoord)
    print('chargetype = ', args.type)
    print('frag = ', args.frag)

    ####################################

    args = parser.parse_args()
    # print('output =', args.output)

    logname = args.incoord
    chgtype = args.type
    frags = args.frag

    ## -- setting end --
    aobj = abmptools.anlfmo()

    natom = aobj.readlog(logname, 'auto')
    chgdf = aobj.getlogchgall(logname, natom, chgtype)
    print(natom)
    print(chgdf)

    # for frag in frags:
        # print(chgdf[chgdf['Frag'] == int(frag)])

    head = os.path.splitext(logname)[0]
    chgfiltdf = chgdf[chgdf['Frag'].isin(frags)]
    chgfiltdf.to_csv(head + '-' + chgtype + '.csv')

