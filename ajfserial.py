import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import abmptools as ampt
import argparse

# Koji Okuwaki: Update 2020/03/22
# moveintocell and assignmolname

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                prog='pdbmodify.py', # program name
                usage='''e.g.)
                python ajfserial.py -i template.ajf -t 1 100 1 -str xxx
                ''', # program usage
                description='description',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--input',
                        help='input template ajf',
                        nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-t', '--time',
                        help='time info',
                        nargs=3,
                        )

    parser.add_argument('-str', '--string',
                        help='rename data',
                        nargs='*',
                        )

    args = parser.parse_args()
    print('input =', args.input)
    print('mode =', args.mode)
    print('move =', args.move)
    print('move pos =', args.pos)
    print('move mol =', args.mol)
    print('add chain =', args.addchain)
    print('intoflag =', args.into)
    print('refreshresid =', args.refreshresid)
    print('assignresname =', args.assignresname)
    print('refreshresid =', args.refreshresid)
    print('refreshatmtype =', args.refreshatmtype)

    ## -- user setting --
    # read info
    mode = args.mode #rfile, resnum
    assignresname = args.assignresname
    refreshatmtype = args.refreshatmtype
    refreshresid = args.refreshresid

    # move info
    moveflag = args.move
    if args.mol == None:
        movemode = 'pos' # pos or mol
    else:
        movemode = 'mol' # pos or mol
        # --- mol mode
        tgtmol = args.mol


    addchain = False
    if args.addchain != None:
        addchain = True
        addres_start = int(args.addchain[0])
        addres_end = int(args.addchain[1])
        chainlab = args.addchain[2]

        addres = [i for i in range(addres_start, addres_end+1)]

    # --- pos mode
    tgtpos = args.pos
    intoflag = args.into


    for i in range(len(args.input)):
        infile = args.input[i]
        head, ext = os.path.splitext(infile)
        print(head, ext)

        if ext != '.ajf':
            out = head.split('.pdb')[0] + ext + '-renamed.pdb'
        else:
            out = head + '-renamed.pdb'
        print(out)

        lines = open(infile, 'r').readlines()
        # print(lines)


        replacelists = []
        for j in range(0, len(args.string), 2):
            rlist = [args.string[j], args.string[j+1]]
            replacelists.append(rlist)
        print(replacelists)
        outf = open(out, 'w')
        # for repinfo in replacelists:
            # print('{0:>3}'.format(repinfo[0]) + '     ' + repinfo[1] + ' ')
        for line in lines:
            for repinfo in replacelists:
                if len(repinfo) == 3:
                    line = line.replace(' {0:>3}'.format(repinfo[0]) + '     ' + repinfo[1] + ' ', '{0:>3}'.format(repinfo[2]) + '     ' + repinfo[1] + ' ')
                if len(repinfo) == 2:
                    line = line.replace(' {0:>3}'.format(repinfo[0]) + ' ', ' {0:>3}'.format(repinfo[1]) + ' ')

            print(line[:-1], file=outf)

        print(out, 'was created.')

