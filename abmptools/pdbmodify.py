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
                python pdbmodify.py -addc 307 312 C -i *.pdb
                python pdbmodify.py -move -p 10 10 10 -aresname -reid -reatm -addc 307 312 C -i *pdb
                python pdbmodify.py -move -mol 3 --into -reres -reatm -addc 307 312 C -i *pdb
                python pdbmodify.py -mode rename -str 001 MRT 002 CD7 003 NA 004 WAT -i *pdb
                ''', # program usage
                description='description',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--input',
                        help='input pdb info',
                        nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-move', '--move',
                        help='move',
                        action='store_true'
                        )

    parser.add_argument('-p', '--pos',
                        help='move position',
                        nargs=3,
                        )

    parser.add_argument('-addc', '--addchain',
                        help='add chain',
                        nargs=3,
                        )

    parser.add_argument('-mol', '--mol',
                        type=int,
                        help='move mol',
                        )

    parser.add_argument('-into', '--into',
                        help='move',
                        action='store_true'
                        )

    parser.add_argument('-aresname', '--assignresname',
                        help='assignresname',
                        action='store_true'
                        )

    parser.add_argument('-reid', '--refreshresid',
                        help='refreshresid',
                        action='store_true'
                        )

    parser.add_argument('-reatm', '--refreshatmtype',
                        help='refreshatmtype',
                        action='store_true'
                        )

    parser.add_argument('-mode', '--mode',
                        help='[rfile], [resnum], [rename]',
                        default='resnum',)

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
    ## -- setting end --

    # main
    # argvs = sys.argv
    # fname = str(argvs[1])

#     for arg in argvs:
#         if arg == '--move':
#             moveflag == True
#         if arg == '--nomove':
#             moveflag == False

    if args.mode == 'rename':
        for i in range(len(args.input)):
            infile = args.input[i]
            head, ext = os.path.splitext(infile)
            print(head, ext)

            if ext != '.pdb':
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

    else:
        for i in range(len(args.input)):
            fname = args.input[i]
            oname, ext = os.path.splitext(fname)
            if ext != '.pdb':
                oname = oname.split('.pdb')[0] + ext.split('.')[1] + '-mod'
            else:
                oname = oname + '-mod'

            aobj = ampt.setfmo()
            aobj.getmode = mode
            aobj.assignresname = assignresname
            aobj.refreshatmtype = refreshatmtype
            aobj.refreshresid = refreshresid

            print('infile:', fname)
            print('oname:', oname)
            #print('centered-molid:', tgtmol - 1)

            # get pdbinfo
            # totalMol, atomnameMol, molnames, nlabmols, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfo(fname)

            aobj = aobj.getpdbinfowrap(fname)
            totalMol = aobj.totalRes
            atomnameMol = aobj.atmtypeRes
            molnames = aobj.resnames
            nlabmols = aobj.gatmlabRes
            posMol = aobj.posRes
            heads = aobj.headRes
            labs = aobj.labRes
            chains = aobj.chainRes
            resnums = aobj.resnumRes
            codes = aobj.codeRes
            occs = aobj.occRes
            temps = aobj.tempRes
            amarks = aobj.amarkRes
            charges = aobj.chargeRes

            mollist = [i for i in range(totalMol)]
    #         cellsize = aobj.getpdbcell(fname)
    #         aobj.cellsize = cellsize

            print('totalMol:', totalMol)

            # print(resnums)
            if addchain == True:
                print(addres)
                alreadys = []
                for i in range(len(resnums)):
                    # print(resnums[i][0])
                    tgt = int(resnums[i][0])
                    if tgt in alreadys:
                        continue
                    if tgt in addres:
                        print(resnums[i])
                        # print(nlabmols[i])
                        for j in range(len(chains[i])):
                            aobj.chainRes[i][j] = chainlab
                        alreadys.append(tgt)

    #         if len(aobj.cellsize) == 0:
    #             aobj.cellsize = 0
    #             print('cellinfo: None')
    #         else:
    #             print('cellsize:', aobj.cellsize)


            if moveflag == True:
                # get center of solute
                if movemode == 'mol':
                    coctgt = aobj.getCenter(posMol[tgtmol-1])
                elif movemode == 'pos':
                    coctgt = tgtpos
                transVec = np.array(-coctgt, dtype='float64')
                # print(transVec)

                # move
                posmoveMol = []
                for i in range(totalMol):
                    posmove = aobj.movemoltranspdb(posMol[i], transVec)
                    posmoveMol.append(posmove)

                if intoflag == True:
                    posintoMol = aobj.moveintocellpdb(posmoveMol, totalMol, cellsize)

                else:
                    posintoMol = copy.deepcopy(posmoveMol)

            else:
                if intoflag == True:
                    posintoMol = aobj.moveintocellpdb(posMol, totalMol, cellsize)

                else:
                    posintoMol = copy.deepcopy(posMol)

            # write
            aobj.posRes = posintoMol
            aobj.exportardpdbfull(oname, mollist)



