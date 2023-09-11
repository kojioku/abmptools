import sys
import os
import pandas as pd
import itertools
import copy
import csv
import abmptools as ampt
import argparse
# Author: Koji Okuwaki


def get_args():
    parser = argparse.ArgumentParser(
                prog='getifiepieda',
                usage='''e.g)
    # python -m abmptools.getifiepieda --frag 1-10 101-200 -i xxx.log
    # python -m abmptools.getifiepieda --frag 10 -d 0.8 -i xxx.log
    # python -m abmptools.getifiepieda --frag 10 --molname WAT -i xxx.log
    # python -m abmptools.getifiepieda --frag 10 --molname WAT -d 0.8 -i xxx.log
    # python -m abmptools.getifiepieda --mol 1-10 -i xxx.log
    # python -m abmptools.getifiepieda --fraginmol 1-10 2 000 1 -i xxx.log
    # python -m abmptools.getifiepieda --ffmatrix 1-100 101-200 -i xxx.log
    # python -m abmptools.getifiepieda --multi 1-100 101-200 -t 100 3100 1000 -i '["6lu7orig_md040j8_163neu", -hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4
    # python -m abmptools.getifiepieda --multi 1-100 101-200 -t 100 3100 1000 -ix 6lu7orig_md040j8_163neuxxx-hopt-ps-mod_forabmp_192n-2p-24t.log --exclude 102 -np 4

    # python -m abmptools.getifiepieda --multi 10 -d 8.0 -t 100 3100 1000 -i '["6lu7orig_md040j8_163neu", "-hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4
    # python -m abmptools.getifiepieda --multi 20 --molname WAT -d 8.0 -t 100 3100 1000 -i '["6lu7orig_md040j8_163neu", "-hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4
    # python -m abmptools.getifiepieda --tfmatrix 1-100 101-200 -t 100 3100 1000 -i '["6lu7orig_md040j8_163neu", "-hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4''',
                description='Analysis script for ABINIT-MP log',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-f', '--frag',
                        help='tgt fragid info',
                        nargs='*'
                        )

    parser.add_argument('-m', '--mol',
                        help='tgt molid info')

    parser.add_argument('-mn', '--molname',
                        help='tgtmolname')

    parser.add_argument('-fi', '--fraginmol',
                        help='tgtmolid, tgt1_localfrag, tgt2molname, tgt2_localfrag',
                        nargs=4,)

    parser.add_argument('-ff', '--ffmatrix',
                        help='generate ffmatrix',
                        nargs=2
                        )

    parser.add_argument('-mul', '--multi',
                        help='use solvation',
                        nargs='*'
                        )

    parser.add_argument('-z', '-zp', '--zp'
                        help='zeropadding',
                        default = 0,
                        type=int
                        )

    parser.add_argument('-tf', '--tfmatrix',
                        help='generate time-frag matrix',
                        nargs=2
                        )

    parser.add_argument('-t', '--time',
                        help='start end interval',
                        nargs=3
                        )

    parser.add_argument('-d', '--dist',
                        type=float,
                        help='screen dist',
                        )

    parser.add_argument('-np', '--pynp',
                        type=int,
                        help='python np',
                        default=1)

    parser.add_argument('-ex', '--exclude',
                        help='assign exclude',
                        nargs='*',
                        type=int,
                        default=[],)

    parser.add_argument('-i', '--input',
                        default=None,
                        help='input file',
                        )

    # Input CPF name, required
    parser.add_argument("-ix", "--inputx",
                        default=None,
                        help="Input log name (e.g., file-xxx-bbb.log)")

    parser.add_argument('-dimeres', '--dimeres',
                        help='get dimer-es info',
                        action='store_true')

    parser.add_argument('-nof90', '--nof90so',
                        help='use f90',
                        action='store_false',
                        default=True)

    parser.add_argument('-nores', '--noresinfo',
                        help='assign resinfo',
                        action='store_false',
                        default=True)

    parser.add_argument('-dimene', '--dimene',
                        help='dimer energy label',
                        nargs=2,
                        type=int,
                        default=[2, 1],)

    parser.add_argument('-momene', '--momene',
                        help='monomer energy label',
                        nargs=1,
                        type=int,
                        default=[1])

    parser.add_argument('-imd', '--is_momdim',
                        help='get monomer dimer energy ',
                        action='store_true',
                        default=False)

    args = parser.parse_args()
    print('## setup info')
    print('tgtfrag =', args.frag)
    print('tgtmolid =', args.mol)
    print('tgt2molname =', args.molname)
    print('fraginmol info =', args.fraginmol)
    print('ffmatrixinfo =', args.ffmatrix)
    print('multiinfo =', args.multi)
    print('tfmatrixinfo =', args.tfmatrix)
    print('time =', args.time)
    print('dist =', args.dist)
    print('process =', args.pynp)
    print('excludefrag =', args.exclude)

    # input setup
    if args.input:
        if args.multi is None and args.tfmatrix is None:
            logname = args.input
        else:
            aobj.ilog_head = eval(args.input)[0]
            aobj.ilog_tail = eval(args.input)[1]
        print('inputlog =', args.input)

    elif args.inputx:
        print('inputlog =', args.inputx)
        args.ilog_head = args.inputx.split('xxx')[0]
        args.ilog_tail = args.inputx.split('xxx')[1]

    else:
        print('Error! No input log name')
        sys.exit()

    print('f90soflag =', args.nof90so)
    print('addresinfo =', args.noresinfo)
    print('zero padding =', args.zp)
    print('monomer energy label =', args.momene)
    print('dimer energy label =', args.dimene)
    print('get monomer dimer energy =', args.is_momdim)

    return args


def setupmode():
    '''setupmode
    Args:
        args
    Returns:
        tgtfrag1, tgtfrag2
    '''

    tgtfrag1 = []
    tgtfrag2 = []

    # setup frag mode
    if args.frag is not None:
        aobj.anlmode = 'frag'
        if len(args.frag) == 2:
            aobj.tgt2type = 'frag'
            tgtfrag1 = args.frag[0]
            tgtfrag2 = args.frag[1]
        if len(args.frag) == 1:
            tgtfrag1 = args.frag[0]

    # check fraginmol mode
    if args.fraginmol is not None:
        aobj.anlmode = 'fraginmol'
        try:
            aobj.tgtmolid = int(args.fraginmol[0])
        except:
            aobj.tgtmolid = args.fraginmol[0]
        aobj.tgt1_lofrag = int(args.fraginmol[1])
        aobj.tgt2molname = args.fraginmol[2]
        aobj.tgt2_lofrag = int(args.fraginmol[3])

    # check molmode
    if args.mol is not None:
        aobj.anlmode = 'mol'
        aobj.tgtmolid = args.mol

    # check dist mode
    if args.dist is not None:
        if args.molname is None:
            aobj.tgt2type = 'dist'
        aobj.dist = args.dist

    # check molname mode
    if args.molname is not None:
        aobj.tgt2type = 'molname'
        aobj.tgt2molname = args.molname

    # check ffmatrix mode
    if args.ffmatrix is not None:
        aobj.anlmode = 'frag'
        aobj.tgt2type = 'frag'
        aobj.matrixtype = 'frags-frags'
        tgtfrag1 = args.ffmatrix[0]
        tgtfrag2 = args.ffmatrix[1]

    # check multimode
    if args.multi is not None:
        aobj.anlmode = 'multi'
        aobj.start = int(args.time[0])
        aobj.end = int(args.time[1])
        aobj.interval = int(args.time[2])
        if len(args.multi) == 1:
            tgtfrag1 = args.multi[0]
        if len(args.multi) == 2:
            tgtfrag1 = args.multi[0]
            tgtfrag2 = args.multi[1]

    # check tfmatrix mode
    if args.tfmatrix != None:
        aobj.anlmode = 'multi'
        aobj.tgt2type='frag'
        aobj.matrixtype='times-frags'
        tgtfrag1 = args.tfmatrix[0]
        tgtfrag2 = args.tfmatrix[1]
        aobj.start = int(args.time[0])
        aobj.end = int(args.time[1])
        aobj.interval = int(args.time[2])

    # check dimeres mode (get or not dimer-es approxymation region)
    if args.dimeres:
        aobj.tgt2type = 'dimer-es'

    aobj.exceptfrag = args.exclude
    aobj.f90soflag = args.nof90so
    aobj.pynp = args.pynp
    aobj.addresinfo = args.noresinfo
    aobj.dimfrag1 = args.dimene[0]
    aobj.dimfrag2 = args.dimene[1]
    aobj.momfrag = args.momene[0]

    # check getting momdimene frag
    aobj.is_momdimene = False
    if args.is_momdim is True:
        aobj.is_momdimene = True

    return tgtfrag1, tgtfrag2


if __name__ == '__main__':
    args = get_args()
    aobj = ampt.anlfmo()

    # setup mode
    tgtfrag1, tgtfrag2 = setupmode()

    # set zero padding
    aobj.zp = args.zp

    # read section
    # multi(read and filter)
    if aobj.anlmode == 'multi':
        if aobj.tgt2type in ['dist', 'molname', 'dimer-es']:
            # multi-fd(frag-dist), f-mname
            aobj = aobj.readifiewrap(tgtfrag1)
        else:
            # multi-ff(frag-frag), tfmatrix
            aobj = aobj.readifiewrap(tgtfrag1, tgtfrag2)

    # single-file(read)
    else:
        # frag-dist
        if aobj.anlmode == 'frag' and aobj.tgt2type == 'dist':
            aobj = aobj.readifiewrap(logname, tgtfrag1)
        # ffmatrix, fragids
        if aobj.anlmode == 'frag' and aobj.tgt2type == 'frag':
            aobj = aobj.readifiewrap(logname, tgtfrag1, tgtfrag2)
        # fraginmol
        if aobj.anlmode == 'fraginmol' or aobj.anlmode == 'mol':
            aobj = aobj.readifiewrap(logname)

    #  filter(for single mode) and write section
    #  with pb
    if aobj.matrixtype == 'frags-frags' and aobj.pbflag:
        aobj = aobj.filterifiewrap()
        aobj.writecsvwrap(word='gas')
        aobj = aobj.filterifiewrap(myifdf=aobj.pbifdf, mypidf=aobj.pbpidf, is_pb=True)
        aobj.writecsvwrap(word='pb', pbwrite=True)

    # only-gas
    else:
        aobj = aobj.filterifiewrap()
        aobj.writecsvwrap()
