import abmptools as ampt
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                prog='pdb2fmo.py', # program name
                usage='python generateajf.py -i xxx.pdb', # program usage
                description='generate ABNITMP input file (ajf)',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--incoord',
                        help='coordinate file (pdb)',
                        # nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-pb', '--solvation',
                        help='parameter file',
                        action='store_true')

    parser.add_argument('-th', '--thicnv',
                        help='threshold',
                        type=float,
                        default=1.0E-5
                        )

    parser.add_argument('-arad', '--atmrad',
                        help='atmrad',
                        default='delphi'
                        )

    parser.add_argument('-ajfv', '--ajfversion',
                        help='ajf version',
                        default='rev23'
                        )

    parser.add_argument('-np', '--npro',
                        help='ajf version',
                        type=int,
                        default=4
                        )

    parser.add_argument('-nopieda', '--nopieda',
                        help='pieda',
                        action='store_false'
                        )

    parser.add_argument('-cmm', '--cmm',
                        help='cmm',
                        action='store_true'
                        )

    parser.add_argument('-nocpf', '--nocpf',
                        help='cpfflag',
                        action='store_false'
                        )

    parser.add_argument('-cpfv', '--cpfver',
                        help='cpf version',
                        default='10'
                        )

    parser.add_argument('-basis', '--basisset',
                        help='basis',
                        default='6-31G*',
                        )

    parser.add_argument('-m', '--method',
                        help='method',
                        default='MP2',
                        )

    parser.add_argument('-ml', '--mldat',
                        help='mldat flag',
                        action='store_true'
                        )

    parser.add_argument('-mll', '--mllimit',
                        help='mldat fraglimit',
                        type=int,
                        default=0
                        )

    parser.add_argument('-disp', '--disp',
                        help='flag disp',
                        action='store_true')

    # WriteMLdata='wstr-1E08_HIS_ES.new2.cmm5.mldat'
    # MLfraglimit=921

    parser.add_argument('-dg', '--dgemm',
                        help='dgemm',
                        action='store_true',
                        )

    parser.add_argument('-rp', '--resp',
                        help='resp',
                        action='store_true',
                        )

    parser.add_argument('-nonbo', '--nonbo',
                        help='nonbo',
                        action='store_false',
                        )

    parser.add_argument('-mem', '--memory',
                        help='memory',
                        default='3000',
                        )

    parser.add_argument('-lc', '--ligandcharge',
                        help='ligand charge',
                        nargs=2,
                        action='append',
                        # default='',
                        )

    parser.add_argument('-rs', '--rsolv',
                        help='rsolv',
                        nargs=2,
                        action='append',
                        # default='',
                        )

    parser.add_argument('-ma', '--manual',
                        help='manual table',
                        default=None,
                        )

    parser.add_argument('-bsse', '--bsse',
                        help='bsse',
                        action='store_true',
                        )

    # get args
    args = parser.parse_args()

    print('coord(pdb) =', args.incoord)
    print('solv = ', args.solvation)
    print('pbcnv = ', args.thicnv)
    print('atmrad =', args.atmrad)
    print('ajfversion =', args.ajfversion)
    print('np =', args.npro)
    print('pieda =', args.nopieda)
    print('cmm =', args.cmm)
    print('cpf =', args.nocpf)
    print('cpfver =', args.cpfver)
    print('basis =', args.basisset)
    print('method =', args.method)
    print('dgemm', args.dgemm)
    print('resp', args.resp)
    print('nbo', args.nonbo)
    print('memory', args.memory)
    print('ligand charge', args.ligandcharge)
    print('rsolv', args.rsolv)
    print('manual', args.manual)
    print('bsse', args.bsse)
    print('mldat', args.mldat)
    print('disp', args.disp)


    aobj = ampt.setfmo()

    # gen ajf file
    aobj.ajf_method = args.method
    aobj.ajf_basis_set = args.basisset
    aobj.abinit_ver = args.ajfversion
    aobj.autofrag = True
    aobj.piedaflag = args.nopieda
    aobj.cpfflag = args.nocpf
    aobj.cpfver = args.cpfver
    aobj.cmmflag = args.cmm
    aobj.npro = args.npro
    aobj.pbmolrad = args.atmrad
    aobj.readgeom = args.incoord
    aobj.solv_flag = args.solvation
    aobj.pbcnv = args.thicnv
    aobj.dgemm = args.dgemm
    aobj.resp = args.resp
    aobj.nbo = args.nonbo
    aobj.memory = args.memory
    aobj.ligchg = args.ligandcharge
    aobj.rsolv = args.rsolv
    aobj.bsseflag = args.bsse
    aobj.disp = args.disp

    if args.mldat:
        aobj.mldatfrag = args.mldat
        aobj.mllimit = args.mllimit

    # aobj.writegeom = os.path.splitext(aobj.readgeom)[0] + '-' + aobj.ajf_method + '-' + aobj.ajf_basis_set.replace('*', 'd') + ".cpf'"

    if args.manual:
        aobj.autofrag = False
        print('manual mode')
        fragfile = args.manual
        aobj.getfragdict([fragfile], 'segment_data.dat')


        # fatomnums, fchgs, fbaas, fatminfos, connects = aobj.getfragtable(tgtmolsets, atomnumsets, nameidMol)
        head, ext = os.path.splitext(fragfile)
        ftemp = head.split('/')[-1]
        print(ftemp)
        aobj.getfragtable([ftemp])
        # print (frag_atoms, frag_charges)

        aobj.saveajf()


    addstr = ''

    if aobj.solv_flag:
        addstr += '-pbcnv' + str(aobj.pbcnv) + '-' + aobj.pbmolrad

    if aobj.bsseflag:
        addstr += '-bsse'

    if aobj.piedaflag is False:
        addstr += '-nopieda'

    if aobj.cpfflag is False:
        addstr += '-nocpf'

    if aobj.nbo is True:
        addstr += '-nbo'

    if aobj.resp is True:
        addstr += '-resp'

    fname = os.path.splitext(aobj.readgeom)[0] + '-' + aobj.ajf_method + '-' + aobj.ajf_basis_set.replace('*', 'd') + addstr
    ajfname = fname + '.ajf'
    aobj.writegeom = "'" + fname + ".cpf'"
    if aobj.mllimit == 0:
        aobj.mldatname = "'" + fname + ".mldat'"
    else:
        aobj.mldatname = "'" + fname + "limit" + str(aobj.mllimit) + ".mldat'"
    aobj.saveajf(ajfname)


    '''
    pb
    bsse
    nopieda
    nbo
    esp
    '''
