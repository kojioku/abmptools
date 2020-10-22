import abmptools as ampt
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                prog='pdb2fmo.py', # program name
                usage='python generateajf.py -c xxx.pdb', # program usage
                description='generate ABNITMP input file (ajf)',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-c', '--coord',
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

    parser.add_argument('-cpf', '--cpf',
                        help='cpf',
                        action='store_true'
                        )

    parser.add_argument('-basis', '--basisset',
                        help='basis',
                        default='6-31G*',
                        )

    parser.add_argument('-m', '--method',
                        help='method',
                        default='MP2',
                        )

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

    parser.add_argument('-ma', '--manual',
                        help='manual table',
                        default=None,
                        )

    # get args
    args = parser.parse_args()

    print('coord(pdb) =', args.coord)
    print('solv = ', args.solvation)
    print('pbcnv = ', args.thicnv)
    print('atmrad =', args.atmrad)
    print('ajfversion =', args.ajfversion)
    print('np =', args.npro)
    print('pieda =', args.nopieda)
    print('cmm =', args.cmm)
    print('cpf =', args.cpf)
    print('basis =', args.basisset)
    print('method =', args.method)
    print('dgemm', args.dgemm)
    print('resp', args.resp)
    print('nbo', args.nonbo)
    print('memory', args.memory)
    print('ligand charge', args.ligandcharge)
    print('manual', args.manual)

    aobj = ampt.setfmo()

    # gen ajf file
    aobj.ajf_method = args.method
    aobj.ajf_basis_set = args.basisset
    aobj.abinit_ver = args.ajfversion
    aobj.autofrag = True
    aobj.piedaflag = args.nopieda
    aobj.cpfflag = args.cpf
    aobj.cmmflag = args.cmm
    aobj.npro = args.npro
    aobj.pbmolrad = args.atmrad
    aobj.readgeom = args.coord
    aobj.solv_flag = args.solvation
    aobj.pbcnv = args.thicnv
    aobj.dgemm = args.dgemm
    aobj.resp = args.resp
    aobj.nbo = args.nonbo
    aobj.memory = args.memory
    aobj.ligchg = args.ligandcharge

    # aobj.writegeom = os.path.splitext(aobj.readgeom)[0] + '-' + aobj.ajf_method + '-' + aobj.ajf_basis_set.replace('*', 'd') + ".cpf'"

    if args.manual:
        print('manual mode')
        fragfile = args.manual
        aobj.getfragdict([fragfile], 'segment_data.dat')


        # fatomnums, fchgs, fbaas, fatminfos, connects = aobj.getfragtable(tgtmolsets, atomnumsets, nameidMol)
        head, ext = os.path.splitext(fragfile)
        ftemp = head.split('/')[-1]
        print(ftemp)
        aobj.getfragtable(ftemp)
        # print (frag_atoms, frag_charges)

        aobj.saveajf()


    if aobj.solv_flag == True:
        fname = os.path.splitext(aobj.readgeom)[0] + '-' + aobj.ajf_method + '-' + aobj.ajf_basis_set.replace('*', 'd') + '-pbcnv' + str(aobj.pbcnv) + '-' + aobj.pbmolrad + '.ajf'
        aobj.saveajf(fname)
    else:
        aobj.saveajf()
