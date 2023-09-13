import sys
import os
import pandas as pd
import itertools
import copy
import csv
import abmptools
import argparse
# Author: Koji Okuwaki


def get_args():
    parser = argparse.ArgumentParser(
                prog='log2cpf',
                usage='',
                description='Analysis script for ABINIT-MP log',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-f', '--frag',
                        help='tgt fragid info',
                        nargs='*'
                        )

    parser.add_argument('-np', '--pynp',
                        type=int,
                        help='python np',
                        default=1)

    parser.add_argument('-i', '--input',
                        default=None,
                        required=True,
                        help='input file',
                        )

    parser.add_argument('-nof90', '--nof90so',
                        help='use f90',
                        action='store_false',
                        default=True)

    args = parser.parse_args()
    print('## setup info')
    print('tgtfrag =', args.frag)
    print('process =', args.pynp)

    return args


if __name__ == '__main__':
    args = get_args()
    aobj = abmptools.LOGManager()

    # aobj.readifiewrap(args.input)
    # print(aobj.ifdf)
    # print(aobj.pidf)

    # def read_fraginfo(self, fname):
    aobj.parse(args.input)
#     frags = aobj.readfraginfo(args.input)
#     print(frags)
#     nf = aobj.getlognf(args.input, 'auto')
#     print(nf)

    # required info
    # self.cpfver = cpfver
    # self.atominfo = pd.DataFrame(atominfo)
    # self.fraginfo = fraginfo
    # self.condition = condition
    # self.static_data = static_data
    # self.mominfo = pd.DataFrame(mominfo)
    # self.diminfo = pd.DataFrame(diminfo)
    # self.labels = labels

    # frag-dist
#     if aobj.anlmode == 'frag' and aobj.tgt2type == 'dist':
#         aobj = aobj.readifiewrap(aobj.ilog_head, tgtfrag1)
#     # ffmatrix, fragids
#     if aobj.anlmode == 'frag' and aobj.tgt2type == 'frag':
#         aobj = aobj.readifiewrap(aobj.ilog_head, tgtfrag1, tgtfrag2)
#     # fraginmol
#     if aobj.anlmode == 'fraginmol' or aobj.anlmode == 'mol':
#         aobj = aobj.readifiewrap(aobj.ilog_head)

#     #  filter(for single mode) and write section
#     #  with pb
#     if aobj.matrixtype == 'frags-frags' and aobj.pbflag:
#         aobj = aobj.filterifiewrap()
#         aobj.writecsvwrap(word='gas')
#         aobj = aobj.filterifiewrap(
#             myifdf=aobj.pbifdf, mypidf=aobj.pbpidf, is_pb=True)
#         aobj.writecsvwrap(word='pb', pbwrite=True)
#
#     # only-gas
#     else:
#         aobj = aobj.filterifiewrap()
#         aobj.writecsvwrap()
