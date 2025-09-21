import os
import abmptools
import argparse

# Author: Koji Okuwaki
# Date: 2025/09/21


def get_args():
    parser = argparse.ArgumentParser(
                prog='log2config',
                usage='',
                description='Convert log to config dict file',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-np', '--pynp',
                        type=int,
                        help='python np',
                        default=1)

    parser.add_argument('-i', '--input',
                        default=None,
                        required=True,
                        help='input file',
                        )

    parser.add_argument('-o', '--output',
                        help='output config file',
                        default='segment_data.dat')

    args = parser.parse_args()
    print('## setup info')
    print('process =', args.pynp)

    return args


if __name__ == '__main__':
    args = get_args()
    aobj = abmptools.LOGManager()

    # parse
    aobj.parse(args.input)

    # extract
    natoms = aobj.fraginfo['natoms']
    baas = aobj.fraginfo['baas']
    connects = aobj.fraginfo['connects']
    atomlabel = aobj.fraginfo['atomlabel']
    chgs = aobj.fraginfo['chgs']

    # write
    f = open(args.output, 'w')
    print("seg_data = [", file=f)
    head, ext = os.path.splitext(args.input)
    print("    {", file=f)
    print("    'name': '" + head.split('/')[-1] + "',", file=f)
    print("    'atom': ", natoms, ',', file=f)
    print("    'charge': ", chgs, ',', file=f)
    print("    'connect_num': ", baas, ',', file=f)
    print("    'seg_info': ", atomlabel, ',', file=f)
    print("    'connect': ", connects, ',', file=f)
    end =  """    'nummol_seg': [1],
    'repeat': [1],
    'pair_file': [],
    'multi_xyz': 'none'
    },"""
    print(end, file=f)
    print(']', file=f)
    f.close()
    print(f'## output file: {args.output}')
