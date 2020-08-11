import sys
import math
import os
import abmptools as ampt
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                prog='generate frag config py file from ABNITMP fragment info (segment_data) file', # program name
                usage='Demonstration of argparser', # program usage
                description='description',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--input',
                        help='input frag info',
                        nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-o', '--output',
                        help='output config file',
                        default='segment_data.dat')


    # get args
    args = parser.parse_args()
    print('input =', args.input)
    print('output =', args.output)


    obj = ampt.abinit_io()
    obj.getfragdict(args.input, args.output)





