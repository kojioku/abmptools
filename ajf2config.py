import sys
import math
import os
import ampt.abinit_io as abio

if __name__ == '__main__':
    argvs = sys.argv
    fnames = argvs[1:]
    print(fnames)

    ofile = 'segment_data.dat'

    obj = abio.abinit_io()
    obj.getfragdict(fnames, ofile)





