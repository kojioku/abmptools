import sys
import math
import os
import abmptools as ampt

if __name__ == '__main__':
    argvs = sys.argv
    fnames = argvs[1:]
    print(fnames)

    ofile = 'segment_data.dat'

    obj = ampt.abinit_io()
    obj.getfragdict(fnames, ofile)





