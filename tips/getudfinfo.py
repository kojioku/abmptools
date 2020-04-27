import numpy as np
from UDFManager import *
import sys
import os
import math
import subprocess
import ampt.udf_io as mbudf

if __name__ == "__main__":

    argvs = sys.argv
    fname = str(argvs[1])

    _udf_ = UDFManager(fname)

    udfio = mbudf.udf_io()
    udfio.getudfinfowrap(_udf_)

