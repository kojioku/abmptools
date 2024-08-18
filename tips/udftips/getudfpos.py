from UDFManager import *
import sys
import abmptools

fname = sys.argv[1]
record = 1
mol = 0
print (fname, record, mol)

uobj = UDFManager(fname)
cobj = abmptools.udf_io()
pos = cobj.getposmolrec(uobj, mol, record)

print (pos)
