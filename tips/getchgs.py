import abmptools as ampt
import sys

fname = sys.argv[1]

aobj = ampt.anlfmo()
natom = aobj.getlognatom(fname)
chgs = aobj.getlogchg(fname, natom)

print('natom', natom)
print('chgs(resp)', chgs)
