import abmptools as ampt

pdb = ampt.pdb_io()
pdb = pdb.readpdb('NMR_1L2Y_10_HET_fixed.pdb')


print (pdb.solutes)
print (pdb.getmode)
print (pdb.assignresname)
print (pdb.refreshatmtype)
print (pdb.refreshresid)
print (pdb.cellsize)
print (pdb.totalRes)
print (pdb.atmtypeRes)
print (pdb.resnames)
print (pdb.gatmlabRes)
print (pdb.posRes)
print (pdb.headRes)
print (pdb.labRes)
print (pdb.chainRes)
print (pdb.resnumRes)
print (pdb.codeRes)
print (pdb.occRes)
print (pdb.tempRes)
print (pdb.amarkRes)
print (pdb.chargeRes)

