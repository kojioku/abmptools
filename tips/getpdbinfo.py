import abmptools as ampt
import os

pdb = ampt.pdb_io()
pdb.readpdb('NMR_1L2Y_10_HET_fixed.pdb')

print ('getmode', pdb.getmode)
print ('assignresname', pdb.assignresname)
print ('refreshatmtype', pdb.refreshatmtype)
print ('refreshresid', pdb.refreshresid)
print ('cellsize', pdb.cellsize)
print ('totalRes', pdb.totalRes)
print ('atomtypelist\n', pdb.atmtypeRes)
print ('Residuenames\n', pdb.resnames)
print ('global atom label\n', pdb.gatmlabRes)
print ('position\n', pdb.posRes)
print ('head\n', pdb.headRes)
print ('labRes\n', pdb.labRes)
print ('chain\n', pdb.chainRes)
print ('residue number\n', pdb.resnumRes)
print ('code\n', pdb.codeRes)
print ('occupancy\n', pdb.occRes)
print ('temperature factor\n', pdb.tempRes)
print ('amark\n', pdb.amarkRes)
print ('charges\n', pdb.chargeRes)
print ('rescount\n', pdb.rescount)

