import abmptools as ampt
import os

aobj = ampt.setfmo()
aobj = aobj.readpdb('NMR_1L2Y_10_HET_fixed.pdb')
# aobj = aobj.


print (aobj.solutes)
print (aobj.getmode)
print (aobj.assignresname)
print (aobj.refreshatmtype)
print (aobj.refreshresid)
print (aobj.cellsize)
print (aobj.totalRes)
print (aobj.atmtypeRes)
print (aobj.resnames)
print (aobj.gatmlabRes)
print (aobj.posRes)
print (aobj.headRes)
print (aobj.labRes)
print (aobj.chainRes)
print (aobj.resnumRes)
print (aobj.codeRes)
print (aobj.occRes)
print (aobj.tempRes)
print (aobj.amarkRes)
print (aobj.chargeRes)


# gen ajf file
aobj.ajf_method = "MP2"
aobj.ajf_basis_set = "6-31G*"
aobj.abinit_ver = 'rev15'
aobj.autofrag = True
aobj.piedaflag = False
aobj.cpfflag = False
aobj.cmmflag = True
aobj.npro = 1
aobj.para_job = 1
aobj.readgeom = 'NMR_1L2Y_10_HET_fixed.pdb'
# aobj.writegeom = os.path.splitext(aobj.readgeom)[0] + '-' + aobj.ajf_method + '-' + aobj.ajf_basis_set.replace('*', 'd') + ".cpf'"

aobj.saveajf()
