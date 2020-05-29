import abmptools as ampt
import os

aobj = ampt.setfmo()

# gen ajf file
aobj.ajf_method = "MP2"
aobj.ajf_basis_set = "6-31G*"
aobj.abinit_ver = 'rev15'
aobj.autofrag = True
aobj.piedaflag = True
aobj.cpfflag = True
aobj.cmmflag = True
aobj.npro = 1
aobj.para_job = 1
aobj.readgeom = 'NMR_1L2Y_10_HET_fixed.pdb'
# aobj.writegeom = os.path.splitext(aobj.readgeom)[0] + '-' + aobj.ajf_method + '-' + aobj.ajf_basis_set.replace('*', 'd') + ".cpf'"

aobj.saveajf()
