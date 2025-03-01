#! /usr/bin/python
#
#####################################################################
## NOTE: This code is PRELIMINARY, it requires MANUAL VERIFICATION ##
##                    ***********              ******************* ##
#####################################################################
#

import re
import sys

def wsadminToList(inStr):
    outList=[]
    if (len(inStr)>0 and inStr[0]=='[' and inStr[-1]==']'):
        tmpList = inStr[1:-1].split() #splits space-separated lists,
    else:
        tmpList = inStr.split("\n")   #splits for Windows or Linux
    for item in tmpList:
        item = item.rstrip();         #removes any Windows "\r"
        if (len(item)>0):
            outList.append(item)
    return outList
#endDef

#!/usr/bin/tclsh

###########################################################################
##########                                                       ##########
##########          SCRIPT FOR MAKING INPUT OF ABINITMP          ##########
##########                   ( CREST-VERSION )                   ##########
##########                                                       ##########
##########              ..... version 20060512 .....             ##########
##########                                                       ##########
###########################################################################
#
#  Ver.4.2 => Ver.7.0
#  Move namelist 'MDCNTRL' above namelist 'XYZ'
#

#######################################################
##########    original_to_abinitmp_crest     ##########
#######################################################

def original_to_abinitmp_crest (  ):

    global nline_in
    global line_in
    global iline

    global n_cntrl, cntrl
    n_cntrl = 0
    global n_fmocntrl, fmocntrl
    n_fmocntrl = 0
    global n_scf, scf
    n_scf = 0
    global n_basis, basis
    n_basis = 0
    global n_optcntrl, optcntrl
    n_optcntrl = 0
    global n_mfmo, mfmo
    n_mfmo = 0
    global n_xuff, xuff
    n_xuff = 0
    global n_sczv, sczv
    n_sczv = 0
    global n_mp2, mp2
    n_mp2 = 0
    global n_mp2dns, mp2dns
    n_mp2dns = 0
    global n_mp2grd, mp2grd
    n_mp2grd = 0
    global n_mp3, mp3
    n_mp3 = 0
    global n_lmp2, lmp2
    n_lmp2 = 0
    global n_lrd, lrd
    n_lrd = 0
    global n_dft, dft
    n_dft = 0
    global n_analysis, analysis
    n_analysis = 0
    global n_bsse, bsse
    n_bsse = 0
    global n_solvation, solvation
    n_solvation = 0
    global n_pbeq, pbeq
    n_pbeq = 0
    global n_pop, pop
    n_pop = 0
    global n_gridcntrl, gridcntrl
    n_gridcntrl = 0

    global n_mcp, mcp
    n_mcp = 0
    global n_cis, cis
    n_cis = 0
    global n_cisgrd, cisgrd
    n_cisgrd = 0
    global n_cafi, cafi
    n_cafi = 0
    global n_pol, pol
    n_pol = 0
    global n_gf2, gf2
    n_gf2 = 0
    global n_ccpt, ccpt
    n_ccpt = 0

    # MD-RELATED ADDED 2013.06.14
    global n_mdcntrl, mdcntrl
    n_mdcntrl = 0

    global buffer

    global xyz_line1
    xyz_line1 = 1
    global xyz_line2
    xyz_line2 = 0
    global chk_coord
    global fragment_line1
    fragment_line1 = 1
    global fragment_line2
    fragment_line2 = 0
    global fragpair_line1
    fragpair_line1 = 1
    global fragpair_line2
    fragpair_line2 = 0

    # MD-RELATED ADDED 2013.06.14
    global vel_line1
    vel_line1 = 1
    global vel_line2
    vel_line2 = 0
    global nhc_line1
    nhc_line1 = 1
    global nhc_line2
    nhc_line2 = 0
    global typfrag_line1
    typfrag_line1 = 1
    global typfrag_line2
    typfrag_line2 = 0

    # EXPROP, COUPLING ADDED 2023.08.10
    global n_exprop, exprop
    n_exprop = 0
    global n_coupling, coupling
    n_coupling = 0

    global atom
    namelists = 1
    values = 2

    i = 0  #forStart
    while ( i < nline_in ):  #forTest

        if re.search('&CNTRL$', line_in[i], re.IGNORECASE) is not None:
            n_cntrl, cntrl = read_item(namelists, values, i)

        if re.search('&FMOCNTRL$', line_in[i], re.IGNORECASE) is not None:
            n_fmocntrl, fmocntrl = read_item(namelists, values, i)

        if re.search('&MFMO$', line_in[i], re.IGNORECASE) is not None:
            n_mfmo, mfmo = read_item(namelists, values, i)

        if re.search('&EXPROP$', line_in[i], re.IGNORECASE) is not None:
            n_exprop, exprop = read_item(namelists, values, i)

        if re.search('&BASIS$', line_in[i], re.IGNORECASE) is not None:
            n_basis, basis = read_item(namelists, values, i)

        if re.search('&MCP$', line_in[i], re.IGNORECASE) is not None:
            n_mcp, mcp = read_item(namelists, values, i)

        if re.search('&SCF$', line_in[i], re.IGNORECASE) is not None:
            n_scf, scf = read_item(namelists, values, i)

        if re.search('&MP2$', line_in[i], re.IGNORECASE) is not None or \
           re.search('&MP2D$', line_in[i], re.IGNORECASE) is not None or \
           re.search('&MP2G$', line_in[i], re.IGNORECASE) is not None:
            n_mp2, mp2 = read_item(namelists, values, i)

        if re.search('&XUFF$', line_in[i], re.IGNORECASE) is not None:
            n_xuff, xuff = read_item(namelists, values, i)

        if re.search('&SCZV$', line_in[i], re.IGNORECASE) is not None:
            n_sczv, sczv = read_item(namelists, values, i)

        if re.search('&LMP2$', line_in[i], re.IGNORECASE) is not None:
            n_lmp2, lmp2 = read_item(namelists, values, i)

        if re.search('&LRD$', line_in[i], re.IGNORECASE) is not None:
            n_lrd, lrd = read_item(namelists, values, i)

        if re.search('&CIS$', line_in[i], re.IGNORECASE) is not None or \
                re.search('&CISG$', line_in[i], re.IGNORECASE) is not None:
            n_cis, cis = read_item(namelists, values, i)

        if re.search('&CISGRD$', line_in[i], re.IGNORECASE) is not None:
            n_cisgrd, cisgrd = read_item(namelists, values, i)

        if re.search('&CAFI$', line_in[i], re.IGNORECASE) is not None:
            n_cafi, cafi = read_item(namelists, values, i)

        if re.search('&POL$', line_in[i], re.IGNORECASE) is not None:
            n_pol, pol = read_item(namelists, values, i)

        if re.search('&OPTCNTRL$', line_in[i], re.IGNORECASE) is not None:
            n_optcntrl, optcntrl = read_item(namelists, values, i)

        if re.search('&GRIDCNTRL$', line_in[i], re.IGNORECASE) is not None:
            n_gridcntrl, gridcntrl = read_item(namelists, values, i)

        if re.search('&GF2$', line_in[i], re.IGNORECASE) is not None:
            n_gf2, gf2 = read_item(namelists, values, i)

        if re.search('&MP3$', line_in[i], re.IGNORECASE) is not None:
            n_mp3, mp3 = read_item(namelists, values, i)

        if re.search('&MP2DNS$', line_in[i], re.IGNORECASE) is not None:
            n_mp2dns, mp2dns = read_item(namelists, values, i)

        if re.search('&MP2GRD$', line_in[i], re.IGNORECASE) is not None:
            n_mp2grd, mp2grd = read_item(namelists, values, i)

        if re.search('&CCPT$', line_in[i], re.IGNORECASE) is not None:
            n_ccpt, ccpt = read_item(namelists, values, i)

        if re.search('&DFT$', line_in[i], re.IGNORECASE) is not None:
            n_dft, dft = read_item(namelists, values, i)

        if re.search('&ANALYSIS$', line_in[i], re.IGNORECASE) is not None:
            n_analysis, analysis = read_item(namelists, values, i)

        if re.search('&BSSE$', line_in[i], re.IGNORECASE) is not None:
            n_bsse, bsse = read_item(namelists, values, i)

        if re.search('&SOLVATION$', line_in[i], re.IGNORECASE) is not None:
            n_solvation, solvation = read_item(namelists, values, i)

        if re.search('&PBE$', line_in[i], re.IGNORECASE) is not None:
            n_pbe, pbe = read_item(namelists, values, i)

        if re.search('&POP$', line_in[i], re.IGNORECASE) is not None:
            n_pop, pop = read_item(namelists, values, i)

        if re.search('&COUPLING$', line_in[i], re.IGNORECASE) is not None:
            n_coupling, coupling = read_item(namelists, values, i)

        if re.search('&MDCNTRL$', line_in[i], re.IGNORECASE) is not None:
            n_mdcntrl, mdcntrl = read_item(namelists, values, i)

        if re.search('&VEL$', line_in[i], re.IGNORECASE) is not None:
            vel_line1, vel_line2 = get_linenumber(i)

        if re.search('&NHC$', line_in[i], re.IGNORECASE) is not None:
            nhc_line1, nhc_line2 = get_linenumber(i)

        if re.search('&TYPFRAG$', line_in[i], re.IGNORECASE) is not None:
            typfrag_line1, typfrag_line2 = get_linenumber(i)

        ###############
        ##### XYZ #####
        ###############
        if re.search('&XYZ', line_in[i], re.IGNORECASE) is not None:
            if (len(line_in[i]) > 1):

                ### GAMESS OUTPUT TYPE COORDINATE
                if re.search('GAMESS-OUT', line_in[i], re.IGNORECASE) is not None:
                    chk_coord = 1

                ### GAUSSIAN OUTPUT TYPE COORDINATE
                if re.search('GAUSSIAN-OUT', line_in[i], re.IGNORECASE) is not None:
                    chk_coord = 2

            xyz_line1 = i + 1
            xyz_line2 = get_end(i) - 1
            iline = xyz_line2

        if re.search('&FRAGMENT$', line_in[i], re.IGNORECASE) is not None:
            fragment_line1, fragment_line2 = get_linenumber(i)

        if re.search('&FRAGPAIR$', line_in[i], re.IGNORECASE) is not None:
            fragpair_line1, fragpair_line2 = get_linenumber(i)

        i += 1  #forNext

####################################
##########    get_end     ##########
####################################

def get_end ( iline ):

    global nline_in, line_in

    i = iline  #forStart
    while ( i < nline_in ):  #forTest
        _J2J_bracket_ = line_in[i].find("&END")
        _J2J_bracket2_ = line_in[i].find("&End")
        if (_J2J_bracket_ >= 0  or _J2J_bracket2_ >= 0  or line_in[i].find("&end") >= 0):
            return i
        #endIf
        if (line_in[i].find("/") >= 0):
            if (line_in[i].rstrip()[-1] == "/"):
                return i
            #endIf
        #endIf
        i += 1  #forNext
    #endWhile  (#endFor)
#endDef

########################################
##########    get_namelist    ##########
########################################

def get_namelist ( line1, line2 ):

    global line_in, buffer

    n_name_list = 0

    namelists = 1
    values = 2

    buffer = [["-" for i in range(200)] for j in range(200)]

    i = line1  #forStart
    while ( i <= line2 ):  #forTest

        leng = len(line_in[i])

        j = 0  #forStart
        while ( j < leng ):  #forTest
            if (line_in[i][j] == "="):

                ### get namelist ###
                tmp1 = line_in[i][0:j].rstrip()
                #tmp2 = tmp1[(len(tmp1) -1)]
                #?PROBLEM? (jacl 734) previous LINDEX may need variable.split("xx")
                namelist = tmp1

                ### get value ###
                tmp1 = line_in[i][(j+1):leng + 1].lstrip()
                #tmp2 = tmp1[0]
                #?PROBLEM? (jacl 739) previous LINDEX may need variable.split("xx")
                value = tmp1
                first_char = value[0]
                last_char = value[(len(value) -1)]
                if (first_char == "\'" and last_char != "\'"):
                    tmp1_leng = len(tmp1)
                    k = 2  #forStart
                    while ( k < tmp1_leng ):  #forTest
                        if (tmp1[k] == "\'"):
                            pos = k
                            break
                        #endIf
                        k += 1  #forNext
                    #endWhile  (#endFor)
                    value = tmp1[:pos + 1]
                #endIf
                if (first_char == "\"" and last_char != "\""):
                    tmp1_leng = len(tmp1)
                    k = 2  #forStart
                    while ( k <= tmp1_leng ):  #forTest
                        if (tmp1[k] == "\""):
                            pos = k
                            break
                        #endIf
                        k += 1  #forNext
                    #endWhile  (#endFor)
                    value = tmp1[:pos + 1]
                #endIf
                ### set namelist and value ###
                if (namelist.find("=") < 0):
                    namelist = namelist.lstrip()
                    value = value.rstrip()
                    n_name_list = (n_name_list+1)
                    buffer[n_name_list][namelists] =  namelist
                    buffer[n_name_list][values] = value
                #endIf
            #endIf
            j += 1  #forNext
        #endWhile  (#endFor)
        i += 1  #forNext
    #endWhile  (#endFor)
    return n_name_list
#endDef

def gamess_out_to_abinitmp ( gamess_in ):

    global atom, atom_count

    atomic_number = int(gamess_in[1] + 0.1)
    #?PROBLEM? (jacl 779) previous LINDEX may need variable.split("xx")
    coord_x = gamess_in[2]
    #?PROBLEM? (jacl 780) previous LINDEX may need variable.split("xx")
    coord_y = gamess_in[3]
    #?PROBLEM? (jacl 781) previous LINDEX may need variable.split("xx")
    coord_z = gamess_in[4]
    #?PROBLEM? (jacl 782) previous LINDEX may need variable.split("xx")
    atomic_symbol = atom[atomic_number][symbol]

    atom_count = atom_count + 1

    abinit_out = "%6s %6s %6s %14.8f %14.8f %14.8f %2s" % (atom_count, atomic_symbol, "1", coord_x, coord_y, coord_z, "1")

    return abinit_out
#endDef

def gaussian_out_to_abinitmp ( gaussian_in ):

    global atom, atom_count

    atomic_number = gaussian_in[1]
    #?PROBLEM? (jacl 799) previous LINDEX may need variable.split("xx")
    coord_x = gaussian_in[3]
    #?PROBLEM? (jacl 800) previous LINDEX may need variable.split("xx")
    coord_y = gaussian_in[4]
    #?PROBLEM? (jacl 801) previous LINDEX may need variable.split("xx")
    coord_z = gaussian_in[5]
    #?PROBLEM? (jacl 802) previous LINDEX may need variable.split("xx")
    atomic_symbol = atom[atomic_number][symbol]

    atom_count = atom_count + 1

    abinit_out = "%6s %6s %6s %14.8f %14.8f %14.8f %2s" % (atom_count, atomic_symbol, "1", coord_x, coord_y, coord_z, "1")

    return abinit_out
#endDef

def mk_atomic_table():
    symbol = 3

    global atom

    ## initialize !!
    i = 1  #forStart
    #while ( i <= 111 ):  #forTest
            #atom[i][symbol] = "-"
            #i += 1  #forNext
    #endWhile  (#endFor)

    atom = [["-" for i in range(4)] for j in range(112)]

    atom[1][symbol] = "H"
    atom[2][symbol] = "He"
    atom[3][symbol] = "Li"
    atom[4][symbol] = "Be"
    atom[5][symbol] = "B"
    atom[6][symbol] = "C"
    atom[7][symbol] = "N"
    atom[8][symbol] = "O"
    atom[9][symbol] = "F"
    atom[10][symbol] = "Ne"
    atom[11][symbol] = "Na"
    atom[12][symbol] = "Mg"
    atom[13][symbol] = "Al"
    atom[14][symbol] = "Si"
    atom[15][symbol] = "P"
    atom[16][symbol] = "S"
    atom[17][symbol] = "Cl"
    atom[18][symbol] = "Ar"
    atom[19][symbol] = "K"
    atom[20][symbol] = "Ca"
    atom[21][symbol] = "Sc"
    atom[22][symbol] = "Ti"
    atom[23][symbol] = "V"
    atom[24][symbol] = "Cr"
    atom[25][symbol] = "Mn"
    atom[26][symbol] = "Fe"
    atom[27][symbol] = "Co"
    atom[28][symbol] = "Ni"
    atom[29][symbol] = "Cu"
    atom[30][symbol] = "Zn"
    atom[31][symbol] = "Ga"
    atom[32][symbol] = "Ge"
    atom[33][symbol] = "As"
    atom[34][symbol] = "Se"
    atom[35][symbol] = "Br"
    atom[36][symbol] = "Kr"
    atom[37][symbol] = "Rb"
    atom[38][symbol] = "Sr"
    atom[39][symbol] = "Y"
    atom[40][symbol] = "Zr"
    atom[41][symbol] = "Nb"
    atom[42][symbol] = "Mo"
    atom[43][symbol] = "Tc"
    atom[44][symbol] = "Ru"
    atom[45][symbol] = "Rh"
    atom[46][symbol] = "Pd"
    atom[47][symbol] = "Ag"
    atom[48][symbol] = "Cd"
    atom[49][symbol] = "In"
    atom[50][symbol] = "Sn"
    atom[51][symbol] = "Sb"
    atom[52][symbol] = "Te"
    atom[53][symbol] = "I"
    atom[54][symbol] = "Xe"
    atom[55][symbol] = "Cs"
    atom[56][symbol] = "Ba"
    atom[57][symbol] = "La"
    atom[58][symbol] = "Ce"
    atom[59][symbol] = "Pr"
    atom[60][symbol] = "Nd"
    atom[61][symbol] = "Pm"
    atom[62][symbol] = "Sm"
    atom[63][symbol] = "Eu"
    atom[64][symbol] = "Gd"
    atom[65][symbol] = "Tb"
    atom[66][symbol] = "Dy"
    atom[67][symbol] = "Ho"
    atom[68][symbol] = "Er"
    atom[69][symbol] = "Tm"
    atom[70][symbol] = "Yb"
    atom[71][symbol] = "Lu"
    atom[72][symbol] = "Hf"
    atom[73][symbol] = "Ta"
    atom[74][symbol] = "W"
    atom[75][symbol] = "Re"
    atom[76][symbol] = "Os"
    atom[77][symbol] = "Ir"
    atom[78][symbol] = "Pt"
    atom[79][symbol] = "Au"
    atom[80][symbol] = "Hg"
    atom[81][symbol] = "Tl"
    atom[82][symbol] = "Pb"
    atom[83][symbol] = "Bi"
    atom[84][symbol] = "Po"
    atom[85][symbol] = "At"
    atom[86][symbol] = "Rn"
    atom[87][symbol] = "Fr"
    atom[88][symbol] = "Ra"
    atom[89][symbol] = "Ac"
    atom[90][symbol] = "Th"
    atom[91][symbol] = "Pa"
    atom[92][symbol] = "U"
    atom[93][symbol] = "Np"
    atom[94][symbol] = "Pu"
    atom[95][symbol] = "Am"
    atom[96][symbol] = "Cm"
    atom[97][symbol] = "Bk"
    atom[98][symbol] = "Cf"
    atom[99][symbol] = "Es"
    atom[100][symbol] = "Fm"
    atom[101][symbol] = "Md"
    atom[102][symbol] = "No"
    atom[103][symbol] = "Lr"
    atom[104][symbol] = "Rf"
    atom[105][symbol] = "Db"
    atom[106][symbol] = "Sg"
    atom[107][symbol] = "Bh"
    atom[108][symbol] = "Hs"
    atom[109][symbol] = "Mt"
    atom[110][symbol] = "Ds"
    atom[111][symbol] = "Rg"

#endDef


def write_item(item, namelists, values):
    """
    global n_cntrl, cntrl
    """
    # n_item = n_ctrl
    n_item = eval('n_' + item)
    if n_item != 0:
        # content = ctrl
        content = eval(item)

    print("&" + item.upper())
    i = 1
    while (i <= n_item):
        print("  %s=%s" % (content[i][namelists], content[i][values]))
        i += 1
    print("/")
    return


def write_item_format2(item, line_in):
    print("&" + item.upper())
    # fragment_line1
    i = eval(item + '_line1')
    # fragment_line2
    line2 = eval(item + '_line2')
    while (i <= line2):
        print(line_in[i].rstrip())
        i += 1
    print("/")
    return


def read_item(namelists, values, i):
    line1 = i
    line2 = get_end(i)
    iline = line2
    n_item = get_namelist(line1, line2)
    j = 1
    item = [["-" for l in range(200)] for m in range(200)]
    while (j <= n_item):
        item[j][namelists] = buffer[j][namelists]
        item[j][values] = buffer[j][values]
        j += 1
    return n_item, item


def get_linenumber(i):
    line1 = i + 1
    line2 = get_end(i) - 1
    iline = line2
    return line1, line2


#################################
##########    MAIN     ##########
#################################
if __name__ == '__main__':
    ###
    ### initialize
    ###
    
    nline_in = 0
    iline = 0
    line_in = []
    
    #saitou+
    namelists = 1
    values = 2
    symbols = 3
    
    ### chk_coord --> 0:ABINIT-MP , 1:GAMESS-OUT , 2:GAUSSIAN-OUT
    chk_coord = 0
    atom_count = 0
    
    ### set atomic table
    mk_atomic_table()
    
    ###
    ### read original input
    ###
    
    lines = sys.stdin.readlines()
    for line in lines:
        temp = line.strip()
        if (len(temp) > 0):
                HeadChar = temp.strip( )[0]
                if (HeadChar != "*" or
                  HeadChar != "!" or
                  HeadChar != "#"):
                        nline_in = (nline_in + 1)
                        line_in.insert(nline_in,line)
                #endIf
        #endIf
    #endWhile
    
    ###
    ### original input to abinitmp ( crest-version ) input
    ###
    
    original_to_abinitmp_crest()
    
    ###
    ### out abinitmp ( crest-version ) ajf-file
    ###
    
    write_item('cntrl', namelists, values)
    write_item('fmocntrl', namelists, values)
    write_item("scf", namelists, values)
    write_item("basis", namelists, values)
    write_item("optcntrl", namelists, values)
    write_item("mfmo", namelists, values)
    write_item("exprop", namelists, values)
    write_item("sczv", namelists, values)
    write_item("xuff", namelists, values)
    write_item("mp2", namelists, values)
    write_item("mp2dns", namelists, values)
    write_item("mp2grd", namelists, values)
    write_item("mp3", namelists, values)
    write_item("lmp2", namelists, values)
    write_item("lrd", namelists, values)
    write_item("dft", namelists, values)
    write_item("analysis", namelists, values)
    write_item("bsse", namelists, values)
    write_item_format2('fragpair', line_in)
    write_item("solvation", namelists, values)
    write_item("pbeq", namelists, values)
    write_item("pop", namelists, values)
    write_item("coupling", namelists, values)
    write_item("gridcntrl", namelists, values)
    write_item("mcp", namelists, values)
    write_item("cis", namelists, values)
    write_item("cisgrd", namelists, values)
    write_item("cafi", namelists, values)
    write_item("pol", namelists, values)
    write_item("gf2", namelists, values)
    write_item("ccpt", namelists, values)
    
    ## XYZ
    print("&XYZ")
    i = xyz_line1  #forStart
    while ( i <= xyz_line2 ):  #forTest
    
        if ( chk_coord == 1 ):
                print(gamess_out_to_abinitmp(line_in[i].rstrip()))
        elif ( chk_coord == 2 ):
                print(gaussian_out_to_abinitmp(line_in[i].rstrip()))
        else:
                print(line_in[i].rstrip())
        #endIf
    
        i += 1  #forNext
    #endWhile  (#endFor)
    print("/")
    
    write_item_format2('fragment', line_in)
    write_item("mdcntrl", namelists, values)
    write_item_format2('vel', line_in)
    write_item_format2('nhc', line_in)
    write_item_format2('typfrag', line_in)
