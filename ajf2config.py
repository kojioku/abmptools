import sys
import math
import os

def flatten(nested_list):
    """2重のリストをフラットにする関数"""
    return [int(e) for inner_list in nested_list for e in inner_list]

# def listtoint(nested_list):
#     """リストの中をintegerに"""
#     return [int(e) for inner_list in nested_list for e in inner_list]

def functor(f, l):
    if isinstance(l,list):
        return [functor(f,i) for i in l]
    else:
        return f(l)

if __name__ == '__main__':
    argvs = sys.argv
    fnames = argvs[1:]
    print(fnames)

    f = open('segment_data.dat', 'w')
    print("seg_data = [", file=f)

    for fname in fnames:
        lines = open(fname, 'r').readlines()
        head, ext = os.path.splitext(fname)

        flag = False
        baseline = 1

        datas = []
        fatomnums = []
        fchgs = []
        fbaas = []
        nfcount = 0
        fatminfo = []
        fatminfos = []
        connects = []
        for i in range(len(lines)):
            itemlist = lines[i].split()
            if i == 0:
                nf = int(itemlist[0])
                if nf > 10:
                    baseline = math.ceil(nf/10)
                    print('baseline', baseline)
            if itemlist[0] == '&FRAGMENT':
                flag = True
                typcount = 0
                lcount = 0
                continue
            # fatom section
            if flag == True and typcount == 0:
                fatomnums.append(itemlist)
                lcount += 1
                if lcount == baseline:
                    typcount += 1
                    lcount = 0
                    fatomnums = flatten(fatomnums)
                    continue
            # chg section
            if flag == True and typcount == 1:
                fchgs.append(itemlist)
                lcount += 1
                if lcount == baseline:
                    typcount += 1
                    lcount = 0
                    fchgs = flatten(fchgs)
                    continue

            # bda section
            if flag == True and typcount == 2:
                fbaas.append(itemlist)
                lcount += 1
                if lcount == baseline:
                    typcount += 1
                    lcount = 0
                    fbaas = flatten(fbaas)
                    continue

            # seginfo section
            if flag == True and typcount == 3:
                # print('typ', typcount)
                # print(nfcount)
                fatomnum = fatomnums[nfcount]
                base2 = math.ceil(int(fatomnum)/10)
                fatminfo.append(itemlist)
                lcount += 1
                if lcount == base2:
                    nfcount += 1
                    lcount = 0
                    fatminfo = flatten(fatminfo)
                    fatminfos.append(fatminfo)
                    fatminfo = []
                    if  nfcount > nf - 1:
                        typcount += 1
                        lcount = 0
                        continue


            # bda-baa atom section
            if flag == True and itemlist[0] == '/':
                break

            if flag == True and typcount == 4:
                connects.append(itemlist)
                connects = functor(int, connects)
                # datas.append(itemlist)

        print("    {", file=f)
        print("    'name': '" + head.split('/')[-1] + "',", file=f)
        print("    'atom': ", fatomnums, ',', file=f)
        print("    'charge': ", fchgs, ',', file=f)
        print("    'connect_num': ", fbaas, ',', file=f)
        print("    'seg_info': ", fatminfos, ',', file=f)
        print("    'connect': ", connects, ',', file=f)
        end =  """    'nummol_seg': [1],
    'repeat': [1],
    'pair_file': [],
    'multi_xyz': 'none'
    },"""
        print(end, file=f)

    print(']', file=f)
    '''
            seg_conf = {'name': seg_name,
                        'atom': [seg_atom],
                        'charge': [0],
                        'connect_num': [0],
                        'connect': [],
                        'seg_info': [[i + 1 for i in range(seg_atom)]],
                        'nummol_seg': [1],
                        'repeat': [1],
                        'pair_file': [],
                        'multi_xyz': 'none'}
    '''

