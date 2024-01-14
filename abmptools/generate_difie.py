import abmptools.cpfmanager
import argparse
import pandas as pd
import os
import datetime
import time
from multiprocessing import Pool
import sys


def get_args():
    parser = argparse.ArgumentParser(description="generate DIFIE")

    # Input CPF name, required
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Input CPF name (e.g., file-xxx-bbb.cpf)")

    # Start, End, Interval, required
    parser.add_argument("-t", "--time",
                        required=True,
                        nargs=3,
                        type=int,
                        help="Start, End, Interval numbers")

    # Zero padding, default is "1"
    parser.add_argument("-z",
                        "--zero-padding",
                        default=1,
                        type=int,
                        help="Number of digits for zero-padding. Default is '1' (no padding).")

    # Representative structure, default is first structure (Start number)
    parser.add_argument("-s",
                        "--structure",
                        type=int,
                        default=0,
                        help="Number of the representative structure for average CPF. Default is the start number.")

    # Target residues, default is removing WAT, HOH at the end of the file
    parser.add_argument("-f",
                        "--fragments",
                        default=0,
                        help="Specify target residues. Default is None")

    # Version, default is "23"
    parser.add_argument("-v",
                        "--version",
                        default=23,
                        help="Specify the output version. Currently, only '23' is available.")

    parser.add_argument("-np",
                        "--np",
                        default=1,
                        type=int,
                        help="Number of processes. Default is '1' (no multiprocessing).")

    args = parser.parse_args()

    print("Arguments received:")
    print("Input CPF:", args.input)
    print("Time:", args.time)
    print("Zero-padding:", args.zero_padding)
    print("Representative structure:", args.structure)
    print("Target fragments:", args.fragments)
    print("Version:", args.version)

    return args


def getcpfobj(intime):
    cpf = abmptools.CPFManager()
    padded = str(intime).zfill(args.zero_padding)
    input_cpf = args.input.replace('xxx', padded)
    print(input_cpf)
    cpf.tgtfrag = cpf.selectfrag(args.fragments)
    cpf = cpf.parse(input_cpf)

    return cpf


def setoutcpfstatic(tgtnum, args, cpfs):
    outcpf = abmptools.CPFManager()
    outcpf.cpfver = args.version
    outcpf.labels = cpfs[tgtnum].labels
    outcpf.atominfo = cpfs[tgtnum].atominfo
    outcpf.fraginfo = cpfs[tgtnum].fraginfo
    outcpf.condition = cpfs[tgtnum].condition
    outcpf.static_data = cpfs[tgtnum].static_data
    outcpf.mominfo = cpfs[tgtnum].mominfo
    # outcpf.diminfo = cpfs[tgtnum].diminfo
    return outcpf


def getavestddf(cpfs):
    # set averaged and stdev diminfo and write
    atominfo = []
    diminfos = []
    for i in range(len(cpfs)):
        atominfo.append(cpfs[i].atominfo)
        diminfos.append(cpfs[i].diminfo)
        # print(cpfs[i].diminfo.head())
    all_atomdfs = pd.concat(atominfo)
    all_dimdfs = pd.concat(diminfos)
    print(all_dimdfs)

    all_chgdfs = all_atomdfs[['alabels'] + outcpf.labels['charge']]
    all_chgdfs.columns = ['M-' + col if col != 'alabels' else col for col in all_chgdfs.columns]
    print(all_chgdfs)

    average_atomdf = all_chgdfs.groupby(['alabels']).mean().reset_index()
    stddev_atomdf = all_chgdfs.groupby(['alabels']).std().reset_index()
    # 'key' カラムを使って df1 と df2 を結合します。
    # result = pd.merge(df1, df2, on='key', how='inner')
    average_dimdf = all_dimdfs.groupby(['fragi', 'fragj']).mean().reset_index()
    stddev_dimdf = all_dimdfs.groupby(['fragi', 'fragj']).std().reset_index()

    print('average_atomdf\n', average_atomdf.head())
    print('stddev_atomdf\n', stddev_atomdf.head())
    sys.exit()
    # print('average_df\n', average_df.head())
    # print('stddev_df\n', stddev_df.head())
    return average_atomdf, stddev_atomdf, average_dimdf, stddev_dimdf


if __name__ == "__main__":
    start_time = time.time()
    args = get_args()
    cpfs = []
    intimes = []
    d_today = str(datetime.datetime.today())
    tgtnum = args.structure
    '''
    個別に残すもの
    - labels
    - atominfo
    - fraginfo
    - condition
    - static_data
    - mominfo
    - fragment dist (diminfo)
    平均、標準偏差を計算するもの
    - 電荷(atominfo)
    - diminfo
    ラベルの追加
    - 電荷カラム(atominfo)
    - IFIEカラム(diminfo)
    '''

    # --- parse (multiprocessing) ---
    with Pool(processes=args.np) as p:
        cpfs = p.map(getcpfobj,
                     [i for i in range(args.time[0],
                                       args.time[1]+1,
                                       args.time[2])])

    # --- set output cpf ---
    difiename = os.path.basename(args.input).split('.')[0] + '-DIFIE.cpf'
    stdevname = os.path.basename(args.input).split('.')[0] + '-STDEV.cpf'

    # get static data for DIFIE cpf
    outcpf = setoutcpfstatic(tgtnum, args, cpfs)

    # get averaged and stdev data for DIFIE cpf
    average_atomdf, stddev_atomdf, average_dimdf, stddev_dimdf \
        = getavestddf(cpfs)

    # outcpf.labels のすべての要素の頭にM-をつける
    charge = ['M-' + i for i in outcpf.labels['charge']]
    charge += ['S-' + i for i in outcpf.labels['charge']]
    dimer = ['M-' + i for i in outcpf.labels['dimer']]
    dimer += ['S-' + i for i in outcpf.labels['dimer']]

    outcpf.labels['charge'] = charge
    outcpf.labels['dimer'] = dimer

    print(outcpf.labels)

    print(average_chgdf)
    sys.exit()

#         labels = {
#             'charge': charge_label,
#             'DPM': DPM_label,
#             'monomer': monomer_label,
#             'dimer': dimer_label,
#             }

    outcpf.diminfo = average_dimdf
    outcpf.write('DIFIE (Generated by ABMPTools ' + d_today + ')', difiename)

    outcpf.diminfo = stddev_dimdf
    outcpf.write('DIFIE-STDEV (Generated by ABMPTools '
                 + d_today + ')', stdevname)

    end_time = time.time()
    print('Elapsed time: ', end_time - start_time, 'seconds')
    print('DIFIE and STDEV CPFs are generated.')

    # --- tips ---
    # cpf = abmptools.CPFManager()
    # cpf.parse('CS4_ligX_md1_12ns_mod_mp2_631gd.cpf')
    # cpf.write('test-from10to23', 'test-i10o23.cpf')

    # print(dir(cpf))
    # print(vars(cpf))

    # print('cpfver\n', cpf.cpfver)
    # print('labels\n', cpf.labels)
    # print('atominfo\n', cpf.atominfo.head())
    # print('fraginfo\n', cpf.fraginfo)
    # print('condition\n', cpf.condition)
    # print('static_data\n', cpf.static_data)
    # print('mominfo\n', cpf.mominfo)
    # print('diminfo\n', cpf.diminfo)
