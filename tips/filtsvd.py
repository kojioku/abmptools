import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import os

argvs = sys.argv
ifile = argvs[1]
df = pd.read_csv(ifile, header=0, index_col=0)

print(df.head())
# get eigvec1 np.array
eigvec1 = df.iloc[:, 0].values
eigvec2 = df.iloc[:, 1].values
eigvec3 = df.iloc[:, 2].values

x = df.index.values
print(x)

datas = [eigvec1, eigvec2, eigvec3]
names = ['eigvec1', 'eigvec2','eigvec3']

for i in range(len(datas)):
    data = datas[i]
    ohead = os.path.splitext(ifile)[0]
    oname = ohead + '-' + names[i] + '-allbar.png'
    ohist = ohead + '-' + names[i] + '-hist.png'

    # plot and save
    plt.figure(figsize = (18, 10))
    plt.rcParams["font.size"] = 28
    plt.title('raw data ' + names[i], fontsize=32)
    plt.ylim(-0.15, 0.15)
    plt.ylabel("eigval", fontsize=28)
    plt.bar(x, data)
    plt.xticks(x, fontsize=8)
    plt.savefig(oname)
    plt.clf()

    plt.title('histgram ' + names[i], fontsize=32)
    plt.xlim(-0.15, 0.15)
    plt.hist(data, bins=50)
    plt.ylabel("num", fontsize=28)
    plt.savefig(ohist)
    plt.clf()

    print('save', oname)
    print('save', ohist)

# filter
maxfil = 0.09
minfil = -0.09
cols = [0, 1, 2]

for i in range(len(cols)):
    col = cols[i]
    # filter
    filtdf = df[(df.iloc[:,col] > maxfil) | (df.iloc[:,col] < minfil)].iloc[:, col]
    # print(filtdf)

    # save csv
    oname = ohead + '-' + names[i] + '-filter-gt' + str(maxfil) + '-lt' + str(minfil)
    filtdf.to_csv(oname + '.csv')
    print('save', oname + '.csv')

    # save fig
    filtvec = filtdf.values
    filtx = filtdf.index.values

    plt.figure(figsize = (18, 10))
    plt.rcParams["font.size"] = 28
    plt.title('pickup data -gt' + str(maxfil) + '-lt' + str(minfil) + ' ' + names[i], fontsize=32)
    plt.ylim(-0.15, 0.15)
    plt.bar(filtx, filtvec)
    plt.ylabel("eigval", fontsize=28)
    plt.xticks(filtx, rotation=90)
    plt.tight_layout()
    plt.savefig(oname + '.png')
    plt.clf()
    print('save', oname + '.png')


