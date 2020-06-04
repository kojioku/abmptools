import numpy as np
import sys
import pandas as pd
import os

if __name__ == "__main__":
    argvs = sys.argv
    fname = argvs[1]

    df = pd.read_csv(fname, header=0, index_col=0, float_precision="high")

    print(df.head())
    # print(df.iloc[:, 1:])

    array = df.as_matrix()
    # array = df.iloc[:, 1:].as_matrix()

    # print(array[0])
    # print(type(array))

    # print(df.columns)
    # print(type(df.columns))

    print('index', len(df.index), 'columns', len(df.columns))
    if len(df.index) <= len(df.columns):
        s_index = df.index
    else:
        s_index = df.columns

    U, s, V = np.linalg.svd(array, full_matrices=True)
#
    df_U = pd.DataFrame(U, columns = df.index, index=df.index)
    df_s = pd.DataFrame(s, index = s_index)
    df_V = pd.DataFrame(V, columns = df.columns, index=df.columns)
# #
#     print(df_U)
#     print(df_s)
#     print(df_V)
# #
    otail = os.path.splitext(fname)[0].split('/')[-1]

    print(otail)
    svd_uname = "svd_U_" + otail + ".csv"
    svd_sname = "svd_s_" + otail + ".csv"
    svd_vname = "svd_VT_" + otail + ".csv"

    df_U.to_csv(svd_uname)
    df_s.to_csv(svd_sname)
    df_V.T.to_csv(svd_vname)

    print(svd_uname, svd_sname, svd_vname, 'are created.')
