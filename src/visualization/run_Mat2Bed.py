"""

run_Mat2Bed.py

Margaret Guo
 inputs a *_abs.bed file and a *matrix file and creats a *matrix.bed file

==> H9D28-B1_5000_abs.bed <==
1	0	5000	1
1	5000	10000	2
1	10000	15000	3
1	15000	20000	4
1	20000	25000	5
1	25000	30000	6
1	30000	35000	7
1	35000	40000	8
1	40000	45000	9
1	45000	50000	10

==> H9D28-B1_5000.matrix <==
3	50010	1
3	509665	1
3	553007	1
4	188	1
4	101332	1
12	305683	1
12	338344	1
13	586882	1
17	278950	1
21	98435	1

==> H9D28-B1_5000.matrix.bed <==
1	10000	15000	2	790000	795000	1
1	10000	15000	17	48115000	48120000	1
1	10000	15000	20	46415000	46420000	1
1	15000	20000	1	935000	940000	1
1	15000	20000	3	14200000	14205000	1
1	55000	60000	8	135595000	135600000	1
1	55000	60000	10	11320000	11325000	1
1	60000	65000	X	53325000	53330000	1
1	80000	85000	8	1930000	1935000	1
1	100000	105000	2	242915000	242920000	1


"""
import pandas as pd
import argparse

import os, glob


def main(**kwargs):
    """
    main wrapper function

    :param kwargs: contains command line arguments
    :return: None
    """
    abs_file, matrix_file = kwargs['abs_file'], kwargs['matrix_file']

    # check tissue matirx
    tissue_abs = abs_file.split("_abs")[0]
    tissue_mat = matrix_file.split('.matrix')[0]

    assert (tissue_abs==tissue_mat)

    output_file = matrix_file+'.bed'
    print('output_file',output_file)

    abs_df = pd.read_csv(abs_file, sep='\t', names = ['chr','start','stop','id'])
    mat_df = pd.read_csv(matrix_file, sep='\t', names = ['l_id','r_id','count'])

    mat_df = mat_df.merge(abs_df, how='left',left_on='l_id',right_on='id')
    mat_df = mat_df.merge(abs_df, how='left',left_on='r_id',right_on='id', suffixes=('_l','_r'))
    mat_df = mat_df[['chr_l','start_l','stop_l','chr_r','start_r','stop_r','count']]
    mat_df.to_csv(output_file, sep='\t',header=None,index=None)
    print('saved' . mat_df.shape)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('abs_file', help='filepath to file with anchor annotations: i.e. H9D28-B1_5000_abs.bed')
    parser.add_argument('matrix_file', help='filepath to file with loop counts based on anchor indices: i.e. H9D28-B1_5000.matrix')
    args = vars(parser.parse_args())

    main(**args)

