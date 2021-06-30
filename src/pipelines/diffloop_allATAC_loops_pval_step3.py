# get pvalue from gwas enrichment shuffling 100 times
# 1/16/2021
# Margaret G

#

import pandas as pd
import os,glob

res_files = glob.glob('result_diffloop_atac/*loops.txt')
# label = ['control']*100+['actual']
diseases = ['anxiety', 'attent', 'autism',
       'bipolar', 'depress', 'ocd', 'panic', 'personality', 'schizo', 'traum', 'alz','diabetes']
diseases_str = '|'.join(diseases)

tissue_dz_pval_dict = dict()
for file in res_files:
    tissue = os.path.basename(file).split('shuffle')[0]
    print(tissue)
    df_orig = pd.read_table(file).fillna(0).drop('File',axis=1)
    df_orig.columns = [x.lower() for x in df_orig.columns]
    df = df_orig.div(df_orig.sum(axis=1), axis=0)
    # subselect columns
    diseases_df = pd.DataFrame()
    diseases_df['dz'] = df.columns
    diseases_df['keep'] = diseases_df.dz.str.contains(diseases_str)
    col_to_keep = diseases_df[diseases_df.keep].dz.to_list()
    df = df[col_to_keep]
    
    dz_to_pval = dict()
    for dz in diseases:
        # print(dz)
        dz_cols = [x for x in col_to_keep if dz in x]
        # print(dz_cols)
        dz_series = df[dz_cols].sum(axis=1)
        pval = pd.Series(dz_series[100]<dz_series[:100]).sum()/100
        dz_to_pval[dz] = pval
    
    tissue_dz_pval_dict[tissue] = dz_to_pval
    print(tissue)
    print(dz_to_pval)

pd.DataFrame.from_dict(tissue_dz_pval_dict)
