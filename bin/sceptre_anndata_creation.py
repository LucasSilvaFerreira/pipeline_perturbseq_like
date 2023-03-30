#!/usr/bin/env python

import numpy as np
import anndata
from statsmodels.stats.multitest import multipletests 
import glob
import sys
import pandas as pd

print (sys.argv[0], 'sys.argv[0]')
print (sys.argv[1], 'sys.argv[1]')
glob_path = f'{sys.argv[1]}/**/results.txt'.replace('//','/')

print (glob_path, 'glob_path')

all_results = glob.glob(glob_path, recursive=True)
print (all_results)
#all_results = !find  '/n/scratch3/users/l/lf114/gasperini_jamboree/gasperini_08' | grep results.txt 
df_all_tests = pd.concat([pd.read_csv(f, sep=' ') for f in all_results ], ignore_index=True)
df_all_tests.head()
df_all_tests['adj_pvalue'] =  multipletests(df_all_tests['p_value'],method='fdr_bh')[1]


set_genes = set(df_all_tests['gene_id'])
set_guide = set(df_all_tests['gRNA_id'])


df_in = pd.DataFrame(np.array([np.nan for x in range (len(set_guide) * len(set_genes))]).reshape(len(set_guide) ,len(set_genes)))
df_in.index = set_guide
df_in.columns = set_genes
df_in.copy()

def populate_df(df_in_mod, var_to_pop):
    for gene_id, grna, value in df_all_tests[['gene_id', 'gRNA_id', var_to_pop]].values:
        df_in_mod.loc[grna][gene_id] = value
    return df_in_mod

df_pvalue = populate_df(df_in.copy(), 'p_value')
adata = anndata.AnnData(df_pvalue)

adata.layers["adj_pvalue"] = populate_df(df_in.copy(), 'adj_pvalue').values
adata.layers["z_value"] = populate_df(df_in.copy(), 'z_value').values
adata.layers["log_fold_change"] = populate_df(df_in.copy(), 'log_fold_change').values
adata.layers["significant"] = adata.layers["adj_pvalue"] < 0.01


adata.write('sceptre_results_ann_data.h5ad')
