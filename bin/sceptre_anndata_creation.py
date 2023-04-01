#!/usr/bin/env python

import numpy as np
import anndata
from statsmodels.stats.multitest import multipletests 
import glob
import sys
import pandas as pd
import mudata as md

import numpy as np
from scipy.stats import combine_pvalues



print (sys.argv[1], 'sys.argv[1]')

glob_path = f'{sys.argv[1]}/**/results.txt'.replace('//','/')

print (glob_path, 'glob_path')

all_results = glob.glob(glob_path, recursive=True)
print (all_results)
#all_results = !find  '/n/scratch3/users/l/lf114/gasperini_jamboree/gasperini_08' | grep results.txt 
df_all_tests = pd.concat([pd.read_csv(f, sep=' ') for f in all_results ], ignore_index=True)
df_all_tests.head()
df_all_tests['adj_pvalue'] =  multipletests(df_all_tests['p_value'],method='fdr_bh')[1]
print (df_all_tests.head())

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


 
mudata_loaded = md.read(sys.argv[2])
mudata_loaded['scRNA'].var.index = mudata_loaded['scRNA'].var['feature_name']
scRNA_use = mudata_loaded['scRNA'].var.drop_duplicates("feature_name")


ann_value = adata
gene_vars = ann_value.var.index.values.tolist()
check_intersection = [x  for x in ann_value.var.index.values.tolist() if x in scRNA_use.index] 
genes_position_info = scRNA_use.loc[check_intersection]
genes_position_info
bool_fiter = [x in genes_position_info.index for x in ann_value.var.index]
ann_data_filtered_genes = ann_value[:, bool_fiter]
guidesann =  mudata_loaded['guides']
check_intersection_guides = [x  for x in ann_value.obs.index.values.tolist() if x in guidesann.var.index] 
ann_value.obs = mudata_loaded['guides'].var.loc[ann_data_filtered_genes.obs.index]
ann_value.var =  mudata_loaded['scRNA'].var.loc[ann_value.var.index]



def pvalue_aggregation(x):
    if x:
        return combine_pvalues(x).pvalue
    else:
        return np.nan

create_target_elements = ann_value
element_p_reconstruction = []


for e in create_target_elements.obs['target_elements']:
    bool_op = create_target_elements.obs['target_elements'] == e 
    recreate_guide_coord = create_target_elements.obs[create_target_elements.obs['target_elements'] == e].iloc[0]
    x_value =  create_target_elements.X[bool_op]
    p_combined = pd.DataFrame(x_value).apply(lambda x : pvalue_aggregation([i  for i in x  if np.isnan(i) == False]) , 0)
    
    element_p_reconstruction.append([e, p_combined, recreate_guide_coord])
    
index_name = [x[0] for x  in element_p_reconstruction]
df_comb_to_ann = pd.DataFrame([x[1] for x  in element_p_reconstruction])

df_comb_to_ann.index = index_name
df_comb_to_ann.columns = create_target_elements.var.index.values
df_comb_to_ann = df_comb_to_ann.drop_duplicates()



ann_combined = anndata.AnnData(df_comb_to_ann)
df_selected_fields = create_target_elements.obs.drop_duplicates('target_elements')
df_selected_fields = df_selected_fields[['target_elements', 'guide_chr', 'guide_start', 'guide_end']]
df_selected_fields.index = df_selected_fields['target_elements'].tolist()
ann_combined.obs = df_selected_fields
ann_combined.obs.columns = 'element element_chr element_start element_end'.split(' ')

ann_combined.var = ann_value.var

# print (ann_combined.obs)
# print (ann_combined.var)
# print (ann_combined)

ann_combined.layers['sig_not_adj'] = ann_combined.X < 0.05


#ann_value.write('sceptre_results_guide_ann_data.h5ad')
#ann_combined.write('sceptre_results_element_ann_data.h5ad')
mdata_out =  md.MuData({"guides":  mudata_loaded['guides'], "scRNA": mudata_loaded['scRNA'] ,'result_guides': ann_value , 'result_elements': ann_combined })
mdata_out.write("mudata_results.h5mu")

