#!/usr/bin/env python
import pandas as pd
import sys

print ('Guides targeting the same element should have the same element name ex: ATP2B4, ATP2B4, ATP2B4. The pipeline will automatically detect them as the same element')
#df_features = pd.read_excel('5p27sgRNA_guide_metainfo_modified.xlsx')
df_features = pd.read_excel(sys.argv[1])
reconstruct_df = []
for k, v in df_features.groupby('Target_name'):
    v_x = v.copy()
    v_x['Target_name'] = [ f'{k}|{n+1}'  for n in range (v_x.shape[0])]
    v_x['pipeline_id'] = v_x.apply(lambda x : f"{x['Target_name']}_sgrna_{x['chr']}:{x['start']}:{x['end']}" ,axis=1)

    reconstruct_df.append(v_x)
df_features_names_changed = pd.concat(reconstruct_df)
df_features_names_changed[[ 'sgRNA_sequences', 'pipeline_id']].to_csv('guide_features.txt', sep='\t', header=None, index=None)

