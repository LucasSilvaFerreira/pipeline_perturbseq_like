#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
import pandas as pd
from GTFProcessing import GTFProcessing
from tqdm import tqdm
import pickle
import sys
import muon as mu
import mudata as md

import argparse

parser = argparse.ArgumentParser(description='Description of your script')

parser.add_argument('-m', '--muon_data', type=str, required=True, help='muon_data with guide and scrna')
parser.add_argument('-f', '--gtf_in', type=str, required=True, help='Description of the third argument')
parser.add_argument('-t', '--in_trans', type=str, required=True, default="FALSE", help='Description of the fourth argument (default: FALSE)')
parser.add_argument('-d', '--distance_from_guide', type=int, default=1000000, help='Description of the fifth argument (default: 1000000)')
parser.add_argument('-g', '--add_gene_names', type=str, required=False, help='gene names separated by comma')


args = parser.parse_args()
MUON_DATA = args.muon_data
GTF_IN = args.gtf_in
IN_TRANS = False if args.in_trans.upper() == 'FALSE' else True
print (f'is in trans? =  {IN_TRANS} given parameter {args.in_trans } ')
DISTANCE_FROM_GUIDE = int(args.distance_from_guide)
ADD_GENE_NAMES = [g for g in args.add_gene_names.split(',')]



NUMBER_FOR_RANDOM_CONTROLS = 10




# IN_GUIDE = '/n/scratch3/users/l/lf114/guillaume_perturb_data/testing_w_option_4/ff/847841e216c2cd64fba0e44c53af2c/results_per_lane/processed_anndata_guides_data.h5ad'
# IN_EXP =   '/n/scratch3/users/l/lf114/guillaume_perturb_data/testing_w_option_4/4b/9c2870f9776840fc09e8de5384ce40/results_per_lane/processed_anndata_transcripts_data.h5ad' 
# GTF_IN =   '/n/scratch3/users/l/lf114/guillaume_perturb_data/testing_w_option_4/cd/d25097d237358a3a48282ece7811da/transcripts.gtf'
# DISTANCE_FROM_GUIDE = 1000000
# NUMBER_FOR_RANDOM_CONTROLS = 10


muon_structure = md.read(MUON_DATA)


ann_exp   =   muon_structure['scRNA'].copy()
ann_guide =   muon_structure['guides'].copy()

print (ann_guide.var.index.values)

df_element_guide = pd.DataFrame([[x.split('|')[0] for x in ann_guide.var.index.values], [x.split('|')[0] for x in ann_guide.var.index.values]]).T
df_element_guide = df_element_guide.drop_duplicates()
df_element_guide.columns = ['Element', 1]
df_element_guide = df_element_guide.set_index('Element')
df_element_guide.values.tolist()

del df_element_guide[1]

def check_guide(guide_group, guide_name):
  if guide_name.split('|')[0]  ==  guide_group:
    return 1
  else:
    return 0

df_to_elements = pd.DataFrame([ [check_guide (guide_group,guide_name ) for guide_name in   ann_guide.var.index.values ]      for guide_group in df_element_guide.index    ])

df_to_elements.index = df_element_guide.index

df_to_elements.columns = ann_guide.var.index.values
# df_to_elements.index.values.tolist()

# df_to_elements

ann_Element_guide =  anndata.AnnData(X=df_to_elements)
#ann_Element_guide.var
dict_capure_guide_coords = {k.split('|')[0] : k.split('_')[-1] for k in df_to_elements.columns}



gtf = GTFProcessing(GTF_IN)
df_gtf_refseq = gtf.get_gtf_df()


df_gtf_refseq_gene =df_gtf_refseq[ df_gtf_refseq['feature'] == 'gene']
set_gene_on_annexp = set(ann_exp.var.index.values)
df_gtf_refseq_gene = df_gtf_refseq_gene[df_gtf_refseq_gene['gene_name'].apply(lambda x : x in set_gene_on_annexp)]

df_gtf_refseq_gene['gene_real_start'] = df_gtf_refseq_gene.apply(lambda x : x['start'] if x['strand'] == "+" else x['end'], axis=1 )


def capture_genes(gene_query,   DISTANCE = 1000000, number_for_random_genes=10, in_trans=False):
  # case coord use a list
#   if 'pos_control_' in gene_query:
#     return set()
#   if 'random' in  gene_query:
#     return set(df_gtf_refseq_gene.sample(number_for_random_genes)['gene_name'].values.tolist())

#   if 'scrambled' in  gene_query:
#     return set(df_gtf_refseq_gene.sample(number_for_random_genes)['gene_name'].values.tolist())
  if in_trans:
    return set(df_gtf_refseq_gene['gene_name'].values.tolist())

  # coord_bool=False
  # if ':' in gene_query:
    #print (gene_query, 'gquery')
    #print (gene_query.split(':')[0].replace('chr', ''))
  coord_bool=True



  #   if 'scrambled' in  gene_query:
#     return set(df_gtf_refseq_gene.sample(number_for_random_genes)['gene_name'].values.tolist())
  if coord_bool:
      gene_query_use = dict_capure_guide_coords[gene_query]
      print (gene_query_use, len(gene_query_use))
      if ':' in gene_query_use:
        chr =  gene_query_use.split(':')[0].replace('chr', '')
        coord =  int(gene_query_use.split(':')[1].split(':')[0])
        
      else:
        chr,  coord = gene_query
      #print (chr, coord)

  # else:
  if len(df_gtf_refseq_gene.query(f'gene_name == "{gene_query}" ')) == 0:
    print (gene_query,'len == 0')
    return set(df_gtf_refseq_gene.sample(number_for_random_genes)['gene_name'].values.tolist())

  
  chr,  coord = df_gtf_refseq_gene.query(f'gene_name == "{gene_query}" ')[['chr', 'gene_real_start']].iloc[0].values.tolist()

  #print (df_gtf_refseq_gene.query(f'chr == "{chr}"  '))
  #should accet raw coordinates as well (to capture enhancer)
  return set(df_gtf_refseq_gene.query(f'chr == "{chr}"  ')[np.abs(df_gtf_refseq_gene.query(f'chr == "{chr}"  ')['gene_real_start'] - coord)< DISTANCE]['gene_name'].values)





def is_gene_included(gene):
  
  if '_TSS' in gene:
    gene = gene.replace('_TSS', '')
  print (gene)
  gene_set = capture_genes(gene, DISTANCE = int(DISTANCE_FROM_GUIDE),
                           number_for_random_genes = NUMBER_FOR_RANDOM_CONTROLS,
                           in_trans=IN_TRANS)
  #print(gene, len(gene_set))
  return gene_set



final_element_binary = []
for element in  tqdm(ann_Element_guide.obs.index):
  print (element)    
  set_of_elements = is_gene_included(element)
  set_of_elements = set_of_elements.union(set(ADD_GENE_NAMES))
  print (set_of_elements, 'set_of_elements')
  presence_gene_e = []
  for e in  ann_exp.var.index.values:
    #print (e)
    if e in set_of_elements:
      presence_gene_e.append(1)
    else:
      presence_gene_e.append(0)
  final_element_binary.append(presence_gene_e)




Elements_x_tested_genes =    pd.DataFrame(final_element_binary, columns = ann_exp.var.index.values, index=ann_Element_guide.obs.index)
ann_Element_x_tested_genes=  anndata.AnnData(X=Elements_x_tested_genes)


def add_test_type(e):
  if '_TSS' in e:
    return 'POSITIVE_CONTROL'
  if 'random' in e or 'scrambled' in e:
    return 'NEGATIVE_CONTROL'
  if 'chr' in e:
    return 'PUTATIVE_ENHANCER'


ann_Element_x_tested_genes.obs['GUIDE_TYPE'] = [add_test_type(e) for e in Elements_x_tested_genes.index.values]

from Perturb_Loader import PERTURB_MANIPULATE
p = PERTURB_MANIPULATE(ann_exp, ann_guide, ann_Element_guide, ann_Element_x_tested_genes)
filename = 'perturbdata.pkl'
outfile = open(filename,'wb')
pickle.dump(p,outfile)
outfile.close()
