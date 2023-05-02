#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import os
import sys


#loading picke file
PICKE_PATH = 'perturbdata.pkl'
DIRECTION = sys.argv[2] #left, right, both
with open(PICKE_PATH, 'rb') as f:
    p = pickle.load(f)




def run_sceprte_by_element(ELEMENT, name_element='positive_control'):
  path_element = f'{ELEMENT}'
  os.system(f'mkdir {path_element}')
  p.genes_to_test_count(ELEMENT) # Genes on a 1MB region  
  p.genes_to_test_count(ELEMENT).T.to_csv(f'{path_element}/gene_exp_one_gene_guide.txt', sep=',')
  df_guides_sceptre = pd.DataFrame([   p.capture_guide(x) for x in p.capture_element_guides(ELEMENT) ])
  df_guides_sceptre.columns  = p.genes_to_test_count(ELEMENT).T.columns.values.tolist()
  df_guides_sceptre.index = p.capture_element_guides(ELEMENT)
  print ('Guide number:', df_guides_sceptre.sum())
  print ('Guide number:', df_guides_sceptre.sum(1))
  dict_guides_to_analyze_in_the_limit = set(df_guides_sceptre.sum(1)[df_guides_sceptre.sum(1) > 30].index.values)
  print('higher than  > guides 30 ', dict_guides_to_analyze_in_the_limit )


  df_guides_sceptre.to_csv(f'{path_element}/guides_one_gene_guide.txt', sep=',')
  df_cov_mod = p.capture_covariates().copy()  # to filter for single factors (case not batch)
  print ('set', set(df_cov_mod['bath_number'].values), len(set(df_cov_mod['bath_number'].values)) == 1 )
  if len(set(df_cov_mod['bath_number'].values)) == 1:
    del df_cov_mod['bath_number']
  df_cov_mod.to_csv(f'{path_element}/covariates.txt', sep=',')


  df_pairs = pd.DataFrame(np.array([ [ [g,e,name_element] for g in p.genes_to_test_count(ELEMENT).columns ] for e in p.capture_element_guides(ELEMENT)   if e in dict_guides_to_analyze_in_the_limit ]  ).reshape(-1,3)     )
  df_pairs.columns  =  'gene_id gRNA_group pair_type'.split(' ')
  df_pairs.to_csv(f'{path_element}/pairs.txt', sep=',', index=False)




  r=f'''
      setwd("{path_element}")
      print(getwd())
      library(sceptre)
      library(tibble)
      library(dplyr)
      exp =  as.matrix(read.table( 'gene_exp_one_gene_guide.txt' , sep=',' ,header=TRUE, row.names=1, check.names=FALSE) )
      #print( tail(colnames(exp), 30) )
      #print( tail(row.names(exp), 30) )

      g_I = as.matrix( read.table('guides_one_gene_guide.txt',     sep=',' , header=TRUE, row.names=1 , check.names=FALSE))
      #g_I[,] = g_I == 1
      #print( tail(colnames(g_I), 30) )
      #print( tail(row.names(g_I), 30) )



      cov = read.table('covariates.txt',                sep=',', header=TRUE, row.names=1)
      if ('bath_number' %in% names(cov) ) {{
        cov$bath_number =  as.factor(cov$bath_number)
      }}



      #print(tail( row.names(cov)  , 30))
      pairs_test = read.table('pairs.txt',              sep=',', header=TRUE , check.names=FALSE)
      pairs_test = as_tibble(pairs_test)
      pairs_test$pair_type = as.factor(pairs_test$pair_type)


      result <- run_sceptre_high_moi(gene_matrix = exp, 
                                            combined_perturbation_matrix = g_I, 
                                            covariate_matrix = cov,
                                            gene_gRNA_group_pairs = pairs_test,

                                            side='{DIRECTION}',
                                            )


      result$p_value
      result$z_value

      write.table(result ,"results.txt" )

      '''

  file_save = open('r_script.r', 'w')
  file_save.write(r)
  file_save.close()
  os.system('Rscript r_script.r')



#p.extract_element_per_type('PUTATIVE_ENHANCER')
elements_sceptre = p.show_elements()
for elements in tqdm(elements_sceptre):
  print ('-'*20, elements, '-'*20)
  print ('Element')
  if not os.path.exists(f'{elements}'):
    print ('-'*20, 'run', elements , '-'*20)
    run_sceprte_by_element(elements,
                           name_element='all_elements_test')
  
  else:
      print ('-'*20, 'not run ', elements, '-'*20)

    
    
    
