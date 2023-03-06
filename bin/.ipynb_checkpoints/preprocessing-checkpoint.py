#!/usr/bin/env python

import pandas as pd
import sys

def preparing_files_initialing(path_list):
    
    print ('preparing files...')
    # create_main_fig_dir = 'results_per_lane'
    # create_case_doesnt_exist(create_main_fig_dir)
    # #capturing files generated (guides and ScRNAseq)
    # #how to capture the file names?
    # assert len(set([f.parent.parent.name for f in path_objects_test])) == len(path_objects_test), 'This seems like multiple runs inside the same directory, use resume to avoid it or change the directory'
    dir_fastqz = [f + '/counts_unfiltered/adata.h5ad' for f in path_list]
    df_files = pd.DataFrame(dir_fastqz, columns=['file_path'])
    df_files['sample'] = df_files['file_path'].apply(lambda x : 'Guide' if 'guide' in x else 'scRNA' )
    df_files['lane'] =   df_files['file_path'].apply(lambda x : x.split('/')[-3].split('_')[1].replace('L', '') )  # check if how it will be handled in the future
    
    print (len(dir_fastqz))
    print (df_files.shape)
    print (df_files)
    df_files.to_csv('initial_preprocessing_file_names.txt', sep='\t', index=None)
    #return df_files

    


path_in_files = sys.argv[1:]
print (sys.argv)
preparing_files_initialing(path_in_files)
