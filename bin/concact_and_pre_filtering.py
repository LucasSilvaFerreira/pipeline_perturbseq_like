#!/usr/bin/env python

import numpy as np
import scanpy as sc
import anndata
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt
import scrublet as scr
import re
from tqdm import tqdm
import os
import pandas as pd
from glob import glob
import pathlib
import gc
import argparse
import seaborn as sns
from anndata import AnnData, read_h5ad



def create_case_doesnt_exist(x):
    if not os.path.exists(x):
       os.makedirs(x)



def analyze_batch(scRNA_ann_FILE,
                  guide_ann_FILE,
                  BATCH_NUMBER,
                  OUT_DIR_FIGURES,
                  EXPECTED_NUMBER_OF_CELLS,
                  SPECIE_MITO,
                  MITO_PERCENTAGE_ALLOWED,
                  UMI_CELL_THRESHOLD
                  ):
    #!mkdir $OUT_DIR_FIGURES
    create_case_doesnt_exist(OUT_DIR_FIGURES)
    # Processing scRNA_seq
    print (scRNA_ann_FILE)
    adata = anndata.read(scRNA_ann_FILE)
    print (adata)
    # Perform SVD
    tsvd = TruncatedSVD(n_components=2)
    tsvd.fit(adata.X)
    X = tsvd.transform(adata.X)
    # Plot the cells in the 2D PCA projection
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(X[:,0], X[:,1], alpha=0.5, c="green")
    #plt.axis('off')
    plt.savefig(OUT_DIR_FIGURES+'/'+ 'svd_batch' + '.png')
    plt.show()
    # Create a plot showing genes detected as a function of UMI counts.
    fig, ax = plt.subplots(figsize=(10, 7))
    x = np.asarray(adata.X.sum(axis=1))[:,0]
    y = np.asarray(np.sum(adata.X>0, axis=1))[:,0]
    ax.scatter(x, y, color="green", alpha=0.25)
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Genes Detected")
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_yscale('log', nonpositive='clip')

    ax.set_xlim((0.5, 4500))
    ax.set_ylim((0.5,2000))
    plt.savefig(OUT_DIR_FIGURES+'/'+ 'saturation_batch' + '.png')
    plt.show()
    #@title Threshold cells according to knee plot { run: "auto", vertical-output: true }
    cutoff = UMI_CELL_THRESHOLD #@param {type:"number" DEFAULT 100}
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    cell_set = np.arange(len(knee))
    num_cells = cell_set[knee > cutoff][::-1][0]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(knee, cell_set, linewidth=5, color="g")
    ax.axvline(x=cutoff, linewidth=3, color="k")
    print (cutoff)
    print (num_cells)
    ax.axhline(y=num_cells, linewidth=3, color="k")
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")
    plt.grid(True, which="both")
    plt.savefig(OUT_DIR_FIGURES+'/'+ 'knee_plot_batch' + '.png')
    plt.show()
    print(f"{num_cells:,.0f} cells passed the {cutoff} UMI threshold")

    sc.pp.filter_cells(adata, min_genes=cutoff)
    sc.pp.filter_cells(adata, min_counts=knee[EXPECTED_NUMBER_OF_CELLS])
    adata.var.index = [x.split('.')[0] for x in adata.var.index]  #removing the dot from ensembl in the ann data
    mito_ensembl_ids = sc.queries.mitochondrial_genes(SPECIE_MITO, attrname="ensembl_gene_id")
    mito_genes = mito_ensembl_ids["ensembl_gene_id"].values
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    sc.pl.scatter(adata, x='n_counts', y='percent_mito')
    plt.savefig(OUT_DIR_FIGURES+'/'+ 'mito_scatter_batch' + '.png')

    #We processed the cDNA UMI count
    #matrix and retained cells with less than 20% mitochondrial reads and at least 850 unique gene
    #UMIs 

    adata = adata[adata.obs.percent_mito < MITO_PERCENTAGE_ALLOWED]  
    sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True)
    plt.savefig(OUT_DIR_FIGURES+'/'+ 'box_plot_batch' + '.png')
    sc.pl.highest_expr_genes(adata, n_top=20, gene_symbols='feature_name')
    plt.savefig(OUT_DIR_FIGURES+'/'+ 'top_genes_batch' + '.png')

    scrub = scr.Scrublet(adata.X)
    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
    scrub.plot_histogram()
    plt.savefig(OUT_DIR_FIGURES+'/'+ 'doublets_batch' + '.png')
    sum(adata.obs['predicted_doublets'])


    adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
    adata = adata[adata.obs['predicted_doublets'] == False] # removing doublets

    # PROCESSING THE GUIDES

    adata_guide = anndata.read(guide_ann_FILE)
    adata_guide.obs['number_of_nonzero_guides'] = [ x.tolist()[0][0] for x in np.sum(adata_guide.X>0, axis=1)]

    sc.pl.violin(adata_guide, ['number_of_nonzero_guides'], jitter=0.4, multi_panel=True)
    plt.savefig(OUT_DIR_FIGURES+'/'+ 'guide_non_zero_batch' +'.png')

    #getting intersection

    cell_bar_intersection = set(adata_guide.obs.index.values).intersection(adata.obs.index.values)


    adata = adata[[ c in cell_bar_intersection for c in adata.obs.index ]] #filtering not shared barcodes (guide and scrna)

    adata_guide = adata_guide[[ c in cell_bar_intersection for c in adata_guide.obs.index ]] #filtering not shared barcodes (guide and scrna)

    adata.obs['batch_number'] = BATCH_NUMBER
    adata_guide.obs['batch_number'] = BATCH_NUMBER
    return adata,adata_guide



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
    
    # print (len(dir_fastqz))
    # print (df_files.shape)
    # print (df_files)
    return df_files



def preparing_files(path):
    create_main_fig_dir = 'results_per_lane'
    create_case_doesnt_exist(create_main_fig_dir)
    #capturing files generated (guides and ScRNAseq)
    #how to capture the file names?
    df_files = pd.read_csv(path, sep='\t')
    # path_objects_test = [f for f in pathlib.Path(path).rglob("*") if 'h5ad' in f.name]
    # assert len(set([f.parent.parent.name for f in path_objects_test])) == len(path_objects_test), 'This seems like multiple runs inside the same directory, use resume to avoid it or change the directory'
    # dir_fastqz = [str(f) for f in path_objects_test]
    # df_files = pd.DataFrame(dir_fastqz, columns=['file_path'])
    # df_files['sample'] = df_files['file_path'].apply(lambda x : 'Guide' if 'guide' in x else 'scRNA' )
    # df_files['lane'] =   df_files['file_path'].apply(lambda x : x.split('/')[-3].split('_')[1].replace('L', '') )  # check if how it will be handled in the future
    return df_files
    

    
    
    
def execute_analysis(df_files, EXPECTED_CELL_NUMBER, MITO_SPECIE, MITO_EXPECTED_PERCENTAGE, UMI_CELL_THRESHOLD ):
    
    ann_scrna_data_to_concat = []
    ann_guide_data_to_concat = []
    
    
    for k, v in df_files.groupby('lane'):
        lane_n = v.query('sample== "Guide" ')['lane'].values[0]
        file_guide = v.query('sample== "Guide" ')['file_path'].values[0]
        file_scRNA = v.query('sample== "scRNA" ')['file_path'].values[0]
        out_name = f'results_per_lane/lane_{lane_n}'
        exp_number_of_cells = EXPECTED_CELL_NUMBER
        specie_mito = MITO_SPECIE
        mito_percent_allowed = MITO_EXPECTED_PERCENTAGE
        tresh_umi_cell = UMI_CELL_THRESHOLD
        print (out_name)
        ann_scrna, ann_guide = analyze_batch(scRNA_ann_FILE = file_scRNA,
                      guide_ann_FILE = file_guide,
                      BATCH_NUMBER = lane_n,
                      OUT_DIR_FIGURES = out_name,
                      EXPECTED_NUMBER_OF_CELLS = exp_number_of_cells,
                      SPECIE_MITO = specie_mito,
                      MITO_PERCENTAGE_ALLOWED = mito_percent_allowed,
                     UMI_CELL_THRESHOLD=  tresh_umi_cell  
                     )
        ann_scrna_data_to_concat.append(ann_scrna)
        ann_guide_data_to_concat.append(ann_guide)
        
        
    return ann_scrna_data_to_concat, ann_guide_data_to_concat



def concact_lanes(ann_scrna_data_to_concat, ann_guide_data_to_concat):
    
    #Concact batch results in a single ANN data file.  
    print ('Concacting lanes...')
    concat_scrna_ann = anndata.concat(ann_scrna_data_to_concat, merge="same")
    concat_scrna_ann.obs_names_make_unique()
    concat_guide_ann = anndata.concat(ann_guide_data_to_concat, merge="same")
    concat_guide_ann.obs_names_make_unique()
    return concat_scrna_ann, concat_guide_ann


def filtering_low_expressed_genes(concat_scrna_ann,PERCENTAGE_OF_CELLS_TO_INCLUDE_TRANSCRIPT = 0.01):
    ### Filtering lower represented transcripts  
    #- Filter transcripts not presented on at least 1% of all cells (all lanes)
    #default  #  1%
    print ('Filtering Low expressed genes...')
    MINIMAL_CELLS_TO_INCLUDE_TRANSCRIPT = int(concat_scrna_ann.shape[0] * PERCENTAGE_OF_CELLS_TO_INCLUDE_TRANSCRIPT)
    sc.pp.filter_genes(concat_scrna_ann, min_cells=MINIMAL_CELLS_TO_INCLUDE_TRANSCRIPT) #remove genes not present more than > 1% of all cells  
    



def saving_anndata_files( concat_scrna_ann,  concat_guide_ann ):
    print ('saving ann data...')
    concat_scrna_ann.write(f'results_per_lane/full_raw_scrna_ann_data.h5ad')
    concat_guide_ann.write(f'results_per_lane/full_raw_guide_ann_data.h5ad')
    
    
  

print ('executing')

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--path', type=str,
                help='The path to the directory having the initialized df (default: current working directory)')

parser.add_argument('--expected_cell_number', type=int, default=8000,
                    help='The expected number of cells in the sample (default: 8000)')
parser.add_argument('--mito_specie', type=str, default='hsapiens',
                    help='The mitochondrial species to use for filtering (default: hsapiens)')
parser.add_argument('--mito_expected_percentage', type=float, default=0.2,
                    help='The expected percentage of mitochondrial reads (default: 0.2)')
parser.add_argument('--percentage_of_cells_to_include_transcript', type=float, default=0.01,
                    help='The percentage of cells to include in transcript counts (default: 0.01)')
parser.add_argument('--transcripts_umi_treshold', type=int, default=100,
                    help='The minimum  UMI count  to consider keep a cell in the analysis (default: 100)')

args = parser.parse_args()

PATH = args.path
EXPECTED_CELL_NUMBER = args.expected_cell_number
MITO_SPECIE = args.mito_specie
MITO_EXPECTED_PERCENTAGE = args.mito_expected_percentage
PERCENTAGE_OF_CELLS_TO_INCLUDE_TRANSCRIPT = args.percentage_of_cells_to_include_transcript
UMI_CELL_THRESHOLD = args.transcripts_umi_treshold

df_processed = preparing_files(PATH)

in_concact_scrna_ann, in_concact_guide_ann = execute_analysis(df_processed, EXPECTED_CELL_NUMBER, MITO_SPECIE, MITO_EXPECTED_PERCENTAGE, UMI_CELL_THRESHOLD )

concat_scrna_ann, concat_guide_ann =  concact_lanes(in_concact_scrna_ann, in_concact_guide_ann)

filtering_low_expressed_genes(concat_scrna_ann) # inplace operation

saving_anndata_files(concat_scrna_ann, concat_guide_ann)

