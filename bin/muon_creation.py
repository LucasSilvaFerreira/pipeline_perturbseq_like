#!/usr/bin/env python

import argparse
import anndata
from GTFProcessing import GTFProcessing
import pandas as pd
import muon as mu

parser = argparse.ArgumentParser()
parser.add_argument('--ann_guide', required=True)
parser.add_argument('--transcript_file', required=True)
parser.add_argument('--gtf_in', required=True)
parser.add_argument('--ann_exp', required=True)
args = parser.parse_args()

ANN_GUIDE = args.ann_guide
TRANSCRIPT_FILE = args.transcript_file
GTF_IN = args.gtf_in
ANN_EXP = args.ann_exp

print(ANN_GUIDE)

print(anndata.read(ANN_GUIDE).var)



def ann_guide_modification():
    ann_guide = anndata.read(ANN_GUIDE)
    ann_guide.var['guide_chr'] = [k.split('|')[1].split(':')[0].split('_')[-1]  for k in ann_guide.var.index]
    ann_guide.var['guide_end'] = [k.split('|')[1].split(':')[-1]  for k in ann_guide.var.index]
    ann_guide.var['guide_start'] = [k.split('|')[1].split(':')[-2]  for k in ann_guide.var.index]
    ann_guide.var['guide_number'] = [k.split('|')[1].split(':')[0].split('_')[0]  for k in ann_guide.var.index]
    ann_guide.var['target_elements'] = [k.split('|')[0]  for k in ann_guide.var.index]
    print ('finishing adding position information to the guides anndata')
    return ann_guide


def get_coords(x_ensg, chr_dict, chr_start, chr_end):
    if x_ensg in chr_dict:
        return chr_dict[x_ensg], chr_start[x_ensg], chr_end[x_ensg]
    else:
        return 'NOT_FOUND', 0, 0

def ann_transcript_modification():
    #capuring genes start sito  to add in the anndata
    #transcripts file format ENST00000456328.2	ENSG00000223972.5	DDX11L1
    ann_exp =   anndata.read(ANN_EXP)
    gtf = GTFProcessing(GTF_IN)
    df_gtf_refseq = gtf.get_gtf_df()
    df_gtf_refseq_gene =df_gtf_refseq[ df_gtf_refseq['feature'] == 'gene']
    transcripts_df = pd.read_csv(TRANSCRIPT_FILE, sep='\t', header=None, names =['transcripts', 'genes', 'gene_name'] )
    transcripts_df['genes'] = transcripts_df['genes'].apply(lambda x : x.split('.')[0] )
    set_genes_gene_name = dict(zip(transcripts_df['genes'].values, transcripts_df['gene_name'].values))
    set_gene_on_annexp = set([set_genes_gene_name[x] for x in ann_exp.var.index.values])
    df_gtf_refseq_gene = df_gtf_refseq_gene[df_gtf_refseq_gene['gene_name'].apply(lambda x : x in set_gene_on_annexp)]
    df_gtf_refseq_gene['gene_real_start'] = df_gtf_refseq_gene.apply(lambda x : x['start'] if x['strand'] == "+" else x['end'], axis=1 )
    df_gtf_refseq_gene['start'] = df_gtf_refseq_gene['gene_real_start'] 
    df_gtf_refseq_gene['end '] = df_gtf_refseq_gene['start']  + 1 
    #Adding coordnates to anndata expression
    chr_dict = dict(zip(df_gtf_refseq_gene['gene_id'].values, df_gtf_refseq_gene['chr'].values))
    chr_start = dict(zip(df_gtf_refseq_gene['gene_id'].values, df_gtf_refseq_gene['start'].values))
    chr_end = dict(zip(df_gtf_refseq_gene['gene_id'].values, df_gtf_refseq_gene['end'].values))
    ann_exp.var['transcript_chr'] = [get_coords(x,chr_dict, chr_start, chr_end)[0] for x in ann_exp.var.index.values]   
    ann_exp.var['transcript_start'] = [get_coords(x, chr_dict, chr_start, chr_end)[1] for x in ann_exp.var.index.values]   
    ann_exp.var['transcript_end'] = [get_coords(x, chr_dict, chr_start, chr_end)[2] for x in ann_exp.var.index.values]  
    print ('finishing adding position information to the  transcripts anndata')

    return ann_exp


mdata =  mu.MuData({"guides":  ann_guide_modification(), "scRNA":  ann_transcript_modification()})
mdata.write("raw_mudata_guide_and_transcripts.h5mu")
