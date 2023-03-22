#!/usr/bin/env python

#Read composition discover
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import subprocess


def get_reads_fasta(fasta):
    if '.gz' in fasta:
        zcat='z'
    else:
        zcat=''
    check_guides = f"{zcat}cat '{fasta}' | head -n 10000 "
    print (check_guides)
    seqs = os.popen(check_guides).read().split('\n')
    #print (seqs)
    return [ x for x in seqs if '@' not in x and '+' not in x and 'F'not in x]
    
def compositional_bias_calculation(fasta_df, name, dir_out) :  
    plt.rcParams["figure.figsize"] = (25,2)
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))

    df_calculate_percentages = pd.DataFrame([  [ a for a in s[:]] for s in fasta_df ])
    (pd.concat([ (df_calculate_percentages == nucleotide).sum()  for nucleotide in ['A', 'C', 'T', 'G']], axis=1).T - df_calculate_percentages.shape[0] ).std().plot()
    plt.xlim(0,df_calculate_percentages.shape[1])
    plt.xticks(rotation=90)
    plt.title(f'{name}:: + Compositional Bias \n (sd between the 4 nucleotides) ')
    plt.tight_layout()
    plt.savefig(f'{dir_out}_composition/{name}_composition.png')
    plt.clf()

    
def plot_compositional_bia(R1,R2, prefix_tag):    
    compositional_bias_calculation(get_reads_fasta(R1), name=f'{prefix_tag}_R1' ,dir_out=prefix_tag)   
    compositional_bias_calculation(get_reads_fasta(R2), name=f'{prefix_tag}_R2', dir_out = prefix_tag )   


read1 = sys.argv[1]
read2 = sys.argv[2]
print ('read1', read1)
#read1 = ' '.join(read1.split(' ')[0])
#read2 = ' '.join(read2.split(' ')[0])
#print (read1.split())
#print (read1)
prefix_tag = sys.argv[-1]
plot_compositional_bia(read1, read2, prefix_tag)
