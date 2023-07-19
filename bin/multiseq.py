#!/usr/bin/env python
import sys
import os


R1_MULTI = sys.argv[1]
R2_MULTI = sys.argv[2]
BARCODES_CELL_LIST_MULTI =     sys.argv[3]
BARCODES_MULTIBAR_LIST_MULTI = sys.argv[4]
BAR_MULTI_0=    sys.argv[5]
BAR_MULTI_1=    sys.argv[6]
UMI_MULTI_0 =     sys.argv[7]
UMI_MULTI_1 =     sys.argv[8]
R2_MULTI_TAG_0 =  sys.argv[9]
R2_MULTI_TAG_1 =  sys.argv[10]
MUON_DATA = sys.argv[11]



print (BARCODES_CELL_LIST_MULTI, BARCODES_MULTIBAR_LIST_MULTI)

code_save = f'''

library(deMULTIplex)
library('ggplot2')

cell.id.vec_df  <- read.table("{BARCODES_CELL_LIST_MULTI}") 
cell.id.vec = as.vector(t(cell.id.vec_df))
bar.ref_load <- read.table("{BARCODES_MULTIBAR_LIST_MULTI}", sep=",")
bar.ref =as.vector(t(bar.ref_load["V1"]))

readTable <- MULTIseq.preProcess(R1 = "{R1_MULTI}",
                                 R2 = "{R2_MULTI}", 
                                 cellIDs = cell.id.vec,
                                 cell=c( {BAR_MULTI_0},   {BAR_MULTI_1}),
                                 umi=c(  {UMI_MULTI_0},   {UMI_MULTI_1}),
                                 tag=c(  {R2_MULTI_TAG_0},{R2_MULTI_TAG_1})) 
 
bar.table <- MULTIseq.align(readTable, cell.id.vec, bar.ref)    

## Visualize barcode space
#print (bar.table)
bar.tsne <- barTSNE(bar.table) #more threads?
## Note: Exclude columns 97:98 (assuming 96 barcodes were used) which provide total barcode UMI counts for each cell. 

pdf("bc.check.pdf")
for (i in 3:ncol(bar.tsne)) {{
    g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
    ##print(g)
}}
dev.off()

## Round 1 -----------------------------------------------------------------------------------------------------
## Perform Quantile Sweep
bar.table.full <- bar.table
write.csv(bar.table ,   "bar_table.csv" ) 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {{
  ##print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
pdf("class_threshrold.check.pdf")
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

dev.off()

## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {{
  ##print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

## Repeat until all no negative cells remain (usually 3 rounds)...
final.calls <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round2.calls),neg.cells)

write.csv(final.calls ,                "final_class.csv" ) 
write.csv(attr(final.calls, 'names') , "final_class_cell_barcode.csv" ) 

'''

with open('multiseq.r', 'w') as file_handle:
    file_handle.write(code_save)
    
os.system('Rscript multiseq.r')


import mudata as md
import pandas as pd
import anndata as ad
from mudata import MuData



muon_read =  md.read(MUON_DATA)
df_count_bar = pd.read_csv('bar_table.csv')
df_count_bar = df_count_bar.set_index('Unnamed: 0')


adata = ad.AnnData(df_count_bar[df_count_bar.columns[:-2]])
#adata.var

df_final_class = pd.read_csv('final_class.csv')
del df_final_class['Unnamed: 0']
df_final_cellbar = pd.read_csv('final_class_cell_barcode.csv')
del df_final_cellbar['Unnamed: 0']
df_final_cellbar.columns = ['cellbar']
df_concat_barcodes  = pd.concat([df_final_class, df_final_cellbar], axis=1)
df_concat_barcodes = df_concat_barcodes.set_index('cellbar')
df_concat_barcodes = df_concat_barcodes.groupby(df_concat_barcodes.index).first()
adata.obs['multiseq_class'] = df_concat_barcodes.loc[adata.obs.index.values.tolist()]
mdata = MuData({"scRNA": muon_read['scRNA'], "guides": muon_read['guides'], "multiseq":adata})
mdata["guides"].obs['multiseq_class'] = df_concat_barcodes.loc[adata.obs.index.values.tolist()]
mdata["scRNA"].obs['multiseq_class'] = df_concat_barcodes.loc[adata.obs.index.values.tolist()]
mdata = mdata[mdata['multiseq'].obs['multiseq_class'].apply(lambda x : x not in {'Double', 'Negative'}), : ]
mdata.write("processed_mudata_guide_and_transcripts_multiseq_filtered.h5mu")

