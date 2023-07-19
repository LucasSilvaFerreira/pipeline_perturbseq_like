#!/usr/bin/env python
import sys
import anndata
import pandas as pd
import mudata as md

MUON_DATA = sys.argv[1]
print ('processing MUON')

ann_read =  md.read(MUON_DATA)['scRNA']
pd.DataFrame(ann_read.obs.index.values).to_csv('cell_barcode_capturing.csv', index=None, header=None)

