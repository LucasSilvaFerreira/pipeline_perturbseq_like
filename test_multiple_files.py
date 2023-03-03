import sys
import pandas as pd



print('MY SYS PYTHON:  ',sys.argv)

for p in sys.argv[1:]:
    print (pd.read_csv(p, header=None, sep='\t'))
