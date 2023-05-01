#!/usr/bin/env python
import sys
import os
import pandas as pd

def format_number(number):
    return f"{number:03d}"

def main():
    input_csv = sys.argv[1]
    df = pd.read_csv(input_csv)
    
    os.makedirs('symlink_dir', exist_ok=True)

    new_names = []
    for index, row in df.iterrows():
        new_name = f"{row['Sample']}_S{format_number(row['Sample_Number'])}_L{format_number(row['Lane_Number'])}_{row['Read_Type']}_{format_number(row['File_Fragment'])}.fastq.gz"
        new_names.append(new_name)
        os.symlink(row['file_dir'], f"symlink_dir/{new_name}")
    
    df['New_Name'] = new_names
    df.to_csv('output.csv', index=False)

if __name__ == "__main__":
    main()
