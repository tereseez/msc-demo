#!/usr/bin/env python
#start mean_CpGsites.py

import pandas as pd
import multiprocessing as mp
import glob
import fnmatch
import re
import os

def mean_methylation(filenames):
    f = filenames[0]
    loop = 0
    
    #merge ctr or test samples within chr
    for fi in filenames:
        #get sample name
        tmp_name = '_'.join(os.path.basename(fi).split('_')[0:2])
        in_df = pd.read_csv(fi, sep='\t', header=None)
        in_df.columns = ["chrom", "start", "end", "percentage",
                         "total_count", "strand"]
        if loop == 0:
            merged_frame = in_df.copy()
        else:
            merged_frame = pd.merge(merged_frame, in_df,
                                    on = ["chrom","start","end","strand"],
                                    how ='outer',
                                    suffixes = ('', '_' + tmp_name)).copy()
        loop =+1
        
    #retrieve all percentage columns and count columns
    percentage_cols = [col for col in merged_frame.columns\
                       if 'percentage' in col]
    count_cols = [col for col in merged_frame.columns\
                  if 'total_count' in col]
                  
    merged_frame.fillna(0, inplace=True)
    merged_frame[percentage_cols] = (merged_frame[percentage_cols]
                                     .astype(float))
    merged_frame[count_cols] = (merged_frame[count_cols]
                                .astype(int))
    
    #compute mean percentage and total counts
    merged_frame['mean_percentage'] = (merged_frame[percentage_cols]
                                       .mean(axis=1))
    merged_frame['mean_total_counts'] = (merged_frame[count_cols]
                                         .mean(axis=1))
                                         
    #select columns for output
    out_df = merged_frame.iloc[:,0:3]
    out_df["percentage"] = merged_frame['mean_percentage'].copy()
    out_df["total_count"] = merged_frame['mean_total_counts'].copy()
    out_df = out_df.join(merged_frame["strand"])
    
    #set output filename
    out_file = 'mean_'
    if fnmatch.fnmatch(f, '*Test*'):
        out_file += 'test_'
    else:
        out_file += 'control_'
        
    chr_value = re.findall(r'(chr([0-9]{1,2}|[a-zA-z]){1})', f)[0][0]
    out_file += chr_value + '.csv'
    
    print(out_file)
    path = 'out_met_count/run2_junbai/CpG_context/mean/'
    output = os.path.join(path, out_file)
    out_df.to_csv(output, sep='\t', index=False, header=None)

#list chr
chr_files = glob.glob('out_met_count/run2_junbai/CpG_context/chr*')

#list test/control sample files for each chromosome
test_files = []
control_files = []
for chr in chr_files:
    test_files.append(glob.glob(f'{chr}/*Test*.csv'))
    control_files.append(glob.glob(f'{chr}/*Ctrl*.csv'))

#do the calculations individually on ctr and test
with mp.Pool() as pool:
    test_results = pool.map(mean_methylation, test_files)
    control_results = pool.map(mean_methylation, control_files)

#end mean_CpGsites.py
