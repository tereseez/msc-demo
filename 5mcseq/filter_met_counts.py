#!/usr/bin/env python
#start filter_met_counts.py

import pandas as pd
import multiprocessing as mp
import glob


def met_filtering(filename):
    """
    input: bismark genome-wide CpG-context methylation report
    - filters file based on:
        - total count (methylated count + unmethylated count) >= 5
        - methylation count  >= 2
    - calculates methylation percentage
    output: bed file
    """
    in_file = filename
    in_df = pd.read_csv(in_file, sep='\t', header = None)
    in_df.columns = ["chrom", "position", "strand", "count_met",
                     "count_unmet", "context", "trinuc_context"]

    filtered_df = (in_df[(in_df["count_met"] >= 2) &
                   ((in_df["count_met"] + in_df["count_unmet"]) >= 5)])

    total_counts = []
    met_percentage = []
    for idx, row in filtered_df.iterrows():
        sum_count = row.count_met + row.count_unmet
        percentage = row.count_met / sum_count
        total_counts.append(sum_count)
        met_percentage.append(percentage)
    filtered_df.reset_index(drop=True)

    chrom = filtered_df.iloc[:,0]
    position = filtered_df.iloc[:,1]
    strand = filtered_df.iloc[:,2]

    df_out = pd.DataFrame(list(zip(chrom, position, position,
                                   met_percentage,
                                   total_counts, strand)))
    out_file = in_file.replace('cov2cytosine.txt.', '')
    out_file = out_file.replace('.CpG_report.txt', '')
    out_file = out_file + '_CpG_methylation.bed' 
    df_out.to_csv(out_file, sep='\t', index=False, header=None)

    return filename

#do the calculcations on all chromosomes in a sample
files = glob.glob('out_met_count/*/coverage2cytosine/*.CpG_report.txt')
with mp.Pool() as pool:
    results = pool.map(met_filtering, files)

#end filter_met_counts.py
