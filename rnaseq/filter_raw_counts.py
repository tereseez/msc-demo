#!/usr/bin/env python
#start filter_raw_counts.py 

import pandas as pd
import csv
import glob

#filter gtf file for protein coding genes only
ref_file = 'ref/genome.gtf'
ref_df = pd.read_csv(ref_file, sep='\t', skiprows=[0,1,2,3,4], header=None)
ref_gene_df = ref_df[ref_df.iloc[:,2]=='gene'].copy()
out_protein_gene_df = (ref_gene_df[ref_gene_df.loc[:,8]
                       .str.contains('protein_coding')].copy())

ref_out_file = ref_file.replace('.gtf', '')
ref_out_file = ref_out_file + '_proteinCoding_only.gtf'
out_protein_gene_df.to_csv(ref_out_file, sep='\t',
                           index=False, header=None,
                           quoting=csv.QUOTE_NONE, quotechar='')

#find gene names in gtf to filter count tables
def find_gene_name(x):
    out_X = ([X for X in x if 'gene_name' in X][0]
             .strip().replace('gene_name','').strip().replace('"',''))
    return out_X

gtf_genename_df = (out_protein_gene_df.iloc[:,8]
                   .apply(lambda x: find_gene_name(x.split(';'))))

#filter raw read counts
in_files = glob.glob('out_counts/*/*/*featureCounts.tsv')
for file in in_files:
    df = pd.read_csv(file, sep='\t', skiprows=[0])
    out_gene_df = df[df['Geneid'].isin(gtf_genename_df)]
    out_file = file.replace('featureCounts.tsv','')
    out_file = out_file + 'ftCounts_filtered.tsv'
    out_gene_df.to_csv(out_file, sep='\t', index=False)

#end filter_raw_counts.py
