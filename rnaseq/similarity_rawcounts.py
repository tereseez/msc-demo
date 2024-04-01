#!/usr/bin/env python
#start similarity_rawcounts.py
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#filtered featureCounts results from HISAT2 and STAR
in_h2_ftCount_file = 'out_counts/hisat2/test/A1_h2_ftCounts_filtered.tsv'
in_h2_ftCount_df = pd.read_csv(in_h2_ftCount_file, sep='\t')
h2_count = in_h2_ftCount_df.iloc[:, np.r_[0, -1]]

in_star_ftCount_file = 'out_counts/star/test/A1_star_ftCounts_filtered.tsv'
in_star_ftCount_df = pd.read_csv(in_star_ftCount_file, sep='\t')
star_count = in_star_ftCount_df.iloc[:, np.r_[0, -1]]

#read count table provided by Novogene
in_novo_file = 'ref/readcount_genename.xls'
in_novo_df = pd.read_csv(in_novo_file, sep='\t')
in_novo_df.rename(columns={'gene_name':'Geneid'}, inplace=True)
novo_count = in_novo_df[['Geneid', 'A1']]

#scatter plot, for example HISAT2 vs STAR
h2_star_merged_df = h2_count.merge(star_count, on='Geneid',how='inner')
f, ax = plt.subplots(figsize=(7,5))
plt.scatter(h2_star_merged_df.iloc[:,1], h2_star_merged_df.iloc[:,2],
            color = '#5b9bd5', alpha = 0.5)
plt.xlabel('HISAT2')
plt.ylabel('STAR')
plt.title('HISAT2 vs STAR: A1 raw read counts')
ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle='--')
plt.show()

#end similarity_rawcounts.py
