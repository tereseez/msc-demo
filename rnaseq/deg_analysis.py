#!/usr/bin/env python
#start deg_analysis.py

import pandas as pd
from scipy import stats
import numpy as np
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import seaborn as sb

#retreive Novogene FPKM values
in_file = 'ref/fpkm_genename.xls'
in_df = pd.read_csv(in_file, sep='\t')

#filter for protein coding transcripts
in_df = in_df.loc[(in_df['gene_biotype']) == 'protein_coding']

a_df = in_df.iloc[:,1:11].join(in_df.iloc[:,31:41])
a_df.index = in_df.gene_id

#log2-transform
log2_fpkm_df = np.log2(a_df)
log2_fpkm_df2 = log2_fpkm_df.replace([np.inf, -np.inf], np.nan).dropna()

#standardize to z-scores
zscore_log_fpkm = log2_fpkm_df2.apply(stats.zscore)
zscore_log_fpkm2 = zscore_log_fpkm.reset_index()

#produce DEG lists with gene names, p-values and t-value
degs = []
pvalues = []
tvalues = []

for idx, row in zscore_log_fpkm2.iterrows():
    group1 = row[1:11]
    g1 = group1.astype(float)
    group2 = row[11:]
    g2 = group2.astype(float)
    statistics, pvalue = stats.ttest_ind(g1,g2)
    if pvalue < 0.05:
        degs.append(row['gene_id'])
        pvalues.append(pvalue)
        tvalues.append(statistics)

df_deg = pd.DataFrame(list(zip(degs, pvalues, tvalues)),
                      columns = ['gene id', 'p-values',
                                 't-statistics'])
out_file_deg = 'degs_p0.05_fpkm.tsv'
df_deg.to_csv(out_file_deg, sep='\t', index=False)

#produce DEGs from subset
in_subset = 'normalization/subset_zscores.tsv'
df_subset = pd.read_csv(in_subset, sep='\t')

degs_sub = []
pvalues_sub = []
tvalues_sub = []

for idx, row in df_subset.iterrows():
    group1 = row[1:6]
    g1 = group1.astype(float)
    group2 = row[6:]
    g2 = group2.astype(float)
    statistics, pvalue = stats.ttest_ind(g1,g2)
    if pvalue < 0.05:
        degs_sub.append(row['Geneid'])
        pvalues_sub.append(pvalue)
        tvalues_sub.append(statistics)

df_subset_deg = pd.DataFrame(list(zip(degs_sub, pvalues_sub, tvalues_sub)), 
        columns = ['gene id', 'p-values', 't-statistics'])
out_file_subset = 'degs_p0.05_subset.tsv'
df_subset_deg.to_csv(out_file_subset, sep='\t', index=False)

#compare DETs (FPKM) to Novogene's results
novo_file = 'ref/Adrenal1vsAdrenal2_DEG_all.xls'
novo_df = pd.read_csv(novo_file, sep='\t')

only_protein_coding = novo_df.loc[((novo_df['gene_biotype'])
                                   == 'protein_coding')]

deg_set = set(df_deg['gene id'])
novo_set = set(only_protein_coding['gene_id'])
total = len(deg_set)

v = venn2([novo_set, deg_set],
      set_labels = ('Novogene`s DETs (10+10)', 'DETs (FPKM, 10+10)'),
      set_colors=('#5b9bd5', 'mediumpurple'),
      subset_label_formatter = lambda x: str(x) +
       "\n(" + f"{(x/total):1.0%}" + ")")

v.get_patch_by_id('01').set_alpha(0.3)
v.get_label_by_id('10').set_text(len(novo_set))
v.get_label_by_id('01').set_text(len(deg_set))
plt.title('Overlap DETs')
plt.show()

#repeat comparison for DEGs (subset) vs Novogene final DEGs
deg_subset_set = set(df_subset_deg['Geneid'])
novo_subset_set = set(only_protein_coding['gene_name'])
total_subset = len(deg_subset_set)

v = venn2([novo_subset_set, deg_subset_set],
      set_labels = ('Novogene`s DEGs (10+10)', 'Own DEGs (5+5)'),
      set_colors=('#5b9bd5', 'mediumpurple'),
      subset_label_formatter = lambda x: str(x) +
       "\n(" + f"{(x/total_subset):1.0%}" + ")")

v.get_patch_by_id('01').set_alpha(0.3)
v.get_label_by_id('10').set_text(len(novo_subset_set))
v.get_label_by_id('01').set_text(len(deg_subset_set))
plt.title('Overlap DEGs')
plt.show()
#end deg_analysis.py
