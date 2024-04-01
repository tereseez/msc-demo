#!/usr/bin/env python
#start normalize_counts.py

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sb
from scipy.stats import zscore

#normalisation functions: quantile and rpkm
def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values, axis=0),
                             index=df.index,
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn = (df.rank(method="min").stack().astype(int)
               .map(df_mean).unstack())
    return(df_qn)

def rpkm(counts, lengths):
    """
    Calculate reads per kilobase transcript per million reads.
    RPKM = (10^9 * C) / (N * L)
    Where:
    C = Number of reads mapped to a gene
    N = Total mapped reads in the experiment
    L = Exon length in base pairs for a gene
    Parameters
    ----------
    counts: array, shape (N_genes, N_samples)
        RNAseq (or similar) count data where columns are
        individual samples and rows are genes.
    lengths: array, shape (N_genes,)
        Gene lengths in base pairs in the same order
        as the rows in counts.
    Returns
    -------
    normed : array, shape (N_genes, N_samples)
        The RPKM normalized counts matrix.
    """
    N = np.sum(counts, axis=0)
    L = lengths
    C = counts
    #add a pesudo count
    normed = np.divide(1e9 *(C+1),
                       N[np.newaxis, :] * L[:, np.newaxis])
    return(normed)

#prepare input data - combine samples
test_files = sorted(
        glob.glob('out_counts/hisat2/test/*h2_ftCounts_filtered.tsv')
        )
ctr_files = sorted(
        glob.glob('out_counts/hisat2/ctr/*h2_ftCounts_filtered.tsv')
        )

merged_df = pd.read_csv(test_files[0], sep='\t')

for file in test_files[1:]:
    df = pd.read_csv(file, sep='\t')
    df = df.iloc[:,[0,-1]]
    merged_df = merged_df.merge(df, on='Geneid',how='inner')

for file in ctr_files:
    df = pd.read_csv(file, sep='\t')
    df = df.iloc[:,[0,-1]]
    merged_df = merged_df.merge(df, on='Geneid',how='inner')

merged_df = merged_df.set_index('Geneid')
merged_df2 = merged_df.iloc[:,5:]

#remove rows with missing values in 50% of columns
n = len(merged_df2.columns)/2
merged_df3 = (merged_df2.replace(0, np.nan)
                        .dropna(axis = 0, thresh = n)
                        .fillna(0).astype(int)
                        )

#quantile-normalize counts
qt_norm_counts = quantile_normalize(merged_df3)
qt_norm_counts.columns = ['A1', 'A2', 'A3', 'A4', 'A5',
                          'A11', 'A12', 'A13', 'A14', 'A15']

out_qt = qt_norm_counts.reset_index()
out_qt_file = 'subset_quant_norm.tsv'
out_qt.to_csv(out_qt_file, sep='\t', index=False)

#compare quant. norm. counts with Novogene's normalised counts
file_novo = 'ref/Adrenal1vsAdrenal2_DEG.xls'
df_novo = pd.read_csv(file_novo, sep='\t')
df_novo.rename(columns={'gene_name':'Geneid'}, inplace=True)

merged_novo_df = qt_norm_counts.merge(df_novo, on='Geneid', how='inner')

#scatterplot of A1 counts
f, ax = plt.subplots(figsize=(8,5))
plt.scatter(merged_novo_df['A1_x'],
            merged_novo_df['A1_y'],
            color = '#5b9bd5', alpha = 0.5)
plt.xlabel('A1: quantile normalised')
plt.ylabel('A1: normalised by Novogene')
ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle='--')
plt.show()

##make numpy array for calculating RPKM values
#get gene lengths
merged_df_mod = merged_df[merged_df.index.isin(merged_df3.index)]
lengths = merged_df_mod.Length.to_numpy()
counts = qt_norm_counts.to_numpy()

#RPKM
rpkm_norm_counts = rpkm(counts, lengths)
rpkm_norm_counts_df = pd.DataFrame(rpkm_norm_counts,
                      columns=['A1', 'A2', 'A3', 'A4', 'A5',
                               'A11', 'A12', 'A13', 'A14', 'A15'],
                      index=qt_norm_counts.index)

out_rpkm_df = rpkm_norm_counts_df.reset_index()
out_rpkm_file = 'subset_rpkm.tsv'
out_rpkm_df.to_csv(out_rpkm_file, sep='\t', index=False)

#log2 transform rpkm values
log2_rpkm_df = np.log2(rpkm_norm_counts_df)
out_log2_rpkm_df = log2_rpkm_df.reset_index()
out_log2_rpkm = 'subset_logrpkm.tsv'
out_log2_rpkm_df.to_csv(out_log2_rpkm, sep='\t', index=False)

#Zscores
zscore_log_rpkm = log2_rpkm_df.apply(zscore)
out_zscores_df = zscore_log_rpkm.reset_index()
out_zscores = 'subset_zscores.tsv'
out_zscores_df.to_csv(out_zscores, sep='\t', index=False)

#distribution of reads: rpkm, log2(rpkm), z-scores of log2(rpkm)
##remove outliers in rpkm table
rpkm_non_outliers = rpkm_norm_counts_df[rpkm_norm_counts_df < 1e3].dropna()

fig = plt.figure(figsize=(15,4))
axes1 = fig.add_subplot(131)
hist1 = rpkm_non_outliers['A1'].hist(bins=50, ax=axes1, color = '#5b9bd5',
                                   alpha = 0.7)
hist1.set_title('A1: RPKM')

axes2 = fig.add_subplot(132)
hist_log1 = log2_rpkm_df['A1'].hist(bins=15, ax=axes2, color = '#5b9bd5',
                                    alpha = 0.7)
hist_log1.set_title('A1: Log2(RPKM)')

axes3 = fig.add_subplot(133)
hist_zscores = zscore_log_rpkm['A1'].hist(bins=15, ax=axes3,
                                          color = '#5b9bd5', alpha = 0.7)
hist_zscores.set_title('A1: Z-scores of log2(RPKM)')
plt.show()

#end normalize_counts.py
