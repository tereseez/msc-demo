#!/usr/bin/env python
#start filtering_deg.py
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as np

in_adrenal_file = '~/old_project/RNAseq/Novogene data/Adrenal/Adrenal1vsAdrenal2_DEG_all.xls'
in_adrenal_df = pd.read_csv(in_adrenal_file, sep='\t')

#Filter for protein coding genes
only_protein_coding = in_adrenal_df.loc[((in_adrenal_df['gene_biotype'])
                                         == 'protein_coding')]
a_df = only_protein_coding.iloc[:,0:25]

#Filtering
a_df = a_df[a_df['pvalue'] < 0.05]
#a_df_padj = a_df[a_df['padj'] < 0.1]

a_df = a_df.drop_duplicates(subset='gene_name')
a_df = a_df.rename(columns={"A1":"T1", "A2":"T2", "A3":"T3", "A4":"T4", "A5":"T5",
                            "A6":"T6", "A7":"T7", "A8":"T8", "A9":"T9", "A10":"T10"})

a_df = a_df.rename(columns={"A11":"C1", "A12":"C2", "A13":"C3", "A14":"C4", "A15":"C5",
                            "A16":"C6", "A17":"C7", "A18":"C8", "A19":"C9", "A20":"C10"})

#Heatmap of DEGs, repeat for the two lists
a_df.index = a_df.gene_name

cluster = sb.clustermap(a_df.iloc[:,1:21], z_score=0, figsize=(10,8),
                        yticklabels=False, cmap='vlag', vmin=-1.5, vmax=2,
                        cbar_kws={"ticks":[-1,0,1,2]})

cluster.fig.subplots_adjust(top=0.95)
cluster.fig.suptitle('Adrenal DEGs (p < 0.05)', fontsize=14)
cluster.ax_heatmap.set_xticklabels(cluster.ax_heatmap.xaxis.get_majorticklabels(),
                                    fontsize = 12)
cluster.ax_heatmap.set_ylabel("")

cluster.ax_cbar.tick_params(axis='y', labelsize=12)
cluster.ax_cbar.set_ylabel('Z-score', size=12)
cluster.ax_cbar.set_position((0.04, .8, 0.05, 0.15))

#plt.show()
plt.savefig("/Users/teresezylla/Desktop/results.png")


# #Barplots of distribution of DEGs
# a_df_up = len(a_df[a_df['log2FoldChange'] > 0])
# a_df_down = len(a_df[a_df['log2FoldChange'] < 0])
# original = [a_df_up, a_df_down]
#
# a_df_padj_up = len(a_df_padj[a_df_padj['log2FoldChange'] > 0])
# a_df_padj_down = len(a_df_padj[a_df_padj['log2FoldChange'] < 0])
# adjusted = [a_df_padj_up, a_df_padj_down]
#
# degs = ['Up', 'Down']
# X_axis = np.arange(len(degs))
#
# fig, ax = plt.subplots()
# bplt1 = ax.bar(X_axis - 0.2, original, 0.4, color='steelblue',
#                label='p-value < 0.05')
# bplt2 = ax.bar(X_axis + 0.2, adjusted, 0.4, color='#5b9bd5',
#                alpha = 0.5, label='FDR < 0.1')
# plt.xticks(X_axis, degs)
# ax.set_xticks(X_axis)
# ax.set_xticklabels(degs)
# plt.ylabel('Number of genes')
# plt.title('Adrenal DETs')
#
# # Shrink current axis by 20%
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# # Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.show()
#
# #Check for potential duplicates
# duplicates_original = (a_df.duplicated(subset=['gene_name'])
#                        .value_counts()[True])
# duplicates_adjusted = (a_df_padj.duplicated(subset=['gene_name'])
#                        .value_counts()[True])
#
# #Prepare for export, repeat for the two lists
# condition1 = list('Testgroup_' + a_df.columns.values[1:11])
# condition2 = list('Controlgroup_' + a_df.columns.values[11:21])
# merged = [condition1, condition2]
#
# import itertools
# merged_2 = list(itertools.chain.from_iterable(merged))
# merged_2.insert(0, 'gene_id')
# merged_2.insert(21, 'log2FoldChange')
# merged_2.insert(22, 'pvalue')
# merged_2.insert(23, 'padj')
# merged_2.insert(24, 'gene_name')
#
# a_df.columns = merged_2
#
# #Export final DEG list (p<0.05)
# out_file = 'adrenal_deg_filtered_p0.05.tsv'
# a_df.to_csv(out_file, sep='\t', index=False)
#
# #Export unique gene list
# out_file_genenames = 'adrenal_deg_filtered_p0.05_unique_genelist.tsv'
# a_df['gene_name'].to_csv(out_file_genenames, index=False)
#
# #repeat for FDR<0.1
# condition1 = list('Testgroup_' + a_df_padj.columns.values[1:11])
# condition2 = list('Controlgroup_' + a_df_padj.columns.values[11:21])
# merged = [condition1, condition2]
#
# merged_2 = list(itertools.chain.from_iterable(merged))
# merged_2.insert(0, 'gene_id')
# merged_2.insert(21, 'log2FoldChange')
# merged_2.insert(22, 'pvalue')
# merged_2.insert(23, 'padj')
# merged_2.insert(24, 'gene_name')
#
# a_df_padj.columns = merged_2
#
# #Export final DEG list (FDR<0.1)
# out_file = 'adrenal_deg_filtered_fdr0.1.tsv'
# a_df_padj.to_csv(out_file, sep='\t', index=False)
#
# #Export unique gene list
# out_file_genenames = 'adrenal_deg_filtered_fdr0.1_unique_genelist.tsv'
# a_df['gene_name'].to_csv(out_file_genenames, index=False)
#
# # Summary analysis report
# import sys
# path = 'adrenal_deg_results_p0.05_fdr0.1_summary.txt'
# sys.stdout = open(path, 'w')
#
# in_res_file = 'ref/Adrenal1vsAdrenal2_DEG.xls'
# res_df = pd.read_csv(in_res_file, sep='\t')
#
# print('Analysis report: Differentially expressed genes in adrenal glands')
# print('Test group: A1-A10 (Exposed to social stress)')
# print('Control group: A11-A20')
# print('Total number of genes that underwent DEG analysis:', len(res_df))
# print('Number of not protein coding genes removed:',
#         len(in_adrenal_df) - len(a_df))
# print('P-value cutoff: 0.05')
# print('Result:', len(a_df), 'DEGs, =',
#         len(a_df)/len(res_df) * 100, '% of input genes')
# print('Upregulated DEGs:', a_df_up, '=',
#         a_df_up/len(a_df) * 100, '% of initial DEGs')
# print('Downregulated DEGs:', a_df_down, '=',
#         a_df_down/len(a_df) * 100, '% of initial DEGs')
# print('Number of gene name duplicates:', duplicates_original)
# print('')
# print('Adjusted p-value cutoff: 0.1')
# print('Result after adjusted:', len(a_df_padj), 'DEGs, =',
#         len(a_df_padj)/len(res_df) * 100, '% of input genes')
# print('Upregulated adjusted DEGs:', a_df_padj_up, '=',
#         a_df_padj_up/len(a_df) * 100, '% of initial DEGs')
# print('Downregulated adjusted DEGs:', a_df_padj_down, '=',
#         a_df_padj_down/len(a_df) * 100, '% of initial DEGs')
# print('Number of gene name duplicates:', duplicates_adjusted)
# sys.stdout = sys.__stdout__
# #end filtering_deg.py
