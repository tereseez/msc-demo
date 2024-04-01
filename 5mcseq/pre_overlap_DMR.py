#!/usr/bin/env python
#start pre_overlap_DMRs.py

import pandas as pd
from matplotlib import pyplot as plt
import glob
import seaborn as sb
import numpy as np

#HMST-seq-analyzer DMRs for all chromosomes combined to one file
hmst_files = glob.glob('hmst_seq_analyzer/chr*/data'\
                       '/test.chr*.mean_vs_control.chr*.mean_5mC'\
                       '_*_imputedWith_zeros_DMRs_all.csv')
all_chromosomes = pd.DataFrame() #empty dataframe

#access chr, start and end in result files and export to same df
for file in hmst_files:
    df = pd.read_csv(file, sep=',')
    name_comb = (df.iloc[:,0].apply(lambda x: x.split('||')[0])\
                 .apply(lambda x: x.split(':')))
    df['chr'] = name_comb.apply(lambda x:x[0])
    df['start'] = name_comb.apply(lambda x:x[1])
    df['end'] = name_comb.apply(lambda x:x[2])
    df['id'] = name_comb.apply(lambda x:x[-1])
    df = df.drop_duplicates(subset=['start', 'end', 'id'])
    all_chromosomes = pd.concat([all_chromosomes,df], axis=0)

#reformat and sort df
all_chromosomes = all_chromosomes.astype({'start' : int,'end' : int})
all_sorted = all_chromosomes.sort_values(by=['chr','start'],
                                         ignore_index=True)
all_sorted = all_sorted[['chr', 'start', 'end', 'id', 'rratio', 'pvals']]
#exporting DMRs for all chr
out_hmst = 'DMR_overlap/hmstseq_all_DMRs_uniqueID_allchr.bed' 
all_sorted.to_csv(out_hmst, sep='\t', index=False, header=None)

#DMR-analysis MRs
dmr_file = 'dmr_analysis/23_chroms_high_miniPercentChange_gt_0.0001'\
           '_Pcutoff_0.05_isSmooth_2_isModTest_0__all'\
           '_dmrRanking_top_0.56_minLogReg_proba_0.7.bed'
dmr_df = pd.read_csv(dmr_file, sep='\t', header=None,
                     names=["chr", "start", "end", "label", "score"])

#Filter DMR-analysis results based on different cutoffs
dmr_filtered = dmr_df[dmr_df["score"] > 0.5] #0.5, 0.6, 0.7, 0.8, 0.9
dmr_filtered.reset_index(inplace=True, drop=True)
dmr_filtered = dmr_filtered.sort_values(by=['chr', 'start'],
                                        ignore_index=True)
#change output_dmr according to score cutoff
out_dmr = 'dmr_analysis/dmranalysis_filtered_05.bed'
dmr_filtered.to_csv(out_dmr, sep='\t', index=False, header=None)

final_dmr_filtered = dmr_df[dmr_df["score"] > 0.5]
fina_dmr_filtered.reset_index(inplace=True, drop=True)
final_dmr_filtered = dmr_filtered.sort_values(by=['chr', 'start'],
                                        ignore_index=True)

#find lengths and log10-transform data for violin plot
hmst_l = np.log10(all_sorted["end"] - all_sorted["start"])
dmr_l = np.log10(final_dmr_filtered["end"] -
                 final_dmr_filtered["start"])
lengths = [hmst_l, dmr_l]

#Violinplot: differences in DMR lengths
ax = sb.violinplot(data=lengths, palette = ["#ed7d31", "darksalmon"])
ax.set_xticklabels(['HMST-Seq-Analyzer','DMR-analysis'])
ax.set_xlabel("Method")
ax.set_ylabel("Log10 DMRs (bp)")
plt.show()

#end pre_overlap_DMRs.py
