#!/usr/bin/env python
#start overlap_DMRs.py
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import glob
import numpy as np

hmst_file = "DMR_overlap/hmstseq_DMRs_uniqueID_allchr.bed"
dmr_files = glob.glob("dmr_analysis/dmranalysis_filtered*")
dmr_files = sorted(dmr_files)
overlap_files = glob.glob("DMR_overlap/overlap_*.bed")
overlap_files = sorted(overlap_files)

percentages = []
cutoffs = np.arange(0.1, 1.0, 0.1)

for file_overlap, file_dmr, i in zip(overlap_files, dmr_files, cutoffs):
    df_hmst = pd.read_csv(hmst_file, sep='\t', header=None)
    hmst_unique_dmr = len(df_hmst)

    df_dmr = pd.read_csv(file_dmr, sep='\t', header=None)
    dmr_analysis_unique_dmr = len(df_dmr)

    df_overlap = pd.read_csv(file_overlap, sep='\t', header=None)
    overlap = len(df_overlap[df_overlap[11]!=0])

    total = dmr_analysis_unique_dmr
    percentages.append(overlap/total * 100)

    #make venn diagram
    v = venn2(subsets = (hmst_unique_dmr, dmr_analysis_unique_dmr, 
        overlap),
        set_labels = ('HMST-Seq (5+5)', 'DMR analysis (5+5)'),
        set_colors=('#ed7d31', 'darksalmon'), alpha = 0.5,
        subset_label_formatter = lambda x: str(x) +
        "\n(" + f"{(x/total):1.0%}" + ")")
    v.get_label_by_id('10').set_text(len(df_hmst))
    v.get_label_by_id('01').set_text(len(df_dmr))
    plt.title("Overlap DMRs (cutoff " + f"{round(i, 1)}" + ")")
    plt.show()

#plot percentage overlaps vs. cutoff values
plt.plot(cutoffs, percentages, "-o", color="#ed7d31", alpha=0.5)
plt.xlabel("DMR-analysis: DMR probability score cutoffs")
plt.ylabel("Overlap (%)")
plt.title("Overlap DMRs: DMR-analysis vs HMST-Seq")
plt.xticks(round(cutoffs,1), cutoffs)
plt.show()

#end overlap_DMRs.py
