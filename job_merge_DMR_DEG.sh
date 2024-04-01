#!/bin/bash
#SBATCH --job-name=merge_dmr_deg   
#SBATCH --account=nn4605k
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=4G --partition=bigmem
#SBATCH --cpus-per-task=6

#keep unique occurences of TSS gene names
sort -u DMRs/DMRs_TSS_genename.tsv -o DMRs/DMRs_TSS_genename_unique.tsv

#sort and merge 5distance up/down unique genes into one file
sort -u DMRs/DMRs_5dist_down_genename.tsv DMRs/DMRs_5dist_up_genename.tsv \
     > DMRs/DMRs_5dist_int_genename_unique.tsv

#merge DMR and DET
comm -12 <(sort DMRs/DMRs_TSS_genename_unique.tsv) \
         <(sort DEGs/adrenal_deg_filtered_p0.05_unique_genelist.tsv) \
         > overlap_DEG_DMR/overlap_DET_p0.05_DMR_TSS.tsv

comm -12 <(sort DMRs/DMRs_5dist_int_genename_unique.tsv) \
         <(sort DEGs/adrenal_deg_filtered_p0.05_unique_genelist.tsv) \
         > overlap_DEG_DMR/overlap_DET_p0.05_DMR_5dist.tsv

comm -12 <(sort DMRs/DMRs_TSS_genename_unique.tsv) \
         <(sort DEGs/adrenal_deg_filtered_fdr0.1_unique_genelist.tsv) \
         > overlap_DEG_DMR/overlap_DET_fdr0.1_DMR_TSS.tsv

comm -12 <(sort DMRs/DMRs_5dist_int_genename_unique.tsv) \
         <(sort DEGs/adrenal_deg_filtered_fdr0.1_unique_genelist.tsv) \
         > overlap_DEG_DMR/overlap_DET_fdr0.1_DMR_5dist.tsv
