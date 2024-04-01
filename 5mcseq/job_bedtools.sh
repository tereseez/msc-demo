#!/bin/bash
#SBATCH --job-name=bedtools   
#SBATCH --account=nn4605k
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=4G --partition=bigmem
#SBATCH --cpus-per-task=6

cutoff='01'

#overlap hmst to dmr_analysis
bedtools intersect -a dmr_analysis/dmranalysis_filtered_${cutoff}.bed \
	-b DMR_overlap/hmstseq_DMRs_uniqueID_allchr.bed \
	-wao > DMR_overlap/overlap_${cutoff}.bed
#repeat for all filtered dmr_analysis files
