#!/bin/bash
#SBATCH --job-name=dmr_annotation   
#SBATCH --account=nn4605k
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=4G --partition=bigmem
#SBATCH --cpus-per-task=6

hmst_seq_analyzer gene_annotation \
                  -l 2000 -xL 5000 \
                  -X 5000 -Y 1000 \
                  -F annotation/ -hu r -n no -M 5000 \
                  -N 1000000 -r ref/refFlat.txt \
                  -g ref/rn6_sorted.txt

# annotate DMRs to TSS, repeat for 5'distance region (up and down)
bedtools intersect -a ../5mcseq/dmr_analysis/dmranalysis_filtered_06.bed \
                   -b annotation/data/\
                   TSS_Up5000_Down1000_removedShort.bed -wa -wb \
                   > DMRS/DMRs_TSS.bed



