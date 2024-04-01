#!/bin/bash
#SBATCH --job-name=gene_exp_levels
#SBATCH --account=nn4605k
#SBATCH --time=9:00:00
#SBATCH --mem-per-cpu=15G --partition=bigmem
#SBATCH --cpus-per-task=16

sample='A1'
group='test'

#using HISAT2 alignments
featureCounts -a ref/genome.gtf \
	-o out_counts/hisat2/${group}/${sample}_h2_featureCounts.tsv \
	-t exon -g gene_name -p -Q 10 -B \
	-C out_mapped/hisat2/${sample}_h2align.sorted.bam -O

#using STAR alignment
featureCounts -a ref/genome.gtf \
	-o out_counts/star/${group}/${sample}_star_featureCounts.tsv \
	-t exon -g gene_name -p -Q 10 -B \
	-C out_mapped/star/${sample}_star_Aligned.sortedByCoord.out.bam -O
