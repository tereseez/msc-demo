#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH --account=nn4605k
#SBATCH --time=9:00:00
#SBATCH --mem-per-cpu=15G --partition=bigmem
#SBATCH --cpus-per-task=16

sample='A1'
group='test'

#USING HISAT2 FOR MAPPING
#Extracting splice-sites and exons
hisat2_extract_splice_sites.py ref/genome.gtf > ref/genome.ss
hisat2_extract_exons.py ref/genome.gtf > ref/genome.exon

#Building index
hisat2-build --ss ref/genome.ss --exon ref/genome.exon ref/genome.fa \
	ref/index_hisat2

#Mapping the reads to the genome
hisat2 -x ref/index_hisat2 -1 in_clean/${group}/paired/${sample}_1_af.paired.fq \
	-2 in_clean/${group}/paired/${sample}_2_af.paired.fq \
	-S out_mapped/hisat2/${sample}_h2align.sam

#Converting SAM to BAM to sorted BAM
samtools view -b out_mapped/hisat2/${sample}_h2align.sam > \
	out_mapped/hisat2/${sample}_h2align.bam
samtools sort out_mapped/hisat2/${sample}_h2align.bam \
	-o out_mapped/hisat2/${sample}_h2align.sorted.bam


#USING STAR FOR MAPPING
#Building index
STAR --runMode genomeGenerate --genomeDir ref/index_star/ --genomeFastaFiles \
     ref/genome.fa --sjdbGTFfile ref/genome.gtf

#Mapping
STAR --genomeDir ref/index_star/ --readFilesIn \
	in_clean/${group}/paired/${sample}_1_af.paired.fq \
	in_clean/${group}/paired/${sample}_2_af.paired.fq \
	--outFileNamePrefix out_mapped/star/star_ \
	--outSAMtype BAM SortedByCoordinate
