#!/bin/bash
#SBATCH --job-name=quality_control
#SBATCH --account=nn4605k
#SBATCH --time=9:00:00
#SBATCH --mem-per-cpu=15G --partition=bigmem
#SBATCH --cpus-per-task=16

sample='A1'
group='test'

#quality check
fastqc in_raw/${sample}*.fq.gz

#Cutting 10 bases off of the beginning of the reads
cutadapt -u 10 -o in_clean/${group}/trimmed/${sample}_1_trimmed.fq.gz \
	in_raw/${sample}_1.fq.gz
cutadapt -u 10 -o in_clean/${group}/trimmed/${sample}_2_trimmed.fq.gz \
	in_raw/${sample}_2.fq.gz

#Quality check of trimmed reads
fastqc in_clean/${group}/trimmed/${sample}*trimmed.fq.gz

#unzip fastq files prior to FASTx-toolkit
gunzip in_clean/${group}/trimmed/${sample}*trimmed.fq.gz

#Quality filter (using FASTx-toolkit)
fastq_quality_filter -Q 33 -v -q 5 -p 50 \
	-i in_clean/${group}/trimmed/${sample}_1_trimmed.fq.gz \
	-o in_clean/${group}/filter/${sample}_1_qf.fq
fastq_quality_filter -Q 33 -v -q 5 -p 50 \
	-i in_clean/${group}/trimmed/${sample}_2_trimmed.fq.gz \
	-o in_clean/${group}/filter/${sample}_2_qf.fq

#Artifact filter (using FASTX-toolkit)
fastx_artifacts_filter -Q 33 -v \
	-i in_clean/${group}/filter/${sample}_1_qf.fq \
	-o in_clean/${group}/filter/${sample}_1_af.fq
fastx_arfifacts_filter -Q 33 -v \
	-i in_clean/${group}/filter/${sample}_2_qf.fq \
	-o in_clean/${group}/filter/${sample}_2_af.fq

#Pairing of filtered reads
fastq_pair in_clean/${group}/filter/${sample}_1_af.fq \
	in_clean/${group}/filter/${sample}_2_af.fq
#single.fq files are discarded

#Quality check of paired reads
fastqc ${sample}_*_af.fq.paired.fq
