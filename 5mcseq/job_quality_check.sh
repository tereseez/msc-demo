#!/bin/bash
#SBATCH --job-name=quality_control
#SBATCH --account=nn4605k
#SBATCH --time=9:00:00
#SBATCH --mem-per-cpu=15G --partition=bigmem
#SBATCH --cpus-per-task=12

sample='A_19'

#quality check
fastqc raw_data/*/clean_data/${sample}/*fq.gz

#cutting off 10 bases from the end of the reads
cutadapt -u -10 -o clean_data/${sample}/trim_1/${sample}_1_trimmed.fq.gz \
	raw_data/*/clean_data/${sample}/${sample}_1.clean.fq.gz
cutadapt -u -10 -o clean_data/${sample}/trim_1/${sample}_2_trimmed.fq.gz \
	raw_data/*/clean_data/${sample}/${sample}_2_clean.fq.gz

fastqc clean_data/${sample}/trim_1/*

cutadapt -u -10 -o clean_data/${sample}/trim_2/${sample}_1_trimmed2.fq.gz \
	clean_data/${sample}/trim_1/${sample}_1_trimmed.fq.gz
cutadapt -u -10 -o clean_data/${sample}/trim_2/${sample}_2_trimmed2.fq.gz \
	clean_data/${sample}/trim_1/${sample}_2_trimmed.fq.gz

fastqc clean_data/${sample}/trim_2/*
