#!/bin/bash
#SBATCH --job-name=bismark   
#SBATCH --account=nn4605k
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=15G --partition=bigmem
#SBATCH --cpus-per-task=16

#Genome preparation and index building
bismark_genome_preparation --bowtie2 ref/

#Alignment and methylation calling
bismark --temp_dir out_temp/ -o out_aligned/ -q --bowtie2 \
	--parallel 14 --genome ref/ \
	-1 clean_data/A_19/trim_2/A_19_1_trimmed2.fq.gz, \
	clean_data/A_2/trim_2/A_2_1_trimmed2.fq.gz \
	-2 clean_data/A_19/trim_2/A_19_2_trimmed2.fq.gz, \
	clean_data/A_2/trim_2/A_2_2_trimmed2.fq.gz

#Deduplication
deduplicate_bismark -p --bam --output_dir out_dedup \
	--outfile A_19_bismark out_aligned/A_19_1_trimmed2_bismark_bt2_pe.bam
deduplicate_bismark -p --bam --output_dir out_dedup \
	--outfile A_2_bismark out_aligned/A_2_1_trimmed2_bismark_bt2_pe.bam


#Extract chromosomes
samtools view -b -L ref/rn6_1to20_XYM.bed \
	out_dedup/A_19_bismark.deduplicated.bam \
	> out_dedup/A_19_bismark_deduplicated_subset.bam
samtools view -b -L ref/rn6_1to20_XYM.bed \
	out_dedup/A_2_bismark.deduplicated.bam \
	> out_dedup/A_2_bismark_deduplicated_subset.bam

#Extract context of methylation
bismark_methylation_extractor --comprehensive --parallel 14 --ignore_r2 2 \
        --ignore_3prime_r2 2 --bedGraph --counts --gzip -p --no_overlap \
        --report out_dedup/A_19_bismark_deduplicated_subset.bam \
	--output out_met_extr
bismark_methylation_extractor --comprehensive --parallel 14 --ignore_r2 2 \
        --ignore_3prime_r2 2 --bedGraph --counts --gzip -p --no_overlap \
        --report out_dedup/A_2_bismark_deduplicated_subset.bam \
	--output out_met_extr

#coverage2cytosine to recieve genome-wide cytosine report
coverage2cytosine --split_by_chromosome --coverage_threshold 1 --gzip \
        --dir out_met_count/A_19/coverage2cytosine --genome_folder ref \
        -o A19_cov2cytosine.txt \
	out_met_extr/A_19/A_19_bismark_deduplicated_subset.bismark.cov.gz

coverage2cytosine --split_by_chromosome --coverage_threshold 1 --gzip \
        --dir out_met_count/A_2/coverage2cytosine --genome_folder ref \
        -o A2_cov2cytosine.txt \
	out_met_extr/A_2/A_2_bismark_deduplicated_subset.bismark.cov.gz
