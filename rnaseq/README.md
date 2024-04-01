ORDER OF SCRIPTS:  
1. job_quality_control.sh  
- uses fastqc to quality check reads  
- trims off 10 bases from beginning of reads, also applies quality filter, artifact filter and pairing of reads  

2. job_alignments.sh  
- aligns the cleaned reads with 2 tools: HISAT2 and STAR  

3. job_gene_exp_levels.sh  
- uses featureCounts to summarize mapped reads  

4. filter_raw_counts.py  
- filter raw read counts to only keep protein coding instances  

5. similarity_rawcounts.py  
- compares raw (filtered) read counts: Novogene, HISAT2 and STAR. Produces scatterplot  

6. normalize_counts.py  
- normalize HISAT2 summarized read counts: quantile normalisation > RPKM > log2 > z-scores  
- plots in-house quant. norm. vs novogene norm. A1 read counts  
- plots distribution of in-house normalised A1 counts (quant. norm., log2(RPKM), zscores of log2(RPKM))  

7. deg_analysis.py  
- identifies DEGs by using t-test on two inputs:  
	1) Novogene FPKM values (are first log2-transformed and zscore-standarized)  
	2) subset (Z-scores of log2(RPKM))  
- compares 1) and 2) to Novogene final DEGs  

8. filtering_deg.py  
- filters the final DEG list, Novogene's final DEGs, for protein coding genes  
- produces clustermaps of DEGs  
