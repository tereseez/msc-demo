ORDER OF SCRIPTS:  
1. job_quality_check.sh  
- uses fastqc to quality check the raw reads and trims 10 bases off from the end of the reads  

2. job_bismark.sh  
- aligns the reads to the reference genom  
- deduplicates the aligned reads  
- extracts the context of the methylation levels by using bismark  
- produces genome-wide cytosine report  

3. filter_met_counts.py  
- filter CpG methylation levels based on a minimum methylation count of 2 and total count of 5  
- calculates methylation percentage per CpG site  

4. overlap_CpG_sites.py  
- verifies CpG methylation levels computed by comparing identified CpG sites and total counts per CpG site  

5. mean_CpGsites.py  
- compute mean CpG methylation percentage and total count for samples in test group and in ctr group  

6. job_hmst_seq_analyzer.sh  
- runs hmst-seq-analyzer to identify DMRs  

7. pre_overlap_DMRs.py  
- preprocesses hmst-seq results: merge all DMRs across all chromosomes into one file  
- preprocesses dmr_analysis results: filters MRs based on different probability score cutoffs  
- compares lengths of DMRs identified with the two tools in a violin plot  

8. job_bedtools.sh
- use bedtools intersect to inspect DMR overlaps between the two tools


9. overlap_DMRs.py
- visualize percentage of DMR that overlaps between the two tools per dmr_analysis DMR probability score cutoff
