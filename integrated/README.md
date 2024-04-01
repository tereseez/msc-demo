ORDER OF SCRIPTS:  
1. job_DMR_annotation.sh  
- define genomic regions by using hmst-seq-analyzer  
- annotate identified DMRs to genomic regions  

2. retrieve_DMRs_genenames.py  
- get names of genes associated with DMRs (TSS and 5' distance regions)  

3. job_merge_DMR_DEG.sh  
- retrieves overlap instances of DEGs (p<0.05 and fdr<0.1) and gene names associated with DMRs (TSS and 5dist)  

4. functional_analysis.R  
- performs ORA with clusterprofiler, export table and plot results  
