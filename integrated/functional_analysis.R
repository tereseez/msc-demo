#!/usr/bin/env Rscript

#start functional_analysis.R

library(clusterProfiler)
organism = "org.Rn.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
library(dplyr)

#repeat for the four overlap files
gene_list = read.csv("overlap_DEG_DMR/overlap_DET_p0.05_DMR_TSS.tsv",  
		     sep='\t', header = FALSE)

gene <- gene_list$V1

ego <- enrichGO(gene = gene,
                OrgDb = organism,
                keyType = "SYMBOL",
                ont = "ALL",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.2,
                minGSSize = 10,
                maxGSSize = 500,
		readable = FALSE,
                pool = FALSE)

write.table(ego, 
	    file="functional_analysis/overlap_DET_p0.05_TSS_enrichGO.tsv", 
            quote=FALSE, sep='\t', col.names = NA)

#produce dotplot of enrichment results
require(DOSE)
dotplot(ego,  split="ONTOLOGY", showCategory = 10, 
        title="Adrenal DETs (p<0.05) associated with TSS DMRs", 
        font.size=12) + 
  facet_grid(ONTOLOGY ~ ., scales="free")

#end functional_analysis.R
