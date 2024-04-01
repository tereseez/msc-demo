#!/usr/bin/env python
#start retrieve_DMRs_genenames.py 
import pandas as pd

#repeat for all files with annotated DMRs 
in_file = 'DMRs/DMRs_5dist_up.bed'
in_df = pd.read_csv(in_file, sep='\t', header=None)

gene_comb = (in_df.iloc[:,-1].apply(lambda x: x.split('||')[2])\
             .apply(lambda x: x.split(':')))
gene = gene_comb.apply(lambda x:x[0])
out_file = 'DMRs/DMRs_5dist_up_genename.tsv'
gene.to_csv(out_file, index=False, header=None)

#end retrieve_DMRs_genenames.py 
