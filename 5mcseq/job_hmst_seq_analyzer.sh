#!/bin/bash
#SBATCH --job-name=hmst_seq_analyzer
#SBATCH --account=nn4605k
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=15G --partition=bigmem
#SBATCH --cpus-per-task=16

CHR='chr5'
echo ${CHR}

TEST_METHOD='Ttest'
echo ${TEST_METHOD}

hmst_seq_analyzer gene_annotation \
	-X 1000 -Y 1000 \
	-l 2000 -xL 5000 \
	-F hmst_seq_analyzer/${CHR} -hu r -n no \
	-M 10000 -N 100000 \
	-r ref/refFlat.txt -g ref/rn6_sorted.txt
echo gene_annotation-DONE

hmst_seq_analyzer data_preprocessing \
	-F hmst_seq_analyzer/${CHR} -z no -m yes \
	-hu r -n no -c 0.4999 \
	-fko out_met_count/run2_junbai/\
	CpG_context/mean/test.${CHR}.mean.csv \
	-fwt out_met_count/run2_junbai/\
	CpG_context/mean/control.${CHR}.mean.csv \
	-g ref/rn6_sorted.txt
echo data_preprocessing-DONE

hmst_seq_analyzer find_MRs \
	-W yes -a 200 -mc1 3 -mc2 5 -mc3 3 \
	-F hmst_seq_analyzer/${CHR} -p 5 \
	-fko hmst_seq_analyzer/${CHR}/list_mC_hmC_files_KO.txt \
	-fwt hmst_seq_analyzer/${CHR}/list_mC_hmC_files_WT.txt \
	-ref hmst_seq_analyzer/data/refFlat_clean_sorted.bed \
	-reg hmst_seq_analyzer/list_region_files.txt
echo find_MRs-DONE

hmst_seq_analyzer prepare_for_DMR_finding \
	-isST 1 -a 200 -mc1 3 -mc2 5 -mc3 3 \
	-F hmst_seq_analyzer/${CHR} -p 4 \
	-ko hmst_seq_analyzer/${CHR}/\
	list_all_filtered_formatted_MRs_KO.txt \
	-wt hmst_seq_analyzer/${CHR}/\
	list_all_filtered_formatted_MRs_WT.txt
echo prepare_for_DMR_finding-DONE

hmst_seq_analyzer DMR_search \
	-F hmst_seq_analyzer/${CHR} \
	-T ${TEST_METHOD} -p 5 \
	-f hmst_seq_analyzer/${CHR}/\
	list_prepared_for_DMR_finding_imputed.txt
echo DMR_search-DONE

hmst_seq_analyzer prep4plot \
	-F hmst_seq_analyzer/${CHR} -c 0.4999 \
	-ko hmst_seq_analyzer/${CHR}/\
	list_all_filtered_formatted_MRs_KO.txt \
	-wt hmst_seq_analyzer/${CHR}/\
	list_all_filtered_formatted_MRs_WT.txt
echo prep4plot-DONE

hmst_seq_analyzer plot_all \
	-F hmst_seq_analyzer/${CHR} \
	-reg hmst_seq_analyzer/${CHR}/\
	list_count_allMRs_regions_files.txt \
	-sit hmst_seq_analyzer/${CHR}/\
	list_count_sites_files.txt \
	-cmc hmst_seq_analyzer/${CHR}/\
	counts_DMR_hypo_hyper_imputed_${TEST_METHOD}_5mC.csv \
	-chmc hmst_seq_analyzer/${CHR}/\
	counts_DMR_hypo_hyper_imputed_${TEST_METHOD}_5hmC.csv \
	-aMR hmst_seq_analyzer/${CHR}/\
	list_TSS_genebody_TES_enhancer_allMRs.txt \
	-oMR hmst_seq_analyzer/${CHR}/\
	list_overlapping_MRs.txt
echo plot_all-DONE
