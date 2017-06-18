#!/bin/bash

WORKING_DIR='/nfs/nobackup/thornton/dem44/methylation_clock/optimize_RRBS/Figures_paper/SuppFig2A/'
PRECOMP_PATH='/nfs/nobackup/thornton/dem44/methylation_clock/optimize_RRBS/pre_computed_files_human_hg38'

while read line; do
  
	echo $line > ${WORKING_DIR}enzyme_tmp.txt
	
	python obtain_distributions_and_fragments_old.py ${WORKING_DIR}enzyme_tmp.txt ${WORKING_DIR}human_clock_sites_annotation.csv \
	${PRECOMP_PATH} ${WORKING_DIR}

	Rscript Score_NF_tradeoff_plots.R -f fl_distributions_*.txt -c fragments_of_interest_*.txt -l ${WORKING_DIR}human_clock_sites_annotation.csv \
	-p ${WORKING_DIR} -r 75 -z FALSE

	rm ${WORKING_DIR}enzyme_tmp.txt
	rm fl_distributions_*.txt
	rm fragments_of_interest_*.txt

done <combinations_to_check.txt
