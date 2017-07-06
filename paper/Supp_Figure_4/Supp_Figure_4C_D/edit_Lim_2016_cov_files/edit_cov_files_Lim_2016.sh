#!/bin/bash

#FILES=/Users/dem44/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Figure_4/Figure_4new/mouse_adipogenesis/Lim_2016_cov_files/SRR*

FILES=/Users/dem44/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Figure_4/Figure_4new/mouse_adipogenesis/Lim_2016_cov_files/all_FastQ_merged_R1_val_1_GRCm38_bismark_bt2_pe.bismark.cov

for f in $FILES
do

	echo "Processing file $f ..."

	cp $f ${f}_mod

	# Change chromosome names from ENSEMBL to UCSC convention

	while read nom1 nom2; do
  		echo $nom1
		sed -i.bak "s/^$nom1/$nom2/" ${f}_mod
	done <GRCm38_ensembl2UCSC.txt

	rm ${f}_mod.bak


	# Delete unused chromosomes

	while read chrd; do
        	echo $chrd
        	sed -i.bak "/^$chrd/d" ${f}_mod
	done <chr_to_delete.txt

	rm ${f}_mod.bak

done
