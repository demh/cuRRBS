#!/bin/bash

cp SRR2721365_all_merged_trimmed_GRCh38_bismark_bt2.bismark.cov SRR2721365_all_merged_trimmed_GRCh38_bismark_bt2.bismark_mod.cov

# Change chromosome names from ENSEMBL to UCSC sonvention

while read nom1 nom2; do
  	echo $nom1
	sed -i.bak "s/^$nom1/$nom2/" SRR2721365_all_merged_trimmed_GRCh38_bismark_bt2.bismark_mod.cov
done <chr_names_Ensembl2UCSC.txt

rm SRR2721365_all_merged_trimmed_GRCh38_bismark_bt2.bismark_mod.cov.bak


# Delete unused chromosomes

while read chrd; do
        echo $chrd
        sed -i.bak "/^$chrd/d" SRR2721365_all_merged_trimmed_GRCh38_bismark_bt2.bismark_mod.cov
done <chr_to_delete.txt

rm SRR2721365_all_merged_trimmed_GRCh38_bismark_bt2.bismark_mod.cov.bak
