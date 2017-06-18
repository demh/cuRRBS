#!/bin/bash

cp Homo_sapiens.GRCh38.87.gtf Homo_sapiens.GRCh38.87_mod.gtf


# Change chromosome names from ENSEMBL to UCSC sonvention

while read nom1 nom2; do
  	echo $nom1
	sed -i.bak "s/^$nom1/$nom2/" Homo_sapiens.GRCh38.87_mod.gtf
done <chr_names_Ensembl2UCSC.txt

rm Homo_sapiens.GRCh38.87_mod.gtf.bak


# Delete unused contigs

while read chrd; do
	echo $chrd
        sed -i.bak "/^$chrd/d" Homo_sapiens.GRCh38.87_mod.gtf
done <chr_to_delete.txt

rm Homo_sapiens.GRCh38.87_mod.gtf.bak

