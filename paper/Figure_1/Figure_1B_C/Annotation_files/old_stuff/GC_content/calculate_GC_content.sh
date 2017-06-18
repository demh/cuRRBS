#!/bin/bash

python create_windows.py chr_lengths.txt 50

bedtools nuc -fi /nfs/nobackup/thornton/dem44/bismark_genomes/human_hg38/hg38.fa -bed windows_50_bp.gff | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $11, $7, $8, $9}'| awk '{if(NR>1)print}' > GC_content_hg38_50bp.gff
