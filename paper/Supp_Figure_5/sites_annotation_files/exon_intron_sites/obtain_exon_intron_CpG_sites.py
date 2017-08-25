# -*- coding: utf-8 -*-
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             06/04/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Obtain the sites annotation file for the CpG sites that are located in a         ####
##### exon-intron boundary in the human genome (hg38). We define the boundary as the   ####
##### genomic regions within plus or minus 5 bp from the 5' and 3' ends of the introns.####
##### Only exon-intron boundaries which contain canonical 5' and 3' splice sites are   ####
##### taken into account.                                                              ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

#### Dependencies ####

import os
import pybedtools
from pybedtools import BedTool
from pyfaidx import Fasta
import re
import numpy as np
import pandas as pd

os.chdir('/Users/dem44/Documents/Manuscripts/cuRRBS/Figures/Figure_3/3X/exon_intron_sites/')


#### 1. Obtain the coordinates for the 5' end and 3' end of all the introns in the
#       human genome (hg38).

print('Retrieving introns from the human genome ...')

path_to_gencode_ann = "/Users/dem44/Documents/Manuscripts/cuRRBS/Figures/Figure_1/1C/Annotation_files/gencode.v25.basic.annotation.gtf"
gencode_ann = BedTool(path_to_gencode_ann).sort()

protein_coding_genes_ann = gencode_ann.filter(lambda x: x[2] == 'gene').filter(lambda x: 'gene_type "protein_coding"' in x[8]).sort()
exon_protein_coding_ann = gencode_ann.filter(lambda x: x[2] == 'exon').filter(lambda x: 'gene_type "protein_coding"' in x[8]).sort()
intron_protein_coding_ann = protein_coding_genes_ann.subtract(exon_protein_coding_ann, s=True).sort()

five_prime_ends = [x.start for x in intron_protein_coding_ann] # 0-based coordinates
three_prime_ends = [(x.end - 1) for x in intron_protein_coding_ann] # 0-based coordinates
chromosomes = [str(x.chrom) for x in intron_protein_coding_ann]


#### 2. Find the coordinates of all the CpG sites that are found in +- 5 bp of both
#       of the ends of the intron (i.e. close to the exon-intron boundary).

print('Finding the CpG coordinates ...')

## Read the human genome (hg38)

path_to_genome = "/Users/dem44/Desktop/methylation_clock/genomes/human_hg38/hg38.fa"
genome = Fasta(path_to_genome, sequence_always_upper=True, as_raw=True)

## Find the CpGs in the exon-intron boundaries

CpG_chromosomes = []
CpG_coordinates = []

for i in range(0,len(five_prime_ends)):
    
    # Consider only canonical 5' and 3' splice sites    
    
    if (str(genome[chromosomes[i]][five_prime_ends[i]:five_prime_ends[i]+2]) == 'GT') and (str(genome[chromosomes[i]][(three_prime_ends[i] - 1):(three_prime_ends[i]+1)]) == 'AG'):
        
        #print('Consensus in ' + str(i))        
        
        # Find CpGs in 5' splice site
        
        start_five = five_prime_ends[i]-5 # 0-based        
        seq_five = str(genome[chromosomes[i]][start_five:start_five+11])
        find_five_CpG = [m.start() for m in re.finditer('CG', seq_five)]
        
        if len(find_five_CpG) > 0:
            #print("Found CpG in 5'")
            five_CpGs = [(start_five + cg + 1) for cg in find_five_CpG] # 1-based coordinates           
            CpG_coordinates.extend(five_CpGs)
            CpG_chromosomes.extend(list(np.repeat(chromosomes[i],len(find_five_CpG))))
        
        # Find CpGs in 3' splice site

        start_three = three_prime_ends[i] - 5 # 0-based       
        seq_three = str(genome[chromosomes[i]][start_three:start_three+11])
        find_three_CpG = [m.start() for m in re.finditer('CG', seq_three)]
        
        if len(find_three_CpG) > 0:
            #print("Found CpG in 3'")
            three_CpGs = [(start_three + cg + 1) for cg in find_three_CpG] # 1-based coordinates     
            CpG_coordinates.extend(three_CpGs)
            CpG_chromosomes.extend(list(np.repeat(chromosomes[i],len(find_three_CpG))))
            


#### 3. Create the final dataframe to export (sites annotation file for cuRRBS).

print('Creating sites annotation file ...')

## Eliminate possible CpG sites which are duplicated

sites_df = pd.DataFrame(data={'Chr':CpG_chromosomes,
                              'Coordinate':CpG_coordinates})
                              
sites_df = sites_df[['Chr', 'Coordinate']] # Order columns
sites_df_uniq =  sites_df.drop_duplicates()

## Create IDs for the sites

site_ids = []

for i in range(1,sites_df_uniq.shape[0]+1):
    site_ids.append('S' + "%05d" % i)    

## Create final dataframe

final_sites_df = pd.DataFrame(data={'Site_ID':site_ids,
                                    'Chr': list(sites_df_uniq['Chr'].values),
                                    'Coordinate': list(sites_df_uniq['Coordinate'].values),
                                    'Weight': np.repeat(1, sites_df_uniq.shape[0])})
final_sites_df = final_sites_df[['Site_ID','Chr', 'Coordinate','Weight']] # Order columns

final_sites_df.to_csv(path_or_buf='exon_intron_human_hg38_sites_annotation.csv', 
                      index=False)


### End of the script ###