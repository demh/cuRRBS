# -*- coding: utf-8 -*-
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             18/05/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Obtain the sites annotation file for the CpG sites that overlap with CpG islands ####
##### in the human genome (hg38). This will be used afterwards to estimate how well    ####
##### cuRRBS prediction performs in the XmaI-RRBS dataset.                             ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

##### Dependencies #####

import os
import pybedtools
from pybedtools import BedTool
from pyfaidx import Fasta 
import re
import itertools
import pandas as pd
import numpy as np

##### Input arguments #####

os.chdir("/Users/dem44/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Figure_4/Figure_4A/CG_sites/")
path_to_genome = "/Users/dem44/Desktop/methylation_clock/genomes/human_hg38/hg38.fa"
path_to_CGI_ann = "CGI_mod.gtf"
output_path = "CpG_in_CGI_annotation.csv"


##### 1. Find the coordinates for all the CpGs that are in CpG Islands in the human genome #####

print('Obtaining CpG coordinates ...')

# CGI annotation.

CGI_ann = BedTool(path_to_CGI_ann).sort()

# Human genome.

genome = Fasta(path_to_genome, sequence_always_upper=True, as_raw=True)

# Initilise variables

CG_chromosomes = []
CG_coordinates = []  # 1-based


# Find the CpGs in each CGI. Report only the coordinate for the C (1-based)

for CGI in CGI_ann:
    
    # Extract raw sequence for the CGI.    
    
    current_chr = str(CGI[0])
    start = int(CGI[3]) - 1 # 0-based
    end = int(CGI[4])
    
    CGI_seq = str(genome[current_chr][start:end])
    
    # Obtain the positions in the sequence where the CpG sites are located
    
    pos = [m.start() for m in re.finditer('CG', CGI_seq)]
    
    # Obtain the 1-based coordinates for the CpGs in the genome.
    
    final_coords = [coord + start + 1 for coord in pos]
    
    # Store the results
    
    CG_chromosomes.extend(list(itertools.repeat(current_chr, len(final_coords))))
    CG_coordinates.extend(final_coords)


##### 2. Create the sites annotation file with the correspondent sites #####

print('Creating sites annotation file ...')

## Eliminate possible CpG sites which are duplicated

sites_df = pd.DataFrame(data={'Chr':CG_chromosomes,
                              'Coordinate':CG_coordinates})
                              
sites_df = sites_df[['Chr', 'Coordinate']] # Order columns
sites_df_uniq =  sites_df.drop_duplicates()

## Create IDs for the sites

site_ids = []

for i in range(1,sites_df_uniq.shape[0]+1):
    site_ids.append('S' + "%07d" % i) 
    
## Create final dataframe

final_sites_df = pd.DataFrame(data={'Site_ID':site_ids,
                                    'Chr': list(sites_df_uniq['Chr'].values),
                                    'Coordinate': list(sites_df_uniq['Coordinate'].values),
                                    'Weight': np.repeat(1, sites_df_uniq.shape[0])})
final_sites_df = final_sites_df[['Site_ID','Chr', 'Coordinate','Weight']] # Order columns

final_sites_df.to_csv(path_or_buf=output_path, 
                      index=False)
                      

##### End of the script #####


