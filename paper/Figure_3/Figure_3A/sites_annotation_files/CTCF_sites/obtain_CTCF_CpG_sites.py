# -*- coding: utf-8 -*-
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             04/04/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Obtain the sites annotation file for the CpG sites that overlap with in vivo CTCF ###
##### binding sites in human (hg38). The CpGs selected have the following properties:  ####
#####      - They overlap with an experimental ChIP-seq peak (see upregulated and      ####
#####        reactivated sites in Maurao et al., 2015).                                ####
#####      - They are found in a highly-scored sequence according to the JASPAR model  ####
#####        for CTCF.                                                                 ####
##### Based on http://biopython-cn.readthedocs.io/zh_CN/latest/en/chr14.html           ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

#### Dependencies ####

import os
import Bio
from Bio import motifs
from Bio.Seq import Seq 
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from pyfaidx import Fasta
import numpy as np
import pandas as pd

os.chdir('/Users/dem44/Documents/Manuscripts/cuRRBS/Figures/Figure_3/3X/CTCF_sites/')


### 1. Model for CTCF.

print('Creating model for CTCF ...')

## Load Jaspar PFM for CTCF in human. 
#http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0139.1.pfm

ctcf_pfm = motifs.read(open("MA0139.1_edited.pfm"),"pfm")

## Create the Position-Specific Scoring matrix
# Take into account that the median GC content of human genome is 41 %

ctcf_pwm = ctcf_pfm.counts.normalize(pseudocounts={'A':0.59, 'C': 0.41, 'G': 0.41, 'T': 0.59})

background = {'A':0.295,'C':0.205,'G':0.205,'T':0.295}
ctcf_pssm = ctcf_pwm.log_odds(background)

distribution = ctcf_pssm.distribution(background=background, precision=10**4)
threshold = distribution.threshold_fnr(0.01) # False negative rate of 1 %



### 2. Find the coordinates for the CpG sites. 

print('Finding coordinates for CpG sites ...')


## Read the human genome (hg38)

path_to_genome = "/Users/dem44/Desktop/methylation_clock/genomes/human_hg38/hg38.fa"
genome = Fasta(path_to_genome, sequence_always_upper=True, as_raw=True)


## Initialise useful variables

CpG_chromosomes = []
CpG_coordinates = [] # 1-based


## Read the ChIP-Seq peaks for CTCF and find the CpGs

peaks_file = open('CTCF_upregulated_and_reactived_sites_hg38.bed', 'r')
all_lines = peaks_file.readlines()

for i in range(0,len(all_lines)):
    
    #print i
    
    peak = all_lines[i]    
    
    # Variables that define the peak    
    
    chr_peak = peak.split('\t')[0]
    start_peak = int(peak.split('\t')[1])
    end_peak = int(peak.split('\t')[2])
    sequence_peak =  Seq(str(genome[chr_peak][start_peak:end_peak]), ctcf_pssm.alphabet) # 0-based coordinates in BED file
    
    if len(sequence_peak) < 19: 
        continue    
    
    # Calculate Score for each one of the positions in the peak sequence
    
    peak_scores = ctcf_pssm.calculate(sequence_peak).tolist()
    
    if not isinstance(peak_scores, list):
        peak_scores = [peak_scores]
    
    # Find the TF binding sites in the current peak (there may be several TFBS under the same peak)
    
    tfbs_starts = [] # 0-based coordinates    
    
    for j in range(0,len(peak_scores)):
        
        if peak_scores[j] > threshold: 
            
            tfbs_starts.append(start_peak + j)
            
    # Store the coordinates if we find CpGs in the TFBS
    
    for tfbs in tfbs_starts:
        
        site_a = str(genome[chr_peak][tfbs+4:tfbs+6]) # Potential first CpG in the CTCF consensus core sequence
        site_b = str(genome[chr_peak][tfbs+14:tfbs+16]) # Potential second CpG in the CTCF consensus core sequence
                
        if site_a == 'CG':
            #print('A')
            CpG_chromosomes.append(chr_peak)
            CpG_coordinates.append(tfbs + 5) # 1-based
        
        if site_b == 'CG':
            #print('B')
            CpG_chromosomes.append(chr_peak)
            CpG_coordinates.append(tfbs + 15) # 1-based
    
peaks_file.close()



### 3. Create the final dataframe to export (sites annotation file for cuRRBS).

print('Creating sites annotation file ...')

## Eliminate possible CpG sites which are duplicated

sites_df = pd.DataFrame(data={'Chr':CpG_chromosomes,
                              'Coordinate':CpG_coordinates})
                              
sites_df = sites_df[['Chr', 'Coordinate']] # Order columns
sites_df_uniq =  sites_df.drop_duplicates()

## Create IDs for the sites

site_ids = []

for i in range(1,sites_df_uniq.shape[0]+1):
    site_ids.append('S' + "%04d" % i)    

## Create final dataframe

final_sites_df = pd.DataFrame(data={'Site_ID':site_ids,
                                    'Chr': list(sites_df_uniq['Chr'].values),
                                    'Coordinate': list(sites_df_uniq['Coordinate'].values),
                                    'Weight': np.repeat(1, sites_df_uniq.shape[0])})
final_sites_df = final_sites_df[['Site_ID','Chr', 'Coordinate','Weight']] # Order columns

final_sites_df.to_csv(path_or_buf='CTCF_human_hg38_sites_annotation.csv', 
                      index=False)


### End of the script ###

