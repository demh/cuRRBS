# -*- coding: utf-8 -*-
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             22/02/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Create a GFF file which contains all the possible non-overlapping windows of a   ####
##### certain window length in a given genome. The output file will be used afterwards ####
##### to calculate the GC content of these windows.                                    ####
###########################################################################################
##### USAGE: python create_windows.py chr_lengths.txt window_length                    ####
#####                                                                                  ####
##### ARGUMENTS:                                                                       ####
#####              chr_lengths.txt: text file containing the chromosome names and      ####
#####                               their lengths. It can be created with the following ###
#####                               command:                                           ####
#####                               faidx ~/Desktop/methylation_clock/genomes/human_hg38/hg38.fa -i chromsizes > chr_lengths.txt
#####              window_length: integer specifying the length of the windows (in bp) ####
###########################################################################################


##### Dependencies #####

import os
import sys
import numpy
from collections import OrderedDict
import pandas


os.chdir(os.path.dirname(os.path.realpath(sys.argv[0])))


##### 1. Read the input from arguments #####

## Chromosome lengths

#with open("chr_lengths.txt") as f:
with open(str(sys.argv[1])) as f:
    raw_chr_l = f.readlines()
raw_chr_l = [x.strip() for x in raw_chr_l] 

chromosomes = [line.split('\t')[0] for line in raw_chr_l]
chr_lengths = [line.split('\t')[1] for line in raw_chr_l]


## Window length

window_length = int(sys.argv[2])
#window_length = 50

 
##### 2. Create GFF file with the genomic intervals for the required windows #####
 
output_f = open('windows_' + str(window_length) + '_bp.gff', 'a')

for i in range(0, len(chromosomes)):
    
    c = chromosomes[i]
    l = int(chr_lengths[i])
    
    print('Creating windows for chromosome ' + c + ' ...')
    
    starts = numpy.arange(1,l+1,window_length)
    ends = starts + window_length - 1
    
    if(ends[-1] > l): ends[-1] = l
    
    nwin = len(starts)
    
    od=OrderedDict([('seqname',[c]*nwin), ('source',['.']*nwin), ('feature',['.']*nwin), 
    ('start',starts), ('end',ends), ('score',[0]*nwin), ('strand',['.']*nwin) , ('frame',['.']*nwin), 
    ('attribute',['.']*nwin)])
    current_df = pandas.DataFrame(data=od)
    
    current_df.to_csv(path_or_buf=output_f, sep='\t', header=False, index=False,
                      mode='a')    
    
output_f.close()

print('The script finished OK.')

### End of the script ###    
