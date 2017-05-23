# -*- coding: utf-8 -*-
###########################################################################################
#########      cuRRBS: customised Reduced Representation Bisulfite Sequencing     #########
###########################################################################################
#
# Created by Daniel E. Martin-Herranz, Antonio J.M. Ribeiro and Thomas M. Stubbs.
#
# Copyright (C) 2016,2017 D.E. Martin-Herranz, A.J.M. Ribeiro, T.M. Stubbs.
#
# This file is part of cuRRBS.
#
# cuRRBS is free software: you can redistribute it and/or modify it under the terms 
# of the GNU General Public License as published by the Free Software Foundation, either 
# version 3 of the License, or (at your option) any later version.
#
# cuRRBS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with cuRRBS. 
# If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################################
#
# DESCRIPTION OF THE SCRIPT:
#
# This script creates the pre-computed digestions for each one of the enzymes in
# enzymes_to_pre_compute.txt. Each one of the precomputed files contains the 
# cleavage sites of the correspoding enzyme for genome.fa. One cleavage site is 
# stored per line as an integer, with the different cleavage sites per chromosome 
# separated by lines starting with '>name_of_chr' (e.g. '>chr1').
#
# e.g.
#
# >chr1
# 1
# 456
# 12789
# >chr2
# 1
# 7895
# 46783
#
# The files with the pre-computed digestions will be used afterwards to speed up
# the rest of the pipeline by removing redundancy in the computations for the 
# different enzyme combinations.
#
###########################################################################################
#
# USAGE:
#
# python pre_compute_digestions.py enzymes_to_pre_compute.txt genome.fa working_directory
#
# COMPULSORY PARAMETERS:
#
#    enzymes_to_pre_compute.txt: absolute path to the text file which contains 
#                                the enzymes for which the digestions will be 
#                                pre-computed. Each enzyme is stored in one line.
#                                
#                                e.g. 
#
#                                AanI
#                                AarI
#                                AasI
#                                Acc36I
#                                Acc65I
#                                
#                                
#    genome.fa: absolute path to the FASTA file which contains the genome to be used
#               for the in silico digestions.
#               
#    working_directory: absolute path to the working directory. The software will 
#                       create a directory called 'pre_computed_files' inside the 
#                       working directory where the pre-computed files will be stored.
#                      
#                       e.g. if the working directory is '~/Desktop/work_dir/' the 
#                           files will be stored in '~/Desktop/work_dir/pre_computed_files/' 
#    
###########################################################################################

#### Dependencies ####

import os

import sys

from pyfaidx import Fasta

import Bio
from Bio import Restriction
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


#### Fix the output directory ####


if(str(sys.argv[3])[-1] == '/'):
    
    work_dir = str(sys.argv[3])

else:
    
    work_dir = str(sys.argv[3]) + '/'

if(os.path.isdir(work_dir) == False):
    
    print('\nThe directory ' + work_dir + ' does not exist. Please provide a correct ' +
    'working directory.\n')
    sys.exit()

else:
    
    output_dir = work_dir + 'pre_computed_files'

    if(os.path.isdir(output_dir)):
        
        print('\nThe output directory ' + output_dir + ' already exists.' + 
        ' Please delete the old "pre_computed_files" directory and restart the software.\n')
        sys.exit()
        
    else:
        
        print('\nReady to pre-compute digestions. The output will be generated in ' +
        output_dir + '/ ...\n')
        os.makedirs(output_dir)
        os.chdir(output_dir)


#### Store the enzymes to index (pre-compute) in a list.

enzymes_to_index = []

file_enzymes = open(str(sys.argv[1]), 'rU')

for line in file_enzymes:
        
    enzymes_to_index.append(line)

file_enzymes.close()

enzymes_to_index = [enzyme.replace('\n', '') for enzyme in enzymes_to_index]
enzymes_to_index = filter(None, enzymes_to_index) # Filter out empty elements derived from additional newlines

print('Enzymes to pre-compute: ' + str(enzymes_to_index) + '\n')
print('This computations can take several hours to finish in an average computer.' + '\n')

#### Load the genome.

path_to_genome = str(sys.argv[2])  # It should be the absolute path

genome = Fasta(path_to_genome, sequence_always_upper=True, as_raw=True)
                          
             
#### Create the index files for the enzymes considered.


for chromosome in genome.keys():
    
    print('Pre-computing ' + chromosome + ' ...')

    sequence = str(genome[chromosome])
    final_sequence = FormattedSeq(Seq(sequence, IUPACAmbiguousDNA()), linear=True)    
    
    for enzyme in enzymes_to_index:
        
        #print('Digesting for enzyme ' + enzyme + ' ...')

        batch = RestrictionBatch([enzyme])
        analysis_digestion = Analysis(batch, final_sequence)
        result = analysis_digestion.with_sites() # Obtain only enzymes that cut at least once
        
        #print('Obtaining restriction sites ...')
        
        restriction_sites_raw = [site for sites in result.values() for site in sites]
        restriction_sites_raw.append(1) # Add beginning
        restriction_sites_raw.append(len(final_sequence)+1) # Add end
        restriction_sites_raw.sort()
        restriction_sites_raw = map(str, restriction_sites_raw) # Convert to characters
        
        #print('Writing file ...')
        
        file_name = enzyme + '_pre_computed.txt'
        
        with open(file_name, "a") as myfile:
            
            myfile.write('>' + chromosome + '\n')
            myfile.write('\n'.join(restriction_sites_raw))
            myfile.write('\n')
            

print('\nThe pre-computing finished correctly.\n')

### End of the script ###
        
        
