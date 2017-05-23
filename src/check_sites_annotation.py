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
# This script checks if the coordinates reported in the site annotation file are correct.
# It performs the following checks:
#       - The first base must be a 'C'. If not, the site ID is reported.
#       - Counts each nucleotide in the second and third base and percentages are reported.
#
###########################################################################################
#
# USAGE:
#
# python check_sites_annotation.py sites_annotation.csv genome.fa 
#
# COMPULSORY PARAMETERS:
#
#    sites_annotation.csv: absolute path to the text file (with header) which contains the necessary information regarding 
#                          the sites of interest, including:
#                          
#                          - Site_ID: unique ID given to the site of interest. It can
#                                     not contain the '_' or the ',' characters.
#                          - Chr: chromosome in which the site of interest is located
#                          - Coordinate: 1-based genomic coordinate in which the site of interest is located
#                          - Weight: in case there is a preference in recovering certain sites
#                                    of interest, their weights should be higher. The weights
#                                    are always positive and are used afterwards to calculate 
#                                    the Score.
#                          
#                          e.g. for 600 sites of interest
#                          
#                          Site_ID,Chr,Coordinate,Weight
#                          S001,chr15,74203013,0.12933661
#                          S002,chr17,40885350,0.002008493
#                          S003,chr10,35605502,0.111945219
#                          S004,chr17,68307004,0.005017857
#                          ...
#                          S600,chr21,32413126,0.869124446
#
#    genome.fa: absolute path to the FASTA file which contains the genome to be used
#               for the in silico digestions.
#
###########################################################################################


### Dependencies ####

import sys
from pyfaidx import Fasta 


### Read the genome of interest ###

path_to_genome = str(sys.argv[2])  # It should be the absolute path
genome = Fasta(path_to_genome, sequence_always_upper=True, as_raw=True)


### Initialise useful variables ###

total_sites = 0
wrong_site_IDs = []

A2 = 0
C2 = 0
T2 = 0
G2 = 0
N2 = 0

A3 = 0
C3 = 0
T3 = 0
G3 = 0
N3 = 0


### Loop over the different sites of interest ###

input_sites = str(sys.argv[1])

with open(input_sites, 'r') as handle:
    
    for line in handle:
        
        # Extract information from current site of interest        
        
        ID = line.split(',')[0]
        
        if(ID=='Site_ID'): # Skip header
            continue
        
        chromosome = line.split(',')[1]
        start = int(line.split(',')[2])   # 1-based coordinate
        seq = str(genome[chromosome][start-1:start+2])
        total_sites += 1

        
        # Check first base
        
        if(list(seq)[0]!='C'):
            
            wrong_site_IDs.append(ID)
            
        # Check second base
            
        if(list(seq)[1]=='A'):
            A2 += 1
        elif(list(seq)[1]=='C'):
            C2 += 1
        elif(list(seq)[1]=='T'):
            T2 +=1
        elif(list(seq)[1]=='G'):
            G2 +=1
        elif(list(seq)[1]=='N'):
            N2 +=1
            
        # Check third base
        
        if(list(seq)[2]=='A'):
            A3 += 1
        elif(list(seq)[2]=='C'):
            C3 += 1
        elif(list(seq)[2]=='T'):
            T3 +=1
        elif(list(seq)[2]=='G'):
            G3 +=1
        elif(list(seq)[2]=='N'):
            N3 +=1
            
          
          
### Create final output ###

print('\nTotal number of sites of interest: ' + str(total_sites))
print('There are ' + str(len(wrong_site_IDs)) + " sites which do not have a C in their first position.")
if(len(wrong_site_IDs) >0):
    print('Their IDs are: ' + ','.join(wrong_site_IDs))
print('\nPercentages for second base:')
print('     A:' + str(A2*100/float(total_sites)))
print('     C:' + str(C2*100/float(total_sites)))
print('     G:' + str(G2*100/float(total_sites)))
print('     T:' + str(T2*100/float(total_sites)))
print('     N:' + str(N2*100/float(total_sites)) + '\n')
print('Percentages for third base:')
print('     A:' + str(A3*100/float(total_sites)))
print('     C:' + str(C3*100/float(total_sites)))
print('     G:' + str(G3*100/float(total_sites)))
print('     T:' + str(T3*100/float(total_sites)))
print('     N:' + str(N3*100/float(total_sites)) + '\n')

        
### End of the script ###        
