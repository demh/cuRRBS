# -*- coding: utf-8 -*-
###########################################################################################
#########      cuRRBS: customized Reduced Representation Bisulfite Sequencing     #########
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
# This script prints the size of the input genome in Mega-basepairs. These values are 
# included in the genome_sizes.txt file for the genomes with pre-computed files available.
#
###########################################################################################
#
# USAGE:
#
# python calculate_genome_size.py genome.fa
#
# COMPULSORY PARAMETERS:
#
#    genome.fa: absolute path to the FASTA file which contains the genome to be checked
#
###########################################################################################

#### Dependencies ####

import sys
from pyfaidx import Fasta

#### Read the genome of interest ####

path_to_genome = str(sys.argv[1])  # It should be the absolute path
genome = Fasta(path_to_genome, sequence_always_upper=True, as_raw=True)


#### Calculate the genome size in Mega-basepairs ####

genome_size_bp = 0

for chrom in genome.keys():
    genome_size_bp += len(genome[chrom])

genome_size_Mbp = genome_size_bp / 1000000.0 
   
print('Size of the input genome (in Mbp): ' + str(genome_size_Mbp))


#### End of the script ####