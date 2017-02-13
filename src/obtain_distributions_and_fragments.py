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
# This script calculates the fragment length distributions (fl_distributions_x_enzymes.txt) and 
# the location of the fragments that contain the sites of interest (fragments_of_interest_x_enzymes.txt)
# for a list of combinations of 'x' enzymes (x=1; x=2). This information will be used by the cuRRBS
#  pipeline afterwards to calculate the Score, NF/1000 and EV values.
#
# The output files contain the following information:
#
# -  fl_distributions_x_enzymes.txt. It stores the theoretical distribution of fragment lengths
#    for the in silico digestion of a given enzyme / enzyme combination. The distributions for
#    the different enzymes are separated by lines starting with '>'. The first column stores the length
#    of the fragment and the second column the number of fragments with that length.
#   
#    e.g. for x=2
#   
#    >AanI_AarI
#    1,222
#    2,148
#    3,188
#    4,191
#    5,148
#    ...
#    3070220,1
#    3102079,1
#    3108192,1
#    >AanI_AasI
#    1,82
#    6,1433
#    7,1577
#    8,1576   
#
# -  fragments_of_interest_x_enzymes.txt. It stores the information necessary to locate
#    the fragments in which the sites of interest are located. Each line contains the 
#    IDs for the sites of interest (as specified in the annotation file and separated 
#    by '_') that are located in a given fragment, followed by the start coordinate and 
#    the length of the fragment.
#   
#    e.g. for x=1
#   
#    >AanI  
#    S001_S002_S003,72169288,3496
#    S004,72172805,2513
#    S005_S006,72199187,5415
#    ...
#    S545,47350313,20194 
#    >AarI
#    S001,72169345,784
#    ...
#
###########################################################################################
#
# USAGE:
#
# python obtain_distributions_and_fragments.py list_of_enzymes.txt sites_annotation.csv path_to_precomputed_files working_directory 
#
# COMPULSARY PARAMETERS:
#
#    list_of_enzymes.txt: absolute path to the text file which contains the enzymes (or enzyme combinations) 
#                         for which the fragment length distributions and the information
#                         for the fragments will be computed. Each enzyme (or enzyme
#                         combination) is stored in one line.
#                                
#                         e.g. for x=1 (digestions of individual enzymes)
#                         
#                         AanI
#                         AarI
#                         AasI
#                         Acc36I
#                         Acc65I
#                                
#                         e.g. for x=2 (digestions for combinations of 2 enzymes)
#                        
#                         AanI,AarI
#                         AanI,AasI
#                         AanI,Acc36I
#
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
#                          
#    path_to_precomputed_files: absolute path to the folder where the files with 
#                               the pre-computed digestions are located.
#                               
#    working_directory: absolute path to the working directory.
#                          
###########################################################################################


#### Dependencies ####

import os
import sys
import csv
import numpy
import bisect
import collections

# working directory
os.chdir(str(sys.argv[4]))

# path to precomputed digestions
index_files_path = str(sys.argv[3])

def find_fragment(cleaving_sites, coordinate):    
    index_coord = bisect.bisect(cleaving_sites, coordinate)- 1 #!!!WHY IS THIS NECESSARY
    if index_coord == len(cleaving_sites):
        index_coord -= 1    # this does not give the corrent end and lenght    
    start = cleaving_sites[index_coord]
    end = cleaving_sites[index_coord + 1] - 1
    length = end - start + 1 
    return (start,end,length)
        
#### Read the input data ####

## Store the enzymes (or enzyme combinations) in list, where each element is a 
#  different enzyme or enzyme combination.

with open(sys.argv[1], 'rU') as enzyme_file:
    enzyme_file_lines = [ line.rstrip() for line in enzyme_file if line.rstrip()]

## Annotation fron sites of interest.
csv_reader = csv.reader(open(sys.argv[2]))
next(csv_reader, None)
site_chr_coord_tuples = [(row[0],row[1],int(row[2])) for row in csv_reader]

individual_enzymes = set([enzyme for line in enzyme_file_lines 
                      for enzyme in line.split(",")])

#read all necessary enzyme files
# return dict: all_data = {'BsaWI': [arr_chr1, arr_chr2]}
def read_precomputed(individual_enzymes):
    all_data = {}
    for enzyme in individual_enzymes:
        sys.stdout.flush()
        print("        Reading precomputed file for {} ...".format(enzyme))
        precomputed_filename = index_files_path + enzyme + '_pre_computed.txt'
        all_data[enzyme] = []
        chr_names = []
        with open(precomputed_filename, 'r') as precomputed_file:
            precomputed_file_str = precomputed_file.read().replace('\n',' ')
            string_by_chr = precomputed_file_str.split('>chr')
            for ls in string_by_chr:
                if ls:
                    ls_list = ls.split()
                    chr_name = "chr{}".format(ls_list[0])
                    all_data[enzyme].append(
                        numpy.array(ls_list[1:], dtype=int))
                    chr_names.append(chr_name)
    return chr_names, all_data

chr_names, all_enzyme_restriction_sites = read_precomputed(individual_enzymes)

for no, enzyme_file_line in enumerate(enzyme_file_lines):
    sys.stdout.flush()
    print('        Processing enzyme (or enzyme combination) {}/{} ...'.format(
        no+1, len(enzyme_file_lines)))        
    chosen_enzymes =  enzyme_file_line.split(',')
    
    all_fragment_sizes = []
    fragment_output_lines = []
    
    fragments_sites_dict = collections.OrderedDict()
    for no_chr, chromosome in enumerate(chr_names):
        chr_restriction_sites = numpy.sort(numpy.concatenate(
            [all_enzyme_restriction_sites[enzyme][no_chr] for
             enzyme in chosen_enzymes]))
        fragment_sizes = numpy.diff(chr_restriction_sites)
        all_fragment_sizes.append(fragment_sizes)
        
        for site_chr_coord in site_chr_coord_tuples:
            if site_chr_coord[1] == chromosome:
                fragment_info = find_fragment(chr_restriction_sites,
                                              site_chr_coord[2])
                if fragment_info in fragments_sites_dict:
                    fragments_sites_dict[fragment_info].append(site_chr_coord[0])
                else:
                    fragments_sites_dict[fragment_info] = [site_chr_coord[0]]
                
    all_fragment_sizes_final = numpy.sort(numpy.concatenate(all_fragment_sizes))
    #counts_fragment_sizes, edges = numpy.histogram(
    #    numpy.array(all_fragment_sizes_final), bins=10, range=(20,30))
    counts_fragment_sizes = numpy.bincount(all_fragment_sizes_final)
    #print(counts_fragment_sizes)

    output_name = 'fl_distributions_{}_enzymes.txt'.format(len(chosen_enzymes))
    fsd_file = open(output_name, 'a') 
    fsd_file.write('>' + '_'.join(map(str,chosen_enzymes)) + '\n')
    writer = csv.writer(fsd_file)
    #for key, value in counts_fragment_sizes.items():
    fg_min = 0
    fg_max = 1000
    for no, count in enumerate(counts_fragment_sizes[fg_min:fg_max+1]):
        writer.writerow([fg_min+no, count])
    fsd_file.close()
    
    output_name_sites = 'fragments_of_interest_{}_enzymes.txt'.format(len(chosen_enzymes))
    
    sites_file = open(output_name_sites, 'a') 
    sites_file.write('>' + '_'.join(map(str,chosen_enzymes)) + '\n')
    for fragment_sites in fragments_sites_dict.items():
        sites_file.write('_'.join(fragment_sites[1]) + ",{0[0]},{0[2]}\n".format(fragment_sites[0]))
    sites_file.close()
    
