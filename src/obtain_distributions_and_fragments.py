# -*- coding: utf-8 -*-

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
        print("Reading precomputed file for {}".format(enzyme))
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
    fg_min = 20
    fg_max = 800
    for no, count in enumerate(counts_fragment_sizes[fg_min:fg_max+1]):
        writer.writerow([fg_min+no, count])
    fsd_file.close()
    
    output_name_sites = 'fragments_of_interest_{}_enzymes.txt'.format(len(chosen_enzymes))
    
    sites_file = open(output_name_sites, 'a') 
    sites_file.write('>' + '_'.join(map(str,chosen_enzymes)) + '\n')
    for fragment_sites in fragments_sites_dict.items():
        sites_file.write('_'.join(fragment_sites[1]) + ",{0[0]},{0[2]}\n".format(fragment_sites[0]))
    sites_file.close()
    
