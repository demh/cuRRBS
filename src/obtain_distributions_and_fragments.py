# -*- coding: utf-8 -*-
"""
-----------------------------
Daniel Elias Martin Herranz
22/09/16
EMBL-EBI
Thornton group
-----------------------------

This script calculates the fragment length distributions (fl_distributions_x_enzymes.txt) and 
the location of the fragments that contain the sites of interest (fragments_of_interest_x_enzymes.txt)
for a list of combinations of 'x' enzymes (x=1; x=2; for x=3 we would need long computational
times due to the huge amount of 3-enzyme combinations). This information will be used 
by the pipeline afterwards to calculate the (Score/NF) ratios for different size ranges
(see score_NF_calculations.R).

The output files contain the following information:

-  fl_distributions_x_enzymes.txt. It stores the theoretical distribution of fragment lengths
   for the in silico digestion of a given enzyme / enzyme combination. The distributions for
   the different enzymes are separated by lines starting with '>'. The first column stores the length
   of the fragment and the second column the number of fragments with that length.
   
   e.g. for x=2
   
   >AanI_AarI
   1,222
   2,148
   3,188
   4,191
   5,148
   ...
   3070220,1
   3102079,1
   3108192,1
   >AanI_AasI
   1,82
   6,1433
   7,1577
   8,1576   


-  fragments_of_interest_x_enzymes.txt. It stores the information necessary to locate
   the fragments in which the sites of interest are located. Each line contains the 
   IDs for the sites of interest (as specified in the annotation file and separated 
   by '_') that are located in a given fragment, followed by the start coordinate and 
   the length of the fragment.
   
   e.g. for x=1
   
   >AanI  
   S001_S002_S003,72169288,3496
   S004,72172805,2513
   S005_S006,72199187,5415
   ...
   S545,47350313,20194 
   >AarI
   S001,72169345,784
   ...
   
-----------------------------

Usage of the script:

python obtain_distributions_and_fragments.py list_of_enzymes.txt sites_annotation.csv path_to_precomputed_files working_directory 

-----------------------------

Parameters:

    list_of_enzymes.txt: absolute path to the text file which contains the enzymes (or enzyme combinations) 
                         for which the fragment length distributions and the information
                         for the fragments will be computed. Each enzyme (or enzyme
                         combination) is stored in one line.
                                
                         e.g. for x=1 (digestions of individual enzymes)
                         
                         AanI
                         AarI
                         AasI
                         Acc36I
                         Acc65I
                                
                         e.g. for x=2 (digestions for combinations of 2 enzymes)
                        
                         AanI,AarI
                         AanI,AasI
                         AanI,Acc36I

                                
    sites_annotation.csv: absolute path to the text file (with header) which contains the necessary information regarding 
                          the sites of interest, including:
                          
                          - Site_ID: unique ID given to the site of interest. It can
                                     not contain the '_' character.
                          - Chr: chromosome in which the site of interest is located
                          - Coordinate: coordinate in which the site of interest is located
                          - Weight: in case there is a preference in recovering certain sites
                                    of interest, their weights should be higher. The weights
                                    are always positive and are used afterwards to calculate 
                                    the Score.
                          
                          e.g. for 600 sites of interest
                          
                          Site_ID,Chr,Coordinate,Weight
                          S001,chr15,74203013,0.12933661
                          S002,chr17,40885350,0.002008493
                          S003,chr10,35605502,0.111945219
                          S004,chr17,68307004,0.005017857
                          ...
                          S600,chr21,32413126,0.869124446
                          
                          
    path_to_precomputed_files: absolute path to the folder where the files with 
                               the pre-computed digestions are located. It has to end in '/'.
                               
    working_directory: absolute path to the working directory.
                          
                          
    
----------------------------- 
"""


#### Dependencies ####

import os

import sys

import csv

import numpy

from collections import Counter
from collections import OrderedDict

import pandas as pd


#### Fix the appropiate paths ####

# Fix the working directory

os.chdir(str(sys.argv[4]))

# Select the path where the files with the pre-computed (indexed) digestions are.

index_files_path = str(sys.argv[3])


#### Some useful functions ####

## Function: given a list with all the restriction cleaving sites for a chromosome and the
# coordinate of the site of interest (located in that chromosome), return the start, end (base inclusive) and 
# length of the fragment that contains the site of interest.
    
    
def find_fragment(cleaving_sites, coordinate):    
   
    ## Create copy of cleaving sites

    cs = list(cleaving_sites[:])
    
    ## Store the coordinate in the cleaving sites list and order it.
    
    cs.append(coordinate)
    cs.sort()
    
    ## Find the index of the place where coord has been inserted.
    
    index_coord = cs.index(coordinate)
    
    ## Extract the start, end and length of fragment
    
    if cs[index_coord] != cs[index_coord+1]: # If the site of interest does not overlap with a cleaving site
    
        index_coord -= 1        
        
    start = cleaving_sites[index_coord]
    end = cleaving_sites[index_coord + 1] - 1
    length = end - start + 1 
    
    result = [start,end,length]
    
    ## Return the values obtained.
    
    return(result)
    
        
## Function: given the chromosome(e.g. 'chr1') and the coordinate (e.g. 78283402), 
#        return the ID for the site of interest (e.g. 'S045') from the annotation 
#        dataframe.
        
def find_ID(chr_site, coord_site, annotation):
    
    chr_cond = annotation['Chr'] == chr_site
    coord_cond = annotation['Coordinate'] == coord_site
    selected_row = annotation[chr_cond & coord_cond]
    result = list(selected_row['Site_ID'])[0]
    return(result)


## Create a Counter class that remembers the order in which the elements are added.
    
class OrderedCounter(Counter, OrderedDict):
    'Counter that remembers the order elements are first encountered'

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)



#### Read the input data ####

## Store the enzymes (or enzyme combinations) in list, where each element is a 
#  different enzyme or enzyme combination.

all_enzymes = []

file_enzymes = open(str(sys.argv[1]), 'rU')
#file_enzymes = open('list_of_enzymes_1.txt', 'rU')

for line in file_enzymes:
    
    all_enzymes.append(line)

file_enzymes.close()

all_enzymes = [line.replace('\n', '') for line in all_enzymes]
all_enzymes = filter(None, all_enzymes) # Filter out empty elements derived from additional newlines



## Annotation fron sites of interest.

sites_annotation = pd.read_csv(str(sys.argv[2]))
#sites_annotation = pd.read_csv('human_clock_sites_annotation.csv')

chr_of_interest = list(sites_annotation.ix[:,1])
sites_of_interest = list(sites_annotation.ix[:,2])



#### We will loop over each one of the elements of all_enzymes (i.e. the 
#    different enzymes or enzyme combinations) ####

print_count = 1

for element in all_enzymes:

    print('Processing enzyme (or enzyme combination) ' + str(print_count) + '/'+ 
    str(len(all_enzymes)) + ' ...')        
    
    chosen_enzymes =  element.split(',')
    number_enzymes = len(chosen_enzymes)
    
    #print('Generating output files for enzyme(s): ' + str(element))
    #print('\n')
    
    
    ##### 1. a) Create a list which contains all the sizes from the different fragments
    #     generated by the digestion with the selected restriction enzymes. 
    #
    #     b) Obtain the start, end and length of the fragment where each one of the 
    #     sites of interest are located.             
    
    
    ### Read the data from the index files.
    
    #print('Reading data from pre-computed digestions files ...')
    
    data_all_enzymes = []  # List of lists containing the lines in the index file for each enzyme considered.
    
    for enzyme in chosen_enzymes:
        
        index_file = open(index_files_path + enzyme + '_pre_computed.txt', 'r')
        data_all_enzymes.append(index_file.read().splitlines()) # Avoid newline characters
        index_file.close()
    
    
    # For each list (enzyme) in data_all_enzymes, find where the delimiters for chromosomes are.
    
    index_for_delimiters = []  # List of lists
    
    for enzyme in range(0, len(chosen_enzymes)):
        
        index_for_delimiters.append([i for i, s in enumerate(data_all_enzymes[enzyme]) if '>' in s])
    
    
    # Extract the chromosome names.
    
    chr_names = []
    
    for i in index_for_delimiters[0]:
        
        chr_names.append(data_all_enzymes[0][i].replace('>', ''))
    
    
    
    ### For a): Initialize a list which will contain all the fragment sizes.
    
    all_fragment_sizes = []
    
    
    ### For b): 
    
    # Initialize lists for the ordered sites of interest, start, end and length of the fragments where the
    #  sites of interest are located for b)
    
    ordered_site_interest = []  # Ordered according to for loop that iterates over chromosomes
    start_interest = []
    end_interest = []
    length_interest = []
    
    
    # Pair up the chromosomes and sites of interest and order them according to
    # chromosome.
    
    pairs_chr_sites = zip(chr_of_interest, sites_of_interest)
    
    
    
    ### Iterate over the different chromosomes to do a) and b).
    
    chr_index = 0
    
    for chromosome in chr_names:
        
        #print('Extracting fragment sizes for ' + chromosome + ' ...')
    
        
        # Initialize the list to store al the restriction sites for all enzymes in the current chromosome.
        
        restriction_sites_raw = []
        
        
        # Extract all the restriction sites in the current chromosome for all enzymes    
        
        for enzyme in range(0,len(chosen_enzymes)):
            
            start_chunk = index_for_delimiters[enzyme][chr_index] + 1 
            
            if chromosome == chr_names[len(chr_names)-1]:  # If we are looking at the last chromosome        
                
                end_chunk = len(data_all_enzymes[enzyme])
            
            else:
                
                end_chunk = index_for_delimiters[enzyme][chr_index+1]
            
            
            restriction_sites_raw.append(data_all_enzymes[enzyme][start_chunk:end_chunk])
            
        
        # Format the restriction sites.
        
        restriction_sites_final = list(set([int(site) for enzyme_sites in restriction_sites_raw for site in enzyme_sites]))
        restriction_sites_final.sort()
               
        
        
        # For a): Calculate the lengths of the fragments generated by the cuts.
    
        #print('Calculating fragment length distribution ...')
    
        fragment_sizes = list(numpy.diff(restriction_sites_final))
        
        # For a): Append the fragment sizes in the initialized list.
        
        all_fragment_sizes.append(fragment_sizes)
        
        #print('Number of fragments found for ' + chromosome + 
        #': ' + str(len(fragment_sizes)))
        
        
        
        
        # For b): Find those sites of interest located in the current chromosome
        
        #print('Calculating fragment lengths for the sites of interest ...')
        #print('\n')
        
        for pair in pairs_chr_sites:
            
            if pair[0] == chromosome:  # For the found sites, obtain start, end, length
            
                ordered_site_interest.append(pair)            
                fragment_info = find_fragment(restriction_sites_final, pair[1])
                start_interest.append(fragment_info[0])
                end_interest.append(fragment_info[1])
                length_interest.append(fragment_info[2])
                
        
        # Update the chromosome index
    
        chr_index += 1
            
    
    ## Avoid nesting in the final list with fragment sizes.
    
    all_fragment_sizes_final = [f for fs in all_fragment_sizes for f in fs]
    
    
    
    ##### 2. Create a dictionary which has the different fragment sizes as keys and
    #        the number of fragments of that size as values. Export it as a CSV file,
    #        with 'key,value' in each line.
    
    
    #print('Creating fragment size distribution ...')
    #print('\n')
    
    ## Order the final list with the fragment sizes and perform the counting,
    #  obtaining the fragment size distribution.
    
    all_fragment_sizes_final.sort()
    counts_fragment_sizes = OrderedCounter(all_fragment_sizes_final)
    
    
    ## Append in output file ()
    
    output_name = 'fl_distributions_' + str(number_enzymes) + '_enzymes.txt'
    
    
    fsd_file = open(output_name, 'a') 
    fsd_file.write('>' + '_'.join(map(str,chosen_enzymes)) + '\n')
    writer = csv.writer(fsd_file)
    for key, value in counts_fragment_sizes.items():
        writer.writerow([key, value])
    fsd_file.close()
        
    
    
    ##### 3. Create a dataframe that contains the following columns: chromosome,
    #        coordinate, fragment start, fragment end, fragment size. Export this 
    #        information as a CSV file.
    
    
    #print('Extracting information from sites of interest ...')
    #print('\n')
    
    ## Create the dataframe.
    
    dataframe_interest = pd.DataFrame(0, index=range(1, len(ordered_site_interest)+1), 
                                      columns=['Site_IDs', 'Chromosome', 'Coordinate', 'Fragment_start', 'Fragment_end', 'Fragment_length'])
    
    
    # Extract information from chromosomes and coordinates.    
    
    extract_chromosomes = [site[0] for site in ordered_site_interest]
    extract_coordinates = [site[1] for site in ordered_site_interest]
    
    
    # Find Site IDs.
    
    
    extract_site_IDs = []

    for i in range(0,len(extract_chromosomes)):
        
        site_ID = find_ID(extract_chromosomes[i], extract_coordinates[i], sites_annotation)
        extract_site_IDs.append(site_ID)
    
    
    dataframe_interest['Site_IDs'] = extract_site_IDs    
    dataframe_interest['Chromosome'] = extract_chromosomes
    dataframe_interest['Coordinate'] = extract_coordinates
    dataframe_interest['Fragment_start'] = start_interest
    dataframe_interest['Fragment_end'] = end_interest
    dataframe_interest['Fragment_length'] = length_interest
    
    
    # Collapse the dataframe for those site IDs which are located in the same fragment.
    
    grouped_dataframe = dataframe_interest.groupby(['Chromosome', 'Fragment_start', 'Fragment_end'], sort=False)
    
    collapsed_dataframe =   pd.DataFrame(0, index=range(1, len(grouped_dataframe)+1), 
                                      columns=['Site_IDs', 'Fragment_start', 'Fragment_length'])
                                                                      
    group_count = 1
    
    for key,group in grouped_dataframe:
        
        collapsed_dataframe.ix[group_count, 0] = '_'.join(list(group['Site_IDs'])) # Store IDs
        collapsed_dataframe.ix[group_count, 1] = key[1] # Store fragment start
        collapsed_dataframe.ix[group_count, 2] = key[2] - key[1] + 1  # Store fragment length
        group_count += 1
        
        
    ## Export the collapsed dataframe.
    
    output_name_sites = 'fragments_of_interest_' + str(number_enzymes) + '_enzymes.txt'
    
    sites_file = open(output_name_sites, 'a') 
    sites_file.write('>' + '_'.join(map(str,chosen_enzymes)) + '\n')
    sites_file.close()
    collapsed_dataframe.to_csv(output_name_sites, mode='a', header=False, index=False)
    
    print_count += 1

#print('The script finished correctly')

### End of the script ###

