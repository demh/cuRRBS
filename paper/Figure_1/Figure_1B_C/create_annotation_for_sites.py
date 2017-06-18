# -*- coding: utf-8 -*-
###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             01/03/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Create a CSV file which contains genomic information about the sites where the   ####
##### different restriction enzymes cut in a given genome. This information will be used ##
##### afterwards to create several figures (e.g. Fig. 1C).                             ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################


##### Dependencies #####

import os
import pybedtools
from pybedtools import BedTool
import numpy

##### Input arguments #####

#os.chdir("/Users/dem44/Documents/Manuscripts/cuRRBS/Figures/Figure_1/1C/")
#
#path_to_pre_computed_folder = "/Users/dem44/Documents/Manuscripts/cuRRBS/Figures/Figure_1/1C/Annotation_files/fake_precomputed"
#path_to_chr_lengths = "Annotation_files/chr_lengths_sorted.txt"
#path_to_genome = "/Users/dem44/Desktop/methylation_clock/genomes/human_hg38/hg38.fa"
#path_to_gencode_ann = "Annotation_files/gencode.v25.basic.annotation.gtf"
#path_to_CGI_ann = "Annotation_files/CGI_mod.gtf"
#output_path = "Annotation_files/genomic_annotation_for_cleavage_sites.csv"

os.chdir("/nfs/nobackup/thornton/dem44/methylation_clock/optimize_RRBS/Figures_paper/1C/")
pybedtools.set_tempdir('/nfs/nobackup/thornton/dem44/methylation_clock/optimize_RRBS/Figures_paper/1C/temp_files/')

path_to_pre_computed_folder = "/nfs/nobackup/thornton/dem44/methylation_clock/optimize_RRBS/pre_computed_files_human_hg38/"
path_to_chr_lengths = "Annotation_files/chr_lengths_sorted.txt"
path_to_genome = "/nfs/nobackup/thornton/dem44/bismark_genomes/human_hg38/hg38.fa"
path_to_gencode_ann = "Annotation_files/gencode.v25.basic.annotation.gtf"
path_to_CGI_ann = "Annotation_files/CGI_mod.gtf"
output_path = "/nfs/nobackup/thornton/dem44/methylation_clock/optimize_RRBS/Figures_paper/1C/genomic_annotation_for_cleavage_sites.csv"

##### Functions #####

#### Create a BedTool object (in GFF format) from a pre-computed file for a given restriction enzyme.

# path_to_pre_computed: absolute path to the corresponding pre-computed file

def pre_computed_to_gff(path_to_pre_computed):
    
    gff_string = ""
    
    # Read all the input    
    
    with open(path_to_pre_computed) as f:
        
        raw_in = f.readlines()
    
    raw_in = [line.strip() for line in raw_in]
    
    
    # Store information in string
    
    for line in raw_in:
        
        if '>' in line:
            
            current_chr = line.replace('>','')
            
        else:
            
            current_str = current_chr + '\t' + 'na' + '\t' + 'na' + '\t' + line + '\t' + line + '\t' + '.' + '\t' + '.' + '\t' + '.' +  '\t' + '.' + '\n'
            gff_string = gff_string + current_str
            
    
    # Create BedTool object and return it
    
    final_gff = BedTool(gff_string, from_string=True)
    
    return(final_gff)        
        
        

#### Calculate the percentage of cleavage sites that overlap with a certain genomic feature.

# enz_gff: BedTool object containing the cleavage sites for the current restriction enzyme
# n: number (integer) of total cleavage sites for the current restriction enzyme
# genomic_feature: (filtered) BedTool object containing the genomic intervals with the current genomic feature

def percentage_sites_in_genomic_feature(enz_gff, n, genomic_feature):
    
        non_overlapping_result = enz_gff.intersect(genomic_feature, v=True)
        n_non_overlapping = 0
        for i in non_overlapping_result: 
            n_non_overlapping +=1
            
        n_overlapping = n - n_non_overlapping
        
        return((float(n_overlapping)/float(n))*100)    



##### Calculate the genomic feature annotation for all the restriction enzymes #####

### 1. Initialize variables.

## Main variables.

enzymes_to_check = []
n_sites_final = []

mean_perc_GC_content = []
mean_perc_CpG_content = []

perc_protein_coding_genes_sites = []
perc_exon_protein_coding_sites = []
perc_intron_protein_coding_sites = []
perc_non_coding_RNA_genes_sites = []
perc_intragenic_sites = []
perc_intergenic_sites = []

perc_CGI_sites = []
perc_shore_sites = []
perc_shelf_sites = []

perc_promoter_CGI_sites = []
perc_promoter_non_CGI_sites = []


## Annotation files.

print('Calculating annotations ...')

gencode_ann = BedTool(path_to_gencode_ann).sort()

protein_coding_genes_ann = gencode_ann.filter(lambda x: x[2] == 'gene').filter(lambda x: 'gene_type "protein_coding"' in x[8]).sort()
exon_protein_coding_ann = gencode_ann.filter(lambda x: x[2] == 'exon').filter(lambda x: 'gene_type "protein_coding"' in x[8]).sort()
intron_protein_coding_ann = protein_coding_genes_ann.subtract(exon_protein_coding_ann, s=True).sort()
non_coding_RNA_genes_ann = gencode_ann.filter(lambda x: x[2] == 'gene').filter(lambda x: ('gene_type "Mt_rRNA"' in x[8]) or
                                                                                         ('gene_type "Mt_tRNA"' in x[8]) or
                                                                                         ('gene_type "miRNA"' in x[8]) or
                                                                                         ('gene_type "misc_RNA"' in x[8]) or
                                                                                         ('gene_type "rRNA"' in x[8]) or
                                                                                         ('gene_type "scRNA"' in x[8]) or
                                                                                         ('gene_type "snRNA"' in x[8]) or
                                                                                         ('gene_type "snoRNA"' in x[8]) or
                                                                                         ('gene_type "ribozyme"' in x[8]) or
                                                                                         ('gene_type "sRNA"' in x[8]) or
                                                                                         ('gene_type "scaRNA"' in x[8]) or 
                                                                                         ('gene_type "lincRNA"' in x[8])).sort()
intragenic_ann = protein_coding_genes_ann.slop(g=path_to_chr_lengths, b=2500).merge().sort()
intergenic_ann = intragenic_ann.complement(g=path_to_chr_lengths).sort()

CGI_ann = BedTool(path_to_CGI_ann).sort()
shore_ann = CGI_ann.flank(g=path_to_chr_lengths, b=2000).sort()
shelf_ann = CGI_ann.flank(g=path_to_chr_lengths, b=4000).subtract(shore_ann).sort()

protein_coding_transcripts_ann = gencode_ann.filter(lambda x: x[2] == 'transcript').filter(lambda x: 'transcript_type "protein_coding"' in x[8]).sort()
promoter_all_ann = protein_coding_transcripts_ann.flank(g=path_to_chr_lengths, l=2500, r=0, s=True).slop(g=path_to_chr_lengths, l=0, r=500, s=True).sort() # There are duplicates
promoter_CGI_ann = promoter_all_ann.intersect(CGI_ann, u=True).sort()
promoter_non_CGI_ann = promoter_all_ann.intersect(CGI_ann, v=True).sort()


## Paths to all the pre-computed files.

pre_computed_files = os.listdir(path_to_pre_computed_folder)
if '.DS_Store' in pre_computed_files: pre_computed_files.remove('.DS_Store')
pre_computed_files = [path_to_pre_computed_folder + '/' + f for f in pre_computed_files]


### 2. Loop over the different enzymes and calculate the values for the different genomic features.

for enz_path in pre_computed_files:

    ## Store the current enzyme
    
    enz = enz_path.split('/')[-1].split('_')[0]
    enzymes_to_check.append(enz)
    
    print('Processing enzyme ' + enz + ' ...')

    ## Read the pre-computed file for the current enzyme and convert it to GFF.

    enzyme_gff = pre_computed_to_gff(enz_path)

    ## A. Calculate the mean GC content of the windows overlapping with the cleavage sites.
    
    print('A ...')    
    
    GC_windows = enzyme_gff.slop(g=path_to_chr_lengths, b=25)
    GC_content = GC_windows.nucleotide_content(fi=path_to_genome)
    
    GC_content_values = []
    n_sites = 0
    
    for i in GC_content: 
        GC_content_values.append(float(i[10]))
        n_sites += 1
    
    mean_perc_GC_content.append(numpy.mean(GC_content_values)*100.0)
    n_sites_final.append(n_sites)
    
    
    ## B. Calculate the mean CpG content of the windows overlapping with the cleavage sites.
    
    print('B ...')    
    
    CpG_windows = enzyme_gff.slop(g=path_to_chr_lengths, b=500)
    CpG_content = CpG_windows.nucleotide_content(fi=path_to_genome, pattern='CG')
    
    CpG_content_values = []
    
    for i in CpG_content:
        CpG_content_values.append(int(i[18]))
    
    mean_perc_CpG_content.append(numpy.mean(CpG_content_values)*100/1000.0)    

    
    ## C. Calculate the % cleavage sites that overlap with different genomic features
    
    print('C ...')    
    
    perc_protein_coding_genes_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, protein_coding_genes_ann))
    perc_exon_protein_coding_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, exon_protein_coding_ann))
    perc_intron_protein_coding_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, intron_protein_coding_ann))
    perc_non_coding_RNA_genes_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, non_coding_RNA_genes_ann))
    perc_intragenic_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, intragenic_ann))
    perc_intergenic_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, intergenic_ann))

    perc_CGI_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, CGI_ann))
    perc_shore_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, shore_ann))
    perc_shelf_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, shelf_ann))
    
    perc_promoter_CGI_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, promoter_CGI_ann))
    perc_promoter_non_CGI_sites.append(percentage_sites_in_genomic_feature(enzyme_gff, n_sites, promoter_non_CGI_ann))


    
### 3. Create the output file, which will be a CSV file with the following columns ():
#      enzyme, n_sites, mean_%_GC_content, mean_%_CpG_content, %_sites_protein_coding_genes,   
#      %_sites_exons, %_sites_introns, %_sites_non_coding_RNA_genes, %_sites_intragenic,
#      %_sites_intergenic, %_sites_CGI, %_sites_shores, %_sites_shelves, %_sites_CGI_promoters,
#      %_sites_non_CGI_promoters    
    
print('Writing output file ...')

f_output = open(output_path, 'w')
    
f_output.write("enzyme,n_sites,mean_%_GC_content,mean_%_CpG_content,%_sites_protein_coding_genes," + 
"%_sites_exons,%_sites_introns,%_sites_non_coding_RNA_genes,%_sites_intragenic," +
"%_sites_intergenic,%_sites_CGI,%_sites_shores,%_sites_shelves,%_sites_CGI_promoters," +
"%_sites_non_CGI_promoters" + '\n')

for i in range(0, len(enzymes_to_check)):
    
    data_in_row = [enzymes_to_check[i], n_sites_final[i], mean_perc_GC_content[i], mean_perc_CpG_content[i], 
                   perc_protein_coding_genes_sites[i], perc_exon_protein_coding_sites[i],
                   perc_intron_protein_coding_sites[i], perc_non_coding_RNA_genes_sites[i],
                   perc_intragenic_sites[i], perc_intergenic_sites[i], perc_CGI_sites[i], 
                   perc_shore_sites[i], perc_shelf_sites[i], perc_promoter_CGI_sites[i],
                   perc_promoter_non_CGI_sites[i]]    
    data_in_row_str = [str(dat) for dat in data_in_row]
    row_string = ','.join(data_in_row_str)
    f_output.write(row_string + '\n')

f_output.close()


print('The script finished correctly.')


### End of the script ###
