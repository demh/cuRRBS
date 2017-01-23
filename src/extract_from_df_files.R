###########################################################################################
#########      cuRRBS: customized Reduced Representation Bisulfite Sequencing     #########
###########################################################################################
#
# Created by Daniel E. Martin-Herranz and Thomas Stubbs.
#
# Copyright (C) 2016,2017 D.E. Martin-Herranz, T. Stubbs.
# 
# This program is free software: you can redistribute it and/or modify it under the terms 
# of the GNU General Public License as published by the Free Software Foundation, either 
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. 
# If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################################
#
# DESCRIPTION OF THE SCRIPT: 
#
# Extract the desired enzymes (or enzyme combinations) from fl_distributions_x_enzymes.txt
# and fragments_of_interest_x_enzymes.txt files. The enzymes (or enzyme combinations) to 
# extract are provided in a list as input. The output is created in the same folder where 
# the original files are located.
#
###########################################################################################
#
# USAGE: Rscript extract_from_df_files.R --help
#
###########################################################################################


###########################################################
##################### Dependencies ########################
###########################################################

suppressWarnings(suppressMessages(library(optparse)));


###########################################################
############## Command-line arguments #####################
###########################################################

## Create arguments and help in a Pythonic style.

option_list <-  list(
  
  make_option(c('-f', '--fld_file'), type='character', default=NULL, 
              help="Absolute path to the fl_distribution_x_enzymes.txt file, which is created by \
                obtain_distributions_and_fragments.py. COMPULSORY",
              metavar='character'),
  
  make_option(c('-c', '--fragments_file'), type='character', default=NULL, 
              help="Absolute path to the fragments_of_interest_x_enzymes.txt, which is created by \
                obtain_distributions_and_fragments.py. COMPULSORY",
              metavar='character'),
  
  make_option(c('-e', '--enzyme_file'), type='character', default=NULL, 
              help="Absolute path to the file containing the enzymes (or enzyme combinations) to \
                extract. COMPULSORY",
              metavar='character'),
  
  make_option(c('-d', '--output_ID'), type='character', default=NULL, 
              help="Suffix used in the output file to allow parallelization. COMPULSORY",
              metavar='character')

);

## Prepare option list and general description of the script.

opt_parser <- OptionParser(option_list=option_list,
                           description="\nExtract the desired enzymes (or enzyme combinations) from fl_distributions_x_enzymes.txt \
and fragments_of_interest_x_enzymes.txt files. The enzymes (or enzyme combinations) \
to extract are provided in a list as input. The output is created in the same \
folder where the original files are located.");  


opt <- parse_args(opt_parser);


## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$fld_file) | is.null(opt$fragments_file) |
   is.null(opt$enzyme_file) | is.null(opt$output_ID)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Input arguments ##

path_to_enzyme_list <- as.character(opt$enzyme_file);
path_to_fld_file <- as.character(opt$fld_file);
path_to_fragments_file <- as.character(opt$fragments_file);
suffix <- as.character(opt$output_ID);


################################################################
################## Running the pipeline ########################
################################################################


#### 1. Read the input data. ####

## Enzymes (or enzyme combinations) to extract.

enz_to_extract <- readLines(path_to_enzyme_list);

## fld and fragment files

raw_fld <- readLines(path_to_fld_file);
raw_fragments <- readLines(path_to_fragments_file);


#### 2. Obtain the information regarding the different data chunks ####

fld_chunk_sep <- grep('>', raw_fld);
fragments_chunk_sep <- grep('>', raw_fragments);

fld_headers <- raw_fld[fld_chunk_sep];
fragments_headers <- raw_fragments[fragments_chunk_sep];

if(length(fld_headers) != length(fragments_headers)){
  
  stop("The input files do not have the same number of headers.");
  
}

fld_enz_combs <- gsub('>', '', fld_headers); # Enzyme combinations
fragments_enz_combs <- gsub('>', '', fragments_headers); # Enzyme combinations

if(!all(fld_enz_combs == fragments_enz_combs)){
  
  stop('The enzyme combinations are not the same for the fl_distribution_x_enzymes.txt and the fragments_of_interest_x_enzymes.txt files.');
  
}



#### 3. Create output files ####


for(enz in enz_to_extract){
  
  ## fld file
  
  index_fld <- which(fld_enz_combs == enz);
  start_fld <- fld_chunk_sep[index_fld];
  
  if(index_fld == length(fld_enz_combs)){
    
    end_fld <- length(raw_fld);
    
  }else{
    
    end_fld <- fld_chunk_sep[index_fld + 1] - 1;
    
  }
  
  fld_chunk <- raw_fld[start_fld:end_fld];
  
  write(fld_chunk, file=paste0(path_to_fld_file, '_', suffix), sep='\n', append=T);
  
  
  
  ## fragment file
  
  index_fragments <- which(fragments_enz_combs == enz);
  start_fragments <- fragments_chunk_sep[index_fragments];
  
  if(index_fragments == length(fragments_enz_combs)){
    
    end_fragments <- length(raw_fragments);
    
  }else{
    
    end_fragments <- fragments_chunk_sep[index_fragments + 1] - 1;
    
  }
  
  fragments_chunk <- raw_fragments[start_fragments:end_fragments];
  
  write(fragments_chunk, file=paste0(path_to_fragments_file, '_', suffix), sep='\n', append=T);

}

#########################################################
########## End of the script ############################
#########################################################
