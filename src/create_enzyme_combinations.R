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
# Create a list with the 2-enzyme combinations that will be checked in the main cuRRBS 
# pipeline.
#
###########################################################################################
#
# USAGE: Rscript create_enzyme_combinations.R --help
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
  
  make_option(c('-e', '--enzymes_to_check'), type='character', default=NULL, 
              help="Absolute path to the enzymes_to_check.txt file. COMPULSORY.",
              metavar='character'),
  
  make_option(c('-s', '--sne_files'), type='character', default=NULL, 
              help="Absolute path to the SNE file for the individual enzymes. COMPULSORY.",
              metavar='character'),
  
  make_option(c('-m', '--mode'), type='character', default=NULL, 
              help="Search mode ('e': exhaustive, 'h': heuristic). COMPULSORY.",
              metavar='character'),
  
  make_option(c('-p', '--output_path'), type='character', default=NULL,
              help="Absolute path for the output directory.  COMPULSORY",
              metavar='character')
);


opt_parser <- OptionParser(option_list=option_list,
                           description="\nCreate a list with the 2-enzyme combinations that will be checked in the main \
cuRRBS pipeline.");

opt <- parse_args(opt_parser);


## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$enzymes_to_check) | is.null(opt$sne_files) |
   is.null(opt$mode) | is.null(opt$output_path)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Input arguments.

enz_to_check_path <- as.character(opt$enzymes_to_check);
SNE_file_path <- as.character(opt$sne_files); # Only needed for heuristic mode
search_mode <- as.character(opt$mode); # Heuristic or exhaustive
output_path <- as.character(opt$output_path);  # Path where the output list will be created 

setwd(output_path);



################################################################
################## Running the pipeline ########################
################################################################


#### 1. Read in all the enzymes to check ####

all_enzymes <- readLines(enz_to_check_path);

if(sum(all_enzymes == "")>0){ # Remove empty lines
  
  all_enzymes <- all_enzymes[-(which(all_enzymes == ""))];
  
}


#### 2a. Create the 2-enzyme combinations for the exhaustive search mode. ####

if(search_mode == 'e'){
  
  if(length(all_enzymes) > 1){
    
    enz_comb <- combn(all_enzymes, 2, function(x){
      
      return(paste0(x[1], ',', x[2]));
      
    });
    
  }else{
    
    enz_comb <- c();
    
  }
}


#### 2b. Create the 2-enzymes combinations for the heuristic search mode. ####

if(search_mode == 'h'){
  
  ## Obtain the enzymes that will serve as seeds (those with Score > (C_Score * max_Score / h))
  
  raw_seed <- read.csv(SNE_file_path, header=F);
  seed_enzymes <- unique(as.character(raw_seed[,1]));
  
  ## Combine all enzymes to check with seeds to form combinations.
  
  all_comb <- expand.grid(seed_enzymes, all_enzymes);
  
  # Remove those combinations with a repeated enzyme.
  
  index_rm <- which(as.character(all_comb[,1]) == as.character(all_comb[,2]));
  all_comb <- all_comb[-(index_rm),];
  
  # Order alphabetically the enzymes in the combinations.
  
  ordered_comb <- apply(all_comb, 1, function(x){
    
    charc <- as.character(x);
    charc_ord <- charc[order(charc)];
    return(paste0(charc_ord[1], ',', charc_ord[2]));
    
  })
  
  # Remove duplicated enzyme combinations.
  
  enz_comb <- unique(ordered_comb);
  
}


#### 3. Generate the output. ####

if(length(enz_comb) > 0){
  
  cat(paste0('        There are ', length(enz_comb), ' enzyme combinations to check.'), sep='\n');
  write(enz_comb, file=paste0(output_path,'/combinations_to_check.txt'), sep='\n');
  
}



#########################################################
########## End of the script ############################
#########################################################


