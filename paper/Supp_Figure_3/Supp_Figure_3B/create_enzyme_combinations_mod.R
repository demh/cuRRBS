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
# Create a list with the enzymes and the 2-enzyme combinations that will be checked in the main cuRRBS 
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
  
  make_option(c('-p', '--output_path'), type='character', default=NULL,
              help="Absolute path for the output directory.  COMPULSORY",
              metavar='character')
);


opt_parser <- OptionParser(option_list=option_list,
                           description="\nCreate a list with the enzymes and the 2-enzyme combinations that will be checked in the main \
cuRRBS pipeline.");

opt <- parse_args(opt_parser);


## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$enzymes_to_check) | is.null(opt$output_path)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Input arguments.

enz_to_check_path <- as.character(opt$enzymes_to_check);
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


#### 2 Create the 2-enzyme combinations ####

if(length(all_enzymes) > 1){
    
  enz_comb <- combn(all_enzymes, 2, function(x){
      
    return(paste0(x[1], ',', x[2]));
      
  });

}else{
    
  enz_comb <- c();
    
}


#### 3. Generate the output. ####

enz_comb <- c(all_enzymes,enz_comb);  # Append individual enzymes and 2-enzyme combinations

if(length(enz_comb) > 0){
  
  cat(paste0('        There are ', length(enz_comb), ' enzyme(s) (combinations to check).'), sep='\n');
  write(enz_comb, file=paste0(output_path,'/combinations_to_check.txt'), sep='\n');
  
}


#########################################################
########## End of the script ############################
#########################################################


