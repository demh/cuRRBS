###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             28/12/2016                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################
###########################################################################################
#####  Optimization of a new RRBS assay using a new combination of restriction enzymes ####
###########################################################################################
##### Given the annotation file with the isoschizomer and methylation-sensitivity      ####
##### information for the comercially available enzymes (e.g. isoschizomers_CpG_annotation.csv),
##### create the "enzymes_to_pre_compute.txt" file that will be used by the            ####
##### "pre_compute_digestions.py" script.                                              ####
##### Only the first enzyme of an isoschizomer family that contains at least one       #### 
##### methylation-insensitive enzyme (0) will be stored in the file. This enzyme       ####
##### will represent all the enzymes of the family to avoid redundant computations.    ####
###########################################################################################
##### USAGE: Rscript iso_annotation_2_pre_computed_enzymes.R --help                    ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

if(suppressWarnings(suppressMessages(require(optparse) == FALSE))){
  
  cat("\nThe 'optparse' package is missing. I will try to install it ...\n");
  cat("\n");
  install.packages('optparse', repos = 'https://mirrors.ebi.ac.uk/CRAN/');
  cat("\n");
  
  if(suppressWarnings(suppressMessages(require(optparse) == FALSE))){
    
    stop("I could not install the 'optparse' package. Please, fix this and rerun the script.");
    
  }
  
  suppressWarnings(suppressMessages(library(optparse)));
  

}else{
  
  suppressWarnings(suppressMessages(library(optparse)));
  
}


###########################################################
############## Command-line arguments #####################
###########################################################

## Create arguments and help in a Pythonic style.

option_list <-  list(
  
  make_option(c('-i', '--iso_file'), type='character', default=NULL, 
              help="Absolute path to the annotation file with the isoschizomer family \
                and methylation-sensitivity information for the comercially available \
                enzymes (e.g. isoschizomers_CpG_annotation.csv).  COMPULSORY", 
              metavar='character'),
  
  make_option(c('-o', '--output_path'), type='character', default=NULL,
              help="Absolute path for the output directory.  COMPULSORY",
              metavar='character')
);



## Prepare option list and general description of the script.

opt_parser <- OptionParser(option_list=option_list,
                           description='\nGiven the annotation file with the isoschizomer and methylation-sensitivity information for the \
comercially available enzymes (e.g. isoschizomers_CpG_annotation.csv), create the "enzymes_to_pre_compute.txt" \
file that will be used by the "pre_compute_digestions.py" script. \
\
Only the first enzyme of an isoschizomer family that contains at least one methylation-insensitive \
enzyme (value of 0) will be stored in the file. This enzyme will represent all the enzymes of the family to \
avoid redundant computations.');        


opt <- parse_args(opt_parser);

## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$iso_file) | is.null(opt$output_path)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Input arguments.

annotation_path <- as.character(opt$iso_file);

if(!file.exists(annotation_path)){
  
  cat('\n');
  stop(paste('The annotation file with the isoschizomer and methylation-sensitivity information could not be found.',
              'Please provide a valid path for the annotation file.', sep='\n'));
  
}
  
  
output_path <- as.character(opt$output_path);

if(!file.exists(output_path)){
  
  cat('\n');
  stop('The output path could not be found. Please provide a valid output path.');
  
}

if(substr(output_path, nchar(output_path), nchar(output_path)) != '/'){
  
  output_path <- paste0(output_path, '/');

}


################################################################
################## Running the code ############################
################################################################


### Read the annotation file ###

raw_annotation <- read.csv(annotation_path);


### Extract the names of the enzymes to store ###

family_names <- as.character(unique(raw_annotation[,1]));
final_enzymes <- c();

for(family in family_names){
  
  family_members <- raw_annotation[(as.character(raw_annotation[,1]) == family),];
  
  # If the family contains a methylation-insensitive enzyme, store the family name in 
  # the final enzymes to store.
  
  if(0 %in% family_members[,3]){ 
    
    final_enzymes <- c(final_enzymes, family);
    
  }
}


### Create the output file (enzymes_to_pre_compute.txt) ###

# Check if there is another enzymes_to_pre_compute.txt file in the output directory.

if(file.exists(paste0(output_path, 'enzymes_to_pre_compute.txt'))){
  
  cat('\n');
  stop(paste0('The "enzymes_to_pre_compute.txt" file already exists in the output directory. ',
              'Please delete the old file and rerun this script.'));

}else{
  
  write(x=final_enzymes, file=paste0(output_path, 'enzymes_to_pre_compute.txt'),
        sep='\n');
  
  
}

cat('\n');
cat('The "enzymes_to_pre_compute.txt" file was created correctly.');
cat('\n');
cat('\n');

#########################################################
########## End of the script ############################
#########################################################