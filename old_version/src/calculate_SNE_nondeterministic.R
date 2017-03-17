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
# Given the fl_distributions_x_enzymes.txt and fragments_of_interest_x_enzymes.txt files 
# obtained with the obtain_distributions_and_fragments.py script, this script finds the best 
# size range (the one that minimizes the Enrichment Value) for some given sets of NF/1000 and 
# Score thresholds and for each one of the enzyme(s)(combinations) in the input files.
# The search is performed non-deterministically, using an implementation of the Differential
# Evolution Optimization algorithm, so the result is not necessarily the best possible 
# solution. This information is stored in the SNE_nondeterministic.csv_XXX file (where XXX is
# an index used for parallelization purposes).
#
###########################################################################################
#
# USAGE: Rscript calculate_SNE_nondeterministic.R --help           
#
###########################################################################################


###########################################################
##################### Dependencies ########################
###########################################################

suppressWarnings(suppressMessages(library(optparse)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(DEoptimR)));


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
  
  make_option(c('-l', '--annotation_file'), type='character', default=NULL, 
              help="Absolute path to the annotation file, which contains the information regarding \
                the sites of interest. COMPULSORY",
              metavar='character'),
  
  make_option(c('-p', '--output_path'), type='character', default=NULL,
              help="Absolute path for the output directory.  COMPULSORY",
              metavar='character'),
  
  make_option(c('-r', '--read_length'), type='integer', default=NULL,
              help="Read length (in bp) for the RRBS experiment. This determines whether a CpG is 'seen' in \
                a size-selected fragment after the sequencing (i.e. only if it is close to one of the ends of \
                the fragment).  COMPULSORY",
              metavar='integer'),
  
  make_option(c('-t', '--NF_thresholds'), type='character', default='1', 
              help="Float (named NT) or vector (e.g. '0.5,1,40') used to calculate the different NF/1000 thresholds \
              which will filter the final output. NT takes values in the interval (0, +Inf). \
              Only those size ranges with NF/1000 <= NT * ref_NF/1000 are kept.  \ 
              By default, ref_NF/1000 is the NF/1000 value for MspI in a theoretical size range of 20-800 bp. \
              i.e. 1343.468. DEFAULT='1'",
              metavar='character'),
  
  make_option(c('-x', '--Score_thresholds'), type='character', default='0.1', 
              help="Float (named ST) or vector e.g. ('0.2,0.5,0.8') used to calculate the different Score thresholds \
              which will filter the final output. ST takes values in the interval (0,1]. 
              Only those size ranges with Score > ST * max_score are kept (max_score is the Score that would be \ 
              obtained if all the sites of interest were captured in the sequencing). DEFAULT='0.1'",
              metavar='character'),
  
  make_option(c('-a', '--min_size'), type='integer', default=1, 
              help="Minimum size of the fragments (in bp) to be considered when screening for \
                different size ranges. DEFAULT=1",
              metavar='integer'),
  
  make_option(c('-b', '--max_size'), type='integer', default=1000, 
              help="Maximum size of the fragments (in bp) to be considered when screening for \
                different size ranges. DEFAULT=1000",
              metavar='integer'),
  
  make_option(c('-s', '--max_iterations'), type='integer', default=150, 
              help="Maximum number of times that the optimization procedure is attempted for \
                each one of the enzyme(s) (combinations). This is necessary to handle those \
                cases that will never find an optimal value due to the NF and Score thresholds. \
                The higher this number, the more computational time is needed (but it is also \
                more likely to find the correct answer). DEFAULT=150",
              metavar='integer'),
  
  make_option(c('-d', '--output_ID'), type='character', default='000', 
              help="Suffix used in the output file to allow parallelization. DEFAULT='000'",
              metavar='character')
  
);



## Prepare option list and general description of the script.

opt_parser <- OptionParser(option_list=option_list,
                           description="\nGiven the fl_distributions_x_enzymes.txt and fragments_of_interest_x_enzymes.txt \
files obtained with the obtain_distributions_and_fragments.py script, this script  \ 
finds the best size range (the one that minimizes the Enrichment Value) for some \
given sets of NF/1000 and Score thresholds and for each one of the enzyme(s)     \
(combinations) in the input files.                                              \
 \
The search is performed non-deterministically, using an implementation of the Differential \
Evolution Optimization algorithm, so the result is not necessarily the best \
possible solution. \  
\
This information is stored in the SNE_nondeterministic.csv_XXX file (where XXX is \ 
an index used for parallelization purposes).");        


opt <- parse_args(opt_parser);

## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$fld_file) | is.null(opt$fragments_file) |
   is.null(opt$annotation_file) | is.null(opt$output_path) | 
   is.null(opt$read_length)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Input arguments.


# Output path / working directory.

output_path <- as.character(opt$output_path);
setwd(output_path);


# Minimum and maximum lengths of fragments to be considered.

min_size <- as.numeric(opt$min_size);
max_size <- as.numeric(opt$max_size);

  
# Other parameters.

max_counter <- as.numeric(opt$max_iterations);  
read_length <- as.numeric(opt$read_length);
ST <- as.numeric(strsplit(opt$Score_thresholds, ',')[[1]]);
NT <- as.numeric(strsplit(opt$NF_thresholds, ',')[[1]]);
output_ID <- as.character(opt$output_ID); # ID at the end of the file to avoid problems derived from parallelization

# Input files.

fld_file_path <- as.character(opt$fld_file);
fragments_file_path <- as.character(opt$fragments_file);
sites_annotation <- read.table(as.character(opt$annotation_file), sep=',', header=T);
sites_annotation[,1] <- as.character(sites_annotation[,1]); 


# Calculate absolute Score thresholds.

max_score <- sum(as.numeric(sites_annotation[,4]));
Score_thr <- ST * max_score;

# Calculate absolute NF/1000 thresholds.

ref_NF_1000 <- 1343.468;  # It can be modified to compare with other protocol.
NF_thr <- NT * ref_NF_1000;


################################################################
################## Functions ###################################
################################################################

## Function: convert the information from the input in fragments_of_interest_x_enzymes.txt from 
#  a collapsed format into an expanded format. The outputted dataframe can be used by the rest of the functions in 
#  the script.
#  
#  Collapsed format: Site_IDs Fragment_start Fragment_length
#  Expanded format: Row_number   Chromosome      Coordinate    Fragment_start    Fragment_end    Fragment_length
#
# input_fragments_df: input dataframe with the same format as the fragments_of_interest_x_enzymes.txt file (collapsed format). i.e. current_ncs object
# sann:  annotation for the sites of interest, as displayed in the annotation file.
#        Site_ID   Chr   Coordinate  Weight


expand_fragments_df <- function(input_fragments_df, sann){
  
  # Expand the input dataframe for those rows with several sites of interest.
  
  expanded_input <- input_fragments_df[!grepl('_', input_fragments_df[,1]),]; # Rows with only one site
  if(length(expanded_input) <= 3){ 
    expanded_input <- matrix(expanded_input, byrow=T, ncol=3);
  }
  colnames(expanded_input) <- c('Site_ID', 'Fragment_start', 'Fragment_length');
  
  for(multi_row in grep('_', input_fragments_df[,1])){
    
    sites <- strsplit(input_fragments_df[multi_row,1], '_')[[1]];
    expanded_sites <- cbind(sites, rep(input_fragments_df[multi_row,2], length(sites)), 
                            rep(input_fragments_df[multi_row,3], length(sites)));
    expanded_input <- rbind(expanded_input, expanded_sites);
    
  }
  
  # Merge the information from the expanded input with the annotation information.

  expanded_input <- as.data.frame(expanded_input);
  expanded_input[,1] <- as.character(expanded_input[,1]);
  merged_expanded <- full_join(expanded_input, as.data.frame(sann), by='Site_ID');
  merged_expanded <- arrange(merged_expanded, Site_ID);
  
  # Create the final dataframe.
  
  starts <- as.numeric(as.character(merged_expanded$Fragment_start));
  lengths <- as.numeric(as.character(merged_expanded$Fragment_length));
  
  final_df <- data.frame(Row_number=1:nrow(merged_expanded), 
                         Chromosome=merged_expanded$Chr, 
                         Coordinate=merged_expanded$Coordinate,
                         Fragment_start=starts, 
                         Fragment_end=starts+lengths-1, 
                         Fragment_length=lengths);
  
  return(final_df);
  
}


## Function: given the fragment size distribution of an enzyme combination and the size range of
#  interest, calculate the appropiate NF/1000 value. NF represents the number of fragments that are
#  captured in the selected size range. NF is divided by 1000 for visualization purposes.

# fld: numeric matrix which contains the fragment length distribution. First column: fragment length. Second column: number of fragments 
#      with that size.
# size_range: range of fragment lengths selected in the RRBS protocol (e.g. 40_220). If NA, then the function returns NA.


NF_calculation <- function(fld, size_range){
  
  if(is.na(size_range)){
    
    return(NA);
    
  }else{
    
  
  
    min_size_f <- as.numeric(strsplit(size_range, '_')[[1]][1]);
    max_size_f <- as.numeric(strsplit(size_range, '_')[[1]][2]);
    selected_fragments <- fld[which(fld[,1] >= min_size_f & fld[,1] <= max_size_f),];
    
    if(length(selected_fragments)>2){
      
      NF_1000 <- sum(as.numeric(selected_fragments[,2])) / 1000;
      
    }
    
    if(length(selected_fragments)==2){
      
      NF_1000 <- as.numeric(selected_fragments[2]) / 1000;
      
    }
    
    if(length(selected_fragments)<2){
      
      NF_1000 <- 0;
      
    }
    
    return(NF_1000);
    
  }
}


## Function: given the information for the sites of interest and the size range selected, 
#  calculate the Score for an enzyme combination. It must take into account
#  how far the site of interest is from the closest end of the fragment. 
#  The Score is calculated as the weighted sum of the sites of interest that are captured 
#  in the selected size range. 
#  The output of the function is a list with the Score and, if ids=TRUE, the 'Site_ID's 
#  (character: e.g. S001,S003) for the sites of interest used to calculate it.

# sinfo: information about the fragments where the sites of interest are located. The format of
#       the dataframe is the expanded one (i.e. after using the function expand_fragments_df on the
#       fragments_of_interest_x_enzymes.txt input file). The last four columns must be 'numeric'.
#       Row_number   Chromosome      Coordinate    Fragment_start    Fragment_end    Fragment_length
# sann:  annotation for the sites of interest, as displayed in the new_sites_annotation.csv file.
#        Site_ID   Chr   Coordinate  Weight
# size_range: range of fragment lengths selected in the RRBS protocol (e.g. 40_220). If NA, then the function returns NA.
# rl: maximum read length for the sequencing library. A site of interest which is located too far 
#     from the end of the fragment won't be sequenced and, therefore, shouldn't be taken into account
#     for the Score calculation.
# ids: if TRUE, the function reports also the Site IDs for the sites of interest.


Score_calculation <- function(sinfo, sann, size_range, rl, ids=FALSE){
  
  if(is.na(size_range)){
    
    return(NA);
    
  }else{

  
    # Subset the selected sites (i.e. sites in fragments inside the size range)
    
    min_size_f <- as.numeric(strsplit(size_range, '_')[[1]][1]);
    max_size_f <- as.numeric(strsplit(size_range, '_')[[1]][2]);
    selected_sites <- sinfo[which(sinfo[,6] >= min_size_f & sinfo[,6] <= max_size_f),];
    
    if(nrow(selected_sites) == 0){  # If there are no selected sites in our size range, exit the function
      
      if(ids){
        
        return(list(Score=0, Site_IDs = NA));
        
        
      }else{
        
        return(list(Score=0));
        
      }
    }
    
    # From the selected sites, choose only those that are in a 'sequenceable' location 
    # (i.e. close to one of the extremes).
    
    distance_to_ends <- cbind((selected_sites[,3] - selected_sites[,4] + 1), 
                              (selected_sites[,5] - selected_sites[,3] + 1));
    distance_to_closest_end <- apply(distance_to_ends, 1, min);
    sequenceable_sites <- selected_sites[distance_to_closest_end <= rl,];
    
    if(nrow(sequenceable_sites) == 0){  # If there are no selected sites in our size range, exit the function
      
      if(ids){
        
        return(list(Score=0, Site_IDs = NA));
        
        
      }else{
        
        return(list(Score=0));
        
      }
      
    }
    
    # If we don't need the Site_IDs and the weights are all the same (i.e. 1): shortcut
    
    if(ids==FALSE & all(as.numeric(sann[,4]) == 1)){
      
      score <- nrow(sequenceable_sites);
      
    }else{ # Otherwise
      
      # For each sequenceable site, obtain its Site_ID and weight.
      
      seq_sites_chr_coord <- sequenceable_sites[,c(2,3)];
      seq_sites_chr_coord[,1] <- as.character(seq_sites_chr_coord[,1]);
      seq_sites_chr_coord[,2] <- as.numeric(seq_sites_chr_coord[,2]);
      seq_sites_info <- apply(seq_sites_chr_coord, 1, function(x){
        
        x <- unlist(x);
        index_site <- which(as.character(sann[,2]) == x[1] 
                            & as.numeric(sann[,3]) == as.numeric(x[2]));
        
        site_weight <- as.numeric(sann[index_site,4]);
        site_id <- as.character(sann[index_site,1]);
        
        return(list(site_id, site_weight));
        
      });
      
      
      final_site_IDs <- unlist(lapply(seq_sites_info, function(x){x[[1]]}));
      order_index <- order(final_site_IDs);
      final_site_IDs <- final_site_IDs[order_index];
      
      final_weights <- unlist(lapply(seq_sites_info, function(x){x[[2]]}));
      final_weights <- final_weights[order_index];
      
      if(length(final_site_IDs)==0 & length(final_weights) == 0){ # If all the values were NAs
        
        return(list(Score=0, Site_IDs = NA));
        
      }
      
      
      # Calculate the Score (sum of the weights of the selected sites).
      
      score <- sum(final_weights);
      
    }
    
    # Return the Score and, if needed, the Site_IDs.
    
    if(ids){
      
      return(list(Score=score, Site_IDs = paste0(final_site_IDs, collapse = ',')));
      
    }else{
      
      return(list(Score=score));
      
    }
  }
}



## Function: given the lower limit of the size range and the window size, calculate the
#  enrichment value (EV) for the corresponding enzyme / enzyme combination. This function
#  will be the one used during the derivative-free optimization process.
#  If the NF/1000 or Score thresholds are not satisfied, the function returns Inf.
#  For those cases where the upper limit of the size range is greater than the maximum
#  value allowed (i.e. lower_limit + window_size - 1 > max_size), the function returns Inf.

# var_vector: vector containing the variables to optimize i.e. c(lower_limit, window_size), which are:

#       - lower_limit: length (in bp) of the smaller fragments that will be captured in the corresponding size range.
#       - window_size: width of the size range. e.g. for a size range of 40-220 bp, the window_size = 181 bp

EV_calculation <- function(var_vector){
  
  lower_limit <- floor(var_vector[1]);
  window_size <- floor(var_vector[2]);
  
  ## Obtain the size range being tested.
  
  upper_limit <- lower_limit + window_size - 1;
  
  if(upper_limit > max_size){ # Return Inf if the upper limit of the size range exceeds the maximum value
    
    return(Inf);
  }
  
  size_range_tested <- paste0(lower_limit, '_', upper_limit);
  
  ## Calculate the NF/1000 value.
  
  NF_1000_value <- NF_calculation(fld=current_fld, size_range=size_range_tested);
  
  if(NF_1000_value > NF_thr_i){ # Return Inf if the NF/1000 threshold is not satisfied
    
      return(Inf);
  }
  
  ## Calculate the Score value.
  
  Score_value <- as.numeric(Score_calculation(sinfo=current_fragments_df, sann=sites_annotation, 
                                              size_range=size_range_tested, rl=read_length,
                                              ids=FALSE));
  
  if(Score_value <= Score_thr_i){ # Return Inf if the Score threshold is not satisfied
    
      return(Inf);
  }
  
  
  ## Calculate the EV value.
  
  EV_value <- - log10((Score_value / (NF_1000_value * 1000)) *
                        (n_sites/max_score));
  
  return(EV_value);
  
}


################################################################
################## Running the pipeline ########################
################################################################


#### 1. For each enzyme or enzyme combination and set of thresholds, we will obtain the best size range
#    (i.e. the one that minimizes the EV) that the optimization algorithm can find. 
#    The output CSV file will contain the following information:
#    Enzyme_combination, NF_and_Score_thr, Size_range, NF_1000, Score, EV


### 1. Calculate the number of enzymes / enzyme combinations in the input files.

#print('Calculating the number of enzyme(s) (combinations) available ...');

fld_headers <- system(paste0("grep -n '>' ", fld_file_path), intern=T);
fragments_headers <- system(paste0("grep -n '>' ", fragments_file_path), intern=T);

if(length(fld_headers) != length(fragments_headers)){
  
  stop("The input files do not have the same number of headers.");
  
}

ncombs <- length(fld_headers);

#print(paste0('There are ', ncombs, ' enzyme(s) (combinations) available in the input files.'));


## Other useful variables.

n_sites <- nrow(sites_annotation);
max_score <- sum(as.numeric(sites_annotation[,4]));



### 2. Read the information for each individual enzyme / enzyme combination.

## Create vectors which contain the line numbers of the chunk separators (e.g. >AanI_AarI). Extract also the enzyme combinations.

fld_chunk_sep <- as.numeric(sapply(strsplit(fld_headers, ':'), function(x){x[1]})); # Line numbers
fld_enz_combs <- sub('>', '', sapply(strsplit(fld_headers, ':'), function(x){x[2]})); # Enzyme combinations

fragments_chunk_sep <- as.numeric(sapply(strsplit(fragments_headers, ':'), function(x){x[1]})); # Line numbers
fragments_enz_combs <- sub('>', '', sapply(strsplit(fragments_headers, ':'), function(x){x[2]})); # Enzyme combinations

## Check that the enzyme combinations are in the same order for both the fld and the fragments files. Otherwise rise an error.

if(!all(fld_enz_combs == fragments_enz_combs)){
  
  stop('The enzyme combinations are not the same for the fl_distribution_x_enzymes.txt and the fragments_of_interest_x_enzymes.txt files.');
  
}


## Obtain the information from the input files and format it for a given enzyme / enzyme combination.

fld_file <- file(fld_file_path);
open(fld_file);

fragments_file <- file(fragments_file_path);
open(fragments_file);

for(i in 1:length(fld_chunk_sep)){
  
  cat(paste0('        Non-deterministic optimisation of enzyme(s) (combinations) ', i, '/', length(fld_chunk_sep),' ...'), sep='\n');
  
  # Read the fld and the fragment chunks.
  
  if(i == length(fld_chunk_sep)){ # For the last combination
    
    current_fld <- readLines(fld_file, n=-1);
    current_fragments <- readLines(fragments_file, n=-1);  
    
    
  }else{
    
    current_fld <- readLines(fld_file, n=fld_chunk_sep[i+1] - fld_chunk_sep[i]); 
    current_fragments <- readLines(fragments_file, n=fragments_chunk_sep[i+1] - fragments_chunk_sep[i]); 
    
  }
  
  
  # Format the fld chunk of information.
  
  current_fld <- current_fld[-1];
  current_fld <- matrix(as.numeric(unlist(strsplit(current_fld, ','))), ncol=2, byrow=T);
  
  # Format the fragment chunk of information.
  
  current_fragments <- current_fragments[-1];
  current_fragments <- matrix(unlist(strsplit(current_fragments, ',')), ncol=3, byrow=T);
  current_fragments_df <- expand_fragments_df(current_fragments, sites_annotation);
  

### 3. Obtain the optimized size ranges for the current enzyme / enzyme combination and each one of the 
#   Score / NF threshold pairs.  
  
  ## Run the Differential Evolution Optimization algorithm for the current enzyme / enzyme combination 
  ## and each one of the Score / NF threshold pairs.
  
  # Initialize the matrix to store the results for the current enzyme / enzyme combination.
  
  results_enz <- data.frame(NF_thresholds=NF_thr, 
                            Score_thresholds=Score_thr, 
                            Lower_limit=rep(NA, length(NF_thr)), 
                            Upper_limit=rep(NA, length(NF_thr)), 
                            EV_values=rep(NA, length(NF_thr)));
  
  
  # Iterate over the different NF and Score thresholds.
  
  for(j in 1:length(NF_thr)){
    
    # Starting conditions.
    
    if(exists('opt_result')){
      
      rm(opt_result);
      
    }
    
    counter <- 1; 
    NF_thr_i <- NF_thr[j];
    Score_thr_i <- Score_thr[j];

    # Run optimization algorithm.
    
    while(!exists('opt_result') & counter < max_counter){
      
      # The errors arising from getting only Inf values for EV at the beginning
      # are supressed. Instead, the optimization algorithm is restarted until
      # some useful result is found.
      
      # The warnings arising from lack of convergence of the algorithm are also supressed.
      
      suppressWarnings(try(opt_result <- JDEoptim(lower=c(min_size,1), upper=c(max_size-1, max_size-min_size+1), fn=EV_calculation,
                                 maxiter=100, NP=20, tol=0.001), silent=T));
      
      counter <- counter + 1; 
      
    }
    
    # Store the results if there is any size range that satisfies the thresholds.
    
    if(exists('opt_result')){
      
      results_enz$Lower_limit[j] <- floor(opt_result$par[1]);
      results_enz$Upper_limit[j] <- floor(opt_result$par[1]) + floor(opt_result$par[2]) - 1;
      results_enz$EV_values[j] <- opt_result$value;
      
    }
  }
  
  
### 4. Append in the output file the best size range found (the one that minimizes EV),
#      for each one of the set of thresholds, in case there was some size range that 
#      satisfied the constraints imposed by the Score and NF thresholds.
  
  # Format: Enzyme_combination, NF_and_Score_thr, Size_range, NF_1000, Score, EV
  
  for(set in 1:nrow(results_enz)){
    
    # Check if there is any size range that satisfied the thresholds
    
    if(sum(is.na(results_enz[set,])) == 0){
      
      # Obtain size range
      
      size_range_set <- paste0(results_enz[set,3], '_', results_enz[set,4]);
      
      # Calculate NF/1000 value
      
      NF_1000_set <- NF_calculation(fld=current_fld, size_range_set);
      
      # Calculate Score value
      
      Score_set <- Score_calculation(sinfo=current_fragments_df, sann=sites_annotation, 
                                     size_range=size_range_set, rl=read_length,
                                     ids=FALSE);
      
      # Append in output file
      
      output_line <- paste(fld_enz_combs[i],
                           paste0(round(results_enz[set,1]/ref_NF_1000, 2), '_', round(results_enz[set,2]/max_score,2)),
                           size_range_set,
                           NF_1000_set,
                           Score_set,
                           results_enz[set,5],
                           sep=',');
      
      output_name <- paste0('SNE_nondeterministic.csv', '_', output_ID);
      write(x=output_line, file=output_name, append=T);
      
    }
  }
}

## Close the input file connections.

close(fld_file);
close(fragments_file);


#print('The script finished correctly.');


#########################################################
########## End of the script ############################
#########################################################
