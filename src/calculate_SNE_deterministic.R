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
# size range (the one that minimizes the Enrichment Value) for some given sets of NF/1000 
# and Score thresholds and for each one of the enzyme(s)(combinations) in the input files.
# The search is performed deterministically, so it can take a while (proportional to the 
# number of size ranges and enzyme combinations to check). This information is stored in 
# the SNE_deterministic.csv_XXX file (where XXX is an index used for parallelization purposes).
#
###########################################################################################
#
# USAGE: Rscript calculate_SNE_deterministic.R --help
#                             
###########################################################################################


###########################################################
##################### Dependencies ########################
###########################################################

suppressWarnings(suppressMessages(library(optparse)));
suppressWarnings(suppressMessages(library(dplyr)));



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
                By default, ref_NF/1000 is the NF/1000 value for MspI in a theoretical size range of 20-800 bp \
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
  
  make_option(c('-s', '--window_step'), type='integer', default=1, 
              help="Step (in bp) taken every time the window slides. DEFAULT=1",
              metavar='integer'),
  
  make_option(c('-w', '--window_width_step'), type='integer', default=1, 
              help="Step taken when sampling different widths for the sliding window (in bp).   \
                DEFAULT=1",
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
The search is performed deterministically, so it can take a while (proportional to   \
the number of size ranges and enzyme combinations to check).                   \

This information is stored in the SNE_deterministic.csv_XXX file (where XXX is \ 
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

# Steps to be taken by the sliding window.

window_step <- as.numeric(opt$window_step);

# Parameters regarding the changing window size.

min_window_size <- 1;
max_window_size <- max_size - min_size + 1;
window_size_step <- as.numeric(opt$window_width_step);
  
# Other parameters.
  
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


################################################################
################## Running the pipeline ########################
################################################################

#### 1. Size ranges matrix: size ranges to test. Create a matrix whose rows are the different minimum sizes
#  (i.e. lower limits) considered in the size ranges and the columns are the different window sizes. The top left half of
#  the matrix will be filled in with the appropiate numbers and the right bottom half will be NA values.

#print('Creating size ranges matrix ...');

# Create a vector with the lower limits for the size ranges.

lower_limits_vector <- seq(from=min_size, to=max_size, by=window_step);

# Create vector with the window sizes to check.

window_sizes_vector <- seq(from=min_window_size, to=max_window_size, by=window_size_step);

# Create the matrix with the size ranges to test.

grid_combinations <- expand.grid(lower_limits_vector, window_sizes_vector);
size_ranges_vector <- apply(grid_combinations, 1, function(x){
  
  min_s <- x[1];
  max_s <- x[1] + x[2] - 1
  
  if(max_s > max_size){ # Make NAs in the appropiate places (i.e. bottom right)
    
    return(NA);
    
  }else{
    
    sr <- paste0(min_s, '_', max_s);
    return(sr);
  }
})

size_ranges_to_test <- matrix(size_ranges_vector,
                              nrow=length(lower_limits_vector), 
                              ncol=length(window_sizes_vector));



#### 2. For each enzyme or enzyme combination and set of thresholds, we will obtain the best size range
#    (i.e. the one that minimizes the EV). The output CSV file will contain the following information:
#     Enzyme_combination, NF_and_Score_thr, Size_range, NF_1000, Score, EV


### 2.1. Calculate the number of enzymes / enzyme combinations in the input files.

#print('Calculating the number of enzyme(s) (combinations) available ...');

fld_headers <- system(paste0("grep -n '>' ", fld_file_path), intern=T);
fragments_headers <- system(paste0("grep -n '>' ", fragments_file_path), intern=T);

if(length(fld_headers) != length(fragments_headers)){
  
  stop("The input files do not have the same number of headers.");
  
}

ncombs <- length(fld_headers);

#print(paste0('There are ', ncombs, ' enzyme(s) (combinations) available in the input files.'));

### 2.2. Read the information for each individual enzyme / enzyme combination.

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
  
  cat(paste0('        Deterministic optimisation of enzyme(s) (combinations) ', i, '/', length(fld_chunk_sep),' ...'), sep='\n');
  
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
  
  

### 2.3. Create three matrices with the same size as 'size_ranges_to_test', that contain,
#    respectively, the NF/1000, Score and EV values for the different size ranges assessed.

  ## Calculate the NF/1000 matrix. 
  
  NF_1000_matrix <- matrix(mapply(NF_calculation, size_ranges_to_test, 
                                  MoreArgs = list(fld=current_fld)), 
                           nrow=nrow(size_ranges_to_test),
                           ncol=ncol(size_ranges_to_test));
  rownames(NF_1000_matrix) <- lower_limits_vector;
  colnames(NF_1000_matrix) <- window_sizes_vector;
  NF_1000_matrix[which(NF_1000_matrix==0)] <- NF_1000_matrix[which(NF_1000_matrix==0)] + 0.001; # Avoid values of 0 for posterior calculations
  
  
  ## Calculate the Score matrix.
  
  Score_matrix <- matrix(unlist(mapply(Score_calculation, size_ranges_to_test, 
                                       MoreArgs = list(sinfo=current_fragments_df, sann=sites_annotation, rl=read_length, ids=FALSE))), 
                         nrow=nrow(size_ranges_to_test),
                         ncol=ncol(size_ranges_to_test));
  rownames(Score_matrix) <- lower_limits_vector;
  colnames(Score_matrix) <- window_sizes_vector;
  
  
  ## Calculate the EV (enrichment value) matrix. 
  
  n_sites <- nrow(sites_annotation);
  max_score <- sum(as.numeric(sites_annotation[,4]));
  EV_matrix <- - log10((Score_matrix / (NF_1000_matrix * 1000)) *
                         (n_sites/max_score));
  rownames(EV_matrix) <- lower_limits_vector;
  colnames(EV_matrix) <- window_sizes_vector;
  
  # Transform the EV matrix so it is optimizable. 
  # Those EV values = Inf are converted into the maximum EV value found in the matrix.
  
  if(length(EV_matrix[is.finite(EV_matrix)]) > 0){ 
    
    max_EV_value <- max(EV_matrix[is.finite(EV_matrix)], na.rm=T);
    
  }else{ # When all the elements in the EV_matrix are Inf (for really bad enzymes, when all the Scores are 0)
    
    max_EV_value <- 9999; # Fix arbitrary high value to allow the optimization to proceed
    
  }
  
  
  EV_matrix_t <- EV_matrix;
  EV_matrix_t[is.infinite(EV_matrix_t)] <- max_EV_value; 
  
  
#### 2.4. Filter the matrices according to the NF/1000 and Score thresholds.
# Those elements that do not satisfy both thresholds will be NAs in all the 3 matrices.
# Store the results for the different thresholds in a list.
  
  # If the thresholds were not specified correctly: quit the script
  
  if(length(NF_thr) != length(Score_thr)){
    
    stop('The thresholds for NF/1000 and Score are not properly specified. Quit the script !!');
    
  }
  
  # Initialize the lists which will contain the filtered matrices for the different thresholds.
  
  NF_filtered_list <- list();
  Score_filtered_list <- list();
  EV_filtered_list <- list();
  
  # Iterate over the different thresholds.
  
  for(j in 1:length(Score_thr)){
    
    # Calculate a matrix that will be used for filtering.
    
    filter_matrix <- (NF_1000_matrix <= NF_thr[j]) & (Score_matrix > Score_thr[j]);
    filter_matrix[filter_matrix==FALSE] <- NA;
    n_rows_f <- dim(filter_matrix)[1];
    n_cols_f <- dim(filter_matrix)[2];
    
    # Filter the 3 matrices.
    
    NF_filtered <- matrix(NF_1000_matrix[filter_matrix], nrow=n_rows_f, ncol=n_cols_f);
    rownames(NF_filtered) <- rownames(NF_1000_matrix);
    colnames(NF_filtered) <- colnames(NF_1000_matrix);
    
    Score_filtered <- matrix(Score_matrix[filter_matrix], nrow=n_rows_f, ncol=n_cols_f);
    rownames(Score_filtered) <- rownames(Score_matrix);
    colnames(Score_filtered) <- colnames(Score_matrix);
    
    EV_filtered <- matrix(EV_matrix_t[filter_matrix], nrow=n_rows_f, ncol=n_cols_f);
    rownames(EV_filtered) <- rownames(EV_matrix_t);
    colnames(EV_filtered) <- colnames(EV_matrix_t);
    
    # Store the matrices in the lists. If all the elements in the matrix are NAs, store NA in the list.
    
    if(sum(is.na(filter_matrix)) == (n_rows_f * n_cols_f)){
      
      NF_filtered_list[[j]] <- NA;
      Score_filtered_list[[j]] <- NA;
      EV_filtered_list[[j]] <- NA;
      
    }else{
      
      NF_filtered_list[[j]] <- NF_filtered;
      Score_filtered_list[[j]] <- Score_filtered;
      EV_filtered_list[[j]] <- EV_filtered;
      
    }
  }
  
  
#### 2.5. Append in the output file the best size range found (the one that minimizes EV),
#        for each one of the set of thresholds.
  
  # Format: Enzyme_combination, NT_and_ST, Size_range, NF_1000, Score, EV
  
  for(j in 1:length(Score_thr)){
    
    if(is.matrix(EV_filtered_list[[j]])){
      
      index_top <- which.min(EV_filtered_list[[j]]);
      output_line <- paste(fld_enz_combs[i],
                           paste0(NT[j], '_', ST[j]),
                           size_ranges_to_test[index_top],
                           NF_filtered_list[[j]][index_top],
                           Score_filtered_list[[j]][index_top],
                           EV_filtered_list[[j]][index_top],
                           sep=',');
      
      output_name <- paste0('SNE_deterministic.csv', '_', output_ID);
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
