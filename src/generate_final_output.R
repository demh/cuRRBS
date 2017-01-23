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
# This script takes the SNE files available and generates the final output in a CSV file, 
# which contains the following columns:                                
#                                                                                  
#       1. Enzyme / enzyme combination                                             
#       2. Optimal experimental size range                                         
#       3. Optimal theoretical size range                                          
#       4. Score value                                                             
#       5. % of maximum Score achieved                                             
#       6. NF/1000 value                                                           
#       7. Cost reduction factor (CRF) = (NF/1000 for original MspI / NF/1000 value)  
#       8. Enrichment Value (EV)                                                   
#       9. Robustness measure                                                      
#       10. Thresholds used (C_Score | C_NF/1000)                                   
#       (11.) Number of sites of interest that will be sequenced theoretically     
#       (12.) IDs of the sites of interest that will be sequenced theoretically    
#                                                                                  
# The enzymes are ordered by the EV value (from minimum to maximum) i.e. from the best one 
# to the worse one (considering only Score and NF/1000 values).
###########################################################################################
#
# USAGE: Rscript generate_final_output.R --help     
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
  
  make_option(c('-s', '--sne_1'), type='character', default=NULL,
              help="Absolute path to the SNE file for individual enzymes. COMPULSORY",
              metavar='character'),
  
  make_option(c('-t', '--sne_2'), type='character', default=NULL,
              help="Absolute path to the SNE file for the 2-enzyme combinations. COMPULSORY",
              metavar='character'),
  
  make_option(c('-f', '--df_files'), type='character', default=NULL, 
              help="Absolute path to the folder containing the fl_distribution_x_enzymes.txt  \
                and the fragments_of_interest_x_enzymes.txt files, which are created by \
                obtain_distributions_and_fragments.py. COMPULSORY",
              metavar='character'),
  
  make_option(c('-e', '--eann'), type='character', default=NULL, 
              help="Absolute path to the enzymes annotation file, which contains the information regarding \
                the different isoschizomer families and their methylation sensitivity. COMPULSORY",
              metavar='character'),
  
  make_option(c('-l', '--sann'), type='character', default=NULL, 
              help="Absolute path to the sites annotation file, which contains the information regarding \
                the sites of interest. COMPULSORY",
              metavar='character'),
  
  make_option(c('-p', '--output_path'), type='character', default=NULL,
              help="Absolute path for the output directory.  COMPULSORY",
              metavar='character'),
  
  make_option(c('-a', '--adapters_size'), type='integer', default=NULL,
              help="Total size (in bp) of the adapters used for the RRBS experiment. This will be \
                used to calculate the experimental size range. COMPULSORY",
              metavar='integer'),
  
  make_option(c('-r', '--read_length'), type='integer', default=NULL,
              help="Read length (in bp) for the RRBS experiment. This determines whether a CpG is 'seen' in \
                a size-selected fragment after the sequencing (i.e. only if it is close to one of the ends of \
                the fragment).  COMPULSORY",
              metavar='integer'),
  
  make_option(c('-d', '--exp_error'), type='integer', default=NULL,
              help="Experimental error assumed when performing the size range selection (in base pairs). COMPULSORY",
              metavar='integer'),
  
  make_option(c('-b', '--top_enzymes'), type='integer', default=NULL,
              help="Maximum number of enzymes / enzyme combinations to display in the final output file. COMPULSORY",
              metavar='integer'),
  
  make_option(c('-i', '--sites_info'), type='character', default=NULL,
              help="If 'Y', the number of sites of interest that will be theoretically be sequenced and \
                their respective IDs are included in the output.  COMPULSORY",
              metavar='character')
  
); 
  

opt_parser <- OptionParser(option_list=option_list,
                           description="\nThis script takes the SNE files available and generates the final output in a \
CSV file, with the following columns: \
\
       1. Enzyme / enzyme combination \
       2. Optimal experimental size range \
       3. Optimal theoretical size range \
       4. Score value \
       5. % of maximum Score achieved \
       6. NF/1000 value \
       7. Cost reduction factor (CRF) = (NF/1000 for original MspI / NF/1000 value) \
       8. Enrichment Value (EV) \
       9. Robustness measure \
       10. Thresholds used (C_Score | C_NF/1000) \
       (11.) Number of sites of interest that will be sequenced theoretically \
       (12.) IDs of the sites of interest that will be sequenced theoretically \
\
The enzymes are ordered by the EV value (from minimum to maximum) i.e. from the \
best one to the worse one (considering only Score and NF/1000 values).");

opt <- parse_args(opt_parser);  

## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$sne_1) | is.null(opt$sne_2) |
   is.null(opt$df_files) | is.null(opt$eann) | 
   is.null(opt$sann) | is.null(opt$output_path) |
   is.null(opt$adapters_size) | is.null(opt$read_length) |
   is.null(opt$exp_error) | is.null(opt$top_enzymes) |
   is.null(opt$sites_info)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Input arguments.

sne_1_path <- as.character(opt$sne_1);
sne_2_path <- as.character(opt$sne_2);
df_files_path <- as.character(opt$df_files);
enz_annotation_path <- as.character(opt$eann);
sites_annotation_path <- as.character(opt$sann);
output_path <- as.character(opt$output_path);

adapters_size <- as.numeric(opt$adapters_size); 
read_length <- as.numeric(opt$read_length); 
experimental_error <- as.numeric(opt$exp_error); 

top_enzymes <- as.numeric(opt$top_enzymes); 
extra_site_information <- (as.character(opt$sites_info) =='Y');



################################################################
################## Functions ###################################
################################################################

## Function: return a string which contains the isoschizomer choices that are not methylation-sensitive. 
#            The output format is the following:
#
#               - Individual enzyme. Different methylation-insensitive isoschizomers separated by OR. e.g. BsiSI OR MspI
#               - 2-enzyme combination. The different methylation-insensitive options for each one of the 2 enzymes are 
#                                        enclosed in parenthesis. e.g. (BsiSI OR MspI) AND (BspQI OR LguI OR SapI).
#
# original_enz_string: string containing the enzyme(s) found in the SNE file (i.e. the names of the isoschizomer families).
#                      e.g. 'BsiSI'  e.g.2. 'BsiSI_BspQI'
# enz_annotation: dataframe with the isoschizomer annotation read from the appropiate file. Columns: family_name, enzyme, methylation_sensitive


find_isoschizomers <- function(original_enz_string, enz_annotation){
  
  n_enzymes <- length(strsplit(original_enz_string, '_')[[1]]);
  
  # For an individual enzyme
  
  if(n_enzymes == 1){
    
    iso_family <- enz_annotation[which(enz_annotation[,1] == original_enz_string),];
    non_ms_iso_enz <- as.character(iso_family[which(iso_family[,3] == 0),2]);
    final_enz_string <- paste(non_ms_iso_enz, collapse=' OR ');
    
  }
  
  # For a 2-enzyme combination
  
  if(n_enzymes == 2){
    
    iso_family_1 <- enz_annotation[which(enz_annotation[,1] == strsplit(original_enz_string, '_')[[1]][1]),];
    non_ms_iso_enz_1 <- as.character(iso_family_1[which(iso_family_1[,3] == 0),2]);
    string_1 <- paste0('(', paste(non_ms_iso_enz_1, collapse=' OR '), ')');
    
    iso_family_2 <- enz_annotation[which(enz_annotation[,1] == strsplit(original_enz_string, '_')[[1]][2]),];
    non_ms_iso_enz_2 <- as.character(iso_family_2[which(iso_family_2[,3] == 0),2]);
    string_2 <- paste0('(', paste(non_ms_iso_enz_2, collapse=' OR '), ')');
    
    final_enz_string <- paste0(string_1, ' AND ', string_2);
    
  }
  
  return(final_enz_string);
  
}


## Function: calculate the experimental size range.
#
# tsr: theoretical size range (in bp), provided as a string with the format 'min_max' (e.g. '40_220')
# adapt_l: length of the adapters, as seen in a gel

calculate_experimental_sr <- function(tsr, adapt_l){
  
  min_e <- as.numeric(strsplit(tsr, '_')[[1]][1]) + adapt_l;
  max_e <- as.numeric(strsplit(tsr, '_')[[1]][2]) + adapt_l;
  
  return(paste0(min_e, '_', max_e));
  
}


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
  
  merged_expanded <- full_join(as.data.frame(expanded_input), as.data.frame(sann), by='Site_ID');
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


## Function: return the robustness measurement as a float. See documentation for details. If the 
#            experimental error is too big and some of the tested size ranges do not make sense,
#            the function returns NA.
#
# fld_info: vector which contains the fragments length distribution for a given enzyme.
#           The first element of the vector has the format ">Enzyme(s)" and the rest "length,number_of_fragments".
# fragments_info: vector with the information regarding the fragments containing the sites of interest.
#                 The first element of the vector has the format ">Enzyme(s)" and the rest "Site_ID(s),start_fragment,length_fragment"
#                 (i.e. collapsed format).
# optimal_sr: optimal theoretical size range (e.g. "40_220")
# rl: read length (in bp) of the potential RRBS experiment. e.g. 100
# exp_error: experimental error (in bp) expected during the size range selection.
# sann: annotation for the sites of interest, as displayed in the sites annotation file. Format:
#        Site_ID   Chr   Coordinate  Weight

calculate_robustness <- function(fld_info, fragments_info, optimal_sr, rl, exp_error, sann){
  
  ## Obtain the different size ranges used to calculate the robustness.
  
  ol <- as.numeric(strsplit(optimal_sr, '_')[[1]]);
  lower_limits <- c(ol[1] - exp_error, ol[1], ol[1] + exp_error);
  upper_limits <- c(ol[2] - exp_error, ol[2], ol[2] + exp_error);
  
  surrounding_sr <- expand.grid(lower_limits, upper_limits);
  surrounding_sr <- apply(surrounding_sr, 1, function(x){
    
    return(paste0(x[1], '_', x[2]));
  });
  
  # Check that all the size ranges make sense. Otherwise return NA.
  
  for(sr in surrounding_sr){
    
    l <- as.numeric(strsplit(sr, '_')[[1]][1]);
    u <- as.numeric(strsplit(sr, '_')[[1]][2]);
    
    if(l < 1 | u < 1 | l >= u){
      
      return(NA);
    }
  }
  
  ## Put input in correct format.
  
  fld_info_f <- fld_info[-1];
  fld_info_f <- matrix(as.numeric(unlist(strsplit(fld_info_f, ','))), ncol=2, byrow=T);
  
  fragments_info_f <- fragments_info[-1];
  fragments_info_f <- matrix(unlist(strsplit(fragments_info_f, ',')), ncol=3, byrow=T);
  fragments_info_f <- expand_fragments_df(fragments_info_f, sann);
  
  ## Calculate Score, NF/1000 and EV values for the different size ranges.
  
  Score_values <- as.numeric(sapply(surrounding_sr, Score_calculation, 
                                    sinfo=fragments_info_f, 
                                    sann=sann, rl=rl, ids=FALSE));
  NF_1000_values <- as.numeric(sapply(surrounding_sr, NF_calculation, fld=fld_info_f));
  
  nsites <- nrow(sann);
  max_Score <- sum(as.numeric(sann[,4]));
  EV_values <- -log10((Score_values * nsites) / (NF_1000_values * 1000 * max_Score));
  
  
  ## Calculate the robustness of the enzyme / enzyme combination
  
  theta <- sum(abs(EV_values - EV_values[5])) / (EV_values[5]);
  R <- exp(-theta);
  
  return(R);
  
}



################################################################
################## Running the pipeline ########################
################################################################


#### 1. Read the input files ####

## Sites annotation

sites_ann <- read.csv(sites_annotation_path);

## SNE files and df files

raw_input <- data.frame();

if(length(Sys.glob(sne_1_path)) && file.exists(Sys.glob(sne_1_path))){
  
  # SNE files    
  
  raw_1 <- read.csv(file=Sys.glob(sne_1_path), header=F);
  raw_input <- rbind(raw_input, raw_1);
  
  # df files
  
  raw_fld_1 <- readLines(paste0(df_files_path, '/fl_distributions_1_enzymes.txt'));
  raw_fragments_1 <- readLines(paste0(df_files_path, '/fragments_of_interest_1_enzymes.txt'));
  
  # df files useful info
  
  fld_chunk_sep_1 <- grep('>', raw_fld_1);
  fragments_chunk_sep_1 <- grep('>', raw_fragments_1);
  fld_headers_1 <- raw_fld_1[fld_chunk_sep_1];
  fragments_headers_1 <- raw_fragments_1[fragments_chunk_sep_1];
  fld_enz_combs_1 <- gsub('>', '', fld_headers_1); 
  fragments_enz_combs_1 <- gsub('>', '', fragments_headers_1);
  
  
  if((length(fld_headers_1) != length(fragments_headers_1))){
    
    stop("The input files for 1 enzyme do not have the same number of headers.");
    
  }
  
  if((!all(fld_enz_combs_1 == fragments_enz_combs_1))){
    
    stop('The enzyme combinations are not the same for the fl_distribution_1_enzymes.txt and the fragments_of_interest_1_enzymes.txt files.');
    
  }
  
}

if(length(Sys.glob(sne_2_path)) && file.exists(Sys.glob(sne_2_path))){

  # SNE file
  
  raw_2 <- read.csv(file=Sys.glob(sne_2_path), header=F);
  raw_input <- rbind(raw_input, raw_2);
  
  # df files
  
  raw_fld_2 <- readLines(paste0(df_files_path, '/fl_distributions_2_enzymes.txt'));
  raw_fragments_2 <- readLines(paste0(df_files_path, '/fragments_of_interest_2_enzymes.txt'));
  
  # df files useful info
  
  fld_chunk_sep_2 <- grep('>', raw_fld_2);
  fragments_chunk_sep_2 <- grep('>', raw_fragments_2);
  fld_headers_2 <- raw_fld_2[fld_chunk_sep_2];
  fragments_headers_2 <- raw_fragments_2[fragments_chunk_sep_2];
  fld_enz_combs_2 <- gsub('>', '', fld_headers_2); 
  fragments_enz_combs_2 <- gsub('>', '', fragments_headers_2);
  
  if((length(fld_headers_2) != length(fragments_headers_2))){
    
    stop("The input files for 2 enzymes do not have the same number of headers.");
    
  }
  
  if((!all(fld_enz_combs_2 == fragments_enz_combs_2))){
    
    stop('The enzyme combinations are not the same for the fl_distribution_2_enzymes.txt and the fragments_of_interest_2_enzymes.txt files.');
    
  }
  
}

if(nrow(raw_input) == 0){
  
  stop("The SNE files could not be found. Quitting ...")
}



#### 2. Obtain the final output dataframe ####

### 2.1. Order and select top enzymes.

top_input <- raw_input[order(raw_input[,6]),];

if(nrow(top_input) > top_enzymes){
  
  top_input <- top_input[1:top_enzymes,];
  
}

### 2.2. Initialize final dataframe.

final_output <- data.frame(matrix(nrow=nrow(top_input), 
                                  ncol=ifelse(extra_site_information, 12, 10)));

colnames(final_output)[1:10] <- c('Enzyme(s)', 'Experimental_size_range', 'Theoretical_size_range', 'Score',
                                   '%_max_Score', 'NF/1000', 'Cost_Reduction_Factor', 'Enrichment_Value',
                                   'Robustness', 'C_Score | C_NF/1000');

if(extra_site_information){
  
  colnames(final_output)[11:12] <- c('Number_of_sites', 'Site_IDs');
  
}



### 2.3. Fill in the output dataframe: straight-forward variables.

## Select the isoschizomers that are not methylation-sensitive.

enz_ann <- read.csv(enz_annotation_path); # Enzyme annotation

final_output[,1] <- as.character(sapply(as.character(top_input[,1]), find_isoschizomers, enz_annotation=enz_ann));

## Size ranges

theoretical_sr <- as.character(top_input[,3]);
final_output[,3] <- theoretical_sr;

final_output[,2] <- as.character(sapply(theoretical_sr, calculate_experimental_sr, adapt_l=adapters_size));

## Score, NF/1000 and related variables

final_output[,4] <- top_input[,5]; # Score

max_score <- sum(as.numeric(sites_ann[,4]));
final_output[,5] <- round((final_output[,4] / max_score) * 100, digits=2);  # % of max Score

final_output[,6] <- top_input[,4]; # NF/1000

ref_NF_1000 <- 655.085; # NF/1000 value for MspI in original protocol
final_output[,7] <- round(ref_NF_1000 / final_output[,6], digits=2); # Cost Reduction Factor

final_output[,8] <- top_input[,6];  # EV

C_Score <- sapply(strsplit(as.character(top_input[,2]), '_'), function(x){return(x[2])});
C_NF_1000 <- sapply(strsplit(as.character(top_input[,2]), '_'), function(x){return(x[1])});
final_output[,10] <- paste(C_Score, C_NF_1000, sep=' | '); # C_Score | C_NF_1000



### 2.4. Fill in the output dataframe: robustness measure, (number of sites), (site IDs)

## Loop over the enzymes combinations

all_enz <- as.character(top_input[,1]);

for(enz in all_enz){
  
  i <- which(all_enz == enz);
    
## Extract fld and fragments information for the current enzyme / enzyme combination
  
  # Individual enzymes
  
  if(length(strsplit(enz, '_')[[1]]) == 1){
    
    index_enz <- which(fld_enz_combs_1 == enz);
    
    # fld chunk
    
    start_fld <- fld_chunk_sep_1[index_enz];
    
    if(index_enz == length(fld_enz_combs_1)){
      
      end_fld <- length(raw_fld_1);
      
    }else{
      
      end_fld <- fld_chunk_sep_1[index_enz + 1] - 1;
      
    }
    
    fld_chunk <- raw_fld_1[start_fld:end_fld];
    
    # fragment chunk
    
    start_fragments <- fragments_chunk_sep_1[index_enz];
    
    if(index_enz == length(fragments_enz_combs_1)){
      
      end_fragments <- length(raw_fragments_1);
      
    }else{
      
      end_fragments <- fragments_chunk_sep_1[index_enz + 1] - 1;
      
    }
    
    fragments_chunk <- raw_fragments_1[start_fragments:end_fragments];
  }
  
  # 2-enzyme combinations
  
  if(length(strsplit(enz, '_')[[1]]) == 2){
    
    index_enz <- which(fld_enz_combs_2 == enz);
    
    # fld chunk
    
    start_fld <- fld_chunk_sep_2[index_enz];
    
    if(index_enz == length(fld_enz_combs_2)){
      
      end_fld <- length(raw_fld_2);
      
    }else{
      
      end_fld <- fld_chunk_sep_2[index_enz + 1] - 1;
      
    }
    
    fld_chunk <- raw_fld_2[start_fld:end_fld];
    
    # fragment chunk
    
    start_fragments <- fragments_chunk_sep_2[index_enz];
    
    if(index_enz == length(fragments_enz_combs_2)){
      
      end_fragments <- length(raw_fragments_2);
      
    }else{
      
      end_fragments <- fragments_chunk_sep_2[index_enz + 1] - 1;
      
    }
    
    fragments_chunk <- raw_fragments_2[start_fragments:end_fragments];
    
  }
  

## Calculate the Robustness for the current enzyme / enzyme combination and store it in the final output.

  R <- calculate_robustness(fld_chunk, fragments_chunk, as.character(top_input[i,3]), 
                            read_length, experimental_error, sites_ann);
  
  final_output[i,9] <- R;

  
## Calculate, if necessary, the number of sites of interest and their IDs
  
  if(extra_site_information){
    
    fragments_info_f <- fragments_chunk[-1];
    fragments_info_f <- matrix(unlist(strsplit(fragments_info_f, ',')), ncol=3, byrow=T);
    fragments_info_f <- expand_fragments_df(fragments_info_f, sites_ann);
    
    Score_and_site_IDs <- Score_calculation(sinfo=fragments_info_f, sann=sites_ann, 
                                            size_range=as.character(top_input[i,3]), 
                                            rl=read_length, ids=TRUE);
    site_IDs <- gsub(',', ';', Score_and_site_IDs[[2]]);
    final_output[i,12] <- site_IDs;
    
    nsites <- length(strsplit(site_IDs, ';')[[1]]);
    final_output[i,11] <- nsites;
    
  }
}



### 2.5. Create the final output file

write.table(x=final_output, file=paste0(output_path, '/final_cuRRBS_output.csv'), quote=F,
            sep=',', row.names=F);


#########################################################
########## End of the script ############################
#########################################################
