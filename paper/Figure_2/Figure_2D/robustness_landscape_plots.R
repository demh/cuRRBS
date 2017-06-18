###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             20/03/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Given the fl_distributions_x_enzymes.txt and fragments_of_interest_x_enzymes.txt ####  
##### files obtained with obtain_distributions_and_fragments_old.py, this script       ####
##### generates contour plots with the NF/1000, Score and EV values as a function of   ####
##### the size range. (i.e. landscapes). The EV values surrounding the optimal are also ###
##### shown to illustrate the calculation of the robustness measure.                   ####
###########################################################################################
##### USAGE: Rscript robustness_landscape_plots.R --help                               ####
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
              help="Absolute path to the fl_distribution_x_y.txt file, which is created by \
              obtain_distributions_and_fragments_old.py. COMPULSORY",
              metavar='character'),
  
  make_option(c('-c', '--fragments_file'), type='character', default=NULL, 
              help="Absolute path to the fragments_of_interest_x_y.txt, which is created by \
              obtain_distributions_and_fragments_old.py. COMPULSORY",
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
  
  make_option(c('-w', '--optimal_sr'), type='character', default=NULL,
              help="Optimal size range (in bp) used to calculate the robustness (e.g. '70-160'). COMPULSORY",
              metavar='character'),
  
  make_option(c('-t', '--NF_thresholds'), type='character', default='0.2', 
              help="Float (named C_NF/1000) or vector (e.g. '0.5,1,40') used to calculate the different NF/1000 thresholds \
              which will filter the final output. C_NF/1000 takes values in the interval (0, +Inf). \
              Only those size ranges with NF/1000 <= C_NF/1000 * ref_NF/1000 are kept.  \ 
              By default, ref_NF/1000 is 41177.15 (estimated NF/1000 in WGBS for human genome). DEFAULT='0.2'",
              metavar='character'),
  
  make_option(c('-x', '--Score_thresholds'), type='character', default='0.25', 
              help="Float (named C_Score) or vector e.g. ('0.2,0.5,0.8') used to calculate the different Score thresholds \
              which will filter the final output. C_Score takes values in the interval (0,1]. 
              Only those size ranges with Score > C_Score * max_Score are kept (max_Score is the Score that would be \ 
              obtained if all the sites of interest were captured in the sequencing). DEFAULT='0.25'",
              metavar='character'),
  
  make_option(c('-a', '--min_size'), type='integer', default=0, 
              help="Minimum size of the fragments (in bp) to be considered when screening for \
              different size ranges. DEFAULT=0",
              metavar='integer'),
  
  make_option(c('-b', '--max_size'), type='integer', default=1000, 
              help="Maximum size of the fragments (in bp) to be considered when screening for \
              different size ranges. DEFAULT=1000",
              metavar='integer'),
  
  make_option(c('-d', '--exp_error'), type='integer', default=20,
              help="Experimental error assumed when performing the size range selection (in base pairs). DEFAULT=20",
              metavar='integer')
  
  
  # make_option(c('-z', '--plot_filtered'), type='logical', default=FALSE, 
  #             help="Generate the contour plots for the NF/1000, Score and EV values as a function \
  #               of the size ranges only for those values that satisfy the thresholds. DEFAULT=FALSE",
  #             metavar='logical')
  
  );



## Prepare option list and general description of the script.

opt_parser <- OptionParser(option_list=option_list,
                           description="\nGiven the fl_distributions_x_enzymes.txt and fragments_of_interest_x_enzymes.txt \
                           files obtained with obtain_distributions_and_fragments_old.py, this script \
                           generates contour plots with the NF/1000, Score and EV values as a function of \
                           the size range. (i.e. landscapes). The EV values surrounding the optimal are also \
                           shown to illustrate the calculation of the robustness measure.");


opt <- parse_args(opt_parser);

## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$fld_file) | is.null(opt$fragments_file) |
   is.null(opt$annotation_file) | is.null(opt$output_path) | 
   is.null(opt$read_length) | is.null(opt$optimal_sr)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Parameters.

print('Reading parameters ...');

# Output path / working directory.

output_path <- as.character(opt$output_path);
#output_path <- '~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Figure_2/Figure_2D/';
setwd(output_path);


# Minimum and maximum lengths of fragments to be considered.

min_size <- as.numeric(opt$min_size);
#min_size <- 0;
max_size <- as.numeric(opt$max_size);
#max_size <- 1000;

# Steps to be taken by the sliding window.

window_step <- as.numeric(opt$exp_error);
#window_step <- 20;

# Parameters regarding the changing window size.

min_window_size <- 1;
max_window_size <- max_size - min_size + 1;
window_size_step <- as.numeric(opt$exp_error);
#window_size_step <- 20;

# Other parameters.

read_length <- as.numeric(opt$read_length);
#read_length <- 75;
C_Score <- as.numeric(strsplit(opt$Score_thresholds, ',')[[1]]);
#C_Score <- 0.25;
C_NF_1000 <- as.numeric(strsplit(opt$NF_thresholds, ',')[[1]]);
#C_NF_1000 <- 0.2;
#activate_plotting_filter <- as.logical(opt$plot_filtered); # Create contour plots or not for the filtered landscapes
#activate_plotting_filter <- TRUE;

# Input files.

fsd_file_path <- as.character(opt$fld_file);
#fsd_file_path <- 'fl_distributions_MspI_BspQI.txt';
si_file_path <- as.character(opt$fragments_file);
#si_file_path <- 'fragments_of_interest_MspI_BspQI.txt';
sites_annotation <- read.table(as.character(opt$annotation_file), sep=',', header=T);
#sites_annotation <- read.table('PlacentalDMRs_hg38_coordinates.CG_positions_cuRRBS_edited.csv', sep=',', header=T);
sites_annotation[,1] <- as.character(sites_annotation[,1]); 

# Final thresholds

max_score <- sum(as.numeric(sites_annotation[,4]));
Score_thr <- C_Score * max_score;

ref_NF_1000 <- 41177.15;  # It can be modified to compare with other protocol.
NF_thr <- C_NF_1000 * ref_NF_1000;

# Optimal size range

opt_sr_input <- as.character(opt$optimal_sr);
#opt_sr_input <- '60-540';

xopt <- as.numeric(strsplit(opt_sr_input, '-')[[1]][1]); # Lower limit in optimal size range
yopt <- as.numeric(strsplit(opt_sr_input, '-')[[1]][2]) - as.numeric(strsplit(opt_sr_input, '-')[[1]][1]) + 1; # Window size in optimal size range
opt_size_range <- paste0(xopt, '-', xopt+yopt-1, ' bp');
exp_error <- as.numeric(opt$exp_error); # Experimental error used to calculate the robustness
#exp_error <- 20;



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

# fsd: numeric matrix which contains the fragment size distribution. First column: fragment size. Second column: number of fragments 
#      with that size.
# size_range: range of fragment sizes selected in the RRBS protocol (e.g. 40_220). If NA, then the function returns NA.


NF_calculation <- function(fsd, size_range){
  
  if(is.na(size_range)){
    
    return(NA);
    
  }else{
    
  
  
    min_size_f <- as.numeric(strsplit(size_range, '_')[[1]][1]);
    max_size_f <- as.numeric(strsplit(size_range, '_')[[1]][2]);
    selected_fragments <- fsd[which(fsd[,1] >= min_size_f & fsd[,1] <= max_size_f),];
    
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

# sinfo: information about the fragments where the sites of interest are located. The format is
#       similar to the one found in sites_info_enzymes.csv files. The last four
#       columns must be 'numeric'.
#       Row_number   Chromosome      Coordinate    Fragment_start    Fragment_end    Fragment_length
# sann:  annotation for the sites of interest, as displayed in the new_sites_annotation.csv file.
#        Site_ID   Chr   Coordinate  Weight
# size_range: range of fragment sizes selected in the RRBS protocol (e.g. 40_220). If NA, then the function returns NA.
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
  NF_1000_values <- as.numeric(sapply(surrounding_sr, NF_calculation, fsd=fld_info_f));
  
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

#### 1. Size ranges matrix: size ranges to test. Create a matrix whose rows are the different minimum sizes
#  (i.e. lower limits) considered in the size ranges and the columns are the different window sizes. The top left half of
#  the matrix will be filled in with the appropiate numbers and the right bottom half will be NA values.

print('Creating size ranges matrix ...');

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


#### 2. Create three matrices with the same size as 'size_ranges_to_test', that contain,
#    respectively, the NF/1000, Score and EV values for the different size ranges assessed.


## Read the data from the input files and format it.

print('Reading input data ...');

current_fsd <- readLines(fsd_file_path, n=-1);
current_fsd_robustness <- current_fsd;
current_fsd <- current_fsd[-1];
current_fsd <- matrix(as.numeric(unlist(strsplit(current_fsd, ','))), ncol=2, byrow=T);

current_si <- readLines(si_file_path, n=-1);
current_si_robustness <- current_si;
current_si <- current_si[-1];
current_si <- matrix(unlist(strsplit(current_si, ',')), ncol=3, byrow=T);
current_si_df <- expand_fragments_df(current_si, sites_annotation);


# Calculate the NF/1000 matrix. 

print('Creating NF/1000 matrix ...');

NF_1000_matrix <- matrix(mapply(NF_calculation, size_ranges_to_test, 
                         MoreArgs = list(fsd=current_fsd)), 
                         nrow=nrow(size_ranges_to_test),
                         ncol=ncol(size_ranges_to_test));
rownames(NF_1000_matrix) <- lower_limits_vector;
colnames(NF_1000_matrix) <- window_sizes_vector;
NF_1000_matrix[which(NF_1000_matrix==0)] <- NF_1000_matrix[which(NF_1000_matrix==0)] + 0.001; # Avoid values of 0 for posterior calculations
  

# Calculate the Score matrix.

print('Creating Score matrix ...');

Score_matrix <- matrix(unlist(mapply(Score_calculation, size_ranges_to_test, 
                              MoreArgs = list(sinfo=current_si_df, sann=sites_annotation, rl=read_length, ids=FALSE))), 
                       nrow=nrow(size_ranges_to_test),
                       ncol=ncol(size_ranges_to_test));
rownames(Score_matrix) <- lower_limits_vector;
colnames(Score_matrix) <- window_sizes_vector;


# Calculate the EV (enrichment value) matrix. 

print('Creating EV matrix ...');

n_sites <- nrow(sites_annotation);
max_score <- sum(as.numeric(sites_annotation[,4]));
EV_matrix <- - log10((Score_matrix / (NF_1000_matrix * 1000)) *
                       (n_sites/max_score));
rownames(EV_matrix) <- lower_limits_vector;
colnames(EV_matrix) <- window_sizes_vector;

# Transform the EV matrix so it is plottable.
# Those EV values = Inf are converted into the maximum EV value found in the matrix.

if(length(EV_matrix[is.finite(EV_matrix)]) > 0){ 
  
  max_EV_value <- max(EV_matrix[is.finite(EV_matrix)], na.rm=T);
  
}else{ # When all the elements in the EV_matrix are Inf (for really bad enzymes, when all the Scores are 0)
  
  max_EV_value <- 9999; # Fix arbitrary high value to allow the plotting
  
}


EV_matrix_t <- EV_matrix;
EV_matrix_t[is.infinite(EV_matrix_t)] <- max_EV_value; 


# Calculate the robustness.

robustness <- round(calculate_robustness(current_fsd_robustness, current_si_robustness, 
                                   paste0(xopt, '_', xopt+yopt-1), read_length, exp_error, 
                                   sites_annotation), 6);



#### 3. Create contour plots for the three complete matrices, showing how the values of NF/1000,
#       Score and EV changes with the size range.

# X: lower limit  Y: window size  Z: matrix values


print('Creating contour plots for the complete landscapes ...');

## Obtain enzymes names.

bits_of_enzyme_name <- strsplit(strsplit(fsd_file_path, '/')[[1]][length(strsplit(fsd_file_path, '/')[[1]])], '[_\\.]')[[1]];
start_bits <- which(bits_of_enzyme_name == 'distributions') + 1;
end_bits <- which(bits_of_enzyme_name == 'txt') - 1;
enzyme_combination <- paste(bits_of_enzyme_name[start_bits:end_bits], collapse='+');
enzyme_combination_2 <- paste(bits_of_enzyme_name[start_bits:end_bits], collapse='_');

## Other useful plotting parameters.

min_x <- min(lower_limits_vector);
max_x <- max(lower_limits_vector);
min_y <- min(window_sizes_vector);
max_y <- max(window_sizes_vector);


## NF/1000 matrix

# Plot.

NF_plot_name <- paste0('NF_plot_all_', enzyme_combination_2, '.pdf');
pdf(NF_plot_name, height=8, width=8);

levels.ramp<-40

filled.contour(x=as.numeric(rownames(NF_1000_matrix)), 
               y=as.numeric(colnames(NF_1000_matrix)), 
               z=NF_1000_matrix, 
               nlevels = levels.ramp,
               color.palette=colorRampPalette(c("hotpink", "lightgrey", "darkgreen"), space = "Lab"),
               xlab='Lower limit in size range (bp)', ylab='Window size (bp)',
               main=paste0('NF/1000 landscape for ', enzyme_combination, 
                           '\n Optimal size range: ', opt_size_range,
                           '\n Robustness (d = ', exp_error, ' bp): ', robustness),
               cex.lab=1.3,
               frame.plot = F, axes=T,
               plot.axes = { axis(1, seq(min_x, max_x+window_step-1, by = 100), cex.axis=0.9);
                 axis(2, seq(min_y-1, max_y+window_size_step-1, by = 100), cex.axis=0.9); 
                 points(xopt,yopt, col='green', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
               });
dev.off();




## Score matrix

# Plot.

Score_plot_name <- paste0('Score_plot_all_', enzyme_combination_2, '.pdf');

pdf(Score_plot_name, height=8, width=8);

levels.ramp<-40

filled.contour(x=as.numeric(rownames(Score_matrix)), 
               y=as.numeric(colnames(Score_matrix)), 
               z=Score_matrix, nlevels = levels.ramp,
               color.palette=colorRampPalette(c("purple4","lightgrey", "saddlebrown"), space = "Lab"),
               xlab='Lower limit in size range (bp)', ylab='Window size (bp)',
               main=paste0('Score landscape for ', enzyme_combination, 
                           '\n Optimal size range: ', opt_size_range,
                           '\n Robustness (d = ', exp_error, ' bp): ', robustness),
               cex.lab=1.3,
               frame.plot = F, axes=T,
               plot.axes = { axis(1, seq(min_x, max_x+window_step-1, by = 100), cex.axis=0.9);
                 axis(2, seq(min_y-1, max_y+window_size_step-1, by = 100), cex.axis=0.9); 
                 points(xopt,yopt, col='green', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
               });

dev.off();


## EV matrix

# Plot

EV_plot_name <- paste0('EV_plot_all_', enzyme_combination_2, '.pdf');

pdf(EV_plot_name, height=8, width=8);

levels.ramp<-40

filled.contour(x=as.numeric(rownames(EV_matrix_t)), 
               y=as.numeric(colnames(EV_matrix_t)), 
               z=EV_matrix_t, nlevels = levels.ramp,
               color.palette = colorRampPalette(c("blue", "lightgrey", "red"), space = "Lab"),
               xlab='Lower limit in size range (bp)', ylab='Window size (bp)',
               main=paste0('EV landscape for ', enzyme_combination, 
                           '\n Optimal size range: ', opt_size_range,
                           '\n Robustness (d = ', exp_error, ' bp): ', robustness),
               cex.lab=1.3,
               frame.plot = F, axes=T,
               plot.axes = { axis(1, seq(min_x, max_x+window_step-1, by = 100), cex.axis=0.9);
                 axis(2, seq(min_y-1, max_y+window_size_step-1, by = 100), cex.axis=0.9); 
                 points(xopt,yopt, col='green', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt-exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
                 points(xopt+exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
               });

dev.off();
  



# #### 4. Filter the matrices according to the NF/1000 and Score thresholds.
# # Those elements that do not satisfy both thresholds will be NAs in all the 3 matrices.
# # Store the results for the different thresholds in a list.
# 
# # If the thresholds were not specified correctly: quit the script
# 
# if(length(NF_thr) != length(Score_thr)){
#   
#   stop('The thresholds for NF/1000 and Score are not properly specified. Quit the script !!');
#   
# }
# 
# 
# # Initialize the lists which will contain the filtered matrices for the different thresholds.
# 
# NF_filtered_list <- list();
# Score_filtered_list <- list();
# EV_filtered_list <- list();
# 
# # Iterate over the different thresholds.
# 
# for(i in 1:length(Score_thr)){
#   
#   print(paste0('Calculating filtered matrices for NF_thr=', NF_thr[i], 
#                ' and Score_thr=', Score_thr[i], ' ...'));
#   
#   # Calculate a matrix that will be used for filtering.
#   
#   filter_matrix <- (NF_1000_matrix <= NF_thr[i]) & (Score_matrix > Score_thr[i]);
#   filter_matrix[filter_matrix==FALSE] <- NA;
#   n_rows_f <- dim(filter_matrix)[1];
#   n_cols_f <- dim(filter_matrix)[2];
#   
#   # Filter the 3 matrices.
#   
#   NF_filtered <- matrix(NF_1000_matrix[filter_matrix], nrow=n_rows_f, ncol=n_cols_f);
#   rownames(NF_filtered) <- rownames(NF_1000_matrix);
#   colnames(NF_filtered) <- colnames(NF_1000_matrix);
#   
#   Score_filtered <- matrix(Score_matrix[filter_matrix], nrow=n_rows_f, ncol=n_cols_f);
#   rownames(Score_filtered) <- rownames(Score_matrix);
#   colnames(Score_filtered) <- colnames(Score_matrix);
#   
#   EV_filtered <- matrix(EV_matrix_t[filter_matrix], nrow=n_rows_f, ncol=n_cols_f);
#   rownames(EV_filtered) <- rownames(EV_matrix_t);
#   colnames(EV_filtered) <- colnames(EV_matrix_t);
#   
#   # Store the matrices in the lists. If all the elements in the matrix are NAs, store NaN in the list.
#   
#   if(sum(is.na(filter_matrix)) == (n_rows_f * n_cols_f)){
#     
#     NF_filtered_list[[i]] <- NA;
#     Score_filtered_list[[i]] <- NA;
#     EV_filtered_list[[i]] <- NA;
#   
#   }else{
#     
#     NF_filtered_list[[i]] <- NF_filtered;
#     Score_filtered_list[[i]] <- Score_filtered;
#     EV_filtered_list[[i]] <- EV_filtered;
#     
#   }
# }
# 
# 
# 
# #### 5. Create contour plots for the filtered matrices, showing how the values of NF/1000,
# #       Score and EV changes with the size range.
# 
# # X: lower limit  Y: window size  Z: matrix values
# 
# 
# if(activate_plotting_filter){
#   
#   print('Creating contour plots for the filtered landscapes ...');
#   
#   ## Obtain enzymes names.
#   
#   bits_of_enzyme_name <- strsplit(strsplit(fsd_file_path, '/')[[1]][length(strsplit(fsd_file_path, '/')[[1]])], '[_\\.]')[[1]];
#   start_bits <- which(bits_of_enzyme_name == 'distributions') + 1;
#   end_bits <- which(bits_of_enzyme_name == 'txt') - 1;
#   enzyme_combination <- paste(bits_of_enzyme_name[start_bits:end_bits], collapse='+');
#   enzyme_combination_2 <- paste(bits_of_enzyme_name[start_bits:end_bits], collapse='_');
#   
#   ## Other useful plotting parameters.
#   
#   min_x <- min(lower_limits_vector);
#   max_x <- max(lower_limits_vector);
#   min_y <- min(window_sizes_vector);
#   max_y <- max(window_sizes_vector);
#   
#   for(i in 1:length(NF_thr)){
#     
#     if(is.matrix(NF_filtered_list[[i]])) {  # Plot only if the matrices were not all NAs
#     
#       ## NF/1000 matrix
#       
#       # Plot.
#       
#       NF_plot_name <- paste0('NF_plot_filtered_', NF_thr[i], '_', Score_thr[i], 
#                              '_', enzyme_combination_2, '.pdf');
#       
#       pdf(NF_plot_name, height=8, width=8);
#       
#       filled.contour(x=as.numeric(rownames(NF_filtered_list[[i]])), 
#                      y=as.numeric(colnames(NF_filtered_list[[i]])), 
#                      z=NF_filtered_list[[i]], 
#                      color.palette = colorRampPalette(c("red", "lightgrey", "blue"), space = "Lab"),
#                      xlab='Lower limit in size range (bp)', ylab='Window size (bp)',
#                      main=paste0('NF/1000 landscape for ', enzyme_combination, '. NF/1000 <= ',
#                                  NF_thr[i], '; Score > ', Score_thr[i], '.'),
#                      cex.lab=1.3,
#                      frame.plot = F, axes=T,
#                      plot.axes = { axis(1, seq(min_x, max_x+window_step-1, by = 100), cex.axis=0.9);
#                        axis(2, seq(min_y-1, max_y+window_size_step-1, by = 100), cex.axis=0.9); 
#                        points(xopt,yopt, col='green', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                      });
#       
#       dev.off();
#       
#       
#       ## Score matrix
#       
#       # Plot.
#       
#       Score_plot_name <- paste0('Score_plot_filtered_', NF_thr[i], '_', Score_thr[i], 
#                                 '_', enzyme_combination_2, '.pdf');
#       
#       pdf(Score_plot_name, height=8, width=8);
#       
#       filled.contour(x=as.numeric(rownames(Score_filtered_list[[i]])), 
#                      y=as.numeric(colnames(Score_filtered_list[[i]])), 
#                      z=Score_filtered_list[[i]], 
#                      color.palette = colorRampPalette(c("blue", "lightgrey", "red"), space = "Lab"),
#                      xlab='Lower limit in size range (bp)', ylab='Window size (bp)',
#                      main=paste0('Score landscape for ', enzyme_combination, '. NF/1000 <= ',
#                                  NF_thr[i], '; Score > ', Score_thr[i], '.'),
#                      cex.lab=1.3,
#                      frame.plot = F, axes=T,
#                      plot.axes = { axis(1, seq(min_x, max_x+window_step-1, by = 100), cex.axis=0.9);
#                        axis(2, seq(min_y-1, max_y+window_size_step-1, by = 100), cex.axis=0.9); 
#                        points(xopt,yopt, col='green', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                      });
#       
#       dev.off();
#       
#       
#       ## EV matrix
#       
#       # Plot
#       
#       EV_plot_name <- paste0('EV_plot_filtered_', NF_thr[i], '_', Score_thr[i], 
#                              '_', enzyme_combination_2, '.pdf');
#       
#       pdf(EV_plot_name, height=8, width=8);
#       
#       filled.contour(x=as.numeric(rownames(EV_filtered_list[[i]])), 
#                      y=as.numeric(colnames(EV_filtered_list[[i]])), 
#                      z=EV_filtered_list[[i]], 
#                      color.palette = colorRampPalette(c("red", "lightgrey", "blue"), space = "Lab"),
#                      xlab='Lower limit in size range (bp)', ylab='Window size (bp)',
#                      main=paste0('EV landscape for ', enzyme_combination, '. NF/1000 <= ',
#                                  NF_thr[i], '; Score > ', Score_thr[i], '.'),
#                      cex.lab=1.3,
#                      frame.plot = F, axes=T,
#                      plot.axes = { axis(1, seq(min_x, max_x+window_step-1, by = 100), cex.axis=0.9);
#                        axis(2, seq(min_y-1, max_y+window_size_step-1, by = 100), cex.axis=0.9); 
#                        points(xopt,yopt, col='green', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt-exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt-exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                        points(xopt+exp_error,yopt+exp_error, col='black', cex=1.5, pch='*', lwd=3);
#                      });
#       
#       dev.off();
#     }  
#   }
# }


print('The script finished correctly.');


#########################################################
########## End of the script ############################
#########################################################
