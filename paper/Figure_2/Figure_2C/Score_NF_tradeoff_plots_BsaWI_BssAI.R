###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             14/02/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Given the fl_distributions_x_enzymes.txt and fragments_of_interest_x_enzymes.txt ####  
##### files obtained with the obtain_distributions_and_fragments_old.py script, calculate #
##### the % max Score and NF/1000 values for different size ranges and generate the plot ##
##### that shows the tradeoff between these 2 variables. Furthermore, append in an     ####
##### output file the Pearson correlation coefficient between these 2 variables.       ####
###########################################################################################
##### USAGE: Rscript calculate_score_NF.R --help                                       ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

suppressWarnings(suppressMessages(library(optparse)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(ggplot2)));


###########################################################
############## Command-line arguments #####################
###########################################################

## Create arguments and help in a Pythonic style.

option_list <-  list(
  
  make_option(c('-f', '--fld_file'), type='character', default=NULL, 
              help="Absolute path to the fl_distributions_x_enzymes.txt file, which is created by \
                obtain_distributions_and_fragments_old.py. It contains information for one enzyme (combination). COMPULSORY",
              metavar='character'),
  
  make_option(c('-c', '--fragments_file'), type='character', default=NULL, 
              help="Absolute path to the fragments_of_interest_x_enzymes.txt, which is created by \
                obtain_distributions_and_fragments.py_old. It contains information for one enzyme (combination). COMPULSORY",
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
  
  make_option(c('-a', '--min_size'), type='integer', default=20, 
              help="Minimum size of the fragments (in bp) to be considered when screening for \
                different size ranges. DEFAULT=20",
              metavar='integer'),
  
  make_option(c('-b', '--max_size'), type='integer', default=1000, 
              help="Maximum size of the fragments (in bp) to be considered when screening for \
                different size ranges. DEFAULT=1000",
              metavar='integer'),
  
  make_option(c('-s', '--window_step'), type='integer', default=20, 
              help="Step (in bp) taken every time the window slides. DEFAULT=20",
              metavar='integer'),
  
  make_option(c('-w', '--window_width_step'), type='integer', default=20, 
              help="Step taken when sampling different widths for the sliding window (in bp).   \
                DEFAULT=20",
              metavar='integer'),
  
  make_option(c('-z', '--create_plot'), type='logical', default=TRUE, 
              help="Should the script create the plot besides calculating the correlation coefficient ? \ 
                DEFAULT=TRUE",
              metavar='logical')
);



## Prepare option list and general description of the script.

opt_parser <- OptionParser(option_list=option_list,
                           description="\nGiven the fl_distributions_x_enzymes.txt and fragments_of_interest_x_enzymes.txt \
files obtained with the obtain_distributions_and_fragments.py script, calculate \
the Scores and NF/1000 values for different size ranges and generate the plots \
that show the tradeoff between these 2 variables. Furthermore, append in an \
output file the Pearson correlation coefficient between these 2 variables.");                                                                         


opt <- parse_args(opt_parser);

## Check that all the compulsory command-line arguments are provided.

if(is.null(opt$fld_file) | is.null(opt$fragments_file) |
   is.null(opt$annotation_file) | is.null(opt$output_path) | 
   is.null(opt$read_length)){
  
  print_help(opt_parser);
  stop("Please, make sure that you have provided all the COMPULSORY arguments.", call.=FALSE);
  
}


## Parameters.

print('Reading parameters ...');

# Output path / working directory.

output_path <- as.character(opt$output_path);
#output_path <- '~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/trial/';
setwd(output_path);


# Minimum and maximum lengths of fragments to be considered.

min_size <- as.numeric(opt$min_size);
#min_size <- 20;
max_size <- as.numeric(opt$max_size);
#max_size <- 1000;

# Steps to be taken by the sliding window.

window_step <- as.numeric(opt$window_step);
#window_step <- 20;

# Parameters regarding the changing window size.

min_window_size <- 1;
max_window_size <- max_size - min_size + 1;
window_size_step <- as.numeric(opt$window_width_step);
#window_size_step <- 20;
  
# Other parameters.
  
read_length <- as.numeric(opt$read_length);
#read_length <- 75;
perc_max_Score_thr <- 25; # Score threshold (it can be changed)
NF_1000_thr <- 150; # NF/1000 threshold (it can be changed). 

# Input files.

fsd_file_path <- as.character(opt$fld_file);
#fsd_file_path <- 'fl_distributions_1_enzymes.txt';
si_file_path <- as.character(opt$fragments_file);
#si_file_path <- 'fragments_of_interest_1_enzymes.txt';
sites_annotation <- read.table(as.character(opt$annotation_file), sep=',', header=T);
#sites_annotation <- read.table('PlacentalDMRs_hg38_coordinates.CG_positions_cuRRBS_edited.csv', sep=',', header=T);
sites_annotation[,1] <- as.character(sites_annotation[,1]);

# Create plot?

create_plot <- as.logical(opt$create_plot);
#create_plot <- FALSE;

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
  
  if(length(expanded_input) == 3){
    
    expanded_input <- matrix(expanded_input, nrow=1, ncol=3);
  }
    
  colnames(expanded_input) <- c('Site_ID', 'Fragment_start', 'Fragment_length');
    
  for(multi_row in grep('_', input_fragments_df[,1])){
    
    sites <- strsplit(input_fragments_df[multi_row,1], '_')[[1]];
    expanded_sites <- cbind(sites, rep(input_fragments_df[multi_row,2], length(sites)), 
                            rep(input_fragments_df[multi_row,3], length(sites)));
    expanded_input <- rbind(expanded_input, expanded_sites);
    
  }
  
  # Merge the information from the expanded input with the annotation information.
  
  expanded_input_df <- as.data.frame(expanded_input);
  expanded_input_df[,1] <- as.character(expanded_input_df[,1]);
  merged_expanded <- full_join(expanded_input_df, sann, by='Site_ID');
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
enzyme_fsd <- gsub('>', '', current_fsd[1]);
current_fsd <- current_fsd[-1];
current_fsd <- matrix(as.numeric(unlist(strsplit(current_fsd, ','))), ncol=2, byrow=T);


current_si <- readLines(si_file_path, n=-1);
enzyme_si <- gsub('>', '', current_si[1]);
current_si <- current_si[-1];
current_si <- matrix(unlist(strsplit(current_si, ',')), ncol=3, byrow=T);
current_si_df <- expand_fragments_df(current_si, sites_annotation);

if(enzyme_fsd != enzyme_si){
  
  stop("The fl_distributions_x_enzymes.txt and fragments_of_interest_x_enzymes.txt contain \
different enzyme (combinations).");
}


# Calculate the NF/1000 matrix. 

print('Creating NF/1000 matrix ...');

NF_1000_matrix <- matrix(mapply(NF_calculation, size_ranges_to_test, 
                         MoreArgs = list(fld=current_fsd)), 
                         nrow=nrow(size_ranges_to_test),
                         ncol=ncol(size_ranges_to_test));
rownames(NF_1000_matrix) <- lower_limits_vector;
colnames(NF_1000_matrix) <- window_sizes_vector;
NF_1000_matrix[which(NF_1000_matrix==0)] <- NF_1000_matrix[which(NF_1000_matrix==0)] + 0.001; # Avoid values of 0 for posterior calculations

# Calculate the % max Score matrix.

print('Creating % max Score matrix ...');

Score_matrix <- matrix(unlist(mapply(Score_calculation, size_ranges_to_test, 
                              MoreArgs = list(sinfo=current_si_df, sann=sites_annotation, rl=read_length, ids=FALSE))), 
                       nrow=nrow(size_ranges_to_test),
                       ncol=ncol(size_ranges_to_test));
rownames(Score_matrix) <- lower_limits_vector;
colnames(Score_matrix) <- window_sizes_vector;

max_Score <- sum(as.numeric(sites_annotation[,4]));
perc_max_Score_matrix <- (Score_matrix / max_Score) * 100;


# Calculate the EV (enrichment value) matrix. 

print('Creating EV matrix ...');

n_sites <- nrow(sites_annotation);
EV_matrix <- - log10((Score_matrix / (NF_1000_matrix * 1000)) *
                       (n_sites/max_Score));
rownames(EV_matrix) <- lower_limits_vector;
colnames(EV_matrix) <- window_sizes_vector; 



#### 3. Create plot which shows the trade-off between the % max Score and the NF/1000 value for
#       different size ranges and a given enzyme (combination).


## Put all the data necessary for the plot in the right format.

final_NF_1000 <- as.numeric(NF_1000_matrix)[!is.na(as.numeric(NF_1000_matrix))];
final_perc_max_Score <- as.numeric(perc_max_Score_matrix)[!is.na(as.numeric(perc_max_Score_matrix))];
final_EV <- as.numeric(EV_matrix)[!is.na(as.numeric(EV_matrix))];
final_size_ranges <- as.character(size_ranges_to_test)[!is.na(as.numeric(EV_matrix))];
final_colours <- ifelse(final_NF_1000 <= NF_1000_thr & final_perc_max_Score > perc_max_Score_thr, 'Passed filtering', 'Did not pass filtering');

if(sum(final_colours == 'Passed filtering') > 0){
  find_optimal_sr <- which(final_EV == min(final_EV[final_colours == 'Passed filtering']));
  final_colours[find_optimal_sr] <- 'Optimal size range';
  optimal_NF <- final_NF_1000[find_optimal_sr];
  optimal_Score <- final_perc_max_Score[find_optimal_sr];
  optimal_size_range <- final_size_ranges[find_optimal_sr];
  arrow_data <- data.frame(x=optimal_NF-24, y=optimal_Score+8, xend=optimal_NF-3, yend=optimal_Score+0.5);
}

final_data <- data.frame(NF_1000=final_NF_1000, Score=final_perc_max_Score, Filtering=final_colours);
corr_coeff <- round(cor(final_data$NF_1000, final_data$Score, method="pearson"), 4);


## Make the plot.

if(create_plot){

  print('Plotting ...');
  
  myColors <- c('cornflowerblue', 'blue', 'darkorange');
  names(myColors) <- c('Did not pass filtering', 'Passed filtering', 'Optimal size range');
  colScale <- scale_colour_manual(name = "Filtering",values = myColors);
  
  plot_tradeoff <- ggplot(final_data, aes(x=NF_1000, y=final_perc_max_Score)) + 
    geom_point(size=0.5, aes(colour=Filtering)) +
    colScale + 
    theme_classic() +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = c(0.50, 0.68),
          legend.margin = margin(),
          legend.box.background = element_rect(colour = "black", size=0.6),
          legend.box.margin = margin(3, 3, 3, 3, "mm"),
          legend.text=element_text(size=10),
          legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    xlab("NF/1000") + ylab("% of maximum Score") + 
    labs(title = gsub('_', ' + ', enzyme_fsd), 
         subtitle = paste0("Pearson's correlation coefficient: ", corr_coeff)) + 
    scale_y_continuous(limit = c(0, 100)) + scale_x_continuous(limit=c(0,180)) +
    geom_hline(aes(yintercept=perc_max_Score_thr), col='blueviolet', linetype='longdash') +
    geom_vline(aes(xintercept=NF_1000_thr), col='red', linetype='longdash');
    
  
  if(sum(final_colours == 'Passed filtering') > 0){
   
    plot_tradeoff <- plot_tradeoff + 
      geom_segment(data=arrow_data, aes(x=x, y=y, xend=xend, yend=yend),
                                                  arrow=arrow(length = unit(0.2, "cm")), col='darkorange', size=0.7) + 
      annotate("text", x=optimal_NF-38, y=optimal_Score+10, label= paste0(gsub('_', '-', optimal_size_range), ' bp'),
               col='darkorange');
  }
  
  ggsave(paste0('Score_NF_tradeoff_plot_', enzyme_fsd, '.pdf'), width = 7, height=7);

}


#### 4. Append in an output file the Pearson correlation coefficient between NF/1000 and Score
#       for the enzyme under study.

print('Storing correlation coefficient ...');

write(paste0(enzyme_fsd, ',', corr_coeff), file='correlation_NF_Score.csv', append=TRUE);

print('The script finished correctly.');


#########################################################
########## End of the script ############################
#########################################################
