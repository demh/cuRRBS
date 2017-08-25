###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             19/05/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Assess how well cuRRBS performs when comparing it with experimental data         ####
##### (XmaI-RRBS).                                                                     ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

#### Dependencies ####

library(dplyr);
library(ggplot2);
library(scales);

setwd('~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Figure_4/Figure_4A_B/');


#### Input arguments ####

print('Reading input data ...');

sites_annotation <- read.table('CpG_in_CGI_annotation.csv', sep=',', header=T);
cuRRBS_output <- readLines('final_cuRRBS_output_90_185_rl200.csv');
bismark_output <- read.table('edit_Bismark_cov_file/SRR2721365_all_merged_trimmed_GRCh38_bismark_bt2.bismark_mod.cov', sep='\t', header=F); # We have deleted the sites in chrM
colnames(bismark_output) <- c('Chr', 'Start', 'Coordinate', 'Perc_meth', 'Reads_meth', 'Reads_no_meth');

max_threshold <- 20; # Maximum depth of coverage threshold to consider


### 1. Obtain two dataframes containing the CpG sites that should be theoretically seen
#   and those which should not in the real data (i.e. predictions). #
# Dataframes format: Site_ID,Chr,Coordinate

print('Calculating theoretical dataframes ...');

theoretical_yes_IDs <- strsplit(strsplit(cuRRBS_output[2], ',')[[1]][12], ';')[[1]];
theoretical_yes_df <- sites_annotation[sites_annotation[,1] %in% theoretical_yes_IDs,1:3];
theoretical_no_df <- sites_annotation[!(sites_annotation[,1] %in% theoretical_yes_IDs),1:3];


### 2. Initialise the dataframes which will contain the information needed to construct the plots #

print('Initialising final dataframes ...');

## Dataframe for plot 1 (barplots).
# Dataframe format: category (TP, FP, FN, TN), threshold (depth of coverage >= threshold), number of sites

plot1_df <- data.frame(category=rep(c('TP','FP','FN','TN'), max_threshold),
                       threshold=rep(1:max_threshold, each=4),
                       number_of_sites=rep(NA, 4*max_threshold));

## Dataframe for plot 2 (sensitivity, specificity)
# Dataframe format: type (sensitivity, specificity), threshold (depth of coverage >= threshold), percentage

plot2_df <- data.frame(type=rep(c('sens', 'spec'), max_threshold),
                       threshold=rep(1:max_threshold, each=2),
                       percentage=rep(NA, 2*max_threshold));


### 3. Filter the Bismark coverage file to keep only those sites that are in the sites annotation file ###

print('Filtering Bismark output ...');

filtered_bismark <- inner_join(bismark_output, sites_annotation, by=c("Chr","Coordinate"));


### 4. Calculate the different variables needed to fill the dataframes for different depth of coverage thresholds. #

depth_of_coverage <- filtered_bismark[,5] + filtered_bismark[,6];

for(thr in 1:max_threshold){
  
  print(paste0('Calculating variables for threshold ', thr, ' ...'));
  
  ## Filter the Bismark output for the current depth of coverage threshold (i.e. keep sites with number of reads >= 1)
  
  bismark_thr <- filtered_bismark[depth_of_coverage >= thr,];
  
  ## Format the dataset that contains the sites that are found experimentally.
  
  experimental_yes_df <- data.frame(Site_ID=bismark_thr$Site_ID,
                                    Chr=bismark_thr$Chr,
                                    Coordinate=bismark_thr$Coordinate);
  
  ## Calculate the number of sites in each category.
  
  TP <- nrow(inner_join(theoretical_yes_df, experimental_yes_df, by=c("Site_ID", "Chr", "Coordinate")));
  FP <- nrow(theoretical_yes_df) - TP;
  FN <- nrow(inner_join(theoretical_no_df, experimental_yes_df, by=c("Site_ID", "Chr", "Coordinate")));
  TN <- nrow(theoretical_no_df) - FN;
  
  ## Make a check that everything went OK.
  
  if(TP+FP+FN+TN != nrow(sites_annotation)){
    stop('Something went wrong when calculating the different subsets.');
  }
  
  ## Calculate sensitivity  and specificity.
  
  sensitivity = (TP/(TP+FN))*100;
  specificity = (TN/(FP+TN))*100;
  
  ## Store the values in the given datasets.
  
  plot1_df[(4*thr-3):(4*thr),3] <- c(TP,FP,FN,TN);
  plot2_df[(2*thr-1):(2*thr),3] <- c(sensitivity,specificity);

}


### 5. Make the plots #.

## Plot1: barplots showing the absolute numbers of TP,FP,FN and TN as a function of depth of coverage thresholds.

plot1 <- ggplot(plot1_df, aes(x=threshold, y=number_of_sites, fill=category)) + 
  geom_bar(stat="identity") + 
  scale_y_continuous(breaks= pretty_breaks()) + 
  scale_fill_manual(name='', 
                    values=c(TP = 'darkgreen', TN = 'blue', FP = 'red', FN = 'orange')) + 
  theme_classic() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold")) + 
  xlab('Depth of coverage threshold') + ylab('Number of sites');

ggsave('positives_negatives_barplots_90_185_rl200.pdf', width=8, height=8);

## Plot2: sensitivity and specificity as a function of depth of coverage thresholds.

plot2 <- ggplot(plot2_df, aes(x=threshold, y=percentage, colour=type)) +
  geom_line() + geom_point() + 
  scale_y_continuous(limits = c(0, 100)) + 
  scale_colour_manual(name='',
                      values=c(sens='green', spec='cyan'),
                      labels=c(expression(Sensitivity == frac(TP,TP+FN) %.% 100),
                               expression(Specificity == frac(TN,FP+TN) %.% 100))) + 
  theme_classic() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.5, 0.3),
        legend.margin = margin(),
        legend.box.background = element_rect(colour = "black", size=0.6),
        legend.box.margin = margin(1, 3, 3, 3, "mm"),
        legend.key = element_rect(size = 5, color='white'),
        legend.key.size = unit(2.5, 'lines')) + 
  xlab('Depth of coverage threshold') + ylab('%');

ggsave('sensitivity_specificity_90_185_rl200.pdf', width=8, height=8);


#### End of the script ####
