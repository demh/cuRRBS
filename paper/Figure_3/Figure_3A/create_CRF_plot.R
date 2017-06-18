###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             26/04/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Create a barplot that displays the Cost Reduction Factor (CRF) values for the    ####
##### different in silico systems under study.                                         ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

setwd('~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Figure_3/Figure_3A/');

###########################################################
##################### Dependencies ########################
###########################################################

suppressWarnings(suppressMessages(library(ggplot2)));


################################################################
################## Running the pipeline ########################
################################################################

#### 1. Create a dataframe that contains the CRF for the best enzyme
#    combination in each one of the in silico systems to plot.

file_names <- dir('cuRRBS_outputs/', pattern =".csv");
n <- length(file_names);
cuRRBS_data <- data.frame(System=rep(NA, n), CRF=rep(NA,n));

for(i in 1:n){
  
  file <- file_names[i];
  raw_file <- read.csv(paste0('cuRRBS_outputs/',file), header=TRUE);
  cuRRBS_data[i,1] <- strsplit(file, '_cuRRBS_output_')[[1]][1];
  cuRRBS_data[i,2] <- raw_file[1,7];
}

cuRRBS_data$System <- c('Arabidopsis CHG sites',
                      'Human epigenetic clock',
                      'Human CTCF sites',
                      'Human exon-intron boundaries',
                      'Human imprinted loci',
                      'Human placental imprinted loci',
                      'Mouse iPSCs demethylated',
                      'Mouse iPSCs maintained',
                      'Mouse NRF1 sites');

# Order the data.

cuRRBS_data_ordered <- cuRRBS_data;
cuRRBS_data_ordered$System <- factor(cuRRBS_data_ordered$System, 
                                     levels = c("Arabidopsis CHG sites",
                                                "Mouse iPSCs demethylated",
                                                "Mouse iPSCs maintained",
                                                "Mouse NRF1 sites",
                                                "Human exon-intron boundaries",
                                                "Human epigenetic clock",
                                                "Human imprinted loci",
                                                "Human CTCF sites",
                                                "Human placental imprinted loci"));

# Prepare to add segment legend

cuRRBS_data_ordered$s <-rep('CRF for traditional RRBS', 9);
colnames(cuRRBS_data_ordered) <- c('System', 'CRF', 'CRF for traditional RRBS');


#### 2. Create the barplot ####

CRF_ref <- 30.65; # CRF for the RRBS protocol of reference (RRBS with MspI in human, size selection with beads 20-800 bp)

plot_CRF <- ggplot(cuRRBS_data_ordered, aes(x=System, y=CRF)) +
  geom_bar(stat="identity", width=0.5, color="black", fill=c('green', rep('blue',3), rep('red',5))) +
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = 45, hjust = 0),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_blank(),
        legend.position = c(0.50, 0.9),
        legend.text=element_text(size=15),
        legend.key.size = unit(4, "cm"),
        plot.margin=unit(c(0,3.2,0.5,0.5),"cm")) +
  xlab('') + ylab('Cost Reduction Factor (CRF)') +
  #geom_segment(aes(x=1.5, y=CRF_ref, xend=9.5, yend=CRF_ref, linetype=`CRF for traditional RRBS`), col='dimgrey', size=1) +
  #scale_linetype_manual(values=c("CRF for traditional RRBS"='dotted')) +
  #guides(linetype = guide_legend(override.aes = list(size = 2))) +
  scale_y_reverse() + scale_x_discrete(position = "top") +
  annotate("rect", xmin = 1.5, xmax = 9.5, ymin = 0, ymax = CRF_ref, fill = "grey", alpha=0.5);

ggsave('CRF_plot.pdf', width = 10, height=10);

### End of the script ###


