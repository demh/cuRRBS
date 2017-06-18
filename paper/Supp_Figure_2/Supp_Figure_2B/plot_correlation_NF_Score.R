###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             21/03/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Create violin plots that show the distribution of Pearson correlation coefficients ##
##### between Score and NF/1000 variables for different size ranges and enzymes.       ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

library(ggplot2);

setwd('~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Supp_Figure_2/Supp_Figure_2B/');

#### Read the input data ####

raw_data <- read.csv('correlation_NF_Score.csv', header=F);

#### Format the data before plotting ####

enzyme_class <- ifelse(grepl('_', raw_data[,1]), '2-enzyme\ncombinations', 'Individual\nenzymes');
format_data <- data.frame(type=enzyme_class, r=raw_data[,2]);
add_all <- format_data;
add_all[,1] <- 'All';
format_data <- rbind(format_data,add_all);
format_data <- format_data[!is.na(format_data[,2]),];

#### Make the violin plots ####

violin_plot <- ggplot(format_data, aes(x=type, y=r)) + 
  geom_violin(fill='blue') + 
  stat_summary(fun.y=median, geom="point", size=2, color="yellow") + # Plot the median as a point
  scale_x_discrete(limits=c("Individual\nenzymes", "2-enzyme\ncombinations", "All")) + 
  scale_y_continuous(limits=c(0,1)) +
  xlab(NULL) + ylab("Pearson's correlation coefficient") + 
  labs(title = "Trade-off between NF/1000 and Score") + 
  theme_classic() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5));

ggsave('NF_Score_correlation_plot.pdf', width=8, height=8);

#### End of the script ####  
