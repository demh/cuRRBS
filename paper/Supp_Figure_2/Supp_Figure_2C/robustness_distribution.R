###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             27/05/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Show the distribution that the values for robustness take in different systems.  ####
##### This can help to differentiate between what can be considered a robust protocol  ####
##### and what not.                                                                    ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

suppressWarnings(suppressMessages(library(ggplot2)));
library(latex2exp);

setwd('~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Supp_Figure_2/Supp_Figure_2C/');


################################################################
################## Running the pipeline ########################
################################################################

#### 1. Read all the cuRRBS output files and extract the robustness values. ####

R_values <- c();

for(file in list.files('cuRRBS_outputs/')){
  
  raw_file <- read.csv(paste0('cuRRBS_outputs/', file));
  R_values <- c(R_values, raw_file$Robustness);
  
};


#### 2. Create a plot with the robustness distribution ####

# Put the data in the appropiate format.

R_df <- as.data.frame(x=R_values);

# Calculate additional summary statistics.

q25 <- quantile(R_values,.25);
q75 <- quantile(R_values,.75);
med_R <- median(R_values);
R.dens <- density(R_values);
df.dens <- data.frame(x = R.dens$x, y = R.dens$y);

# Make the plot

distribution_plot <- ggplot(data=R_df) + 
  geom_density(aes(x=R_values, y=..density..), color='black', size=1) + 
  geom_area(data = subset(df.dens, x >= q25 & x <= q75), 
            aes(x=x,y=y,fill = 'orange')) +
  geom_area(data = subset(df.dens, x < q25), 
            aes(x=x,y=y, fill = 'red')) +
  geom_area(data = subset(df.dens, x > q75), 
            aes(x=x,y=y, fill = 'forestgreen')) +
  geom_segment(x = med_R, xend=med_R, y=-1, yend=max(df.dens)+1, col='blue', linetype='dashed') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.5, 0.5),
        legend.margin = margin(),
        legend.box.background = element_rect(colour = "black", size=0.6),
        legend.box.margin = margin(1, 3, 3, 3, "mm")) +
  scale_fill_manual(values=c('forestgreen', 'orange', 'red'), 
                    labels=c("R > Q3 = 0.9834",
                             "Q1 <= R <= Q3",
                             "R < Q1 = 0.9580"),
                    name='') +
  xlab('Robustness (R)') + ylab('Density');

ggsave(paste0('robustness_distribution.pdf'), width = 7, height=7);


### End of the script ###


