###########################################################################################
#########                                                                         #########
#########                           Antonio J. Ribeiro                            #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             31/05/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Study how different parameters affect the computational time for cuRRBS runs.    ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

#### Dependencies ####

library(ggplot2);

setwd('~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Supp_Figure_5/Supp_Figure_5A_D/');


#### 1. Read the raw data ####

time_n_enzymes <- data.frame(n_enzymes = c(5,10,20,30,40,50,60,70,80,90,100,110,120), 
                             time1 = c(3.5,7.4,39.2,41.6,123,152.4,240.6,317.9,452.8,545,672,749.4,862.8),
                             time2 = c(4.8,12,28,76,137,163,226,376,500,601,658,847,917),
                             time3 = c(1.8,10,26,51,114,143,252,303,478,570,667.7,758,986));

time_exp_error <- data.frame(exp_error = c(5,10,15,20,25,30,35,40,45,50),
                             time1 = c(171,140,112,105,101,103,105,102,103,101),
                             time2 = c(165,123,120,112,116,104,109,116,107,112),
                             time3 = c(171,114,121,114,109,121,106,113,112,117));

time_n_sites <- data.frame(n_sites = c(100000,50000,10000,5000,1000,500,100,50,10,5),
                           time1 = c(331,222,143,117,110,115,114,112,107,109),
                           time2 = c(315,208,127,116,109,109,108,108,106,107),
                           time3 = c(313,209,126,117,110,108,108,107,106,109));

time_genome_size <- data.frame(genome_size = c(4.3,0.15,3.8),
                               time1 = c(141,50,125),
                               time2 = c(145,48,121),
                               time3 = c(133,47,139));


#### 2. Calculate the mean time and the standard deviation between the 3 runs ####

time_n_enzymes$mean_time <- apply(time_n_enzymes[,2:4], 1, mean);
time_n_enzymes$sd <- apply(time_n_enzymes[,2:4], 1, sd);

time_exp_error$mean_time <- apply(time_exp_error[,2:4], 1, mean);
time_exp_error$sd <- apply(time_exp_error[,2:4], 1, sd);

time_n_sites$mean_time <- apply(time_n_sites[,2:4], 1, mean);
time_n_sites$sd <- apply(time_n_sites[,2:4], 1, sd);

time_genome_size$mean_time <- apply(time_genome_size[,2:4], 1, mean);
time_genome_size$sd <- apply(time_genome_size[,2:4], 1, sd);

#### 3. Make the plots ####

p_n_enzymes <- ggplot(time_n_enzymes, aes(x=n_enzymes, y=mean_time)) + 
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=mean_time-sd, ymax=mean_time+sd), width=.1, col='red') +
  theme_classic() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold")) + 
  xlab('Number of enzymes') + ylab('Mean time (s)');

ggsave('plot_time_vs_number_ezymes.pdf', width=8, height=8);

p_exp_error <- ggplot(time_exp_error, aes(x=exp_error, y=mean_time)) + 
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=mean_time-sd, ymax=mean_time+sd), width=.1, col='red') +
  theme_classic() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold")) + 
  xlab('Experimental error (bp)') + ylab('Mean time (s)');

ggsave('plot_time_vs_exp_error.pdf', width=8, height=8);

p_n_sites <- ggplot(time_n_sites, aes(x=n_sites, y=mean_time)) + 
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=mean_time-sd, ymax=mean_time+sd), width=.1, col='red') +
  theme_classic() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold")) + 
  xlab('Number of sites of interest') + ylab('Mean time (s)');

ggsave('plot_time_vs_number_sites.pdf', width=8, height=8);

p_genome_size <- ggplot(time_genome_size, aes(x=genome_size, y=mean_time)) + 
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=mean_time-sd, ymax=mean_time+sd), width=.1, col='red') +
  theme_classic() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold")) + 
  xlab('Genome size (GB of pre-computed files)') + ylab('Mean time (s)');

ggsave('plot_time_vs_genome_size.pdf', width=8, height=8);


### End of the script ###
