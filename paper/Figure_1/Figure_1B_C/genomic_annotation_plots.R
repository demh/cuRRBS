###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             16/03/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Create several plots which show how different restriction enzymes target different ##
##### genomic features in the genome.                                                  ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

#### Dependencies ####

library(ggplot2);
library(tidyr);

#### Input arguments ####

setwd('~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Figure_1/Figure_1B_C/');

## Read the raw data ##

raw_data <- read.csv('genomic_annotation_for_cleavage_sites.csv');


#### Run the pipeline ####

### 1. Scatterplot of % of sites in CGI vs % of sites in promoters ###

raw_data_1 <- raw_data;
raw_data_1$highlight <- ifelse(raw_data_1$enzyme == 'BsiSI', 'highlight', ifelse(raw_data_1$enzyme == 'Bse118I', 'better', 'normal'));
col_1 <- c("highlight" = "red", "better" = "blue", "normal" = "gray42");
arrow_data_1 <- data.frame(x=12.13, y=10.5, xend=12.13, yend=11.35);
arrow_data_1b <- data.frame(x=15.08, y=12, xend=15.08, yend=12.85);
raw_data_1$total_promoters = raw_data_1$X._sites_CGI_promoters + raw_data_1$X._sites_non_CGI_promoters;

plot_1 <- ggplot(raw_data_1, aes(x=X._sites_CGI, y=total_promoters)) +
  geom_point(aes(colour = highlight, size=n_sites), alpha=0.6) +
  scale_size_continuous(range = c(1, 6), name='Total number\n of sites') +
  scale_color_manual(values = col_1) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold")) +
  guides(colour=FALSE) +
  xlab('% sites in CpG islands') + ylab('% of sites in promoters') +
  annotate("text", x=12.13, y=10.3, label = 'MspI', col='red', size=4) +
  geom_segment(data=arrow_data_1, aes(x=x, y=y, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.2, "cm")), col='red', size=0.5) +
  annotate("text", x=14.9, y=11.8, label = 'BssAI', col='blue', size=4) +
  geom_segment(data=arrow_data_1b, aes(x=x, y=y, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.2, "cm")), col='blue', size=0.5);

ggsave("plot_CGI_vs_promoters.pdf", height=7, width=7);


### 2. Scatterplot of % intergenic sites vs % of sites in non-coding RNA genes ###

raw_data_2 <- raw_data;
raw_data_2$highlight <- ifelse(raw_data_2$enzyme == 'BsiSI', 'highlight', ifelse(raw_data_2$enzyme == 'MfeI', 'better', 
                                                                                 ifelse(raw_data_2$enzyme == 'BsmI', 'better2', 'normal')));
col_2 <- c("highlight" = "red", "better" = "blue", "better2" = "blue", "normal" = "gray42");
arrow_data_2 <- data.frame(x=43.2, y=5.65, xend=41.85, yend=5.5);
arrow_data_2b <- data.frame(x=61.8, y=7.76, xend=60.8, yend=7.92);
arrow_data_2c <- data.frame(x=63.3, y=5.71, xend=64.3, yend=5.87);

plot_2 <- ggplot(raw_data_2, aes(x=X._sites_intergenic, y=X._sites_non_coding_RNA_genes)) +
  geom_point(aes(colour = highlight, size=n_sites), alpha=0.6) +
  scale_size_continuous(range = c(1, 6), name='Total number\n of sites') +
  scale_color_manual(values = col_2) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold")) +
  guides(colour=FALSE) +
  xlab('% of intergenic sites') + ylab('% of sites in non-coding RNA genes') + 
  annotate("text", x=62.7, y=7.73, label = 'MfeI', col='blue', size=4) + 
  geom_segment(data=arrow_data_2b, aes(x=x, y=y, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.2, "cm")), col='blue', size=0.5) +
  annotate("text", x=44, y=5.7, label = 'MspI', col='red', size=4) + 
  geom_segment(data=arrow_data_2, aes(x=x, y=y, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.2, "cm")), col='red', size=0.5) + 
  annotate("text", x=62.2, y=5.68, label = 'BsmI', col='blue', size=4) + 
  geom_segment(data=arrow_data_2c, aes(x=x, y=y, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.2, "cm")), col='blue', size=0.5);


ggsave("plot_intergenic_vs_ncRNA_genes.pdf", height=7, width=7);


### 3. Scatterplots matrix. ###

## Put the data in the right format

all_plots_data <- raw_data[,3:ncol(raw_data)];
colnames(all_plots_data) <- c('Mean GC \ncontent (%)', 'Mean CpG \ncontent (%)', '% of sites in \nprotein-coding \ngenes',
                              '% of sites \nin exons', '% of sites \nin introns', '% of sites \nin non-coding \nRNA genes',
                              '% of intragenic \nsites', '% of intergenic \nsites', '% of sites \nin CGI', 
                              '% of sites \nin shores', '% of sites \nin shelves', '% of sites in \nCGI-containing \npromoters',
                              '% of sites in non \nCGI-containing \npromoters');


## Function to plot the Pearson correlation coefficient
# Adapted from https://www.r-bloggers.com/scatter-plot-matrices-in-r/.

panel.cor <- function(x, y, digits = 3, cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr));
  par(usr = c(0, 1, 0, 1));
  r <- round(cor(x, y), digits=digits);
  text(0.5, 0.6, r, cex=1.5, col=);
}


## Plot

cols <- character(nrow(all_plots_data));
cols[] <- "black";
cols[raw_data$enzyme == 'BsiSI'] <- "red";
pdf('scatterplot_matrix.pdf', height = 10, width = 14);
pairs(all_plots_data, upper.panel=panel.cor, cex=0.4, cex.labels=0.9, col=cols, pch=19);
dev.off();

#### End of the script ####