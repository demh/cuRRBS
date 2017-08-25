###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             04/04/2017                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Lift the CTCF sites coordinates from hg19 to hg38.                               ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

#### Dependencies ####

library(GenomicRanges);
library(rtracklayer);

setwd('~/Documents/Manuscripts/cuRRBS/Figures/Figure_3/3X/CTCF_sites/');


#### Code ####

# Read the coordinates for CTCF sites in hg19.

raw_hg19 <- read.table(file='CTCF_upregulated_and_reactived_sites_hg19.bed',
                       sep='\t');

# Create a GRanges object with this information.

grange_hg19 <- GRanges(seqnames=as.character(raw_hg19[,1]),
                       ranges=IRanges(start=as.numeric(raw_hg19[,2]) + 1,
                                      end=as.numeric(raw_hg19[,3])));

# Lift-over the coordinates.

chain_hg <- import.chain('hg19ToHg38.over.chain');
lifted_hg38 <- liftOver(x=grange_hg19, chain=chain_hg);

# Create the final BED file.

final_bed <- data.frame(chr=as.character(seqnames(unlist(lifted_hg38))),
                        start=start(unlist(lifted_hg38)) - 1,
                        end=end(unlist(lifted_hg38)));

write.table(x=final_bed, file="CTCF_upregulated_and_reactived_sites_hg38.bed", quote=F, sep='\t', 
            row.names = F, col.names = F);

### End of the script ###