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
##### Remove the 'G' coordinates from the sites annotation file in the ESC and IAP     ####
##### 'more' mouse systems.                                                            ####
###########################################################################################
##### USAGE: Rscript remove_Gs.R sites_annotation.csv                                  ####
###########################################################################################


args <- commandArgs(trailingOnly=TRUE);

## Read the file ##

raw_file <- read.csv(args[1]);

## Check that all the CpG sites are reported as 'CG' pairs.

all_putative_Cs_index <- seq(1,nrow(raw_file),2);

for(i in all_putative_Cs_index){
  
  coord_C <- as.numeric(raw_file[i,3]);
  coord_G <- as.numeric(raw_file[i+1,3]);
  
  if(coord_G-coord_C != 1){
    
    print('ERROR: CpG site incorrectly reported !');
    
  }
}


## Remove the 'G's ##

no_Gs <- raw_file[all_putative_Cs_index,];
no_Gs[,1] <- 1:nrow(no_Gs); # New indexing

## Export the new sites annotation file ##

output_path <- paste0(strsplit(args[1], '.csv')[[1]][1], '_noGs.csv');
write.csv(no_Gs, file=output_path, quote=F, row.names=F);

### End of the script ###