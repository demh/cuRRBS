###########################################################################################
#########                                                                         #########
#########                           Thomas M. Stubbs                              #########
#########                             03/03/2017                                  #########
#########                         Babraham Institute                              #########
#########                             Reik group                                  #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####                                 cuRRBS paper                                     ####
###########################################################################################
##### Study the phylogenetics and the variability in the different restriction enzyme   ###
##### recognition motifs.                                                              ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

################## analysis of enzyme motif distances ################ 
library(ape)
library(msa)
library(seqinr)
library(ggtree)  
library(RColorBrewer)

setwd("~/Desktop/methylation_clock/optimize_RRBS/cuRRBS_paper/Supp_Figure_1/Supp_Figure_1A_B/")
test<-read.table("recognition_sites_enzymes.csv", header = F, sep = ",")

###### distance calculation #####
DNA<-as.vector(test[,2])
DNAset<-DNAStringSet(DNA,use.names = T)
names(DNAset)<-paste0(test[,2])
output<-msaMuscle(DNAset, type = "dna", gapOpening = 400, gapExtension = 0)
output
output2 <- msaConvert(output, type="seqinr::alignment")

distances <- dist.alignment(output2, "identity")
dist2<-as.matrix(dist.alignment(output2, matrix = "identity"))

##### phylogeny ####
phylo <- njs(distances)

names<-strsplit(phylo$tip.label,split = "")

numberAT<-sapply(names, GC, exact=TRUE)
numberAT[which(numberAT>=0 & numberAT<0.25)] <- brewer.pal(n=5, name = "YlGn")[2]
numberAT[which(numberAT>=0.25 & numberAT<0.5)] <- brewer.pal(n=5, name = "YlGn")[3]
numberAT[which(numberAT>=0.5 & numberAT<0.75)] <- brewer.pal(n=5, name = "YlGn")[4]
numberAT[which(numberAT>=0.75 & numberAT<=1)] <- brewer.pal(n=5, name = "YlGn")[5]

tipColour<-numberAT

####### plotting phylogeny #######
pdf(file = "Phylogeny_restriction_motifs.pdf",width = 20,height = 20,paper = "special")
plot.phylo(phylo, type="r",tip.color =tipColour,
           cex=2,font=2,edge.width=c(4),no.margin=T)
legend("topleft",legend = c("0-25% GC content", "25-50% GC content","50-75% GC content","75-100% GC content"), fill = c(
  brewer.pal(n=5, name = "YlGn")[2],brewer.pal(n=5, name = "YlGn")[3],brewer.pal(n=5, name = "YlGn")[4],brewer.pal(n=5, name = "YlGn")[5])
)
dev.off()


####### euclidean distance converting NAs to a maximal distance of 1 ########
for (i in 1:length(dist2[,1])){
  as.numeric(gsub(NaN,1,dist2[i,]))->dist2[i,]
}

###### defining the PCA #####
pca<-prcomp(dist2)
data_to_plot<-cbind(pca$x[,1],pca$x[,2])
dim(data_to_plot)
#plot(pca$x[,1], pca$x[,2])

####### plot info ######
test_motifs<-strsplit(rownames(data_to_plot),"")

width_to_plot<-sapply(test_motifs,function(x){length(x)})

numberAT<-sapply(test_motifs, GC, exact=TRUE)
numberAT[which(numberAT>=0 & numberAT<0.25)] <- brewer.pal(n=5, name = "YlGn")[2]
numberAT[which(numberAT>=0.25 & numberAT<0.5)] <- brewer.pal(n=5, name = "YlGn")[3]
numberAT[which(numberAT>=0.5 & numberAT<0.75)] <- brewer.pal(n=5, name = "YlGn")[4]
numberAT[which(numberAT>=0.75 & numberAT<=1)] <- brewer.pal(n=5, name = "YlGn")[5]

###### plot function ######
pdf(file = "PCA_restriction_motifs.pdf",width = 10,height = 10,paper = "special")
symbols(x=data_to_plot, circles=width_to_plot, xlim=c(-5,3),ylim = c(-3,3),
        inches=1/5,ann=F, bg=numberAT, fg=NULL)
title(xlab="First component", ylab="Second component")
legend("topleft",legend = c("0-25% GC content", "25-50% GC content","50-75% GC content","75-100% GC content"), fill = c(
  brewer.pal(n=5, name = "YlGn")[2],brewer.pal(n=5, name = "YlGn")[3],brewer.pal(n=5, name = "YlGn")[4],brewer.pal(n=5, name = "YlGn")[5])
)
dev.off()
