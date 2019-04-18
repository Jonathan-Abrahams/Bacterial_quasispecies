args = commandArgs(trailingOnly=TRUE)


#cirlize full
library("circlize")
framey=read.csv(args[1],sep=" ",header=F)
kolp=read.delim("kolp.txt",sep=" ")
framey=framey[which(framey$V1%in%kolp$IS),]
framey$V1="Genome"
#library(DB)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
#colnames(framey)=c("gene","start","end","read")
#Great plot!

jpeg(args[2])
circos.genomicInitialize(framey,major.by = 500000,axis.labels.cex = 2)

for(i in c(1:nrow(framey)))
{
  circos.link("Genome", framey$V2[i], "Genome", framey$V3[i],h.ratio = 0.5)
  
}
    dev.off()

