args = commandArgs(trailingOnly=TRUE)


#cirlize full
library("circlize")

library("stringr")
framey=read.csv(args[1],sep=" ",header=F,stringsAsFactors = F)
kolp=read.delim(args[4],sep=" ",header=F,stringsAsFactors = F)

framey$V1=framey$V1-1
for(i in c(1:nrow(framey)))
{
  lenny=nchar(framey$V1[i])
  
  framey$V1[i]=paste("x",paste(rep("0",c(8-lenny)),collapse=""),c(framey$V1[i]),sep="")
}
framey=framey[framey$V1%in%kolp$V1,]
head(kolp)

#framey=framey[which(framey$V1%in%kolp$V1),]
framey$Category="ZERO"



#Lets categorise breakpoints:
#If its an inversion then prebreakpoint and post breakpoint will be opposite directions:

for(i in c(1:nrow(framey)))
{
  #all_UK54_final.fasta.blast_results/x00000001.fa_blast_table.txt
 lenny=nchar(framey$V1[i]) 
 #cat(args[2],"/","x",rep("0",c(8-lenny)),framey$V1[i],".fa_blast_table.txtkek",sep="") 
 #cat(args[2],"/","x",rep("0",c(8-lenny)),framey$V1[i],".fa_blast_table.txtkek",sep="") 
 
 temp123=read.delim(paste(args[2],"/",framey$V1[i],".fa_blast_table.txtkek",sep=""),sep=" ",header=F)
  
 # temp123=read.delim(paste(args[2],"/",framey$V1[i],".fa_blast_table.txt_kek",sep=""),sep=" ",header=F)
  
  
  
  tempy_kek=vector()
  for(q in c(1:c(nrow(temp123))))
  {
    tempy_kek[q]=temp123$V10[q]-temp123$V9[q]
    
  }
  any_pos=any(tempy_kek>0)
  any_neg=any(tempy_kek<0)
  if(any_pos==T&any_neg==T)
  {
    framey$Category[i]="INV"
  }
}
framey_OG=framey
write.csv(framey,"Reads.csv")


#library(DB)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
#colnames(framey)=c("gene","start","end","read")
framey$V1="Genome"
#Great plot!

jpeg(args[3])
circos.genomicInitialize(framey[,c(1:3)],axis.labels.cex = 0.8)
#TF=c(framey$V2>1260000&framey$V2<1730000)&c(framey$V3>1260000&framey$V3<1750000)
for(i in c(1:nrow(framey)))
{
 
    if(framey$Category[i]=="INV")
    {
      circos.link("Genome", framey$V2[i], "Genome", framey$V3[i],h.ratio = 0.5,col="blue")
      print("kek")
      
    }
   else{
       circos.link("Genome", framey$V2[i], "Genome", framey$V3[i],h.ratio = 0.5,col="red")
       
    
   }
  
  
}
dev.off()

# 
# 
# for(i in c(1:10))
# {
#   
#   if(framey$Category[i]=="INV")
#   {
#     circos.link("Genome", framey$V2[i], "Genome", framey$V3[i],h.ratio = 0.5,col="blue")
#     print("kek")
#     
#   }
#   else{
#     circos.link("Genome", framey$V2[i], "Genome", framey$V3[i],h.ratio = 0.5,col="red")
#     
#     
#   }
#   
#   
# }


