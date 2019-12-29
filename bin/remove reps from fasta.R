#Remove repetitive sequnces from fasta file
args = commandArgs(trailingOnly=TRUE)

library(Biostrings)
library(stringr)
library(IRanges)
#Genome
kek=readDNAStringSet(args[1])

seq123=seq(1,nchar(kek[[1]]),200)
kolp=DNAStringSet(x = rep("A",c(length(seq123)-10)))
for(k in c(1:c(length(seq123)-10)))
{
  kolp[[k]]=kek[[1]][seq123[k]:c(seq123[k]+1000)]
  
 
}
length(kolp)
names(kolp)=c(1:length(kolp))
writeXStringSet(x = kolp,filepath = paste(args[1],"_1kb_200_step.fa",sep=""))
print("first over with")
#query=readDNAStringSet("UK54_1kb_200_step.fa")
#head(query)
#names(query)=c(1:length(query))
#writeXStringSet(x = query,filepath = "UK54_1kb_200_step.fa")

#Run this command on server to find rep genes:
system(paste("blastn -task megablast -query",paste(args[1],"_1kb_200_step.fa",sep=""), "-db", args[1]," -outfmt 6  -qcov_hsp_perc 50 -out", paste(args[1],"_blast_redund",sep="")))

UK54_reps=read.delim(paste(args[1],"_blast_redund",sep=""),row.names = NULL,stringsAsFactors = F,header=F)
dups=names(table(UK54_reps$V1)[as.numeric(which(table(UK54_reps$V1)>1))])
dups_1=data.frame(Start_rep=UK54_reps$V9[UK54_reps$V1%in%dups],
                  End_rep=UK54_reps$V10[UK54_reps$V1%in%dups])
dups_2=unique(dups_1)
#Put that onto CD hit website and get this file back

seq=rep(1,nchar(kek[[1]]))
for(k in c(1:nrow(dups_2)))
{
  if(k%%1000==0)
  {
    print(k)
  }
  starty=dups_2$Start_rep[k]
    endy=dups_2$End_rep[k]
    seq[starty:endy]=0
    
    
}
kek[[1]]=kek[[1]][which(seq==1)]
writeXStringSet(x = kek,filepath = paste(args[1],"_blast_redund_final_data.fa",sep=""))
data111=data.frame(All=seq)
#write.table(data111,"UK54_depleted_dictionary_genome.txt")
