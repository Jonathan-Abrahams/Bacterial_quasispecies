#Remove repetitive sequnces from fasta file
args = commandArgs(trailingOnly=TRUE)

library(Biostrings)
library(stringr)
library(IRanges)
kek=readDNAStringSet("Bp_UK54_R941.flipflop.porechop.ctg.lay.fa")
kek=readDNAStringSet(args[1])

seq123=seq(1,nchar(kek[[1]]),200)
kolp=DNAStringSet(x = rep("A",length(seq123)-1))

qoa=lapply(seq123,function(x)
  {
  # print(x)
  if(x<c(length(kek[[1]])-2000))
  {
    return(kek[[1]][x:c(x+1000)])
    
  }
})
listy123=DNAStringSet(unlist(qoa,recursive = T))
names(listy123)=c(1:length(listy123))
# for(k in c(1:c(length(seq123)-1)))
# {
#   kolp[[k]]=kek[[1]][seq123[k]:c(seq123[k]+1000)]
#   
#   
# }
 writeXStringSet(x = listy123,filepath = paste(args[1],"_non_redunc.fasta",sep=""))


#Run this command on server to find rep genes:
system(paste("makeblastdb -in", args[1], "-parse_seqids -dbtype nucl"))
system(paste("blastn -task megablast -query", paste(args[1],"_non_redunc.fasta",sep=""), "-db",args[1], "-outfmt 6  -out",paste(args[1],"_non_redunc.fasta_blast",sep="")))

UK54_reps=read.delim(paste(args[1],"_non_redunc.fasta_blast",sep=""),row.names = NULL,stringsAsFactors = F,header=F)
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
writeXStringSet(x = kek,filepath = paste(args[1],"depleted.fasta",sep="_"))
data111=data.frame(All=seq)
write.table(data111,"UK54_depleted_dictionary_genome.txt")
