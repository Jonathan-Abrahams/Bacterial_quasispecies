#Remove repetitive sequnces from fasta file
args = commandArgs(trailingOnly=TRUE)

library(Biostrings)
library(stringr)
library(IRanges)
kek=readDNAStringSet("Bp_UK54_R941.flipflop.porechop.ctg.lay.fa")
kek=readDNAStringSet(args[1])

seq123=seq(1,nchar(kek[[1]]),200)
kolp=DNAStringSet(x = rep("A",length(seq123)-1))
for(k in c(17592:length(seq123)))
{
  kolp[[k]]=kek[[1]][seq123[k]:c(seq123[k]+1000)]
  
  
}
writeXStringSet(x = kolp[1:20406],filepath = "UK54_1kb_200_step.fa")


#Run this command on server to find rep genes:
system("blastn -task megablast -query UK54_1kb_200_step_rename.fa -db Bp_UK54_R941.flipflop.porechop.ctg.lay.fa -outfmt 6  -out UK54_rep_blast_1kb_200_step_no_redunc")

UK54_reps=read.delim("UK54_rep_blast_1kb_200_step_no_redunc",row.names = NULL,stringsAsFactors = F,header=F)
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
writeXStringSet(x = kek,filepath = "UK54_depleted7.fa")
data111=data.frame(All=seq)
write.table(data111,"UK54_depleted_dictionary_genome.txt")
