#compare the filtered blast results with the real blast resutls and answer the question:
library(IRanges)
args = commandArgs(trailingOnly=TRUE)

#"Is there an IS element at the junction?"
#gff=read.delim("sequence (4).gff3",header=F,stringsAsFactors = F)
#args1 is the full blast results folder
#args2 is the hits blast results_folder
is_IS=function(read.name)
{
  Full_blast=read.delim(paste(args[1],read.name,".fa_blast_table.txt",sep=""),header=F)
  Hit_blast=read.delim(paste(args[2],read.name,".fa_blast_table.txtkek",sep=""),header=F,sep=" ")
  boomlies=read.delim(paste(args[2],read.name,".fa_blast_table.txtlast_frame",sep=""),header=F,sep="")
  
  pre_break=Hit_blast$V8[which(Hit_blast$V10%in%boomlies$V2)]
  post_break=Hit_blast$V7[which(Hit_blast$V10%in%boomlies$V2)+1]
  tabley1=table(Full_blast$V7)
  tabley2=table(Full_blast$V8)
  reppy_seq=tabley1[which(as.numeric(tabley1)>=30)]
  
  which(Hit_blast$V10%in%boomlies$V2)
  
  kek1=IRanges(pre_break-1500,pre_break)
  kek2=IRanges(post_break,post_break+1500)
  Read_IS=IRanges(as.numeric(names(reppy_seq)),as.numeric(names(reppy_seq)))
  
  o1=findOverlaps(Read_IS,kek1)
  o2=findOverlaps(Read_IS,kek2)
  
  if(length(o1)>=1|length(o2)>=1)
  {
    return(read.name)
  }
}
ready=read.delim(args[3],header=F)
results_vector=vector()
for(i in c(1:nrow(ready)))
{
  print(paste("Args[3]: ",i))
  kol=is_IS(ready$V1[i])
  print(kol)
  if(length(kol>=1))
  {
    results_vector[i]=kol
  }
  
  
}
print("All_done")
kolp=data.frame(IS=results_vector)
kolp=kolp[!is.na(kolp$IS),]
write.table(kolp,"kolp.txt",col.names = F,row.names = F)
