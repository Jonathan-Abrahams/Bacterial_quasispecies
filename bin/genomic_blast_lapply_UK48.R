
args = commandArgs(trailingOnly=TRUE)

#genomic blast lapply
library(tictoc)
#debug
library(BiocGenerics)
#Can we detect rearangements that conflict with the prevailing rearangment?

#In technical terms, can we detect orders of genes which are discordant with the order of genes that is normal?




#gff=read.delim("sequence (4).gff3",header=F,stringsAsFactors = F)

#trans_gff=gff[which(gff$V10==1),]

library(IRanges)
#P_blast_table=kek$Frame
find_max_overlappies=function(P_blast_table)
{
  results_vector=vector()
  for( i in c(1:nrow(P_blast_table)))
  {
    
    Rangey=IRanges(P_blast_table$`Start_of_alignment(Q.)`[i],P_blast_table$`End_of_alignment(Q.)`[i])
    Rest_rows=c(1:nrow(P_blast_table))[!c(1:nrow(P_blast_table))%in%i]
    Rest_range=IRanges(P_blast_table$`Start_of_alignment(Q.)`[Rest_rows],P_blast_table$`End_of_alignment(Q.)`[Rest_rows])
    overlaps123=findOverlaps(Rangey,Rest_range)
    hit_rows=Rest_rows[subjectHits(overlaps123)]
    maxxy=max(P_blast_table$Bit_score[c(hit_rows,i)])
    top_hit=c(hit_rows,i)[which(P_blast_table$Bit_score[c(hit_rows,i)]%in%maxxy)]
    results_vector[i]=top_hit[1]
  }
  testy_final=P_blast_table[unique(results_vector),]
  
  return(testy_final)
}
#Trans_range=IRanges(start=trans_gff$V4,end=trans_gff$V5)

files=list.files(args[1])
#files=vector()
#files[1]="100008.fa_blast_table.txt_"
#blast_table=blast_results
Test_rangement=function(blast_table)
{
  column_names=c("Query_ID","Subject_ID","Percentage_base_match","Alignment_Length",
                 "Number_of_mismatches","Gaps","Start_of_alignment(Q.)","End_of_alignment(Q.)",
                 "Start_of_alignment(S.)","End_of_alignment(S.)","E_value","Bit_score")
  colnames(blast_table)=column_names
  
  
  
  blast_table=blast_table[which(blast_table$`Alignment_Length`>=2500),]
  if(nrow(blast_table)>=1)
  {
    #print(percentOverlap)
    testy131=find_max_overlappies(blast_table)
    #So we have made it! now to do the actual test
    
    genome_gaps=abs(testy131$`Start_of_alignment(S.)`[-1]-testy131$`End_of_alignment(S.)`[-length(testy131$`Start_of_alignment(S.)`)])
    read_gaps=abs(testy131$`Start_of_alignment(Q.)`[-1]-testy131$`End_of_alignment(Q.)`[-length(testy131$`Start_of_alignment(Q.)`)])
    
    tingy=which(abs(genome_gaps-read_gaps)>=30000)
    if(any(which(abs(genome_gaps-read_gaps)>=30000)))
    {
      testy131=find_max_overlappies(testy131)
      
      results_frame=list(Frame=testy131,tingy=tingy)
      return(results_frame)
    }
    
    #however this messes up the order. The order of this variable is ntohing to do with the order of the blast table.
    
  }
 
}

#print(files[4000:5000])

#x=files[1]
funky=function(x)
{
  
  #print(x)
  blast_results=read.delim(paste(args[1],"/",x,sep=""),header=F)
  column_names=c("Query_ID","Subject_ID","Percentage_base_match","Alignment_Length",
                 "Number_of_mismatches","Gaps","Start_of_alignment(Q.)","End_of_alignment(Q.)",
                 "Start_of_alignment(S.)","End_of_alignment(S.)","E_value","Bit_score")
  colnames(blast_results)=column_names
  blast_results$UID=c(1:nrow(blast_results))
  #blast_table=blast_results
  
  
  kek=Test_rangement(blast_results)
  
  if(length(kek)>=2&is.null(nrow(kek$Frame))==F)
  {
    kek$Frame=kek$Frame[order(kek$Frame[,7]),]
    
    if(length(unique(kek$Frame[,12]))==nrow(kek$Frame))
    {
      genome_gaps=abs(kek$Frame$`Start_of_alignment(S.)`[-1]-kek$Frame$`End_of_alignment(S.)`[-length(kek$Frame$`Start_of_alignment(S.)`)])
      read_gaps=abs(kek$Frame$`Start_of_alignment(Q.)`[-1]-kek$Frame$`End_of_alignment(Q.)`[-length(kek$Frame$`Start_of_alignment(Q.)`)])
      
      tingy=which(abs(genome_gaps-read_gaps)>=30000)
      more_than=which(abs(genome_gaps-read_gaps)>=30000)
      less_than=which(abs(genome_gaps-read_gaps)<=3400000)
      if(length(which(more_than==less_than))>=1)
      {
        #kek$Frame=find_max_overlappies(kek$Frame)
        
        # results_frame=list(Frame=kek$Frame,tingy=tingy)
        print(x)
        write.table(kek$Frame,paste(args[2],"/",x,"kek",sep=""),row.names = F,col.names = F)
        
        #Need to write the before,after and read name
        before=kek$Frame[which(more_than==less_than),10]
        after=kek$Frame[which(more_than==less_than)+1,9]
        read=kek$Frame[1,1]
        last_frame=data.frame(read,before,after)
        write.table(last_frame,paste(args[2],"/",x,"last_frame",sep=""),row.names = F,col.names = F)
        
        return(kek)
      }
  
    #}
  }
  

  }
}
library(lme4)
tic()
lol=mclapply(files, funky,mc.cores = 8)
#lol=lapply(files, funky)


#lol=lapply(files, funky)

toc()
