#fastq-dump $1

#bash  ./bin/blast_prep_data.sh $1

#mkdir .$1.blast_results

#makeblastdb -in $2 -parse_seqids -dbtype nucl

ls -1 $1.reads >all_$1.reads.txt
cat all_$1.reads.txt|xargs -d '\n' -P 8 -n 1 -I file blastn -task megablast -query $1.reads/file -db $2 -outfmt 6 -out ./$1.blast_results/file_blast_table.txt
cd $1.blast_results

find . -size  0 -print0 |xargs -0 rm --

cd ..
mkdir $1.blast_results_final

Rscript ./bin/genomic_blast_lapply_UK48.R $1.blast_results $1.blast_results_final
#For now we change this as we are already in bin. CHANGE THIS BACK FOR NEXT RUN

#Rscript genomic_blast_lapply_UK48.R $1.blast_results $1.blast_results_final

cd $1.blast_results_final
find . -iname "*kek"|sed 's/\.\///g'|sed 's/.fa_blast_table.txtkek//g' >all_$1_reads

cd ..

#For now we change this as we are already in bin. CHANGE THIS BACK FOR NEXT RUN

Rscript ./bin/IS\ there\ an\ IS.R $1.blast_results/ $1.blast_results_final/ ./$1.blast_results_final/all_$1_reads

#Rscript .//IS\ there\ an\ IS.R $1.blast_results/ $1.blast_results_final/ ./$1.blast_results_final/all_$1_reads


cd $1.blast_results_final


cat *last_frame >all_$1_last_frames

cd .. 

#Rscript temp_circlize.R all_UK54_last_frames UK54_new_processed_blast_results UK54_colour2.jpg kolp.txt

Rscript Circlize_basketball.R ./$1.blast_results_final/all_$1_last_frames $1.blast_results_final/ $1_colour.jpg kolp.txt
