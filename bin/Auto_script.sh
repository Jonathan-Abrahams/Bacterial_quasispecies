fastq-dump $1

bash  blast_prep_data.sh $1.fastq

mkdir $1.fastq_blast_results

makeblastdb -in $2 -parse_seqids -dbtype nucl

ls -1 $1.fastq_reads >all_$1.fastq_reads.txt
cat all_$1.fastq_reads.txt|xargs -d '\n' -P 8 -n 1 -I file blastn -task megablast -query $1.fastq_reads/file -db $2 -outfmt 6 -out ./$1.fastq_blast_results/file_blast_table.txt
cd $1.fastq_blast_results

find . -size  0 -print0 |xargs -0 rm --

cd ..
mkdir $1.fastq_blast_results_final

Rscript ./bin/genomic_blast_lapply_UK48.R $1.fastq_blast_results $1.fastq_blast_results_final
cd $1.fastq_blast_results_final
find . -iname "*kek"|sed 's/\.\///g'|sed 's/.fa_blast_table.txtkek//g' >all_$1_reads

cd ..

Rscript ./bin/IS\ there\ an\ IS.R $1.fastq_blast_results/ $1.fastq_blast_results_final/ ./$1.fastq_blast_results_final/all_$1_reads
cd $1.fastq_blast_results_final


cat *last_frame >all_$1_last_frames

cd .. 

Rscript Circlize\ generic.R ./$1.fastq_blast_results_final/all_$1_last_frames all_$1_last_frames.jpg

