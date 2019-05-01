#Script that takes a supplied fastq file and blasts it against a supplied consensus sequnece, one read at a time.

#$1 is a fastq

mkdir $1_reads

#turn to fasta
sed -n '1~4s/^@/>/p;2~4p' $1 >$1.fasta

#rename the reads
awk '/^>/{print ">" ++i; next}{print}' < $1.fasta >$1.rename.fasta

#split into individual reads
cp $1.rename.fasta ./$1_reads
cd $1_reads
cat $1.rename.fasta|awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}print $0 > filename}'

rm $1.rename.fasta

cd ..

mkdir $1_blast_results





