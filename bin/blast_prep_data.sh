#Script that takes a supplied fastq file and blasts it against a supplied consensus sequnece, one read at a time.

#$1 is a fastq or fasta 

mkdir $1.reads
if [[ ${file: -5} == ".fastq" ]]
then
        echo Fastq file detected. Turning to fasta.
	sed -n '1~4s/^@/>/p;2~4p' $1 >$1.fasta

else
        echo Fasta file detected! adding another .fasta to the end to ease the process along!
	cat $1 >$1.fasta
fi

#turn to fasta
#sed -n '1~4s/^@/>/p;2~4p' $1 >$1.fasta

#rename the reads
awk '/^>/{print ">" ++i; next}{print}' < $1.fasta >$1.rename.fasta

#split into individual reads
cp $1.rename.fasta ./$1.reads
cd $1.reads
split --verbose -a 8 -d --additional-suffix=.fa -l 2 $1.rename.fasta

rm $1.rename.fasta

cd ..

mkdir $1.blast_results
