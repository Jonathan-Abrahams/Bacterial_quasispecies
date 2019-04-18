# Bacterial_quasispecies
Finds read which are discordant to the consensus sequence

Download fastq file of long reads

```bash
fastq-dump SRR5851457
```

Prepare the data by splitting each read to its own file and making folders etc.

```bash
bash  blast_prep_data.sh SRR5851457.fastq
```

Create a blastdb

```bash
makeblastdb -in GCF_001987595.1.fasta -parse_seqids -dbtype nucl
```




Blast each read against the consensus

```bash
 cat all_SRR5826165_reads.txt|xargs -d '\n' -P 8 -n 1 -I file blastn -task megablast -query SRR5826165_reads/file -db GCF_001985065.1.fasta -outfmt 6 -out ./SRR5826165_blast_results/file_blast_table.txt
```

Remove blast result files with size 0:

```
find . -size  0 -print0 |xargs -0 rm --
```

### Analyse all blast hits for discordant reads


Run the script
```bash
Rscript genomic_blast_lapply_UK48.R SRR5826165_blast_results SRR5826165_blast_results_final
```
Generate all the hit names

```bash
 find . -iname "*kek"|sed 's/\.\///g'|sed 's/.fa_blast_table.txtkek//g' >all_reads
```

Note any reads which have repetitive sequences at the breakpoint

```bash
Rscript IS\ there\ an\ IS.R SRR5851457.fastq_blast_results/ SRR5851457.fastq_blast_results_final/
```

Generate circos plots
```bash
 Rscript Circlize\ generic.R all_SRR5851457_last_frames all_SRR5851457_last_frames.pdf
```
