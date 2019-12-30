# Bacterial_quasispecies
Finds read which are discordant to the consensus sequence

Blasting gets very complicated when there are multiple regions of homology in the reference sequence. To simplify the process we therefore deplete the reference genome of any high homology regions.
(UNTESTED ON SERVER)

```
Rscript remove_reps_from_fasta.R example.fasta
```


Prepare the data by splitting each read to its own file and making folders etc.

```bash
bash  blast_prep_data.sh UK54.fastq
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


If the blast has been done against a reference genome without any repeat elements in it, the reads which were flagged up must be tested
for IS content at the breakpoint.

Reads therefore must be extracted and blasted against the reference with all repeats still in it.
Note any reads which have repetitive sequences at the breakpoint

```bash
ls -1 all_UK54_final.fasta.rename.fasta.blast_results_final_new_script_final|grep "kek"|sed 's/_blast_table.txtkek//g' >all_reads_depleted_but_not_chopped_to_reblast.txt
mkdir all_reads_depleted_but_not_chopped_to_reblast
mkdir all_reads_depleted_but_not_chopped_to_reblast_results
mkdir all_reads_depleted_but_not_chopped_to_reblast_results_final
cat all_reads_depleted_but_not_chopped_to_reblast.txt|xargs -d '\n' -P 8 -n 1 -I file cp ./all_UK54_final.fasta.rename.fasta.reads/file all_reads_depleted_but_not_chopped_to_reblast/
Rscript ./bin/genomic_blast_lapply_UK48.R all_reads_depleted_but_not_chopped_to_reblast_results all_reads_depleted_but_not_chopped_to_re
blast_results_final/
Rscript ./bin/IS\ there\ an\ IS.R all_reads_depleted_but_not_chopped_to_reblast_results/ all_reads_depleted_but_not_chopped_to_reblast_results_final/  all_reads_depleted_but_not_chopped_to_reblast_results_final/all_reads

```

```bash
Rscript IS\ there\ an\ IS.R SRR5851457.fastq_blast_results/ SRR5851457.fastq_blast_results_final/
```

Generate circos plots
```bash
Rscript temp_circlize.R all_UK54_last_frames UK54_new_processed_blast_results UK54_colour2.jpg kolp.txt
```
