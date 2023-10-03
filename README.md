[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)

# Document of Team 7

## Setting up the environment
We activate the conda environment 'day1', that includes all of the softwares that we need for this analyzis.

```
conda activate day1
```

## Input data

We create virtual links to fasta files in our working directory:
```
for fq in $(ls ./*.fq.gz); do
ln -s ~/course/data/day2/fastq/PRI-TJPGK-CATN-96-98.fq.gz .
done
```


## Pre-processing
  ### Check if trimmed
We looked at the size length of the reads of the different files and they were already trimmed to more than 30 bp.

```
for fq in $(ls ./*.fq.gz); do
echo ${fq/.fq.gz/}
zcat $fq | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'
done
```

### Deduplication

Here we remove duplicated sequences with vsearch.
```
for fq in $(ls ./*.fq.gz); do
echo ${fq/.fq.gz/}
vsearch --fastx_uniques $fq --fastqout ${fq/.fq.gz/}.vs.fq --strand both
echo “rm duplicates of ${fq/.fq.gz/}”
gzip ${fq/.fq.gz/}.vs.fq
done
```

## Mapping
  ### Choice of parameters

  We use bowtie2 for mapping against the aegenomics.db database.
  We used 5 threads to speed the analyzis, with parameter k sat to 100 which means that any reads that have more than  alignments is not reported, this is high enough to keep most of the reads that mapped are kept. (to be discussed together)
 
  
  For practical reasons we run the mapping through tmux to have our job running in the background.

```
tmux new -s name
#come back to your screen
tmux a -t name #link back
#to kill the session
tmux kill-session -t name
#to list number of opened sessions
mux ls #list
#exit the screen
ctrl+B    +    D 
```
  Mapping:
```
for fq in $(ls ./*.vs.fq.gz); do
bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U $fq --no-unal | samtools view -bS - > ${fq/.vs.fq.gz/}.bam
echo “map ${fq/.vs.fq.gz/}”
done
```
   Sorting the bam files by read name:
 ```
for fq in $(ls ./*.vs.fq.gz); do
        samtools sort -n ${fq/.vs.fq.gz/}.bam -@ 5 > ${fq/.vs.fq.gz/}.sorted.bam
done
 ```

## Classification and damage mapping:

Damage patterns are specific to ancient DNA or very degraded such as DNA preserved in formaldehyde for instance. The extremisties could have cytosine deaminations or methylation of the Cytosines in CpG islands resulting in TpG.
MetaDMG is a program specifically designed to analyze damage patterns from metagenomes. It is also computationnally efficient because it combines both the taxonomic classification and the damage analysis by using directly the bam files.
The taxonomic classification is using a LCA algorythm standing for last common ancestor (here it is ngsLCA). => https://www.geeksforgeeks.org/lowest-common-ancestor-binary-tree-set-1/
The damage pattern algorythm compares the reads aligned in the sample to the reference and measure the difference between the reference and the sequence. It also includes background noise estimation (ie. sequencing errors).
metaDMG provides a damage pattern estimation at the lowest taxonomic resolution possible.


We included in the command line directly the "custom database" parameter with --custom-database and sat cores per sample.

```
metaDMG config *.sorted.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp --custom-database --cores-pe--cores-per-sample 4
```
Then we set the differents parameters, for instance the amount of mismatches that we allow, here between 0 and 5 mismatches.
```
vim config.yaml
```

Computing the statistics:
  ```
metaDMG compute config.yaml
  ```

 Visualtion of the damage patterns using GUI:
 ```
metaDMG dashboard config.yaml

 ```
