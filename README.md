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
