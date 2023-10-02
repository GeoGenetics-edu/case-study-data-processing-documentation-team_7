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
  ### Trimming
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
  ### Choice of parameter

  We use bowtie2 for mapping against the aegenomics.db database.
  
```
for fq in $(ls ./*.fq.gz); do
bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U ${fq/.fq.gz/}.vs.fq.gz --no-unal | samtools view -bS - > ${fq/.fq.gz/}.bam
echo “map ${fq/.fq.gz/}”
done
```
