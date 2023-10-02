[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)

# Document of Team 7

## Setting up the environment

'''
conda activate day1
'''

## Input data

Create virtual links to fasta files:
virtual_links.sh:
'''
ln -s ~/course/data/day2/fastq/PRI-TJPGK-CATN-96-98.fq.gz .
'''
Launch the script:
'''
'''

## Pre-processing
  ### Trimming
  We trim sequences under 30 base pairs.
  ### Deduplication
  Then we remove duplicate sequences.
Both steps are done with vsearch:

deduplicate.sh
'''
vsearch --fastx_uniques PRI-TJPGK-CATN-96-98.fq.gz --fastqout PRI-TJPGK-CATN-96-98.vs.fq --minseqlength 30 --strand both
'''

## Mapping
  ### Choice of parameter

  We use bowtie2 for mapping against the aegenomics.db database.
  
