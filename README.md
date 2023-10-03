
# Document of Team 7 - part of the Ancient Environmental Genomics Course 2023

# Julia, Chenyu and Mathilde
[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)
![20231003_140557](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/32941862/6165fe93-e43a-4eb3-b0fa-e5fbb7de17e3)
![analysis](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/db593307-ff54-4810-bf97-58451cc6895f)


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

=> MetaDMG is a program specifically designed to analyze damage patterns from metagenomes. It is also computationnally efficient because it combines both the taxonomic classification and the damage analysis by using directly the bam files.
      
  * The taxonomic classification, included in metaDMG, is using a LCA algorythm standing for last common ancestor (here it is ngsLCA).
    * => https://www.geeksforgeeks.org/lowest-common-ancestor-binary-tree-set-1/
      
  * The damage pattern algorythm compares the reads aligned in the sample to the reference and measure the difference between the reference and the sequence. It also includes background noise estimation (ie. sequencing errors). metaDMG provides a damage pattern estimation at the lowest taxonomic resolution possible.


We included in the command line directly the "custom database" parameter with --custom-database and sat 4 cores per sample.

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
### Results observations:
* Exemple of one sample:
![image](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/6b4b8ea4-bab4-4dc3-a1d1-cb8cfbb545c8)


On this plot you can see on the horizontal axis the damage significance and on the vertical axis the damage frequency. Each dot represents one taxa and the size of the dot is correlated with the number of reads of the taxa.


## Statistical plotting in R:
Activate the R environment:
```
conda activate R
```
We generate a header for the new file:
```
echo "ID,Readlength,Count" > long_table.csv
```
Now let's plot the read lenth from the LCA output with ggplot2:

```
for file in *lca.gz
do
    # Extract the ID from the input filename
    id=$(basename "$file" .lca.gz)

    # Sort and count the unique numbers in the file, skip first line, and append to long_table.csv
    zcat "$file" | cut -d':' -f9 | sort | uniq -c | sed 1d | awk -v id="$id" '{print id "," $2 "," $1}' >> long_table.csv

    # Generate a scatterplot of the data using ggplot2
    Rscript -e "library(ggplot2); data <- read.csv('long_table.csv'); myplot <- ggplot(data, aes(x=Readlength, y=Count)) + geom_point() + facet_wrap(~ID, ncol = 2, scales = 'free_y'); ggsave('readlength_distributionTotal.pdf', myplot, width = 6, height = 4)"
done
```
Let's plot at the superkingdom level:
```
echo "ID,Category,Readlength,Count" > long_table2.csv

for file in *lca.gz
do
    # Extract the ID from the input filename
    id=$(basename "$file" .lca.gz)

    # Sort and count the unique numbers in the file, skip first line, and append to long_table.csv
    zcat "$file"  | grep -e 'Archaea' | cut -d':' -f9 | sort -n | uniq -c | sed 1d | awk -v id="$id" '{print id "," "Archaea" "," $2 "," $1 }' >> Archaea_long_table.csv
    zcat "$file"  | grep -e 'Bacteria' -e 'bacteria' | cut -d':' -f9 | sort -n | uniq -c | sed 1d | awk -v id="$id" '{print id "," "Bacteria" "," $2 "," $1 }' >> Bacteria_long_table.csv
    zcat "$file"  | grep -e 'Viridiplantae' | cut -d':' -f9 | sort -n | uniq -c | sed 1d | awk -v id="$id" '{print id "," "Viridiplantae" "," $2 "," $1 }' >> Viridiplantae_long_table.csv
    zcat "$file"  | grep -e 'Virus' | cut -d':' -f9 | sort -n | uniq -c | sed 1d | awk -v id="$id" '{print id "," "Virus" "," $2 "," $1 }' >> Virus_long_table.csv
    zcat "$file"  | grep -e 'Metazoa' | cut -d':' -f9 | sort -n | uniq -c | sed 1d | awk -v id="$id" '{print id "," "Metazoa" "," $2 "," $1 }' >> Metazoa_long_table.csv
    # Create a long table of all long tables
    cat Archaea_long_table.csv Bacteria_long_table.csv Viridiplantae_long_table.csv Virus_long_table.csv Metazoa_long_table.csv >> long_table2.csv

    # Generate a scatterplot of the data using ggplot2
    Rscript -e "library(ggplot2);library(viridisLite); data <- read.csv('long_table2.csv'); plot <- ggplot(data, aes(x = Readlength, y = Count, fill = Category)) + geom_area(alpha = 0.8) + scale_fill_viridis_d() + facet_wrap(~ID, scales = 'free_y', ncol = 2) + labs(x = 'Read Length', y = 'Count', title = 'Read Length Distribution by Category') + theme(plot.title = element_text(hjust = 0.5)); ggsave('readlength_distributionPerKingdom.pdf', plot, width = 5, height = 7)"
done
```
Let's plot at the genus level but only for plants:
```
echo "ID,Category,Readlength,Count" > long_table2.csv

for file in *lca.gz
do
    # Extract the ID from the input filename
    id=$(basename "$file" .lca.gz)

    # Sort and count the unique numbers in the file, skip first line, and append to long_table.csv
    zcat "$file"  | grep -e 'Viridiplantae' | cut -d':' -f9 | sort -n | uniq -c | sed 1d | awk -v id="$id" '{print id "," "Viridiplantae" "," $2 "," $1 }' >> Viridiplantae_long_table.csv
    

    # Generate a scatterplot of the data using ggplot2
    Rscript -e "library(ggplot2);library(viridisLite); data <- read.csv('Viridiplantae_long_table.csv'); plot <- ggplot(data, aes(x = Readlength, y = Count, fill = Category)) + geom_area(alpha = 0.8) + scale_fill_viridis_d() + facet_wrap(~ID, scales = 'free_y', ncol = 2) + labs(x = 'Read Length', y = 'Count', title = 'Read Length Distribution by Category') + theme(plot.title = element_text(hjust = 0.5)); ggsave('readlength_distributionPerGenus_plants.pdf', plot, width = 5, height = 7)"
done
```
