
# Document of Team 7 - part of the Ancient Environmental Genomics Course 2023

# Julia, Chenyu and Mathilde
[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)
![20231003_140557](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/32941862/6165fe93-e43a-4eb3-b0fa-e5fbb7de17e3)
![analysis](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/db593307-ff54-4810-bf97-58451cc6895f)

Note: code verified a second time until mapping so far, ie corrections of typos etc.

## Setting up the environment
We activate the conda environment 'day1', that includes all of the softwares that we need for this analyzis.

```
conda activate day1
```

## Input data

We create virtual links to fasta files in our working directory:
```
for fq in $(ls ~/course/data/day2/fastq/*.fq.gz); do
ln -s ${fq} .
ls -l
done
```


## Pre-processing
  ### Check if trimmed
We looked at the size length of the reads of the different files and they were already trimmed to more than 30 bp.

```
for fq in $(ls ./*.fq.gz); do
echo ${fq}
zcat ${fq} | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'
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

```
conda activate day1
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
metaDMG config *.sorted.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp --custom-database --cores-per-sample 4
```

Then we set the differents parameters in the [configuration file](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/blob/main/config.yaml), for instance the amount of mismatches that we allow, here between 0 and 5 mismatches.
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

* Overview of the damage results for all the samples, with a filter of damage significance > 2.3 and damage frequency > 3%.
![image_720](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/cda7937d-06b3-4ff0-a71e-5c8e9bb2efd0)

On this plot you can see on the horizontal axis the damage significance and on the vertical axis the damage frequency. Each dot represents one taxa and the size of the dot is correlated with the number of reads of the taxa.
A large part of the taxa sit below 10% estimate.

* On this plot, we aggregated the taxa at the Kingdom level. We can observe that this results in the display of the mean of damage frequency aggregrated. As a results damage frequency appear lower overall. 
![newplot](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/9cdd7d0f-d929-438b-b05a-3443ea7bb81d)

* Here we filtered at 10% of damage frequency, and the same threshold for damage significance:
![newplot__2__720](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/90b4c612-3669-47ff-b4d4-2e12684367da)

Having a more drastic threshold for damage frequency can help to avoid signal from living organisms from the sediment, ie bacteria.
=> However, we will discuss this more in depth in the coming days as those organisms while being not necessarily dead could be also informative in past environment/climate reconstruction (ie. organisms from [Antonio's preprint][Reference_preprint_Antonio], [spores in the marine sediment][ref_notdead_yet]).


## Statistical plotting in R: /!\ to be done
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


## DNA Authentication
We have filtered the damage results further to find the threshold where the results are most reliable. We are filtering for the number of reads, the length of the reads and the significance of the damage.


Filtering for DNA fragments with more than 35 base pairs and a damage significance of 5.
```
filtered_data_viridiplantae <- filtered_data %>% filter(N_reads >= 100, mean_L > MinLength3, MAP_significance  > MapSig3,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", YearsBP))

> MinLength3
[1] 35
> MapSig3
[1] 5

```
What species can be found?
```
> select(filtered_data_viridiplantae, tax_name, MAP_damage, MAP_significance, N_reads, YearsBP)
       tax_name MAP_damage MAP_significance N_reads YearsBP
1 g__Potentilla 0.009030714        5.117909  21725   1200
2 g__Epilobium 0.010413541        6.293386  33841   1200
> unique(filtered_data_viridiplantae$tax_name)
[1] "g__Potentilla" "g__Epilobium"
```
When you filter for a significance of 5 it is too stringent because it only results in two species. 

[aeCourse.DNAdamageLRJitterPlot.pdf](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/files/12801320/aeCourse.DNAdamageLRJitterPlot.pdf)

[aeCourse.DNAdamageModelJitterPlot.pdf](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/files/12801321/aeCourse.DNAdamageModelJitterPlot.pdf)


## EUKA: mapping multiple samples to the pangenome graph (needs descriptions)

creat an output directory
```
mkdir /home/user-5/course/wdir/day1_group/euka
```
input a single sample (which works)
```
~/course/data/vgan/bin/vgan euka -fq1 <(zcat ~/course/wdir/mapping/PRI-TJPGK-CATN-96-98.fq.gz.vs.fq.gz) -o PRI-TJPGK-CATN-96-98 -t 5
```
input multiple samples (which did not work)
```
~/course/data/vgan/bin/vgan euka -fq1 <(zcat ~/course/wdir/day1_group/*.vs.fq.gz) -o all_sample -t 5 --euka_dir /home/user-5/course/wdir/day1_group/euka
```
visualization:

?where is the [tree.nw](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/blob/main/tree.nw)

drag the tree.csv file to https://itol.embl.de/

## Pathphynder: Phylogenetic placement (needs description)
 Create VCF with snp-sites from multiple sequence alignment (MSA) file

snp-sites -v -c -o Aln_mafft_taxa_references.vcf All_ref_only_realigned.fa

Make sure package stringr is installed
R
library(string)

if not install it:
install.packages("stringr")
quit()

fix vcf file
Rscript /path_to_Rscript/fix_vcf.R Aln_mafft_taxa_references.vcf Aln_mafft_taxa_references_output.vcf
fix consensus naming problem in vcf file (replace 1’s with “consensus”)

awk '{ if ($1 == "1") $1="consensus";}1' Aln_mafft_taxa_references_output.vcf | sed 's/ /\t/g' > Aln_mafft_taxa_references_fixed_output.vcf

Install biopython

pip install biopython

Create consensus for MSA.fa and index it

python /path_To_Script/get_consensus.txt All_ref_only_realigned.fa Cons_Aln_mafft_taxa_references.fa
bwa index Cons_Aln_mafft_taxa_references.fa
Create directory in pathPhynder_analysis folder

mkdir map_to_cons
bwa mapping of Ovis reads to consensus: conda activate day1

bwa aln -l 1024 -n 0.001 -t 5 /path_to_folder/Consensus_Mafft_All_references.fa /path_to_folder/library.fq | bwa samse /path_to_folder/Consensus_Mafft_All_references.fa  - /path_to_folder/library.fq | samtools view -F 4 -q 25 -@ 5 -uS - | samtools sort -@ 5 -o library.sort.bam
Count number of mapped reads

samtools view -c file.bam
8. PathPhynder
Install pathPhynder

git clone https://github.com/ruidlpm/pathPhynder.git
make a new .bash_rc

touch ~/.bash_rc
vi ~/.bash_rc
Add manually these lines, using the path where pahPhynder is installed:


alias pathPhynder="Rscript /path_to_pathPhynder_folder/pathPhynder.R"
save file

Esc 
:wq
source it
```
source ~/.bash_profile
```
test it
pathPhynder -h

Create directory in pathPhynder_analysis folder

mkdir pathphynder_results
Assign informative SNPs to tree branches: use 
```
conda activate day2

phynder -B -o ./branches.snp ./Mafft_All_BEAST4.nwk ./Mafft_All_references_fixed_consensus.vcf

conda activate r
```
Install samtools
```
mamba install -c bioconda samtools
```
Prepare data - this will output a bed file for calling variants and tables for phylogenetic placement
```
pathPhynder -s prepare -i ./Mafft_All_BEAST4.nwk -p taxa_pathphynder_tree -f ./branches.snp -r ./Cons_Aln_mafft_taxa_references.fa
```

PathPhynder command1: ### only transversions
```
pathPhynder -s all -t 100 -m transversions -i ./Mafft_All_BEAST4.nwk -p ./taxa_pathphynder_tree -l ./bamlist.txt -r ./Cons_Aln_mafft_taxa_references.fa
```

PathPhynder command2: ###transitions and transvertions--do not specify -m to extract transversions
```
pathPhynder -s all -t 100 -i ./Mafft_All_BEAST4.nwk -p ./taxa_pathphynder_tree -l ./bamlist.txt -r ./Cons_Aln_mafft_taxa_references.fa
```

## Population genomics
Based on pseudo-haploid genotypes:
(1) plot PCA to determine the relationships among populations.
(2) f-statistics is used to place the sample on the phylogenetic tree based on the distance among tips (leaves, samples).

prepare the environment and directory:
```
conda activate day3
mamba install r-tidyverse r-ggrepel
mkdir popgen
cd popgen
cp ~/course/data/day3/popgen/* .
```
To analyze SNP genotype data, PLINK format ('https://www.cog-genomics.org/plink/') is used
and the formats includes .bed, .bim and .fam ('https://www.cog-genomics.org/plink/1.9/formats')
here we have modern_polar_mexican.bed/bim/fam
```
plink --bfile modern_polar_mexican --missing --out modern_polar_mexican


$head modern_polar_mexican.imiss
        FID            IID MISS_PHENO   N_MISS   N_GENO   F_MISS
      Kenai         AAACGG          Y     1163   174734 0.006656
       West         AACTGA          Y     6055   174734  0.03465
       West         AAGACG          Y     1697   174734 0.009712
  Southwest         AAGCTA          Y      862   174734 0.004933
  Southwest         AATATC          Y     1937   174734  0.01109
       East         AATGAG          Y      768   174734 0.004395
...
```
From this amouont of missing genotype rate (F_MISS), assuming we could still be
able to obtain meaning for results.

(1) PCA analysis

make parameter file `modern_all.smartpca.par​​`
```
genotypename:    modern_polar_mexican.bed
snpname:         modern_polar_mexican.bim
indivname:       modern_polar_mexican.fam
evecoutname:     modern_all.evec
evaloutname:     modern_all.eval
familynames:     NO
numoutevec:      20
numthreads:	 2
numoutlieriter:	 0
poplistname:   	 modern_all.pops.txt
lsqproject:  YES
pordercheck:  NO
autoshrink: YES
```
`numoutevec` - the number of princiapl components to be returned

​`poplistname` - name of a file containing the population IDs for samples to be used to infer the principal components. When using PLINK format, population IDs for individuals have to be specified in column 6 of the .fam file. All individuals in the dataset with population ID not included in the file will be projected onto the inferred components

`lsqproject` - PCA projection algorithm appropriate for samples with high amount of missing data.

here `modern_all.pops.txt` specifies the popoulations included for calculating the coordinations in PCA. Samples in the genotype but not in this list will be projected on the map but will not influence the calculation or layout of other populations. This feature can be used to show the influence of outgroup to PCA.
```
smartpca -p modern_all.smartpca.par | tee modern_all.smartpca.log
#check log
modern_all.smartpca.log

#check output coordinates (eigenvalue results)
head modern_all.evec
```
prepare the file with a list of IDs for samples to be highlighted in the plot
```
cat label_inds.txt 
UE1210
UE1212
UE1605
```
visualized with
```
Rscript plot_pca.R modern_all.evec label_inds.txt modern_all.pdf
#modern_all.pdf is the output PCA plot
```

Examples as below:

a. with polar bear in the list for .par

<img width="608" alt="Screenshot 2023-10-06 at 14 20 31" src="https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/53645977/583167e9-bba4-4cf5-b004-b9c9861de535">

b. without
<img width="619" alt="Screenshot 2023-10-06 at 14 19 27" src="https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/53645977/fbd7cb11-aca6-4416-9bf3-d30081fec989">

This shows that polar bear as an outgroup can compress the difference among modern populations of bears. Also, there is an ancient sample slight biased towards polar bear in the PCA with all populations. This might be due to damage or outgroup attraction effect.



## (2) f-statistic

conda activate r
mamba install bioconductor-ggtree r-ape

use the f-statistic framework, implemented in the R package admixtools

then run `f_statistics.R` in interactive mode

For f3 statistic, there are 3 populations: A, B (sources) and C (target).
f3(A,B;C)=⟨(c−a)(c−b)⟩

`f3<0 significantly (f3 < -3)` C might be in between of A and B.

f3-statistic
<img width="498" alt="image" src="https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/53645977/b489cabb-9097-4f47-be04-4ccc179a1f64">

f4 is assuming Population((A,B),(X,D)), here X is query sample.

`f4>0` BD is more similar to each other.

`f4=0` The order is fine.

`f4<0` BC is more similar to each other.

f4-statistic
<img width="642" alt="image" src="https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/53645977/16b9c40c-9c98-4972-9cdd-c958c9a5da9e">

The f4 statistics as above can be explained by the graph (Pedersen, 2021).

<img width="437" alt="image" src="https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/53645977/50e2b786-4906-40b8-9275-c33991cf3432">




## Microbiome
The goal of this part is to tailor specifically the analysis for microbes. Here we focus in particular on bacteria, archea and viruses.

This is the pipeline we thought of prior to the course:

![image](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/4b27a043-8bf5-4196-8e0c-f00b3bab9151)

We were for instance reflecting on how to entangle the ancient signal and the ancient dormant/alive signal.

One of the questions for the next part of this analysis is to be able to see whether or not we can find microbes with damages in our samples.

### Extension of the reads, dereplication and mapping:
Preliminary statistics of our samples:
![image](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/3258a88c-220a-4e97-b34d-9dc452bbd48a)


These first steps are similar to what we have done for the eukaryotes. However here we extend the reads first.

####Description of the extension method:
* This method uses other kmer from the metagenomic data of the samples to extend existing reads. It is very conservative so it doesn't extend on damages and if it finds bubbles in the graph it stops (ie if there are several position on which the read could extend). The consequences of that is that the modern reads will be more extended than ancient ones in general. 
* On other consequence of extending the reads is that is allows to remove more duplicates as we made some reads artificially different at the adapter removal and trimming step.

For the extension we use tadpole loop through our samples:
```
for fq in $(ls /home/user-13/course/data/day2/fastq/*.fq.gz); do echo ${fq}; tadpole.sh -Xmx10G   k=17   in=${fq}   out=$(echo '$(basename '$fq')' | sed ‘s/\.fq\.gz$//’).extended.fq.gz   mode=extend   ibb=f   prefilter=0   el=100 er=100   threads=5   overwrite=true   trimends=9   ecc=f ecco=f   filtermem='${MEM}'   conservative=t   ignorebadquality; done 
```

To visualize the output statistics we use seqkit:
```
seqkit stats -j 5 -T data/*.fq.gz | csvtk -t pretty
```
Here are the statistics of the raw reads:
![image](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/471e3582-0340-41f6-ba65-1e5b33cec879)

Here are the statistics of the extended reads:
![image](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/dec00817-25c7-46a7-8758-b7670bb5b210)


* Dereplication: here we use seqkit which does the same job as vsearch but faster:
  ```
  seqkit rmdup -j 5 -s -o PRI-TJPGK-CATN-160-162.extended.derep.fq.gz PRI-TJPGK-CATN-160-162.extended.fq.gz
  ```
Statistitics of the extended dereplicated reads:
![image](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team_7/assets/111506710/b30d9cef-6f2f-4487-8c49-4e8748a17d0b)

#### Database:
The database includes mitochondrial and chloroplasts genomes from NCBI as proxies for bacterial genomes.
It is important to tailor the database to our target environment (ie. include sequences from Tara Oceans if marine).

Here the particularity is that the database is built with each genome of reference being concatenated in one fasta. It is both easier to sort (lower number of references) and faster for Bowtie2.

#### Mapping:
The mapping is similar here than what we did previously for the eukaryotes.
We map we the original reads fished with the extended deduplicated reads for several reasons including that the impact of extension on the damage has not been tested yet, + you loose the short reads authentification information which will impact the coverage and breadth of coverage.
We also use Bowtie2. Note: the parameters here are sat for teh tutorial but for a real analysis we need to thune the parameters (ie -L is the parameter to set the size of the seed, the shorter, the more sensitive, however the mapping will be slower).
```
bowtie2 -p 5 -k 100 -D 10 -R 2 \
    -N 0 -D 5 -R 1 -L 22 -i S,0,2.50 \
    --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.1" \
    -x data/db/aegenomics.db \
    -q PRI-TJPGK-CATN-160-162.mapping.fastq.gz --no-unal \
    | samtools view -F 4 -b \
    | samtools sort -@ 32 -m 8G -o PRI-TJPGK-CATN-160-162.sorted.bam
```
#### Filtering:
Here we use Samba to remove duplicates:
```
sambamba markdup -r -t 5 -p PRI-TJPGK-CATN-160-162.sorted.bam PRI-TJPGK-CATN-160-162.sorted.rmdup.bam
```
Then we filter the noise:
```
filterBAM --chunk-size 25 \
  --bam PRI-TJPGK-CATN-160-162.sorted.rmdup.bam \
  -N \
  -r data/misc/aegenomics.db.len.map \
  -A 92 \
  -a 94 \
  -n 100 \
  -b 0.75 \
  -B 0.01 \
  -t 5 \
  --sort-memory 8G \
  --include-low-detection \
  --stats PRI-TJPGK-CATN-160-162.stats.tsv.gz \
  --stats-filtered PRI-TJPGK-CATN-160-162.stats-filtered.tsv.gz \
  --bam-filtered PRI-TJPGK-CATN-160-162.filtered.bam
```
A few parameters are really important to have in mind.



Look at the results of filtering:
```
zcat PRI-TJPGK-CATN-160-162.stats-filtered.tsv.gz \
  | csvtk cut -t -T -f "reference,n_reads,read_ani_mean,read_ani_std,coverage_mean,breadth,exp_breadth,breadth_exp_ratio,norm_entropy,norm_gini,cov_evenness,tax_abund_tad" \
  | csvtk grep -r -t -v -f reference -p _plas -p _mito \
  | csvtk sort -t -T -k "n_reads:Nr" \
  | tabview -
```

These are 4 reference examples to explore our output:
```
printf "GCA_002781685.1\nIMGVR_UViG_3300027782_000260\nGCA_014380485.1" > ref-list.txt
getRPercId --bam PRI-TJPGK-CATN-160-162.sorted.rmdup.bam --reference-list ref-list.txt --threads 5 --sort-memory 8G
```
To see the coverage patterns:
```
/home/antonio/opt/conda/envs/day2/bin/bamcov -w 0 -m GCA_002781685.1.bam 
/home/antonio/opt/conda/envs/day2/bin/bamcov -w 0 -m GCA_014380485.1.bam
/home/antonio/opt/conda/envs/day2/bin/bamcov -w 0 -m IMGVR_UViG_3300027782_000260.bam 
/home/antonio/opt/conda/envs/day2/bin/bamcov -w 0 -m data/misc/example-4.bam
```

#### Damage:
The particularity here is that we look at the damages in local mode. The local mode allows to take into account the size of the genome of each reference.
```
metaDMG config --config-file PRI-TJPGK-CATN-160-162.local.yaml \
  --metaDMG-cpp /usr/local/bin/metaDMG-cpp \
  --parallel-samples 1 \
  --cores-per-sample 5 \
  --output-dir PRI-TJPGK-CATN-160-162.local \
  --max-position 35 \
  --min-similarity-score 0.92 \
  --damage-mode local \
  PRI-TJPGK-CATN-160-162.filtered.bam
```
Compute:
```
metaDMG compute PRI-TJPGK-CATN-160-162.local.yaml
```
Export to CSV:
```
metaDMG convert --add-fit-predictions \
  --output PRI-TJPGK-CATN-160-162.csv.gz \
  --results PRI-TJPGK-CATN-160-162.local/results
```

Explore:
```
metaDMG dashboard --results PRI-TJPGK-CATN-160-162.local/results/
```
Compare with the LCA mode:
```
conda activate metaDMG
metaDMG config --config-file PRI-TJPGK-CATN-160-162.lca.yaml \
  --custom-database \
  --names data/taxonomy/names.dmp \
  --nodes data/taxonomy/nodes.dmp \
  --acc2tax data/taxonomy/acc2taxid.map.gz \
  --metaDMG-cpp /usr/local/bin/metaDMG-cpp \
  --parallel-samples 1 \
  --cores-per-sample 5 \
  --output-dir PRI-TJPGK-CATN-160-162.lca \
  --max-position 35 \
  --lca-rank '' \
  --min-similarity-score 0.92 \
  --damage-mode lca \
  --weight-type 1 \
  PRI-TJPGK-CATN-160-162.filtered.bam
```


#### Fonctionnal profiling:

## References
[Reference_preprint_Antonio]: <https://www.biorxiv.org/content/10.1101/2023.06.10.544454v2.abstract> "Fernandez-Guerra et al., 2023, bioRxiv"
[ref_notdead_yet]: <https://bsapubs.onlinelibrary.wiley.com/doi/10.1002/ajb2.1780> "Sanyal et al., 2021, American Journal of Botany"

## Ancient pathogens:
The goal of this part is to apply our dataset to the search of pathogenes. In this case we focus on three different strings of the genus Yesinia:
*
*
*

### Duplicate removal:
Here we remove duplicates after the mapping bc it removes both exact duplicates and reads that starts at the same exact position and keeps the one with the best mapping quality. This strategy comes from Hominin/single species workflows. Samba seams to  work the same way than Picard but does a reverse engeneering than Picard. Samtools does only replicates removal.

The third genome has a lower breadth of coverage.

### Damage

### Coverage distribution:

Count/Coverage plot: way to vizualize how many positions are covered by 1-2-3... reads.
Y. pestis and Y. pseudotuberculosis have a mean of coverage at 8X with similar distribution. Y. pestis is a clone of Y. pseudotubercuosis.
Y. intermedia has most of the positions with low coverage.

Plasmids:
We don't find any reads for Y. pseudotuberculosis. But find find reads from all plasmids of Y. pestis.
Different number of copies of plasmids could come from the fact that they have different number of chromosomes.

Edit distance:
Y. intermedia maps at least with one mismatch.
Mapping of Y. pestis is slightly better than Y. pseudotuberculosis.

Coverage plot for genes of Y. pestis, we don't find genes that map to this particular plasmid.
