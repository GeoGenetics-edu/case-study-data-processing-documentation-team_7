# Population genomics
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

This shows that polar bear as an outgroup can compress the difference among modern populations of bears. Also, there is an ancient sample slight biased towards polar bear in the PCA with all populations. This might be due to ???

damage
?bais



## (2) f-statistic






