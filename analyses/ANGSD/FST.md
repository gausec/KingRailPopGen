
### Aim: 
#### I will be using [ANGSD to estimate F<sub>ST</sub>](http://www.popgen.dk/angsd/index.php/Fst) between populations. [F<sub>ST</sub>](https://www.nature.com/articles/nrg2611) is a measurment of genetic variaton *within* populations reative to genetic variation *between* populations.  F<sub>ST</sub> is based on heterozygosity (a simple measure of genetic variation, 1 - the sum of the allele frequencies squared) and is reported on a scale of 0-1 (0 meaning there is 0% variation between populations and 100% variation within populations.  F<sub>ST</sub>=1  would mean there is 100% variation between and 0% among).
&nbsp;
---
### Steps:  
#### 1. SFS estimation for each population
&nbsp; 1.1 Create a text file list of bam files for each population
- NC_list.txt
- OH_list.txt
- AR_list.txt
- FL_list.txt  

&nbsp; *These can be found in the data subdirectory.*
      
&nbsp;

&nbsp; 1.2 Generate site allele frequency likelihoods in SAF format for each population.
```
../../angsd/angsd -b NC_list.txt -anc CLRAindex/Rallus_crepitans_1.0.fasta -gl 2 -doSaf 1 -doMaf 1 -SNP_pval 0.05 -minind 5 -nthreads 8 -doMajorMinor 1 -out FST/SFS.NC

```

&nbsp; *Repeat for the other 3 populations.*

&nbsp;

---

#### 2. Use [RealSFS](http://www.popgen.dk/angsd/index.php/RealSFS) to calculate 2d sfs for each pair of populations
```
../../../angsd/misc/realSFS AR/SFS.gl2.AR.saf.idx NC/SFS.NC.saf.idx -P 30 > AR_NC.ml;
../../../angsd/misc/realSFS AR/SFS.gl2.AR.saf.idx OH/SFS.OH.saf.idx -P 30 > AR_OH.ml; 
../../../angsd/misc/realSFS AR/SFS.gl2.AR.saf.idx FL/SFS.FL.saf.idx -P 30 > AR_FL.ml; 
../../../angsd/misc/realSFS NC/SFS.NC.saf.idx OH/SFS.OH.saf.idx -P 30 > NC_OH.ml; 
../../../angsd/misc/realSFS NC/SFS.NC.saf.idx FL/SFS.FL.saf.idx -P 30 > NC_FL.ml; 
../../../angsd/misc/realSFS OH/SFS.OH.saf.idx FL/SFS.FL.saf.idx -P 30 > OH_FL.ml

```

---
#### 3.  For computing the pairwise fst, ANGSD reccomends performing pairwise FST calculations for each population separately. Generate separate `.fst` files for each pairwise population comparison
&nbsp; 
```
../../../angsd/misc/realSFS fst index AR/SFS.gl2.AR.saf.idx NC/SFS.NC.saf.idx -sfs AR_NC.ml -fstout FstOut.AR_NC;
../../../angsd/misc/realSFS fst index AR/SFS.gl2.AR.saf.idx OH/SFS.OH.saf.idx -sfs AR_OH.ml -fstout FstOut.AR_OH;
../../../angsd/misc/realSFS fst index AR/SFS.gl2.AR.saf.idx FL/SFS.FL.saf.idx -sfs AR_FL.ml -fstout FstOut.AR_FL;
../../../angsd/misc/realSFS fst index NC/SFS.NC.saf.idx OH/SFS.OH.saf.idx -sfs NC_OH.ml -fstout FstOut.NC_OH;
../../../angsd/misc/realSFS fst index NC/SFS.NC.saf.idx FL/SFS.FL.saf.idx -sfs NC_FL.ml -fstout FstOut.NC_FL;
../../../angsd/misc/realSFS fst index OH/SFS.OH.saf.idx FL/SFS.FL.saf.idx -sfs OH_FL.ml -fstout FstOut.OH_FL

```
- *These files contain the pairwise FST values.*

---
#### 4. Use realSFS to extract the the fst values from the fst binary files
&nbsp;
```
../../../angsd/misc/realSFS fst stats FstOut.AR_NC.fst.idx > FstOut.AR_NC.fst;
../../../angsd/misc/realSFS fst stats FstOut.AR_OH.fst.idx > FstOut.AR_OH.fst;
../../../angsd/misc/realSFS fst stats FstOut.AR_FL.fst.idx > FstOut.AR_FL.fst;
../../../angsd/misc/realSFS fst stats FstOut.NC_OH.fst.idx > FstOut.NC_OH.fst;
../../../angsd/misc/realSFS fst stats FstOut.NC_FL.fst.idx > FstOut.NC_FL.fst;
../../../angsd/misc/realSFS fst stats FstOut.OH_FL.fst.idx > FstOut.OH_FL.fst

```

---
### Visualization in R
#### 5. Load packages
```
library(hierfstat) #FST statistics
library("ggplot2") #plotting
library("reshape2") #plotting
library("adegenet") #FST statistics and data storage
```
#### 6. Read in pairwise *.fst output files
```
fst_AR_FL <- read.table("FstOut.AR_FL.fst", header = FALSE)
fst_AR_NC <- read.table("FstOut.AR_NC.fst", header = FALSE)
fst_AR_OH <- read.table("FstOut.AR_OH.fst", header = FALSE)
fst_NC_FL <- read.table("FstOut.NC_FL.fst", header = FALSE)
fst_NC_OH <- read.table("FstOut.NC_OH.fst", header = FALSE)
fst_OH_FL <- read.table("FstOut.OH_FL.fst", header = FALSE)
```
#### 7. Combine the values from .fst files into a single matrix
```
fst_AR_FL <- fst_AR_FL[, -2]
fst_AR_NC <- fst_AR_NC[, -2]
fst_AR_OH <- fst_AR_OH[, -2]
fst_NC_FL <- fst_NC_FL[, -2]
fst_NC_OH <- fst_NC_OH[, -2]
fst_OH_FL <- fst_OH_FL[, -2]

```

# Create a matrix with appropriate dimensions
```{r}
populations <- c("AR", "FL", "NC", "OH")
combined_matrix <- matrix(0, nrow = length(populations), ncol = length(populations))
rownames(combined_matrix) <- populations
colnames(combined_matrix) <- populations
```
