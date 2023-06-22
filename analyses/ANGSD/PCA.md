&nbsp;
## Principal component analysis
#### I am completing a principal component analysis using [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd) to better understand population structure. 
---
&nbsp;

### Steps in the command line: 
#### 1. Prepare input data by obtaining genotype likelihoods in Beagle format from BAM files
```
../../angsd/angsd -GL 2 -out ../../PCA/AR_genolike -nThreads 10 -ref CLRAindex/Rallus_crepitans_1.0.fasta -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam AR_list.txt;
../../angsd/angsd -GL 2 -out ../../PCA/NC_genolike -nThreads 10 -ref CLRAindex/Rallus_crepitans_1.0.fasta -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam NC_list.txt;
../../angsd/angsd -GL 2 -out ../../PCA/FL_genolike -nThreads 10 -ref CLRAindex/Rallus_crepitans_1.0.fasta -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam FL_list.txt;
../../angsd/angsd -GL 2 -out ../../PCA/OH_genolike -nThreads 10 -ref CLRAindex/Rallus_crepitans_1.0.fasta -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam OH_list.txt
```
#### 2. Run PCAngsd
##### &nbsp; 2.1 Install and build PCAngsd
```
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd
python setup.py build_ext --inplace
pip3 install -e .
```
##### &nbsp; 2.2 Perform principal component analysis on the genoltype liklihood beagle files
```
../angsd/pcangsd -b AR_genolike.beagle.gz -b FL_genolike.beagle.gz -b NC_genolike.beagle.gz -b OH_genolike.beagle.gz -sites_save -t 20 -o output.pcangsd
```
---
&nbsp;

### Steps in R
#### 1. Load packages
```{r}
library(ggplot2)
library(readr)
```

#### 2. Set working directory
```{r}
setwd("C:/Users/CarolPC/Documents")
```

#### 3. Data
```{r}
cov_matrix <- as.matrix(read.table("output.pcangsd.cov", header = TRUE))
```

#### 4. Visualize
```{r}

```

