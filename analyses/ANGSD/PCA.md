&nbsp;
## Principal component analysis
#### I am completing a principal component analysis using [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd) to better understand population structure. 
---
&nbsp;

### Steps in the command line: 
#### 1. Prepare input data by obtaining genotype likelihoods in Beagle format from BAM files
##### &nbsp; 1.1 I have a BAM file list containing samples and their population location. I would like to check that I have the right number for each population.
```
cut -f 1 -d " " PopInfo | sort | uniq -c
```
##### &nbsp; 1.2 Obtain genotype likelihoods in Beagle format
```
../../angsd/angsd -GL 2 -out ../../PCA/Agenolike -nThreads 10 -ref CLRAindex/Rallus_crepitans_1.0.fasta -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam BamFileList;
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
pcangsd -b genolike.beagle.gz --sites_save -t 20 -o output.pcangsd
```
---
&nbsp;

### Steps in R
#### 1. Load packages
```{r}
library(ggplot2)
library(readr)
library(ggcorrplot)
library(FactoMineR)
library(RColorBrewer)
```

#### 2. Set working directory
```{r}
setwd("C:/Users/CarolPC/Documents")
```

#### 3. Data
##### &nbsp; 3.1 Read in PCANGSD output covariance matrix
```{r}
cov_matrix <- as.matrix(read.table("output.pcangsd.cov", header = FALSE))
```
##### &nbsp; 3.2 Read in population info and assign column headers
```{r}
pop<-read.csv("PopInfo.csv", header = FALSE)
colnames(pop) <- c("Sample_ID","Location")
```
#### 4. ID principal components. Compute the eigenvalues from covariance matrix with `eigen`
```{r}
e<-eigen(cov_matrix)
# eigenvalue -> a number representing the amount of variance present in the data for a given direction.
```

#### 5. Plot
```{r}

```

