&nbsp;
## Principal component analysis
#### I am completing a principal component analysis using [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd) to better understand population structure. I am using one [CLRA](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR23269683&display=download) in my analysis.
---
&nbsp;

### Steps in the command line: 
#### 1. Prepare input data by obtaining genotype likelihoods in Beagle format from BAM files
##### &nbsp; 1.1 I have a file containing sample names and their population location. Each field is tab delimited. I would like to check that I have the right number for each population.
```
cut -f 1 -d $'\t' PopInfo.txt | sort | uniq -c
```
##### &nbsp; 1.2 Obtain genotype likelihoods in Beagle format
```
angsd/angsd -bam BamFileList.txt -ref CLRAindex/Rallus_crepitans_1.0.fasta -GL 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -minind 5 -SNP_pval 1e-6 -minMapQ 30 -minQ 20  -doGlf 2 -out input -nThreads 8
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
 pcangsd -b input.beagle.gz --sites_save -t 20 -o output.pcangsd
```
---
&nbsp;

### Steps in R:
#### 1. Load packages
```{r}
library(ggplot2)
library(readr)
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
##### &nbsp; 5.1 Extract the first two eigenvectors
```{r}
eigenvectors <- as.matrix(e$vectors[, 1:2]) 
```

##### &nbsp; 5.2 Transformed data
```{r}
transformed_data <- as.matrix(cov_matrix) %*% eigenvectors
```

##### &nbsp; 5.3 Combine transformed data with population information
```{r}
pca.vectors <- data.frame(pop, V1 = transformed_data[, 1], V2 = transformed_data[, 2])
```
##### &nbsp; 5.4 Plot PCA
```{r}
pca <- ggplot(data = pca.vectors, aes(x = V1, y = V2, colour = Location, label = Sample_ID)) + 
  geom_point() + 
  labs(title = "Individual Allele Frequency", x = "PC1", y = "PC2") + 
  theme(plot.title = element_text(face = "bold"))

plot(pca)
```
##### &nbsp; 5.5 Save plot as a .png file
```{r}
ggsave("pca_plot.png", plot = pca, width = 8, height = 6, dpi = 300)
```

---

### DAPC
> I will be using DAPC to further investigate structure between populations using R package [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html).

&nbsp;
 
Helpful tutorial: https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
&nbsp;

#### 1. Cross validation to ID PCs
```
pc_range <- 40:70
xval_results <- list()

for (num_pcs in pc_range) {
  # Perform cross-validation with num_pcs
  xval_result <- xvalDapc(cov_matrix, grp = pca.vectors$Location, n.pca.max = num_pcs, n.rep = 30)
  xval_results[[as.character(num_pcs)]] <- xval_result
}
```
&nbsp;

#### 2. Extract
```
xval_result$`Number of PCs Achieving Lowest MSE`
```
&nbsp;

---
### Neighbour-joining tree

&nbsp;

#### Construct neighbour-joining tree of samples from estimated covariance matrix estimated based on indivdual allele frequencies.


```
 pcangsd -b input.beagle.gz --tree --tree_sample PopInfo.csv -o output.tree
```
--tree
--tree_samples


