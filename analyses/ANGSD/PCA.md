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
angsd/angsd -bam BamFileList.txt -ref CLRAindex/Rallus_crepitans_1.0.fasta -GL 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -minind 5 -SNP_pval 1e-6 -minMapQ 30 -minQ 20  -doGlf 2 -out input -nThreads 8
```
- `-GL 1`: SAMtools model: GATK assumes sequencing errors are independent, while Samtools believes the second error (unrelated sequencing error that occurs independently of the first one and affects the same genomic position) comes at a higher chance, especially at high depths.

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
 pcangsd -b input.beagle.gz -t 20 -o output.pcangsd
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
##### &nbsp; 5.4 Plot PCA using the [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) package
```{r}
pca <- ggplot(data = pca.vectors, aes(x = V1, y = V2, colour = Location, label = Sample_ID)) + # What data to plot, color-coding, and legend
  geom_point() + # data are represented as points 
  labs ( x = "PC1 (%)", y = "PC2 (%)") + # axis labels 
  theme_minimal() + theme(axis.title.x = element_text(margin = margin(t = 20))) + # Change the x-axis margin spacing 
  theme(axis.title.y = element_text(margin = margin(r = 20))) + # Change y-axis margin spacing  
  theme(text = element_text(size = 20))  # Increase the text size


plot(pca)
```

##### &nbsp; 5.5 Sometimes, I want to label points by their sample ID. I can do that using the [ggrepel](https://ggrepel.slowkow.com/) package.
```
library(ggrepel)
# sometimes I only want to look at specific points to see where they are in relation to other points. I can label just those:
sample_ids_to_label <- c("NC48827.downsampled.bam", "NC43642.downsampled.bam", "43623.downsampled.bam") # example IDs

# subset to include only the rows with the specified samples
samples_to_label <- subset(pca.vectors, Sample_ID %in% sample_ids_to_label)

# add to the PCA plot
pca +  geom_label_repel(data = samples_to_label, aes(label = Sample_ID),
                   box.padding = 0.55, point.padding = 0.05, segment.color = "grey50", 
                   max.overlaps = 10, show.legend = FALSE) # setting show.legend to false gets rid of those annoying "a" letters that ggrepel sometimes puts in your legend!

# Adjust max.overlaps as needed to see the labels you need
```

##### &nbsp; 5.6 Save plot as a .png file
```{r}
ggsave("pca_plot.png", plot = pca, width = 8, height = 6, dpi = 300)

# another option:
# png(filename = "Just_KIRA_newgenomes.png", width = 1000, height = 800)
# Plot(pca)
# dev.off()
```

---

### DAPC
> I will be using DAPC to further investigate structure between populations using R package [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html).

&nbsp;
 
Helpful tutorial: https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
&nbsp;

#### 1. Cross validation to ID PCs

##### &nbsp; 1.1 cross validation
```
pc_range <- 5:30
xval_results <- list()

for (num_pcs in pc_range) {
  # cross-validation with num_pcs
  xval_result <- xvalDapc(cov_matrix, grp = pca.vectors$Source, n.pca.max = num_pcs, n.rep = 1000)
  xval_results[[as.character(num_pcs)]] <- xval_result
}
```
##### &nbsp; 1.2 summarize
```
xval_result$DAPC
````
##### &nbsp; 1.3 extract number of PCs to conserve
```
xval_result$`Number of PCs Achieving Lowest MSE`
```

#### 2. DAPC  
```
dapc <- dapc(cov_matrix, grp = pca.vectors$Source, n.pca = 25, n.da = 5)
# n.pca is the number of PCs achieving lowest MSE
# n.da is the number of distinct groups
```
#### 3. Plot

```
scatter(dapc, scree.da = FALSE, bg = "white", cex = 1.15, cstar = 0, col = custom_colors, clab = 0, leg = TRUE,  pch = 20) +
 stat_ellipse() # confidence elipses

# adding an inset plot for the cumulative variance 
myInset <- function(){
temp <- dapc$pca.eig
temp <- 100* cumsum(temp)/sum(temp)
plot(temp, col=rep(c("black","lightgrey"),
c(dapc$n.pca,1000)), ylim=c(0,100),
xlab="PCA eigenvalues", ylab="Cumulated variance (%)",
cex=1, pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="bottomleft",
inset=c(-0.05,-0.12), ratio=.28,
bg=transp("white"))
```

&nbsp;

---
### Neighbour-joining tree

&nbsp;

#### Construct neighbor-joining tree of samples from estimated covariance matrix estimated based on indivdual allele frequencies.


```
 pcangsd -b input.beagle.gz --tree --tree_sample PopInfo.csv --filterSites sites.txt -t 20 -o output.tree
```

&nbsp;

#### Plot in R

```
library(ape)
```
```
tree <- read.tree("output.tree.tree")
```
```
pdf("tree.pdf")
plot(tree, type = "fan", cex=0.25)
dev.off()
```

