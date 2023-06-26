
# Steps to plot PCANGSD results in R. Let's explore genetic variance! 


# 1. Load packages
library(ggplot2)
library(readr)
library(RColorBrewer)


# 2. Set your working directory. Example:
setwd("C:/Users/PC/Documents")


# 3. Data
# 3.1 Read in PCANGSD output covariance matrix

cov_matrix <- as.matrix(read.table("output.pcangsd.cov", header = FALSE))

# 3.2 Read in population info and assign column headers. Column 1 = your sample names, column 2 = location name
pop<-read.csv("PopInfo.csv", header = FALSE)
colnames(pop) <- c("Sample_ID","Location")

# 4. ID principal components. Compute the eigenvalues from covariance matrix with `eigen`
e<-eigen(cov_matrix)
# eigenvalue -> a number representing the amount of variance present in the data for a given direction.

# 5. Plot
# 5.1 Extract eigenvectors
eigenvectors <- as.matrix(e$vectors[, 1:2])  # Extract the first two eigenvectors


# 5.2 Transformed data
transformed_data <- as.matrix(cov_matrix) %*% eigenvectors

# 5.3 Combine transformed data with population information
pca.vectors <- data.frame(pop, V1 = transformed_data[, 1], V2 = transformed_data[, 2])

# 5.4 Plot PCA
pca <- ggplot(data = pca.vectors, aes(x = V1, y = V2, colour = Location, label = Sample_ID)) + 
  geom_point() + 
  labs(title = "Individual Allele Frequency", x = "PC1", y = "PC2") + 
  theme(plot.title = element_text(face = "bold"))

plot(pca)
