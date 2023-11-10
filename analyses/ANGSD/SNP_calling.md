### SNP calling with ANGSD
> Single nucleotide polymorphisms (SNPs) are variations at a single DNA base pair position. SNP calling is the process of identifying SNPs. I am using ANGSD to [call SNPs based on allele frequencies](http://www.popgen.dk/angsd/index.php/SNP_calling).

&nbsp;
&nbsp;

#### 1. Likelihood ratio test (using a chi-square distribution with one degree of freedom) to ID SNPs

```
../../angsd/angsd -bam BamFileList.txt -out SNPs -doMaf 2 -GL 1 -SNP_pval 1e-6 -fai CLRAindex/Rallus_crepitans_1.0.fasta.fai -doMajorMinor 1
```
&nbsp;

---
&nbsp;

### Estimate the site frequency spectrum 

> SFS = the distribution of allele frequencies across sites in the genome of a species

&nbsp;

```
angsd/angsd -b BamFileList -anc CLRAindex/Rallus_crepitans_1.0.fasta -gl 2 -doSaf 1 -doMaf 1 -SNP_pval 0.05 -minind 5 -nthreads 8 -doMajorMinor 1 -out SFS
```
- saf files are generated using `/angsd -doSaf`. 
- .saf = site allele frequency likelihood file
- `doSaf 1` = Calculate the site allele frequency likelihood based on individual genotype likelihoods assuming HWE

