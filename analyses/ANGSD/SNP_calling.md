## SNP calling with ANGSD
> Single nucleotide polymorphisms (SNPs) are variations at a single DNA base pair position. SNP calling is the process of identifying SNPs. I am using ANGSD to [call SNPs based on allele frequencies](http://www.popgen.dk/angsd/index.php/SNP_calling).

&nbsp;
&nbsp;

#### 1. Likelihood ratio test (using a chi-square distribution with one degree of freedom) to ID SNPs

```
../../angsd/angsd -bam BamFileList.txt -out SNPs -doMaf 2 -GL 1 -SNP_pval 1e-6 -fai CLRAindex/Rallus_crepitans_1.0.fasta.fai -doMajorMinor 1
```
